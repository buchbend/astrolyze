"""The :class:`Collection` facade: open a published corpus and browse it (#57).

A :class:`Collection` is the consumer-facing handle on a published corpus — a directory of
astrolyze-flavoured Zarr stores indexed by one ``catalog.parquet`` at its root. It is **read
only**: it opens stores into lazy Cubes and never writes to the corpus (the bridge to a writable
:class:`~astrolyze.experiment.Experiment` is the origin Metadata each Cube carries, not a shared
mutation path — PRD #56). The catalog-format knowledge lives entirely in
:mod:`astrolyze.collection.catalog`; this layer adds the *presentation* (object-first grouping)
and the *open* path (record -> Cube with origin provenance).

Path handling routes through fsspec from day one (:func:`fsspec.core.url_to_fs`): the same
``Collection.open(...)`` call shape serves a local directory today and an ``s3://`` URL later
with no API change. Today the open path materialises the local store path the Zarr backend
expects; the fsspec-URL-end-to-end work (a store mapper for remote bytes) is #63, and the seam is
:meth:`Record.open` — it is the single place a remote-store opener slots in.

Extension points for the slices that build on this tracer:

- **#60 describe/query** add ``Collection.describe(object)`` / ``query(**filters)`` over the same
  :attr:`records`; the typed :class:`~astrolyze.collection.catalog.CatalogRow` is the filter axis.
- **#62 covering** adds ``Collection.covering(SkyCoord)`` reading the footprint columns
  (``ra_deg`` / ``dec_deg`` / ``radius_deg``) already parsed on each row, then deciding final
  containment from the candidate store's own WCS at open time.
- **#61 scan-builder** adds an alternate constructor (``Collection.from_store_dir``) building a
  :class:`~astrolyze.collection.catalog.Catalog` by scanning stores — it reuses this same facade.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from fsspec.core import url_to_fs

from .catalog import Catalog, CatalogRow, read_catalog


@dataclass(frozen=True)
class ObjectSummary:
    """The object-first overview row for one source (what :meth:`Collection.list` returns).

    One per distinct ``object`` in the catalog, aggregating that source's stores: the surveys and
    species it appears in, how many stores cover it, and its beam range (coarsest to finest major
    axis). This is a *consumer* presentation concern — the catalog stores flat per-store rows and
    the collection groups them (the spec keeps grouping out of storage)."""

    # ``object`` may be None for an unnamed source — a partial row, like the catalog itself.
    object: str | None
    surveys: tuple[str, ...]
    species: tuple[str, ...]
    n_stores: int
    beam_range_arcsec: tuple[float | None, float | None]


@dataclass(frozen=True)
class Record:
    """A single catalog entry bound to its corpus context — a record that can :meth:`open`.

    Wraps a typed :class:`~astrolyze.collection.catalog.CatalogRow` together with everything
    needed to resolve its store: the corpus root URI and the catalog version. The identity
    accessors (``object`` / ``survey`` / ``species`` / ``transition``) read straight off the row,
    so the #60 query layer can filter a list of records without unwrapping. :meth:`open` is the
    seam where the Zarr store becomes a lazy Cube stamped with origin provenance."""

    row: CatalogRow
    root_uri: str
    catalog_version: str

    @property
    def object(self) -> str | None:
        return self.row.object

    @property
    def survey(self) -> str | None:
        return self.row.survey

    @property
    def species(self) -> str | None:
        return self.row.species

    @property
    def transition(self) -> str | None:
        return self.row.transition

    @property
    def store_path(self) -> str:
        """The store's path relative to the corpus root (POSIX), straight from the catalog."""
        return self.row.store_path

    @property
    def store_uri(self) -> str:
        """The store's absolute fsspec URI (root URI + relative ``store_path``).

        The provenance anchor stamped onto an opened cube and the locator a remote opener (#63)
        will hand to fsspec. Relative-to-root joining is the contract that keeps a corpus
        relocatable (move the root, the URIs follow); ``store_path`` is relative by contract, so
        a stray leading slash is stripped rather than escaping the root."""
        return f"{self.root_uri.rstrip('/')}/{self.row.store_path.lstrip('/')}"

    def open(self):
        """Open this record's Zarr store into a lazy, dask-backed :class:`~astrolyze.core.Cube`.

        The store opens through the existing Zarr backend (``Cube.from_zarr``), so the cube is
        dask-backed and context-carrying exactly as a directly-loaded store is (#23). On top of
        that, the cube's :class:`~astrolyze.io.Metadata` is stamped with its **origin** — the
        source store URI and the catalog version — so a product later saved into an experiment
        traces back to the exact corpus snapshot it came from (PRD #56). This stamping is the
        only thing the collection adds over a bare open; it never mutates the corpus.

        The remote-store path (an ``s3://`` URI opened via an fsspec mapper) is #63; this is the
        single seam it slots into. Today the store opens from the resolved local path."""
        from astrolyze.core import Cube

        cube = Cube.from_zarr(_local_store_path(self.store_uri))
        cube.metadata = replace(
            cube.metadata,
            origin_store_uri=self.store_uri,
            origin_catalog_version=self.catalog_version,
        )
        return cube


class Collection:
    """A read-only handle on a published corpus: open it, list it, open its records.

    Built by :meth:`open`; holds the parsed :class:`~astrolyze.collection.catalog.Catalog` and
    the resolved corpus root URI. :attr:`records` is the flat per-store view (a :class:`Record`
    each), and :meth:`list` is the grouped object-first overview. The corpus is never written
    to — lineage flows out through each opened cube's origin Metadata, not back in."""

    def __init__(self, catalog: Catalog, root_uri: str):
        self._catalog = catalog
        self._root_uri = root_uri

    @classmethod
    def open(cls, path) -> "Collection":
        """Open the published corpus at *path* (a local directory today, an ``s3://`` URL later).

        Path handling goes through fsspec (:func:`fsspec.core.url_to_fs`), so the call shape is
        identical for both — the local case resolves to a ``file://`` root, an object-store case
        to its scheme. The catalog is read + schema-validated by
        :func:`~astrolyze.collection.catalog.read_catalog`; an unknown catalog MAJOR raises there
        (:class:`~astrolyze.collection.catalog.CatalogSchemaError`) rather than mis-reading."""
        fs, resolved = url_to_fs(str(path))
        root_uri = fs.unstrip_protocol(resolved)
        catalog = read_catalog(root_uri)
        return cls(catalog, root_uri)

    @property
    def root_uri(self) -> str:
        """The corpus root as an fsspec URI (the prefix every store URI is relative to)."""
        return self._root_uri

    @property
    def catalog_version(self) -> str:
        """The ``catalog_schema_version`` the opened catalog declared."""
        return self._catalog.schema_version

    @property
    def records(self) -> list[Record]:
        """Every published store as a :class:`Record`, in catalog order (the flat per-store view).

        The axis the #60 query layer filters and #62 covering prefilters; each record opens into a
        lazy Cube via :meth:`Record.open`."""
        return [
            Record(
                row=row,
                root_uri=self._root_uri,
                catalog_version=self._catalog.schema_version,
            )
            for row in self._catalog.rows
        ]

    def list(self) -> list[ObjectSummary]:
        """The object-first overview: one :class:`ObjectSummary` per source.

        Groups the flat catalog rows by ``object`` and summarises each group — the surveys and
        species it spans, its store count, and its beam range (coarsest to finest major axis) — so
        a researcher sees at a glance what data exists for each source (PRD #56 user story 2).
        Sorted by object name for a stable rendering."""
        grouped: dict[str, list[CatalogRow]] = {}
        for row in self._catalog.rows:
            grouped.setdefault(row.object, []).append(row)

        summaries = []
        for obj in sorted(grouped, key=lambda name: (name is None, name)):
            rows = grouped[obj]
            summaries.append(
                ObjectSummary(
                    object=obj,
                    surveys=_distinct(r.survey for r in rows),
                    species=_distinct(r.species for r in rows),
                    n_stores=len(rows),
                    beam_range_arcsec=_beam_range(rows),
                )
            )
        return summaries


# -- aggregation helpers ---------------------------------------------------------------
def _distinct(values) -> tuple:
    """Distinct non-null values, sorted — a stable summary of a categorical column."""
    return tuple(sorted({v for v in values if v is not None}))


def _beam_range(rows) -> tuple[float | None, float | None]:
    """The (min, max) beam major axis across *rows* in arcsec (``(None, None)`` if none stated).

    The coarsest-to-finest span a source covers — a single-dish 13" beam and an interferometer
    1.5" beam in the same source read as ``(1.5, 13.0)``."""
    majors = [r.beam_major_arcsec for r in rows if r.beam_major_arcsec is not None]
    if not majors:
        return (None, None)
    return (min(majors), max(majors))


def _local_store_path(store_uri: str) -> str:
    """The local filesystem path a ``file://`` store URI names (the Zarr backend opens a path).

    The narrow open path the tracer needs (local corpora today); a remote (``s3://``) opener that
    hands an fsspec mapper to the Zarr backend is #63, slotting in at :meth:`Record.open`. Uses
    fsspec to strip the protocol so the resolution is the same machinery the open went through."""
    fs, resolved = url_to_fs(store_uri)
    return resolved


__all__ = ["Collection", "Record", "ObjectSummary"]
