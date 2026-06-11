"""The :class:`Collection` facade: open a published corpus and browse it (#57).

A :class:`Collection` is the consumer-facing handle on a published corpus — a directory of
astrolyze-flavoured Zarr stores indexed by one ``catalog.parquet`` at its root. It is **read
only**: it opens stores into lazy Cubes and never writes to the corpus (the bridge to a writable
:class:`~astrolyze.experiment.Experiment` is the origin Metadata each Cube carries, not a shared
mutation path — PRD #56). The catalog-format knowledge lives entirely in
:mod:`astrolyze.collection.catalog`; this layer adds the *presentation* (object-first grouping)
and the *open* path (record -> Cube with origin provenance).

Path handling routes through fsspec from day one (:func:`fsspec.core.url_to_fs`): the same
``Collection.open(...)`` call shape serves a local directory, a ``file://`` URL, and an ``s3://``
URL with no API change. :meth:`Record.open` hands the store's fsspec URI straight to the Zarr
backend, which opens local and remote stores through one seam and stays lazy over a remote store
(#63) — so a corpus on disk and the same corpus on S3 are analysed with identical code.

Extension points for the slices that build on this tracer:

- **#60 describe/query** add ``Collection.describe(object)`` / ``query(**filters)`` over the same
  :attr:`records`; the typed :class:`~astrolyze.collection.catalog.CatalogRow` is the filter axis.
- **#62 covering** adds ``Collection.covering(SkyCoord)`` reading the footprint columns
  (``ra_deg`` / ``dec_deg`` / ``radius_deg``) already parsed on each row, then deciding final
  containment from the candidate store's own WCS at open time.
- **#61 scan-builder** wires :meth:`Collection.open` to fall back to
  :func:`~astrolyze.collection.scan.build_catalog` when the root carries no ``catalog.parquet``,
  building an equivalent :class:`~astrolyze.collection.catalog.Catalog` by scanning the stores —
  it reuses this same facade, the trigger being the index file's presence.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace

from fsspec.core import url_to_fs

from .catalog import CATALOG_FILENAME, Catalog, CatalogRow, read_catalog


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
class StoreDetail:
    """One store's physical detail row (what :meth:`Collection.describe` returns, one per store).

    Built **catalog-row-first**: every field except the three ``*_kms`` velocity-axis fields is
    read straight off the typed :class:`~astrolyze.collection.catalog.CatalogRow`, so describing a
    source costs no store read (cheap on a remote corpus — PRD #56 user story 3). The velocity-axis
    fields the catalog does not carry — native channel width and velocity coverage — are filled
    **only on request** (``describe(object, deep=True)``), which opens each store and reads its own
    spectral axis; on a shallow describe they stay ``None`` (the field is unknown, never guessed —
    ADR-0003). Velocities are in km/s, the analysis-side convention."""

    object: str | None
    survey: str | None
    telescope: str | None
    species: str | None
    transition: str | None
    rest_frequency_hz: float | None
    beam_major_arcsec: float | None
    beam_minor_arcsec: float | None
    beam_pa_deg: float | None
    bunit: str | None
    store_path: str
    store_uri: str
    content_checksum: str | None
    # Deep-only (filled when describe(deep=True) opens the live store; None otherwise):
    channel_width_kms: float | None = None
    velocity_min_kms: float | None = None
    velocity_max_kms: float | None = None


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

        The store URI is handed straight to the Zarr backend, which opens a local path and a
        remote fsspec URL (``file://`` / ``memory://`` / ``s3://``) through the same seam and
        stays lazy over a remote store (#63) — so a corpus on disk and the same corpus on S3 open
        with no code change."""
        from astrolyze.core import Cube

        cube = Cube.from_zarr(self.store_uri)
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
        """Open the corpus at *path* (a local directory today, an ``s3://`` URL later).

        Path handling goes through fsspec (:func:`fsspec.core.url_to_fs`), so the call shape is
        identical for both — the local case resolves to a ``file://`` root, an object-store case
        to its scheme.

        **The trigger is the presence of ``catalog.parquet``.** When the root carries the published
        index it is read + schema-validated by
        :func:`~astrolyze.collection.catalog.read_catalog` (an unknown catalog MAJOR raises there,
        :class:`~astrolyze.collection.catalog.CatalogSchemaError`, rather than mis-reading). When it
        does **not**, the collection transparently falls back to the scan-builder (#61): it
        reconstructs an equivalent catalog by reading each store's attrs and computing its footprint
        from its own WCS (:func:`~astrolyze.collection.scan.build_catalog`), so a bare directory of
        self-describing Zarr stores is browsable with zero extra tooling (PRD #56 user story 29).
        The scan is read-only — it builds an in-memory catalog and never writes one into the corpus
        (persist explicitly with :func:`~astrolyze.collection.scan.write_catalog`)."""
        fs, resolved = url_to_fs(str(path))
        root_uri = fs.unstrip_protocol(resolved)
        if fs.exists(f"{resolved.rstrip('/')}/{CATALOG_FILENAME}"):
            catalog = read_catalog(root_uri)
        elif not fs.exists(resolved):
            # A root that exists-but-lacks-a-catalog is a scan target (below); a root that does not
            # exist at all is a user error (typo'd path, wrong bucket) — fail clearly here rather
            # than letting the scan silently return an empty collection (#61/#63 seam).
            raise FileNotFoundError(
                f"no corpus at {root_uri}: it carries neither a {CATALOG_FILENAME} "
                "nor a directory of Zarr stores to scan"
            )
        else:
            # Catalog-less directory: build one by scanning the stores (the #61 fallback). The scan
            # stays a local-directory path today; the remote-store scan rides the same #63 seam as
            # the remote open path.
            from .scan import build_catalog

            catalog = build_catalog(resolved)
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

    def describe(self, object: str, *, deep: bool = False) -> list[StoreDetail]:
        """Expand one source into per-store physical detail: one :class:`StoreDetail` per store.

        The drill-down behind :meth:`list` (PRD #56 user story 3): where ``list`` summarises a
        source object-first, ``describe`` opens it back up — every store of *object* with its beam
        (major/minor/PA), bunit, rest frequency, transition, and provenance summary
        (survey/telescope/checksum), so a researcher can judge dataset suitability before loading.

        **The answer comes from the catalog row first.** With ``deep=False`` (the default) nothing
        is opened: every field is read off the typed catalog row, so describing a source on a
        remote corpus stays cheap (no store read — the acceptance contract). The velocity-axis
        fields the catalog does not carry (native channel width, velocity coverage) are filled
        **only on request**: ``deep=True`` opens each store's spectral axis (still lazy — only the
        1-D axis is read, not the cube) and reads its channel width and min/max velocity in km/s.
        On a shallow describe those three fields stay ``None`` (unknown, never guessed — ADR-0003).

        An *object* with no store raises a descriptive :class:`KeyError` naming the known objects,
        rather than returning an empty list (an empty result would silently hide a typo)."""
        rows = [row for row in self._catalog.rows if row.object == object]
        if not rows:
            known = sorted(
                {r.object for r in self._catalog.rows if r.object is not None}
            )
            raise KeyError(
                f"no source {object!r} in this collection; known objects: "
                f"{', '.join(known) if known else '(none)'}"
            )
        details = []
        for row in rows:
            record = Record(
                row=row,
                root_uri=self._root_uri,
                catalog_version=self._catalog.schema_version,
            )
            spectral = _spectral_axis_kms(record) if deep else (None, None, None)
            details.append(
                StoreDetail(
                    object=row.object,
                    survey=row.survey,
                    telescope=row.telescope,
                    species=row.species,
                    transition=row.transition,
                    rest_frequency_hz=row.rest_frequency_hz,
                    beam_major_arcsec=row.beam_major_arcsec,
                    beam_minor_arcsec=row.beam_minor_arcsec,
                    beam_pa_deg=row.beam_pa_deg,
                    bunit=row.bunit,
                    store_path=row.store_path,
                    store_uri=record.store_uri,
                    content_checksum=row.content_checksum,
                    channel_width_kms=spectral[0],
                    velocity_min_kms=spectral[1],
                    velocity_max_kms=spectral[2],
                )
            )
        return details

    def query(self, **filters) -> "Collection":
        """Filter the catalog along any metadata axis; return a composable sub-:class:`Collection`.

        Slices the corpus along the catalog's own fields (``species=`` / ``survey=`` /
        ``telescope=`` / ``transition=`` / ``object=`` / …, PRD #56 user story 4) without parsing
        filenames. Multiple filters are **conjunctive** (a store must match all). The result is a
        new ``Collection`` over the same corpus root — so every facade method composes on it
        (``query(...).list()``, ``query(...).describe(obj)``, ``query(...).query(...)``), which is
        the more useful return than a bare ``list[Record]`` (a list would dead-end the chain).

        An **unknown filter key** raises :class:`ValueError` naming the offending key and the valid
        axes — a silent no-op would hide a typo (``surveys=`` vs ``survey=``) as an empty result
        (ADR-0003: refuse, never silently mis-match). The filter axes are exactly the
        :class:`~astrolyze.collection.catalog.CatalogRow` fields."""
        valid = {f.name for f in fields(CatalogRow)}
        unknown = set(filters) - valid
        if unknown:
            raise ValueError(
                f"unknown collection filter key(s): {', '.join(sorted(unknown))}; "
                f"valid filter axes are {', '.join(sorted(valid))}"
            )
        matched = tuple(
            row
            for row in self._catalog.rows
            if all(getattr(row, key) == value for key, value in filters.items())
        )
        # A sub-Collection over the same root: the matched rows keep their relative store_path, so
        # every record still resolves against the unchanged root_uri (covering()/open stay valid).
        return Collection(
            replace(self._catalog, rows=matched),
            self._root_uri,
        )

    # NOTE: covering(SkyCoord) lands here next (#62) — radius/MOC prefilter over the footprint
    # columns, final containment decided by the candidate store's own WCS. Left unimplemented here
    # to keep this slice (#60) additive; the seam is the same `records` axis query() filters.


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


def _spectral_axis_kms(record: "Record"):
    """The (channel_width, v_min, v_max) of *record*'s store in km/s — the deep-describe read.

    Opens the record's store (still lazy — :meth:`Record.open` is dask-backed, so only the 1-D
    spectral axis is materialised, not the cube data) and reads its native channel width and
    velocity coverage off the cube's own spectral axis. These are the fields the catalog does not
    carry, so they are only ever read on a ``describe(..., deep=True)`` request.

    Returns ``(None, None, None)`` when the store has no usable spectral axis (e.g. a 2-D map or an
    axis that does not convert to velocity) — the field stays unknown, never invented (ADR-0003)."""
    import astropy.units as u

    try:
        axis = record.open()._sc.spectral_axis.to(u.km / u.s)
    except Exception:  # noqa: BLE001 — a non-velocity / map store has no velocity coverage to add
        return (None, None, None)
    if len(axis) < 2:
        return (None, None, None)
    width = abs(float((axis[1] - axis[0]).value))
    return (width, float(axis.min().value), float(axis.max().value))


__all__ = ["Collection", "Record", "ObjectSummary", "StoreDetail"]
