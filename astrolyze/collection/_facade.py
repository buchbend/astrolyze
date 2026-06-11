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
- **#62 covering** adds ``Collection.covering(SkyCoord)``: it *prefilters* on the footprint columns
  parsed on each row — the always-present center+radius (``ra_deg`` / ``dec_deg`` / ``radius_deg``)
  and, when mocpy is installed and the row carries one, the exact ``moc`` — then makes the final
  containment decision from the candidate store's own celestial WCS at open time (the catalog is an
  index, the WCS is the authority).
- **#61 scan-builder** wires :meth:`Collection.open` to fall back to
  :func:`~astrolyze.collection.scan.build_catalog` when the root carries no ``catalog.parquet``,
  building an equivalent :class:`~astrolyze.collection.catalog.Catalog` by scanning the stores —
  it reuses this same facade, the trigger being the index file's presence.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from typing import TYPE_CHECKING

from fsspec.core import url_to_fs

from .catalog import CATALOG_FILENAME, Catalog, CatalogRow, read_catalog

if TYPE_CHECKING:
    from astropy.coordinates import SkyCoord


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

    def covering(self, position: "SkyCoord") -> "Collection":
        """Stores whose footprint contains *position*; return a composable sub-:class:`Collection`.

        The sky-coverage search (PRD #56 user stories 5/6/7): which cubes cover a sky position,
        across every survey, without knowing which observed it. The semantics are deliberately
        **two-layered** (catalog-schema spec §2.4/§2.5/§5) — the catalog footprint is a *prefilter
        only*, never the geometry authority:

        1. **Prefilter (cheap, no store opened).** Every record is first filtered by the catalog's
           own footprint columns: the always-present center+radius circle
           (``ra_deg``/``dec_deg``/``radius_deg``) — works on a bare install — and, *when mocpy is
           installed and the row carries a ``moc``*, the exact-coverage MOC (eliminating the
           corner-of-the-circle false positives radius alone admits). See
           :func:`_prefilter_candidates`.
        2. **WCS authority (the final decision).** For every *surviving* candidate the store is
           opened (lazily — :meth:`Record.open` is dask-backed, so only the celestial WCS + grid
           shape are touched, not the cube data) and the position is tested against that store's
           **own** celestial pixel grid (:func:`_position_in_store_wcs`). A candidate the prefilter
           admitted but whose real WCS does **not** contain the position is REJECTED. This guards
           against a stale or coarse index ever yielding a wrong containment answer (user story 6):
           the catalog is an index, the WCS is the authority.

        The result is a new ``Collection`` over the same corpus root (like :meth:`query`), so every
        facade method composes on it (``covering(pos).query(...)``, ``covering(pos).list()``). A
        position covered by nothing returns an **empty** collection (an honest "no coverage", never
        an error)."""
        candidates = _prefilter_candidates(self.records, position)
        covered = tuple(
            record.row
            for record in candidates
            if _position_in_store_wcs(record, position)
        )
        return Collection(replace(self._catalog, rows=covered), self._root_uri)

    def stack(self, position_or_sources, *, size, partial: str = "raise"):
        """Gather sky-coordinate cutouts into a :class:`~astrolyze.collection.stack.Stack` (#64).

        The stage-1 stack builder (PRD #56 user stories 10/11): assemble a multi-survey,
        multi-line view of a target — or a multi-target sample — by taking a
        :meth:`~astrolyze.core.Cube.cutout` of angular *size* from every cube that covers each
        position. The result is a :class:`~astrolyze.collection.stack.Stack`: an aligned-cutout
        container that is **always safe to construct and browse**, however heterogeneous its
        members. Co-addition is a separate explicit stage (#65) — this call never co-adds.

        Parameters
        ----------
        position_or_sources : ~astropy.coordinates.SkyCoord or list
            **Either** a single sky position (a scalar :class:`~astropy.coordinates.SkyCoord`) —
            :meth:`covering` finds every cube containing it and one cutout is taken from each
            (user story 10) — **or** a list of catalog object **names** and/or SkyCoords, building
            a multi-target sample (user story 11). A *name* is resolved against **this collection's
            catalog** (no SIMBAD/Sesame in the library — PRD #56 out of scope); an unknown name
            raises :class:`KeyError`. A non-scalar SkyCoord array is treated as a list of its
            positions.
        size : ~astropy.units.Quantity
            The angular stamp size handed to :meth:`~astrolyze.core.Cube.cutout` (a scalar square
            or a ``(height, width)`` pair). Keyword-only — a stack always declares its stamp size.
        partial : {"raise", "trim"}, default "raise"
            The cutout edge rule for a covering cube whose footprint cannot fit the full stamp
            (passed through to :meth:`~astrolyze.core.Cube.cutout`): ``"raise"`` rejects a partial
            window, ``"trim"`` opts into the trimmed stamp (no silent clipping either way).

        Returns
        -------
        ~astrolyze.collection.stack.Stack
            The gathered members (each a cutout carrying origin provenance + identity), with the
            request's :class:`~astrolyze.collection.stack.Selection` provenance recorded.
        """
        from .stack import gather_stack

        return gather_stack(self, position_or_sources, size=size, partial=partial)


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


# -- covering(): prefilter then WCS authority (#62) -------------------------------------
def _prefilter_candidates(records, position: "SkyCoord") -> list["Record"]:
    """The records whose catalog footprint *might* cover *position* — the cheap prefilter layer.

    The catalog footprint is an **index, not the authority** (catalog-schema spec §2.4/§2.5): a
    candidate that survives here is still decided by its store's own WCS in
    :func:`_position_in_store_wcs`. Two prefilters, in increasing sharpness, never opening a store:

    - **center+radius** (always available, the v1.0 footprint): the position is within
      ``radius_deg`` of (``ra_deg``, ``dec_deg``). A row missing any footprint column cannot be
      prefiltered, so it is kept as a candidate (the WCS authority then decides — a null footprint
      never silently drops a store).
    - **MOC** (only when mocpy is importable AND the row carries a ``moc``): the exact-coverage
      polygon test, which removes the corner-of-the-circle false positives radius alone admits.
      mocpy is an **optional consumer extra** (``astrolyze[coverage]``); absent, this layer is
      simply skipped and center+radius stands alone — the bare install works end to end (PRD #56).

    A candidate must pass *every* prefilter it can be evaluated against (radius AND, where present,
    MOC). The MOC import is lazy and per-call so importing astrolyze never imports mocpy."""
    moc_module = _try_import_moc()
    candidates = []
    for record in records:
        row = record.row
        if not _within_radius(row, position):
            continue
        # MOC is the finer prefilter: applied only when mocpy is present AND this row carries one.
        # A row without a moc (a v1.0 catalog, or a null cell) falls back to the radius result.
        if moc_module is not None and row.moc is not None:
            if not _within_moc(moc_module, row.moc, position):
                continue
        candidates.append(record)
    return candidates


def _within_radius(row: CatalogRow, position: "SkyCoord") -> bool:
    """Whether *position* is within the row's center+radius prefilter circle.

    ``True`` when the angular separation from (``ra_deg``, ``dec_deg``) is ``<= radius_deg``. A row
    with no center/radius footprint cannot be prefiltered on geometry, so it returns ``True`` (kept
    as a candidate for the WCS authority to decide — a missing footprint never drops a store
    silently, ADR-0003)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    if row.ra_deg is None or row.dec_deg is None or row.radius_deg is None:
        return True
    center = SkyCoord(row.ra_deg * u.deg, row.dec_deg * u.deg)
    return bool(center.separation(position).to_value(u.deg) <= row.radius_deg)


def _within_moc(moc_module, moc_string: str, position: "SkyCoord") -> bool:
    """Whether *position* falls inside the row's exact-coverage MOC (the finer prefilter).

    Deserializes the mocpy ASCII MOC (catalog-schema spec §2.6) and tests containment. A MOC that
    cannot be deserialized or queried is treated as *non-discriminating* (returns ``True``): a
    malformed index entry must never silently drop a store — the store's WCS stays the authority
    (ADR-0003)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    try:
        moc = moc_module.from_string(moc_string, format="ascii")
        contains = moc.contains_skycoords(
            SkyCoord(position.icrs.ra * u.deg, position.icrs.dec * u.deg)
        )
    except Exception:  # noqa: BLE001 — a bad MOC is non-discriminating, not a hard failure
        return True
    return bool(_as_scalar_bool(contains))


def _position_in_store_wcs(record: "Record", position: "SkyCoord") -> bool:
    """Whether *position* falls inside the candidate store's OWN celestial pixel grid — the authority.

    This is the final containment decision (catalog-schema spec §5, PRD #56 user story 6): the
    store is opened (lazily — :meth:`Record.open` is dask-backed, so this touches only the celestial
    WCS and the grid shape, never the cube data) and the position is mapped to a pixel on the
    store's real celestial WCS. Containment uses the **same** rule as :meth:`Cube.cutout`'s
    footprint test (issue #58): a world coordinate maps to a pixel centre at integer indices, so the
    valid extent is the half-open ``[-0.5, n-0.5]`` cell range; a NaN pixel (an undefined projection
    at this position, e.g. outside a SIN field of view) is off the footprint too.

    A store that cannot be opened, or that has no celestial WCS, cannot confirm containment, so it
    returns ``False`` — an unverifiable candidate is *not* asserted to cover the position (refuse to
    guess, ADR-0003), rather than trusting the prefilter as the answer."""
    import numpy as np

    try:
        cube = record.open()
        wcs = cube._sc.wcs
        if not getattr(wcs, "has_celestial", False):
            return False
        celestial = wcs.celestial
        ny, nx = cube.shape[-2], cube.shape[-1]
        x, y = celestial.world_to_pixel(position)
    except Exception:  # noqa: BLE001 — an unopenable/footprint-less store cannot confirm coverage
        return False
    on_x = bool(np.isfinite(x)) and -0.5 <= float(x) <= nx - 0.5
    on_y = bool(np.isfinite(y)) and -0.5 <= float(y) <= ny - 0.5
    return on_x and on_y


def _try_import_moc():
    """The ``mocpy.MOC`` class if mocpy is importable, else ``None`` (the optional-extra gate).

    mocpy is an optional consumer extra (``astrolyze[coverage]``): a bare install must run the
    center+radius prefilter + WCS authority end to end without it. The import is lazy and per-call
    so ``import astrolyze`` never pulls mocpy (the cost is one cached import on a covering call)."""
    try:
        from mocpy import MOC
    except ImportError:
        return None
    return MOC


def _as_scalar_bool(value) -> bool:
    """Collapse mocpy's containment result (a 0-d / 1-element array or a bool) to a scalar bool."""
    import numpy as np

    arr = np.asarray(value)
    return bool(arr.reshape(-1)[0]) if arr.size else False


__all__ = ["Collection", "Record", "ObjectSummary", "StoreDetail"]
