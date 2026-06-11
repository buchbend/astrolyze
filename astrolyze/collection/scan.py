"""The scan-builder: a :class:`~astrolyze.collection.catalog.Catalog` from a bare store directory (#61).

A *published corpus* carries a curated ``catalog.parquet`` at its root (the ifm-published index,
read by :mod:`astrolyze.collection.catalog`). But a directory of astrolyze-flavoured Zarr stores
with **no** such index is still browsable: every store is self-describing (it carries the
:class:`~astrolyze.io.Metadata` schema and a verbatim WCS in its group attrs, issue #23), so a
catalog can be *reconstructed* by reading those attrs and computing each store's sky footprint
from its own WCS. This is the catalog-less fallback for :meth:`~astrolyze.collection.Collection.open`
**and** the community path for cataloguing any compatible collection with zero extra tooling
(PRD #56 user story 29).

The builder targets the **same** catalog format the reader validates (the spec's `CatalogRow` /
`CATALOG_COLUMNS`), so a catalog built here is interchangeable with a published one (catalog-schema
spec §6): :func:`build_catalog` returns a :class:`~astrolyze.collection.catalog.Catalog` of typed
:class:`~astrolyze.collection.catalog.CatalogRow` records, and :func:`write_catalog` persists it as
``catalog.parquet`` so the build -> persist -> read round-trip passes the same schema validation as
a published catalog. All catalog-format knowledge stays behind the catalog seam — this module
reuses ``CatalogRow`` / ``CATALOG_COLUMNS`` rather than re-deriving the column set.

**What a scanned row carries, and what it cannot.** The store supplies the physics (beam, bunit,
rest frequency, species, transition, systemic velocity) and the sky footprint (center+radius from
its celestial WCS — the same `calc_footprint` + max-corner-separation circle the ifm publisher
computes per ADR-0001). What a bare directory has *no* source for stays ``null``, never guessed
(the catalog honesty rule, ADR-0003): ``survey`` and ``content_checksum`` are manifest identity
the curation pipeline owns, and ``v_sys_kms`` is curated context — a scan reads a store's recorded
systemic velocity onto the metadata but does not *curate* one, so the catalog's curated-v_sys
column is left null here (no curation source in a bare directory).

**MOC is an optional consumer extra.** When ``mocpy`` is importable the builder also computes each
store's exact-coverage MOC (the v1.1 ``moc`` column) and the built catalog declares
``catalog_schema_version`` ``"1.1"``; without ``mocpy`` the MOC is left null, the catalog declares
``"1.0"``, and center+radius prefiltering still works end-to-end (PRD #56: mocpy a hard *producer*
dependency, only an optional *consumer* one — the bare install must work). The in-memory
:class:`~astrolyze.collection.catalog.CatalogRow` models only the v1.0 footprint columns, so a
built MOC rides alongside in :class:`ScanResult` and is written into the persisted parquet's
trailing ``moc`` column (additive evolution, §3) — a v1.0 reader simply does not see it.

**Incomplete stores are skipped with a warning, never mis-catalogued.** A directory that is not an
astrolyze store, or a store missing the celestial WCS a footprint needs, is reported via a
:class:`ScanWarning` and left out of the catalog rather than entered with fabricated geometry
(PRD #56 acceptance: incomplete attrs -> clear warning, not silent mis-catalogue).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

from .catalog import (
    CATALOG_COLUMNS,
    CATALOG_FILENAME,
    Catalog,
    CatalogRow,
)

# The v3 Zarr group sentinel — a directory carrying this file is a Zarr node (group or array). A
# cube store written by the astrolyze backend (issue #23) is a *group*, but so are its data/coord
# arrays' parents and its noise/validity companion subgroups, so the sentinel alone does not
# identify a cube store — the attrs do (see _is_cube_store).
_ZARR_V3_SENTINEL = "zarr.json"

# The provenance marker the astrolyze Zarr backend writes on a *cube store's* group attrs
# (``_save_zarr`` -> ``{"provenance": {"backend": "zarr"}}``). It is the precise cube-store
# discriminant: a store's data/coord arrays are ``node_type: array`` (no provenance), and its
# noise/validity companion subgroups are groups but carry method/version attrs, not this marker.
# Matching on the backend's own marker means the scan catalogues exactly the cube stores, never a
# phantom inner array or companion (no name-pattern guessing, ADR-0003).
_PROVENANCE_ATTR = "provenance"
_PROVENANCE_BACKEND = "zarr"

# The MOC resolution (HEALPix order) the footprint polygon is rasterized to when mocpy is present —
# matched to the producer's choice (catalog-schema spec §2.6 / ifm ``_MOC_MAX_DEPTH``) so a
# scan-built MOC is interchangeable with a published one. Order 12 is ~3.2" cells, finer than the
# arcminute-to-degree fields a single-dish/interferometer corpus holds.
_MOC_MAX_DEPTH = 12

# The catalog version a scan declares. A scan with MOCs is a v1.1 catalog (it carries the moc
# column); without mocpy it is a plain v1.0 center+radius catalog. Both are valid published-format
# catalogs the reader accepts (a v1.0 reader ignores a v1.1's trailing moc — §3.1).
_SCHEMA_VERSION_WITH_MOC = "1.1"
_SCHEMA_VERSION_BARE = "1.0"


class ScanWarning(UserWarning):
    """Warned (not raised) when a directory under the scan root cannot be catalogued.

    A directory that is not an astrolyze store, or a store missing the celestial WCS a footprint
    needs, is skipped and surfaced through this warning rather than entered with fabricated
    geometry — the scan never silently mis-catalogues a store (PRD #56, ADR-0003). The warning
    names the offending path and why it was skipped so the gap is actionable."""


@dataclass(frozen=True)
class ScanResult:
    """One scanned store: its typed v1.0 :class:`CatalogRow` plus the optional v1.1 ``moc``.

    The :class:`CatalogRow` carries the published-format v1.0 columns (identity, physics,
    footprint) the reader models. The exact-coverage ``moc`` (the v1.1 additive column) cannot live
    on the v1.0 row, so it rides here alongside it — ``None`` when mocpy is absent or no footprint
    MOC could be built. :func:`write_catalog` folds it into the persisted parquet's trailing ``moc``
    column; the in-memory :class:`Catalog` keeps to the v1.0 row the reader returns."""

    row: CatalogRow
    moc: str | None


def build_catalog(root) -> Catalog:
    """Scan the store directory at *root* and build a :class:`Catalog` from each store's attrs + WCS.

    The catalog-less counterpart of :func:`~astrolyze.collection.catalog.read_catalog`: instead of
    reading a curated ``catalog.parquet``, it reconstructs the same typed :class:`Catalog` by
    reading every astrolyze Zarr store under *root*, reading its physical context off the
    :class:`~astrolyze.io.Metadata` in the group attrs, and computing its sky footprint
    (center+radius, and a MOC when mocpy is present) from the store's own celestial WCS.

    *root* is a local directory path (the bare-directory community path; the fsspec-URL form is the
    same shape but the remote store opener is #63). Stores are located by their Zarr group sentinel,
    recursively, with the catalog rows ordered by store path for a stable rendering. A directory
    that is not a readable astrolyze store, or one without the celestial WCS a footprint needs, is
    skipped with a :class:`ScanWarning` rather than catalogued with invented geometry.

    The built catalog declares ``catalog_schema_version`` ``"1.1"`` when mocpy is importable (it
    carries the exact-coverage MOC), else ``"1.0"`` (center+radius only) — both are valid
    published-format catalogs (catalog-schema spec §3.1)."""
    results = scan_directory(root)
    version = _scan_schema_version()
    rows = tuple(result.row for result in results)
    return Catalog(rows=rows, schema_version=version)


def scan_directory(root) -> list[ScanResult]:
    """Scan *root* for astrolyze Zarr stores, returning one :class:`ScanResult` per readable store.

    The low-level walk behind :func:`build_catalog` / :func:`write_catalog`: it discovers the
    stores, reads each one's attrs + WCS, and folds them into ``(row, moc)`` pairs — exposed
    separately so a caller that needs the MOCs alongside the rows (the persistence path) does not
    rebuild them. Unreadable / footprint-less stores are skipped with a :class:`ScanWarning`.
    Ordered by store path (POSIX) for a stable catalog."""
    root = Path(root)
    results: list[ScanResult] = []
    for store in _find_stores(root):
        result = _summarize_store(store, root)
        if result is not None:
            results.append(result)
    return results


def write_catalog(root, *, overwrite: bool = False) -> Path:
    """Scan *root*, persist the built catalog as ``catalog.parquet`` at *root*, and return its path.

    The community catalogue-once path (PRD #56 user story 29): turns a bare directory of stores
    into a *published* corpus by writing the same ``catalog.parquet`` the reader validates, so a
    later :func:`~astrolyze.collection.catalog.read_catalog` (or :meth:`Collection.open`) consumes
    it exactly as it would an ifm-published one — the build -> persist -> read round-trip passes the
    same schema validation. The parquet carries the v1.0 footprint columns always and the v1.1
    ``moc`` column when mocpy built one (additive, stamped ``"1.1"``); ``overwrite`` guards an
    existing index. Raises :class:`FileNotFoundError` when *root* does not exist."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"scan root {root} is not a directory")
    out = root / CATALOG_FILENAME
    if out.exists() and not overwrite:
        raise FileExistsError(
            f"a catalog already exists at {out}; pass overwrite=True to rebuild it "
            "(the scan never silently clobbers a published index)"
        )

    results = scan_directory(root)
    version = _scan_schema_version()
    table = _to_table(pa, results, version=version)
    pq.write_table(table, out)
    return out


# -- store discovery -------------------------------------------------------------------
def _find_stores(root: Path) -> list[Path]:
    """Every astrolyze cube store under *root*, identified by the backend's provenance marker.

    Walks the ``zarr.json`` sentinels and keeps only the directories that are cube stores — a group
    whose attrs carry the astrolyze backend marker (see :func:`_is_cube_store`). This excludes the
    store's own data/coord arrays (``node_type: array``) and its noise/validity companion subgroups
    (groups without the marker), so neither is catalogued as a phantom cube. Sorted by POSIX path
    for a stable catalog order."""
    stores = [
        sentinel.parent
        for sentinel in root.rglob(_ZARR_V3_SENTINEL)
        if _is_cube_store(sentinel)
    ]
    return sorted(stores, key=lambda p: p.relative_to(root).as_posix())


def _is_cube_store(sentinel: Path) -> bool:
    """Whether the Zarr node at *sentinel* (a ``zarr.json``) is an astrolyze cube store.

    A cube store's group attrs carry the backend provenance marker the astrolyze Zarr writer stamps
    (:data:`_PROVENANCE_ATTR` = ``{"backend": "zarr"}``). Reading the marker straight off the
    ``zarr.json`` keeps discovery cheap (no per-candidate store open) and precise — an inner data
    array (``node_type: array``) and a noise/validity companion group (no marker) both fail it. A
    malformed/unreadable ``zarr.json`` is simply not a cube store (it is skipped silently here; a
    genuinely broken candidate store surfaces its error at :func:`_summarize_store`)."""
    import json

    try:
        meta = json.loads(sentinel.read_text())
    except (OSError, ValueError):
        return False
    if meta.get("node_type") != "group":
        return False
    provenance = meta.get("attributes", {}).get(_PROVENANCE_ATTR)
    return (
        isinstance(provenance, dict)
        and provenance.get("backend") == _PROVENANCE_BACKEND
    )


# -- per-store summary -----------------------------------------------------------------
def _summarize_store(store: Path, root: Path) -> ScanResult | None:
    """Read *store*'s attrs + WCS into a :class:`ScanResult`; ``None`` (with a warning) if skipped.

    Reads the physical context off the store's :class:`~astrolyze.io.Metadata` and computes the sky
    footprint from its own celestial WCS — the thin path the ifm publisher uses (ADR-0001), via the
    io seam so nothing materialises (the load is lazy/dask-backed, #23). A store that cannot be read
    as an astrolyze dataset, or that lacks the celestial WCS a footprint needs, is skipped with a
    :class:`ScanWarning` rather than catalogued with invented geometry (ADR-0003)."""
    from astrolyze.io import load

    rel = store.relative_to(root).as_posix()
    try:
        loaded = load(store)
    except Exception as exc:  # noqa: BLE001 — any unreadable store is reported, not fatal
        warnings.warn(
            f"skipping {rel}: not a readable astrolyze store ({type(exc).__name__}: {exc})",
            ScanWarning,
            stacklevel=2,
        )
        return None

    footprint = _footprint(loaded.wcs, loaded.data.shape)
    if footprint is None:
        warnings.warn(
            f"skipping {rel}: no celestial WCS to compute a sky footprint from "
            "(astrolyze never invents a footprint, ADR-0003)",
            ScanWarning,
            stacklevel=2,
        )
        return None
    ra_deg, dec_deg, radius_deg, corner_coords = footprint

    md = loaded.metadata
    beam_major, beam_minor, beam_pa = _beam_columns(md.beam)
    transition = md.lines[0].transition if md.lines else None
    record = {
        # Identity the store records; survey is manifest-only -> null in a bare scan.
        "object": md.object,
        "survey": None,
        "telescope": md.telescope,
        "species": md.species,
        "transition": transition,
        "rest_frequency_hz": _hz(md.rest_frequency),
        # Physics off the store.
        "beam_major_arcsec": beam_major,
        "beam_minor_arcsec": beam_minor,
        "beam_pa_deg": beam_pa,
        "bunit": str(md.bunit) if md.bunit is not None else None,
        # Locator: the store's path relative to the scan root (POSIX), never null.
        "store_path": rel,
        # No curation source in a bare directory -> no content checksum.
        "content_checksum": None,
        # Footprint from the store's own WCS (center + max-corner radius).
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "radius_deg": radius_deg,
        "catalog_schema_version": _scan_schema_version(),
    }
    row = CatalogRow.from_record(record, version=_scan_schema_version())
    moc = _footprint_moc(corner_coords)
    return ScanResult(row=row, moc=moc)


def _footprint(wcs, data_shape):
    """The sky footprint (ra_deg, dec_deg, radius_deg, corner SkyCoords) from a store's WCS.

    Center+radius is the coarse prefilter the published catalog carries (catalog-schema spec §2.4):
    the geometric-center pixel and the radius of the smallest circle about it enclosing the four
    image corners — exactly the ifm publisher's computation (ADR-0001). Returns ``None`` when the
    WCS has no celestial part (no footprint can be computed; the store is skipped, never invented).
    The corner SkyCoords are returned alongside so the optional MOC reuses the same geometry."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    if not getattr(wcs, "has_celestial", False):
        return None
    celestial = wcs.celestial
    # data is in FITS array order (spectral, y, x) for a cube or (y, x) for a map; the celestial
    # grid is the trailing (ny, nx). calc_footprint wants axes=(nx, ny) and returns corner sky
    # coords in degrees; pixel_to_world gives the center of the (nx, ny) grid.
    ny, nx = data_shape[-2], data_shape[-1]
    corners = celestial.calc_footprint(axes=(nx, ny))
    center = celestial.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    corner_coords = SkyCoord(corners[:, 0] * u.deg, corners[:, 1] * u.deg)
    radius_deg = float(center.separation(corner_coords).max().to_value(u.deg))
    return float(center.ra.deg), float(center.dec.deg), radius_deg, corner_coords


def _footprint_moc(corner_coords) -> str | None:
    """The exact-coverage MOC (ascii string) of the footprint polygon — ``None`` without mocpy.

    The v1.1 ``moc`` column (catalog-schema spec §2.5/§2.6): a MOC of the four WCS sky corners,
    serialized to mocpy's portable ascii form, rasterized at :data:`_MOC_MAX_DEPTH` to match the
    producer. mocpy is an **optional** consumer extra — absent, this returns ``None`` and the scan
    falls back to the always-present center+radius prefilter (PRD #56: the bare install must work).
    A degenerate footprint mocpy cannot polygonise is also ``None`` (exact coverage is exact or
    absent, never a radius circle dressed up as a MOC — spec §5)."""
    try:
        from mocpy import MOC
    except ImportError:
        return None
    try:
        moc = MOC.from_polygon_skycoord(corner_coords, max_depth=_MOC_MAX_DEPTH)
    except Exception:  # noqa: BLE001 — a degenerate footprint -> null MOC, not a scan failure
        return None
    return moc.to_string(format="ascii")


# -- metadata -> column helpers --------------------------------------------------------
def _beam_columns(beam):
    """A store beam projected onto the catalog's (major, minor, PA) arcsec/deg columns.

    ``(None, None, None)`` when the store states no beam — a partial row is valid (the catalog never
    invents a beam, the same honesty rule as the store schema)."""
    import astropy.units as u

    if beam is None:
        return None, None, None
    return (
        float(beam.major.to_value(u.arcsec)),
        float(beam.minor.to_value(u.arcsec)),
        float(beam.pa.to_value(u.deg)),
    )


def _hz(rest_frequency):
    """A rest-frequency ``Quantity`` projected onto the catalog's ``rest_frequency_hz`` (Hz float);
    ``None`` passes through (an uncurated rest frequency is null, never guessed)."""
    import astropy.units as u

    if rest_frequency is None:
        return None
    return float(rest_frequency.to_value(u.Hz))


def _scan_schema_version() -> str:
    """The catalog version a scan declares: ``"1.1"`` when mocpy is importable (the scan carries the
    exact-coverage MOC column), else ``"1.0"`` (center+radius only). Both are valid published-format
    versions a reader accepts (catalog-schema spec §3.1)."""
    try:
        import mocpy  # noqa: F401
    except ImportError:
        return _SCHEMA_VERSION_BARE
    return _SCHEMA_VERSION_WITH_MOC


# -- parquet persistence ---------------------------------------------------------------
def _to_table(pa, results: list[ScanResult], *, version: str):
    """Build a pyarrow table of *results* in the contract column order, stamping *version*.

    The v1.0 :data:`CATALOG_COLUMNS` always; the v1.1 ``moc`` column appended (after the footprint,
    before ``catalog_schema_version`` per spec §2) only when this scan produced MOCs — additive
    evolution, so a v1.0 reader reads the file ignoring the trailing ``moc``. The version is written
    into the parquet key-value metadata too, so a reader gates from the header (spec §2.8/§3)."""
    rows = [result.row for result in results]
    columns = {name: [getattr(row, name) for row in rows] for name in CATALOG_COLUMNS}
    has_moc = version == _SCHEMA_VERSION_WITH_MOC
    if has_moc:
        # Insert the v1.1 moc column at its contract position: after radius_deg (the last v1.0
        # footprint column), before the trailing catalog_schema_version.
        moc_values = [result.moc for result in results]
        ordered = {}
        for name, values in columns.items():
            if name == "catalog_schema_version":
                ordered["moc"] = moc_values
            ordered[name] = values
        columns = ordered
    table = pa.table(columns)
    return table.replace_schema_metadata({"catalog_schema_version": version})


__all__ = [
    "build_catalog",
    "scan_directory",
    "write_catalog",
    "ScanResult",
    "ScanWarning",
]
