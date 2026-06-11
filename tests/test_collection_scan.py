"""Tests for the scan-builder — a catalog from a bare directory of Zarr stores (issue #61).

Written first (red/green TDD). The scan-builder reconstructs a
:class:`~astrolyze.collection.catalog.Catalog` from any directory of astrolyze-flavoured Zarr
stores with **no** published ``catalog.parquet``: it reads each store's group attrs (the
:class:`~astrolyze.io.Metadata` schema, issue #23) and computes its sky footprint from the store's
own celestial WCS. This is both the catalog-less fallback for
:meth:`~astrolyze.collection.Collection.open` and the community path for cataloguing arbitrary
collections (PRD #56 user story 29).

The contract under test:

1. **Build -> Catalog round-trips.** :func:`build_catalog` returns typed
   :class:`~astrolyze.collection.catalog.CatalogRow` records; :func:`write_catalog` persists a
   ``catalog.parquet`` that :func:`~astrolyze.collection.catalog.read_catalog` reads back and
   schema-validates exactly as a published one (build -> persist -> read).
2. **Footprints from WCS.** The center+radius footprint columns are computed from a store with a
   known WCS and match a hand-computed expectation.
3. **mocpy is an optional consumer extra.** Without mocpy the build still works end-to-end
   (center+radius only, schema ``"1.0"``); the MOC assertion is guarded behind ``importorskip``.
4. **Incomplete stores are skipped with a warning**, not silently mis-catalogued.
5. **``Collection.open`` on a catalog-less directory** transparently scans and returns a working
   collection whose records open into lazy cubes.

The reference corpus mirrors the published-catalog suite: NGC 3521 in two surveys' worth of stores
plus NGC 0628, so a scan-built catalog and a published one are interchangeable.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from astropy.wcs import WCS

from astrolyze.core import Cube
from astrolyze.io.access import LoadedData
from astrolyze.io.schema import Metadata

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")  # the catalog persistence/read path needs parquet
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz
REST_HI = 1.420405751e9  # HI 21cm Hz

# A known WCS: a SIN-projected grid centred at (RA, Dec) = (170.0, -0.04) deg, 2e-4 deg/pixel.
CRVAL1 = 170.0
CRVAL2 = -0.04
CDELT = 2e-4  # deg/pixel


def _header(*, obj, rest_hz, species, telescope, bmaj_arcsec, ny=4, nx=5):
    """A small schema-complete cube header with a real celestial+spectral WCS (the known grid)."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", CRVAL1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -CDELT, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", CRVAL2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = CDELT, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = obj
    h["TELESCOP"] = telescope
    h["BUNIT"] = "K"
    h["RESTFRQ"] = (rest_hz, "Hz")
    h["HIERARCH ASTROLYZE SPECIES"] = species
    beam = radio_beam.Beam(
        major=bmaj_arcsec * u.arcsec, minor=(bmaj_arcsec - 2) * u.arcsec, pa=30 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    return h, beam


def _loaded(data, header):
    """A LoadedData built straight from an in-memory array + header (no file)."""
    wcs = WCS(header)
    return LoadedData(
        data=data,
        wcs=wcs,
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )


def _write_store(root, rel_dir, *, obj, rest_hz, species, telescope, bmaj_arcsec):
    """Write one synthetic astrolyze Zarr store under ``root/rel_dir`` and return its store path."""
    header, _ = _header(
        obj=obj,
        rest_hz=rest_hz,
        species=species,
        telescope=telescope,
        bmaj_arcsec=bmaj_arcsec,
    )
    data = np.arange(3 * 4 * 5, dtype="float32").reshape(3, 4, 5)
    cube = Cube.from_loaded(_loaded(data, header))
    target = root / rel_dir
    target.mkdir(parents=True, exist_ok=True)
    cube.to_zarr(target)
    return next(target.glob("*.zarr"))


@pytest.fixture
def store_dir(tmp_path):
    """A bare directory of 3 astrolyze Zarr stores (2 objects), with **no** catalog.parquet."""
    root = tmp_path / "bare_corpus"
    root.mkdir()
    specs = [
        dict(
            rel_dir="l1/heracles",
            obj="NGC3521",
            rest_hz=REST_CO21,
            species="CO",
            telescope="IRAM30M",
            bmaj_arcsec=13.0,
        ),
        dict(
            rel_dir="l1/phangs",
            obj="NGC3521",
            rest_hz=REST_CO21,
            species="CO",
            telescope="ALMA",
            bmaj_arcsec=1.5,
        ),
        dict(
            rel_dir="l1/things",
            obj="NGC0628",
            rest_hz=REST_HI,
            species="HI",
            telescope="VLA",
            bmaj_arcsec=11.0,
        ),
    ]
    for spec in specs:
        _write_store(root, spec.pop("rel_dir"), **spec)
    return root


# --------------------------------------------------------------------------------------
# build_catalog: a Catalog from a bare directory, reusing the published-format CatalogRow
# --------------------------------------------------------------------------------------
def test_build_catalog_returns_typed_rows_for_every_store(store_dir):
    from astrolyze.collection.catalog import CatalogRow
    from astrolyze.collection.scan import build_catalog

    catalog = build_catalog(store_dir)
    assert len(catalog.rows) == 3
    assert all(isinstance(r, CatalogRow) for r in catalog.rows)
    # The physics the store records is read onto the row (object/species/beam/bunit/rest freq).
    by_object = {r.object for r in catalog.rows}
    assert by_object == {"NGC3521", "NGC0628"}
    co = next(
        r for r in catalog.rows if r.object == "NGC3521" and r.telescope == "IRAM30M"
    )
    assert co.species == "CO"
    assert co.bunit == "K"
    assert co.rest_frequency_hz == pytest.approx(REST_CO21)
    assert co.beam_major_arcsec == pytest.approx(13.0)
    assert co.store_path.endswith(".zarr")


def test_scanned_rows_leave_uncurated_columns_null(store_dir):
    """A bare directory has no manifest: survey + content_checksum have no source and stay null —
    the scan never invents curation identity (the catalog honesty rule, ADR-0003)."""
    from astrolyze.collection.scan import build_catalog

    catalog = build_catalog(store_dir)
    assert all(r.survey is None for r in catalog.rows)
    assert all(r.content_checksum is None for r in catalog.rows)


# --------------------------------------------------------------------------------------
# Footprint columns: computed from a store's own (known) WCS
# --------------------------------------------------------------------------------------
def test_footprint_matches_hand_computed_wcs(store_dir):
    from astrolyze.collection.scan import build_catalog

    catalog = build_catalog(store_dir)
    row = next(r for r in catalog.rows if r.object == "NGC0628")

    # Reconstruct the same footprint the publisher computes: center pixel + max-corner separation.
    header, _ = _header(
        obj="NGC0628",
        rest_hz=REST_HI,
        species="HI",
        telescope="VLA",
        bmaj_arcsec=11.0,
    )
    celestial = WCS(header).celestial
    ny, nx = 4, 5  # the cube grid (3, 4, 5) -> (ny, nx) = (4, 5)
    corners = celestial.calc_footprint(axes=(nx, ny))
    center = celestial.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    from astropy.coordinates import SkyCoord

    corner_coords = SkyCoord(corners[:, 0] * u.deg, corners[:, 1] * u.deg)
    expected_radius = float(center.separation(corner_coords).max().to_value(u.deg))

    assert row.ra_deg == pytest.approx(center.ra.deg, abs=1e-9)
    assert row.dec_deg == pytest.approx(center.dec.deg, abs=1e-9)
    assert row.radius_deg == pytest.approx(expected_radius, rel=1e-9)
    # Sanity: the radius brackets the half-diagonal of a ~(5x4)-pixel field at 2e-4 deg/pixel.
    assert 0 < row.radius_deg < 0.01


# --------------------------------------------------------------------------------------
# Build -> persist -> read round-trips through the same schema validation as a published catalog
# --------------------------------------------------------------------------------------
def test_write_catalog_round_trips_through_read_catalog(store_dir):
    from astrolyze.collection.catalog import read_catalog
    from astrolyze.collection.scan import write_catalog

    out = write_catalog(store_dir)
    assert out.name == "catalog.parquet"
    assert out.exists()

    # The persisted catalog reads back through the SAME reader a published corpus uses.
    catalog = read_catalog(str(store_dir))
    assert len(catalog.rows) == 3
    assert {r.object for r in catalog.rows} == {"NGC3521", "NGC0628"}
    # Footprint columns survive the parquet round-trip.
    assert all(r.ra_deg is not None and r.radius_deg is not None for r in catalog.rows)


def test_write_catalog_refuses_to_clobber_without_overwrite(store_dir):
    from astrolyze.collection.scan import write_catalog

    write_catalog(store_dir)
    with pytest.raises(FileExistsError):
        write_catalog(store_dir)
    # ...but overwrite=True rebuilds it.
    assert write_catalog(store_dir, overwrite=True).exists()


# --------------------------------------------------------------------------------------
# mocpy is an OPTIONAL consumer extra: the bare path works without it; MOC lights up with it
# --------------------------------------------------------------------------------------
def test_build_works_without_mocpy_center_radius_only(store_dir):
    """The bare install must work end-to-end: a scan without mocpy declares schema 1.0 and carries
    center+radius footprints. This is the default state of the test venv (no mocpy)."""
    from astrolyze.collection import scan

    catalog = scan.build_catalog(store_dir)
    if scan._scan_schema_version() == "1.0":
        # No mocpy in the environment: every row has its center+radius prefilter and the catalog is
        # a plain v1.0 one. No moc column is required.
        assert catalog.schema_version == "1.0"
        assert all(r.radius_deg is not None for r in catalog.rows)
    else:
        # mocpy IS installed: the scan is a v1.1 catalog instead (covered by the MOC test below).
        assert catalog.schema_version == "1.1"


def test_scan_schema_version_reflects_mocpy_presence(store_dir):
    """The declared version is 1.1 iff mocpy is importable; otherwise 1.0 — both valid."""
    from astrolyze.collection import scan

    expected = "1.1" if _mocpy_present() else "1.0"
    assert scan.build_catalog(store_dir).schema_version == expected


def test_moc_column_is_populated_when_mocpy_installed(store_dir):
    pytest.importorskip("mocpy")  # guard the MOC assertion behind mocpy's presence
    from astrolyze.collection.catalog import read_catalog
    from astrolyze.collection.scan import scan_directory, write_catalog

    # The in-memory ScanResults carry the exact-coverage MOC string when mocpy is present.
    results = scan_directory(store_dir)
    assert all(r.moc is not None for r in results)
    from mocpy import MOC

    # ...and it round-trips through mocpy's ascii form (the spec's portable serialization).
    assert MOC.from_string(results[0].moc, format="ascii") is not None

    # Persisted, the catalog declares 1.1 and the trailing moc column is present (additive).
    write_catalog(store_dir)
    catalog = read_catalog(
        str(store_dir)
    )  # a v1.0-aware reader ignores the trailing moc
    assert catalog.schema_version == "1.1"


def _mocpy_present() -> bool:
    try:
        import mocpy  # noqa: F401
    except ImportError:
        return False
    return True


# --------------------------------------------------------------------------------------
# Incomplete / non-store directories are skipped with a warning, not silently mis-catalogued
# --------------------------------------------------------------------------------------
def test_store_without_celestial_wcs_is_skipped_with_warning(store_dir):
    """A store with the astrolyze backend marker but no celestial WCS (incomplete attrs) has no
    footprint source: it is skipped with a clear ScanWarning, not catalogued with invented geometry
    (PRD #56 acceptance: incomplete attrs -> warning, never silent mis-catalogue)."""
    from astrolyze.collection.scan import ScanWarning, build_catalog
    from astrolyze.io import save

    # A 1D "store" carrying the full astrolyze schema (so it is a discovery candidate) but no
    # celestial WCS header — _footprint can't be computed.
    header = fits.Header()
    md = Metadata.from_header(header)
    save(np.arange(5, dtype="float32"), md, store_dir / "l1" / "wcsless", format="zarr")

    with pytest.warns(ScanWarning, match="sky footprint"):
        catalog = build_catalog(store_dir)
    # The real 3 stores are catalogued; the WCS-less one is not (skipped, not entered with garbage).
    assert len(catalog.rows) == 3
    assert all("wcsless" not in r.store_path for r in catalog.rows)


def test_non_astrolyze_group_is_not_discovered_as_a_store(store_dir):
    """A bare Zarr group with no astrolyze backend marker is not a cube store: it is filtered out at
    discovery (not even a load candidate), so it never appears as a phantom row."""
    from astrolyze.collection.scan import build_catalog

    bogus = store_dir / "l1" / "bogus.zarr"
    bogus.mkdir(parents=True)
    (bogus / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')

    catalog = build_catalog(store_dir)
    assert len(catalog.rows) == 3
    assert all("bogus" not in r.store_path for r in catalog.rows)


def test_companion_subgroups_are_not_catalogued_as_stores(store_dir):
    """A store's noise/validity companion subgroups are Zarr groups carrying their own zarr.json but
    *not* the astrolyze cube-store provenance marker; the scan must not enter them as phantom cubes.
    """
    from astrolyze.collection.scan import build_catalog

    # A noise/ companion subgroup as the astrolyze backend writes it: a group with method/version
    # attrs but no {"provenance": {"backend": "zarr"}} cube-store marker.
    store = next((store_dir / "l1" / "things").glob("*.zarr"))
    noise = store / "noise"
    noise.mkdir()
    (noise / "zarr.json").write_text(
        '{"zarr_format": 3, "node_type": "group", '
        '"attributes": {"method": "mad_std", "version": 1}}'
    )

    catalog = build_catalog(store_dir)
    # Still exactly the 3 cube stores — the noise companion is not a fourth row.
    assert len(catalog.rows) == 3
    assert all("/noise" not in r.store_path for r in catalog.rows)


# --------------------------------------------------------------------------------------
# Collection.open on a catalog-less directory transparently scans and works
# --------------------------------------------------------------------------------------
def test_collection_open_falls_back_to_scan_on_catalog_less_directory(store_dir):
    from astrolyze.collection import Collection

    collection = Collection.open(str(store_dir))
    # A working collection built by scanning: object-first list + flat records, no catalog file.
    assert not (store_dir / "catalog.parquet").exists()
    assert len(collection.records) == 3
    by_object = {s.object for s in collection.list()}
    assert by_object == {"NGC3521", "NGC0628"}


def test_scanned_collection_records_open_into_lazy_cubes(store_dir):
    import dask

    from astrolyze.collection import Collection

    collection = Collection.open(str(store_dir))
    record = next(r for r in collection.records if r.object == "NGC0628")
    cube = record.open()
    assert isinstance(cube, Cube)
    # Lazy: the scanned record opens through the same Zarr backend, so it stays dask-backed (#23).
    assert dask.is_dask_collection(cube._sc._data)
    # Origin provenance is stamped from the scanned catalog's version, exactly as a published one.
    assert cube.metadata.origin_store_uri is not None
    assert cube.metadata.origin_store_uri.endswith(".zarr")


def test_published_catalog_still_takes_precedence_over_scan(store_dir):
    """When a catalog.parquet IS present the reader path is used (scan is only the fallback)."""
    from astrolyze.collection import Collection
    from astrolyze.collection.scan import write_catalog

    write_catalog(store_dir)
    collection = Collection.open(str(store_dir))
    # The declared version comes from the persisted file's header, not a re-scan.
    expected = "1.1" if _mocpy_present() else "1.0"
    assert collection.catalog_version == expected
    assert len(collection.records) == 3
