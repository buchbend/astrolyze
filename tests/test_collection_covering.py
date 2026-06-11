"""Tests for sky-coverage search — ``Collection.covering(SkyCoord)`` (issue #62).

Written first (red/green TDD). On top of the #57 facade and the #60 query layer, #62 adds the
sky-coverage operation: ``Collection.covering(SkyCoord)`` returns the sub-collection of stores whose
footprint *contains* a position. The whole point is the **coverage semantics** — and they are
deliberately two-layered (PRD #56 user stories 5/6/7, catalog-schema spec §2.4/§2.5/§5):

1. **The catalog footprint is a PREFILTER ONLY.** Center+radius (``ra_deg`` / ``dec_deg`` /
   ``radius_deg``) is always available and works on a bare install. With mocpy + a ``moc`` column,
   an exact-MOC prefilter narrows the candidates further *before any store is opened*.
2. **The store's own WCS is the FINAL authority.** Every surviving candidate is opened and the
   position is tested against that store's real celestial pixel grid. A candidate that passes the
   radius/MOC prefilter but whose actual WCS does **not** contain the position is REJECTED — the
   catalog is an index, never the geometry authority.

The three mandatory acceptance cases:

- **HIT** — a position genuinely inside a store's WCS footprint is returned.
- **MISS** — a position far outside every footprint returns an empty collection.
- **PREFILTER-FALSE-POSITIVE** — a position inside a store's center+radius *circle* but outside its
  actual (wide, short) WCS *pixel rectangle* must be REJECTED by the WCS authority. This is the
  proof the catalog footprint is not trusted as geometry (user story 6).

Plus: the MOC prefilter path runs when mocpy is present (``importorskip``), the bare center+radius
path returns correct final answers without mocpy, and the CLI ``collection covering`` smoke-tests.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from astrolyze.core import Cube

warnings.filterwarnings("ignore", module="spectral_cube")

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

REST_CO21 = 230.538e9  # CO(2-1) Hz
REST_HI = 1.420405751e9  # HI 21cm Hz

# Two well-separated field centres so a MISS position can sit far from both.
CENTRE_A = (170.0, -0.04)  # the NGC3521-ish field
CENTRE_B = (200.0, 30.0)  # a far-away field (the MISS lives near here, off both grids)


# --------------------------------------------------------------------------------------
# Synthetic corpus helpers (mirroring the #57/#60 patterns: real SIN WCS, lazy Zarr stores)
# --------------------------------------------------------------------------------------
def _header(*, obj, rest_hz, species, telescope, bmaj_arcsec, crval1, crval2, nx, ny):
    """A schema-complete cube header with a real celestial+spectral WCS of the requested shape."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", crval1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", crval2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
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


def _write_store(
    root,
    rel_dir,
    *,
    obj,
    rest_hz,
    species,
    telescope,
    bmaj_arcsec,
    crval1,
    crval2,
    nx,
    ny,
):
    """Write one synthetic astrolyze Zarr store with the given shape; return its store dir."""
    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    header, beam = _header(
        obj=obj,
        rest_hz=rest_hz,
        species=species,
        telescope=telescope,
        bmaj_arcsec=bmaj_arcsec,
        crval1=crval1,
        crval2=crval2,
        nx=nx,
        ny=ny,
    )
    data = np.arange(3 * ny * nx, dtype="float32").reshape(3, ny, nx)
    loaded = LoadedData(
        data=data,
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )
    target = root / rel_dir
    target.mkdir(parents=True, exist_ok=True)
    Cube.from_loaded(loaded).to_zarr(target)
    return next(target.glob("*.zarr")), beam


def _wcs_footprint(crval1, crval2, nx, ny):
    """The (ra_deg, dec_deg, radius_deg, corner_coords) prefilter a store's WCS yields.

    Mirrors :func:`astrolyze.collection.scan._footprint` so the catalog rows the tests build carry
    exactly the prefilter the producer would compute — center pixel + max-corner separation."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", crval1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", crval2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    celestial = WCS(h)
    corners = celestial.calc_footprint(axes=(nx, ny))
    center = celestial.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    corner_coords = SkyCoord(corners[:, 0] * u.deg, corners[:, 1] * u.deg)
    radius_deg = float(center.separation(corner_coords).max().to_value(u.deg))
    return float(center.ra.deg), float(center.dec.deg), radius_deg, corner_coords


def _moc_string(corner_coords):
    """The exact-footprint MOC ascii string (mocpy present) or ``None`` — matches scan._footprint_moc."""
    try:
        from mocpy import MOC
    except ImportError:
        return None
    moc = MOC.from_polygon_skycoord(corner_coords, max_depth=12)
    return moc.to_string(format="ascii")


# Three stores: two square fields at CENTRE_A (a matched multi-survey source) and one square
# field far away at CENTRE_B. Plus a deliberately WIDE-AND-SHORT field at CENTRE_A used only by
# the false-positive test (big corner gaps between its bounding circle and its pixel rectangle).
_SPECS = [
    dict(
        rel_dir="l1/heracles",
        obj="NGC3521",
        rest_hz=REST_CO21,
        species="CO",
        telescope="IRAM30M",
        survey="HERACLES",
        transition="2-1",
        bmaj_arcsec=13.0,
        crval1=CENTRE_A[0],
        crval2=CENTRE_A[1],
        nx=21,
        ny=21,
    ),
    dict(
        rel_dir="l1/phangs",
        obj="NGC3521",
        rest_hz=REST_CO21,
        species="CO",
        telescope="ALMA",
        survey="PHANGS-ALMA",
        transition="2-1",
        bmaj_arcsec=1.5,
        crval1=CENTRE_A[0],
        crval2=CENTRE_A[1],
        nx=21,
        ny=21,
    ),
    dict(
        rel_dir="l1/things",
        obj="NGC0628",
        rest_hz=REST_HI,
        species="HI",
        telescope="VLA",
        survey="THINGS",
        transition="21cm",
        bmaj_arcsec=11.0,
        crval1=CENTRE_B[0],
        crval2=CENTRE_B[1],
        nx=21,
        ny=21,
    ),
]

# The false-positive store: 81x5 px (wide, short) at CENTRE_A. Its bounding circle has a radius set
# by the long x-axis, so a point offset diagonally sits inside the circle but far outside the y-rows.
_WIDE_SPEC = dict(
    rel_dir="l1/wide",
    obj="NGC9000",
    rest_hz=REST_CO21,
    species="CO",
    telescope="SURVEY-W",
    survey="WIDE",
    transition="1-0",
    bmaj_arcsec=8.0,
    crval1=CENTRE_A[0],
    crval2=CENTRE_A[1],
    nx=81,
    ny=5,
)


def _row(spec, store, root, beam, *, with_moc):
    """A catalog row dict for *spec*'s store, carrying its WCS-derived footprint (and optional MOC)."""
    ra, dec, radius, corners = _wcs_footprint(
        spec["crval1"], spec["crval2"], spec["nx"], spec["ny"]
    )
    rel = store.relative_to(root).as_posix()
    row = {
        "object": spec["obj"],
        "survey": spec["survey"],
        "telescope": spec["telescope"],
        "species": spec["species"],
        "transition": spec["transition"],
        "rest_frequency_hz": float(spec["rest_hz"]),
        "beam_major_arcsec": beam.major.to_value(u.arcsec),
        "beam_minor_arcsec": beam.minor.to_value(u.arcsec),
        "beam_pa_deg": beam.pa.to_value(u.deg),
        "bunit": "K",
        "store_path": rel,
        "content_checksum": f"sha256:{spec['rel_dir']}",
        "ra_deg": ra,
        "dec_deg": dec,
        "radius_deg": radius,
        "catalog_schema_version": "1.1" if with_moc else "1.0",
    }
    if with_moc:
        row["moc"] = _moc_string(corners)
    return row


def _catalog_table(rows, *, schema_version):
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    columns = {name: [row.get(name) for row in rows] for name in CATALOG_COLUMNS}
    if schema_version == "1.1":
        # Insert the v1.1 moc column at its contract position: after radius_deg, before the version.
        ordered = {}
        for name, values in columns.items():
            if name == "catalog_schema_version":
                ordered["moc"] = [row.get("moc") for row in rows]
            ordered[name] = values
        columns = ordered
    table = pa.table(columns)
    return table.replace_schema_metadata({"catalog_schema_version": schema_version})


def _build_corpus(tmp_path, *, with_wide, with_moc):
    """Write the synthetic corpus + its catalog.parquet; return the root path."""
    root = tmp_path / "ism_corpus"
    root.mkdir(parents=True)
    specs = list(_SPECS)
    if with_wide:
        specs.append(_WIDE_SPEC)
    rows = []
    for spec in specs:
        store, beam = _write_store(
            root,
            spec["rel_dir"],
            obj=spec["obj"],
            rest_hz=spec["rest_hz"],
            species=spec["species"],
            telescope=spec["telescope"],
            bmaj_arcsec=spec["bmaj_arcsec"],
            crval1=spec["crval1"],
            crval2=spec["crval2"],
            nx=spec["nx"],
            ny=spec["ny"],
        )
        rows.append(_row(spec, store, root, beam, with_moc=with_moc))
    version = "1.1" if with_moc else "1.0"
    pq.write_table(
        _catalog_table(rows, schema_version=version), root / "catalog.parquet"
    )
    return root


def _mocpy_present() -> bool:
    try:
        import mocpy  # noqa: F401
    except ImportError:
        return False
    return True


@pytest.fixture
def corpus(tmp_path):
    """A published corpus of 3 square fields (NGC3521 x2 at A, NGC0628 at B), center+radius only."""
    return _build_corpus(tmp_path, with_wide=False, with_moc=False)


# --------------------------------------------------------------------------------------
# (a) HIT: a position genuinely inside a store's WCS footprint is returned
# --------------------------------------------------------------------------------------
def test_covering_returns_stores_whose_wcs_contains_the_position(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    # CENTRE_A is the field centre of both NGC3521 stores -> both must be covering.
    position = SkyCoord(CENTRE_A[0] * u.deg, CENTRE_A[1] * u.deg)
    covering = collection.covering(position)

    objects = {r.object for r in covering.records}
    assert objects == {"NGC3521"}
    # Both surveys of the matched source are returned (the multi-survey acceptance).
    surveys = {r.survey for r in covering.records}
    assert surveys == {"HERACLES", "PHANGS-ALMA"}
    assert len(covering.records) == 2


def test_covering_returns_a_composable_subcollection(corpus):
    """covering() returns a Collection (like query()): the result composes with the facade."""
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    position = SkyCoord(CENTRE_A[0] * u.deg, CENTRE_A[1] * u.deg)
    covering = collection.covering(position)
    assert isinstance(covering, Collection)
    # Composable: query/list/describe all work on the covering sub-collection.
    assert covering.query(survey="HERACLES").records[0].survey == "HERACLES"
    assert {s.object for s in covering.list()} == {"NGC3521"}


# --------------------------------------------------------------------------------------
# (b) MISS: a position far outside every footprint returns empty
# --------------------------------------------------------------------------------------
def test_covering_far_outside_every_footprint_is_empty(corpus):
    from astrolyze.collection import Collection

    collection = Collection.open(str(corpus))
    # A position degrees away from both field centres -> outside every radius and every WCS.
    far = SkyCoord(10.0 * u.deg, -60.0 * u.deg)
    covering = collection.covering(far)
    assert covering.records == []
    assert covering.list() == []


# --------------------------------------------------------------------------------------
# (c) PREFILTER-FALSE-POSITIVE: inside the radius circle, outside the WCS rectangle -> rejected
# --------------------------------------------------------------------------------------
def test_covering_rejects_radius_false_positive_by_wcs_authority(tmp_path):
    """A wide-short field's bounding circle contains a corner-gap position its pixel grid does NOT.

    The catalog footprint (center+radius) is a coarse circle around a 81x5-px field, so a point
    offset diagonally from the centre sits *inside the circle* but well *above the 5 image rows*.
    The radius prefilter therefore admits it as a candidate — and the store's own WCS must REJECT
    it. This is the proof the catalog is an index, never the geometry authority (user story 6)."""
    from astrolyze.collection import Collection

    root = _build_corpus(tmp_path, with_wide=True, with_moc=False)
    collection = Collection.open(str(root))

    ra, dec, radius, _ = _wcs_footprint(
        _WIDE_SPEC["crval1"], _WIDE_SPEC["crval2"], _WIDE_SPEC["nx"], _WIDE_SPEC["ny"]
    )
    center = SkyCoord(ra * u.deg, dec * u.deg)
    # A point at 0.9*radius along a 45-deg diagonal: inside the circle, off the short y-axis.
    false_positive = center.directional_offset_by(
        position_angle=45 * u.deg, separation=0.9 * radius * u.deg
    )

    # Sanity 1: it IS inside the catalog radius prefilter (so the prefilter would admit it).
    assert center.separation(false_positive).to_value(u.deg) <= radius

    # Sanity 2: it is NOT inside the wide store's real WCS pixel footprint (the authority).
    wide_record = next(r for r in collection.records if r.object == "NGC9000")
    cube = wide_record.open()
    cel = cube._sc.wcs.celestial
    px, py = cel.world_to_pixel(false_positive)
    ny, nx = cube.shape[1], cube.shape[2]
    in_bounds = (-0.5 <= float(px) <= nx - 0.5) and (-0.5 <= float(py) <= ny - 0.5)
    assert not in_bounds

    # The contract: covering() must NOT return the wide store for this radius false positive.
    covering = collection.covering(false_positive)
    assert all(r.object != "NGC9000" for r in covering.records)


def test_covering_prefilter_candidate_count_exceeds_wcs_matches(tmp_path):
    """The radius prefilter admits MORE candidates than the WCS confirms — the two layers differ.

    Directly exercises the seam: for the false-positive position the radius prefilter must yield the
    wide store as a *candidate* while the final WCS authority drops it (candidates > final matches)."""
    from astrolyze.collection import Collection
    from astrolyze.collection import _facade

    root = _build_corpus(tmp_path, with_wide=True, with_moc=False)
    collection = Collection.open(str(root))

    ra, dec, radius, _ = _wcs_footprint(
        _WIDE_SPEC["crval1"], _WIDE_SPEC["crval2"], _WIDE_SPEC["nx"], _WIDE_SPEC["ny"]
    )
    center = SkyCoord(ra * u.deg, dec * u.deg)
    false_positive = center.directional_offset_by(
        position_angle=45 * u.deg, separation=0.9 * radius * u.deg
    )

    candidates = _facade._prefilter_candidates(collection.records, false_positive)
    final = collection.covering(false_positive).records
    candidate_objs = {r.object for r in candidates}
    final_objs = {r.object for r in final}
    # The wide store is a prefilter candidate but NOT a final WCS match.
    assert "NGC9000" in candidate_objs
    assert "NGC9000" not in final_objs


# --------------------------------------------------------------------------------------
# MOC prefilter path (mocpy present) vs bare center+radius path (mocpy absent)
# --------------------------------------------------------------------------------------
def test_bare_center_radius_path_returns_correct_final_answers(tmp_path, monkeypatch):
    """Without mocpy the prefilter is center+radius only — and the final answers are still correct.

    Forces the no-mocpy branch (monkeypatch the lazy import to fail) and re-checks HIT + MISS, so
    the bare-install path is proven end-to-end regardless of whether mocpy is in the test venv."""
    import builtins

    from astrolyze.collection import Collection

    real_import = builtins.__import__

    def _no_mocpy(name, *args, **kwargs):
        if name == "mocpy" or name.startswith("mocpy."):
            raise ImportError("mocpy disabled for the bare-path test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_mocpy)

    root = _build_corpus(tmp_path, with_wide=False, with_moc=False)
    collection = Collection.open(str(root))

    hit = SkyCoord(CENTRE_A[0] * u.deg, CENTRE_A[1] * u.deg)
    assert {r.object for r in collection.covering(hit).records} == {"NGC3521"}

    miss = SkyCoord(10.0 * u.deg, -60.0 * u.deg)
    assert collection.covering(miss).records == []


def test_moc_prefilter_runs_when_mocpy_present(tmp_path):
    """With mocpy + a moc column, the MOC prefilter narrows candidates before any store opens.

    Guarded behind importorskip. The corpus carries the v1.1 moc column; the position is the field
    centre (a HIT). The MOC prefilter must run (the moc column is surfaced on the row) and the final
    answer is identical to the radius path — MOC sharpens the prefilter, the WCS still decides."""
    pytest.importorskip("mocpy")
    from astrolyze.collection import Collection
    from astrolyze.collection import _facade

    root = _build_corpus(tmp_path, with_wide=False, with_moc=True)
    collection = Collection.open(str(root))
    # The reader carried the v1.1 moc column through onto the rows (not dropped as trailing-unknown).
    assert collection.catalog_version == "1.1"
    assert all(r.row.moc is not None for r in collection.records)

    position = SkyCoord(CENTRE_A[0] * u.deg, CENTRE_A[1] * u.deg)
    # The MOC prefilter path is taken (mocpy present + moc column populated).
    candidates = _facade._prefilter_candidates(collection.records, position)
    assert {r.object for r in candidates} == {"NGC3521"}
    assert {r.object for r in collection.covering(position).records} == {"NGC3521"}


def test_moc_prefilter_still_defers_to_wcs_authority(tmp_path):
    """Even with the exact MOC, the wide-field false positive is decided by the WCS, not the index.

    The MOC of a wide-short polygon excludes most of the corner gap, so the false-positive point may
    already be MOC-rejected — but the contract is that *whatever* the prefilter admits, the WCS is
    the final authority. The wide store must never be a final match for the corner-gap position."""
    pytest.importorskip("mocpy")
    from astrolyze.collection import Collection

    root = _build_corpus(tmp_path, with_wide=True, with_moc=True)
    collection = Collection.open(str(root))

    ra, dec, radius, _ = _wcs_footprint(
        _WIDE_SPEC["crval1"], _WIDE_SPEC["crval2"], _WIDE_SPEC["nx"], _WIDE_SPEC["ny"]
    )
    center = SkyCoord(ra * u.deg, dec * u.deg)
    false_positive = center.directional_offset_by(
        position_angle=45 * u.deg, separation=0.9 * radius * u.deg
    )
    covering = collection.covering(false_positive)
    assert all(r.object != "NGC9000" for r in covering.records)


# --------------------------------------------------------------------------------------
# The reader must carry the v1.1 ``moc`` column through (not drop it as a trailing unknown)
# --------------------------------------------------------------------------------------
def test_reader_surfaces_moc_column_on_catalog_row(tmp_path):
    """A v1.1 catalog's moc column reaches CatalogRow.moc — covering's MOC prefilter needs it.

    Before #62 the reader dropped moc as a trailing-unknown column; covering teaches it to carry it
    through so the optional MOC prefilter has a value to deserialize."""
    from astrolyze.collection.catalog import read_catalog

    root = _build_corpus(tmp_path, with_wide=False, with_moc=True)
    catalog = read_catalog(str(root))
    assert catalog.schema_version == "1.1"
    if _mocpy_present():
        assert all(row.moc is not None for row in catalog.rows)
    # A v1.0 catalog has no moc column -> the field is simply None (back-compatible).
    bare_root = _build_corpus(tmp_path / "bare", with_wide=False, with_moc=False)
    bare = read_catalog(str(bare_root))
    assert bare.schema_version == "1.0"
    assert all(row.moc is None for row in bare.rows)
