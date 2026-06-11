"""Tests for the Stack container — ``Collection.stack`` + :class:`Stack` (issue #64, PRD #56 stage 1).

Written first (red/green TDD). On top of #62 ``covering`` and #58 ``Cube.cutout``, #64 adds the
stage-1 stack: ``Collection.stack(position_or_sources, size=...)`` gathers sky-coordinate cutouts
from every covering cube into a :class:`~astrolyze.collection.Stack` — a container that is **always
safe to construct and browse**, however heterogeneous its members. Co-addition is a *separate*
stage (#65) and is deliberately NOT exercised here.

The obligations pinned (the acceptance criteria of #64):

1. ``stack(SkyCoord)`` gathers one member per covering cube, each a correct cutout (user story 10).
2. A source list (catalog object names and/or SkyCoords) builds a multi-target sample; an unknown
   name raises (user story 11). Names resolve against THIS catalog — no SIMBAD.
3. ``filter()`` narrows to a homogeneous subset and preserves selection provenance; an unknown
   criterion raises (user story 13).
4. ``plot_grid()`` renders one house-style panel per member, labelled by identity (user story 12) —
   structure asserted, no image goldens.
5. Members carry origin provenance (origin_store_uri + catalog version on the cutout's Metadata).
6. The stack records selection provenance (what position/sources + size produced it).
7. Browsing a HETEROGENEOUS stack always works — no homogeneity precondition at this stage.
8. The read-only homogeneity report is the #65 coadd seam, not a stage-1 gate.

The synthetic corpus mirrors ``tests/test_collection_covering.py``: NGC3521 in HERACLES (CO 2-1)
and PHANGS-ALMA (CO 2-1) at one field centre (the matched multi-survey source), plus NGC0628 in
THINGS (HI 21cm) at a far field — built with real SIN WCS lazy Zarr stores + a catalog.parquet.
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

CENTRE_A = (170.0, -0.04)  # the NGC3521-ish field (two surveys cover it)
CENTRE_B = (200.0, 30.0)  # a far-away field (NGC0628 / THINGS)

PIX_DEG = 2e-4
PIX_ARCSEC = (
    PIX_DEG * 3600.0
)  # 0.72 arcsec/pixel — a 7*0.72" stamp brackets ~7 px (fits a 21px store)
STAMP = 7 * PIX_ARCSEC * u.arcsec


# --------------------------------------------------------------------------------------
# Synthetic corpus helpers (mirroring test_collection_covering.py: real SIN WCS, lazy Zarr)
# --------------------------------------------------------------------------------------
def _header(*, obj, rest_hz, species, telescope, bmaj_arcsec, crval1, crval2):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", crval1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -PIX_DEG, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", crval2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = PIX_DEG, 1.0, "deg"
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
    root, rel_dir, *, obj, rest_hz, species, telescope, bmaj_arcsec, crval1, crval2
):
    """Write one synthetic astrolyze Zarr store centred at (crval1, crval2); return its store dir."""
    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    nx = ny = 21
    header, beam = _header(
        obj=obj,
        rest_hz=rest_hz,
        species=species,
        telescope=telescope,
        bmaj_arcsec=bmaj_arcsec,
        crval1=crval1,
        crval2=crval2,
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


def _wcs_footprint(crval1, crval2):
    """The (ra_deg, dec_deg, radius_deg) prefilter a 21x21 store's WCS yields (mirrors scan._footprint)."""
    nx = ny = 21
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", crval1
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -PIX_DEG, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", crval2
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = PIX_DEG, 1.0, "deg"
    celestial = WCS(h)
    corners = celestial.calc_footprint(axes=(nx, ny))
    center = celestial.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    corner_coords = SkyCoord(corners[:, 0] * u.deg, corners[:, 1] * u.deg)
    radius_deg = float(center.separation(corner_coords).max().to_value(u.deg))
    return float(center.ra.deg), float(center.dec.deg), radius_deg


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
    ),
]


def _row(spec, store, root, beam):
    ra, dec, radius = _wcs_footprint(spec["crval1"], spec["crval2"])
    return {
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
        "store_path": store.relative_to(root).as_posix(),
        "content_checksum": f"sha256:{spec['rel_dir']}",
        "ra_deg": ra,
        "dec_deg": dec,
        "radius_deg": radius,
        "catalog_schema_version": "1.0",
    }


def _catalog_table(rows):
    from astrolyze.collection.catalog import CATALOG_COLUMNS

    columns = {name: [row.get(name) for row in rows] for name in CATALOG_COLUMNS}
    table = pa.table(columns)
    return table.replace_schema_metadata({"catalog_schema_version": "1.0"})


@pytest.fixture
def corpus(tmp_path):
    """A published corpus: NGC3521 (HERACLES CO + PHANGS CO) at A, NGC0628 (THINGS HI) at B."""
    root = tmp_path / "ism_corpus"
    root.mkdir(parents=True)
    rows = []
    for spec in _SPECS:
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
        )
        rows.append(_row(spec, store, root, beam))
    pq.write_table(_catalog_table(rows), root / "catalog.parquet")
    return root


@pytest.fixture
def collection(corpus):
    from astrolyze.collection import Collection

    return Collection.open(str(corpus))


def _centre_a():
    # The field's FOOTPRINT centre (pixel ~10,10), not CRVAL (which sits at the corner pixel 0,0):
    # a stamp must be centred on the data, and this is also the centre the catalog row records.
    ra, dec, _ = _wcs_footprint(CENTRE_A[0], CENTRE_A[1])
    return SkyCoord(ra * u.deg, dec * u.deg)


def _centre_b():
    ra, dec, _ = _wcs_footprint(CENTRE_B[0], CENTRE_B[1])
    return SkyCoord(ra * u.deg, dec * u.deg)


# --------------------------------------------------------------------------------------
# (1) stack(SkyCoord) gathers one member per covering cube, each a correct cutout
# --------------------------------------------------------------------------------------
def test_stack_position_gathers_one_member_per_covering_cube(collection):
    from astrolyze.collection import Stack

    stack = collection.stack(_centre_a(), size=STAMP)
    assert isinstance(stack, Stack)
    # CENTRE_A is covered by both NGC3521 surveys -> two members; THINGS (far) contributes none.
    assert len(stack) == 2
    surveys = {member.survey for member in stack}
    assert surveys == {"HERACLES", "PHANGS-ALMA"}
    assert {member.object for member in stack} == {"NGC3521"}


def test_stack_members_are_cutouts_with_full_spectral_axis(collection):
    """Each member is a real cutout: shrunk on the sky to the stamp, full spectral axis kept."""
    stack = collection.stack(_centre_a(), size=STAMP)
    for member in stack:
        # The parent store is 3 channels x 21 x 21; the stamp brackets ~7 sky px, all channels.
        assert member.cube.shape[0] == 3
        assert member.cube.shape[1] < 21 and member.cube.shape[2] < 21
        assert member.cube.shape[1] <= 9 and member.cube.shape[2] <= 9


def test_stack_uncovered_position_is_empty_but_valid(collection):
    """A position nothing covers yields an empty (but fully valid, browsable) stack — not an error."""
    far = SkyCoord(10.0 * u.deg, -60.0 * u.deg)
    stack = collection.stack(far, size=STAMP)
    assert len(stack) == 0
    assert list(stack) == []
    # An empty stack is vacuously homogeneous and still carries its selection provenance.
    assert stack.is_homogeneous
    assert stack.selection.kind == "position"


# --------------------------------------------------------------------------------------
# (2) source-list input: names and coords build a multi-target sample; unknown names raise
# --------------------------------------------------------------------------------------
def test_stack_source_list_of_names_builds_multi_target_sample(collection):
    """A list of catalog object NAMES builds a sample spanning both sources' covering cubes."""
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    # NGC3521 -> 2 covering members (HERACLES + PHANGS); NGC0628 -> 1 (THINGS). Total 3.
    assert len(stack) == 3
    assert {member.object for member in stack} == {"NGC3521", "NGC0628"}
    assert {member.survey for member in stack} == {"HERACLES", "PHANGS-ALMA", "THINGS"}


def test_stack_source_list_mixes_names_and_coords(collection):
    """A list may mix a catalog NAME and a bare SkyCoord — both resolve to their covering cubes."""
    stack = collection.stack(["NGC3521", _centre_b()], size=STAMP)
    assert {member.object for member in stack} == {"NGC3521", "NGC0628"}
    assert stack.selection.kind == "sources"


def test_stack_unknown_name_raises_descriptive_keyerror(collection):
    """An object name absent from THIS catalog raises (no SIMBAD; a typo must not be silent)."""
    with pytest.raises(KeyError) as exc:
        collection.stack(["NGC9999"], size=STAMP)
    message = str(exc.value)
    assert "NGC9999" in message
    # The error names the known objects so the typo is obvious.
    assert "NGC3521" in message or "NGC0628" in message


def test_stack_single_name_as_bare_string_is_one_source(collection):
    """A bare string source is treated as one name (not iterated character-by-character)."""
    stack = collection.stack("NGC3521", size=STAMP)
    assert {member.object for member in stack} == {"NGC3521"}
    assert len(stack) == 2


# --------------------------------------------------------------------------------------
# (3) filter() narrows to a homogeneous subset, preserves provenance; unknown key raises
# --------------------------------------------------------------------------------------
def test_filter_narrows_to_a_homogeneous_subset(collection):
    """filter(survey=...) selects one survey's members; the result is a smaller Stack."""
    from astrolyze.collection import Stack

    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    co = stack.filter(species="CO")
    assert isinstance(co, Stack)
    assert {member.species for member in co} == {"CO"}
    # CO is the two NGC3521 surveys; the HI member is dropped.
    assert len(co) == 2
    assert {member.survey for member in co} == {"HERACLES", "PHANGS-ALMA"}


def test_filter_is_conjunctive(collection):
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    hit = stack.filter(species="CO", survey="PHANGS-ALMA")
    assert len(hit) == 1
    assert hit[0].survey == "PHANGS-ALMA"
    # CO AND THINGS matches nothing (THINGS is HI) -> empty, not an error.
    assert len(stack.filter(species="CO", survey="THINGS")) == 0


def test_filter_preserves_selection_provenance(collection):
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    co = stack.filter(species="CO")
    # The narrowed stack still names the request it descended from.
    assert co.selection is stack.selection
    assert co.selection.kind == "sources"
    assert co.selection.size == STAMP


def test_filter_unknown_criterion_raises(collection):
    stack = collection.stack(_centre_a(), size=STAMP)
    with pytest.raises(ValueError) as exc:
        stack.filter(specie="CO")  # typo'd axis
    message = str(exc.value)
    assert "specie" in message
    assert "species" in message  # the valid axis is named


# --------------------------------------------------------------------------------------
# (4) plot_grid() renders one house-style panel per member, labelled by identity
# --------------------------------------------------------------------------------------
def test_plot_grid_one_panel_per_member(collection):
    import matplotlib

    matplotlib.use("Agg")

    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    fig, axes = stack.plot_grid()
    # One panel per member (3 members), structure only — no image golden.
    assert len(axes) == len(stack) == 3
    titles = {ax.get_title() for ax in axes}
    # Each panel is labelled by member identity (object + survey appear in the title).
    assert any("HERACLES" in t for t in titles)
    assert any("THINGS" in t for t in titles)
    assert any("NGC3521" in t for t in titles)


def test_plot_grid_works_on_heterogeneous_stack(collection):
    """Browsing a mixed-species, mixed-unit stack must always work — no homogeneity precondition."""
    import matplotlib

    matplotlib.use("Agg")

    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    assert not stack.is_homogeneous  # CO + HI, distinct transitions
    fig, axes = stack.plot_grid()  # must not raise despite heterogeneity
    assert len(axes) == 3


def test_plot_grid_empty_stack_raises(collection):
    import matplotlib

    matplotlib.use("Agg")

    far = SkyCoord(10.0 * u.deg, -60.0 * u.deg)
    empty = collection.stack(far, size=STAMP)
    with pytest.raises(ValueError):
        empty.plot_grid()


# --------------------------------------------------------------------------------------
# (5) members carry origin provenance
# --------------------------------------------------------------------------------------
def test_members_carry_origin_provenance(collection):
    """Each member cutout is stamped with its source store URI + catalog version (PRD user story 19)."""
    stack = collection.stack(_centre_a(), size=STAMP)
    for member in stack:
        meta = member.cube.metadata
        assert meta.origin_store_uri is not None
        assert meta.origin_store_uri.endswith(".zarr")
        assert meta.origin_catalog_version == collection.catalog_version
        # The member surfaces the same origin via its convenience accessor.
        assert member.origin_store_uri == meta.origin_store_uri


# --------------------------------------------------------------------------------------
# (6) the stack records selection provenance
# --------------------------------------------------------------------------------------
def test_position_stack_records_selection_provenance(collection):
    position = _centre_a()
    stack = collection.stack(position, size=STAMP)
    sel = stack.selection
    assert sel.kind == "position"
    assert sel.targets is position
    assert sel.size == STAMP
    assert sel.catalog_version == collection.catalog_version


def test_sources_stack_records_selection_provenance(collection):
    sources = ["NGC3521", "NGC0628"]
    stack = collection.stack(sources, size=STAMP)
    sel = stack.selection
    assert sel.kind == "sources"
    assert sel.targets == sources
    assert sel.size == STAMP


# --------------------------------------------------------------------------------------
# (7) sequence protocol + (8) homogeneity report (the #65 coadd seam, read-only here)
# --------------------------------------------------------------------------------------
def test_stack_sequence_protocol(collection):
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    assert len(stack) == 3
    assert stack[0] is stack.members[0]
    # Slicing yields a sub-Stack carrying the same provenance (a browse op, not a new gather).
    head = stack[:2]
    from astrolyze.collection import Stack

    assert isinstance(head, Stack)
    assert len(head) == 2
    assert head.selection is stack.selection
    # Iteration yields members in order.
    assert [m.survey for m in stack] == [m.survey for m in list(stack)]


def test_homogeneous_subset_reports_homogeneous(collection):
    """A CO-only filtered subset is homogeneous (one species/transition/bunit) — the coadd gate."""
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    co = stack.filter(species="CO")
    assert co.is_homogeneous
    report = co.homogeneity_report()
    assert report.species == ("CO",)
    assert report.transitions == ("2-1",)
    assert report.bunits == ("K",)
    assert report.conflicts == ()


def test_heterogeneous_stack_report_names_conflicts(collection):
    """The mixed CO+HI stack is NOT homogeneous; the report names species/transition as conflicts."""
    stack = collection.stack(["NGC3521", "NGC0628"], size=STAMP)
    report = stack.homogeneity_report()
    assert not stack.is_homogeneous
    assert "species" in report.conflicts
    assert "transition" in report.conflicts
    # Both surveys share bunit K -> bunit is not a conflict here.
    assert "bunit" not in report.conflicts
    assert set(report.species) == {"CO", "HI"}


# --------------------------------------------------------------------------------------
# broadcasting: a per-member Cube -> Cube op across the stack (the #65 alignment seam)
# --------------------------------------------------------------------------------------
def test_map_broadcasts_per_member_op_and_preserves_identity(collection):
    """map(fn) applies a Cube->Cube op to each member and returns a new Stack with identity kept."""
    from astrolyze.collection import Stack

    stack = collection.stack(_centre_a(), size=STAMP)

    # A trivial per-member op: take a 1-channel sub-cube (still a Cube, the #65 op shape).
    def first_channel(cube):
        return cube[0:1, :, :]

    mapped = stack.map(first_channel)
    assert isinstance(mapped, Stack)
    assert len(mapped) == len(stack)
    # The op was applied (each member now has 1 channel) and identity is preserved.
    for original, new in zip(stack, mapped):
        assert new.cube.shape[0] == 1
        assert new.object == original.object
        assert new.survey == original.survey
        assert new.species == original.species
    # Selection provenance carried through the broadcast.
    assert mapped.selection is stack.selection


def test_map_requires_fn_to_return_a_cube(collection):
    """A broadcast that does not return a Cube fails loudly (no malformed stack — ADR-0003)."""
    stack = collection.stack(_centre_a(), size=STAMP)
    with pytest.raises(TypeError):
        stack.map(lambda cube: "not a cube")


# --------------------------------------------------------------------------------------
# no co-addition / alignment at this stage (the seam belongs to #65)
# --------------------------------------------------------------------------------------
def test_stack_has_no_coadd_or_alignment_methods_yet(collection):
    """Stage 1 is the container only — coadd/alignment are #65 and must NOT exist here."""
    stack = collection.stack(_centre_a(), size=STAMP)
    for forbidden in ("coadd", "to_common_beam", "to_velocity_grid", "shift_to_rest"):
        assert not hasattr(stack, forbidden), (
            f"{forbidden} is the #65 alignment/coadd stage and must not exist on the stage-1 Stack"
        )
