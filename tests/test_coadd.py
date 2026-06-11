"""Tests for Stack stage-2: alignment + homogeneity-gated noise-weighted coadd (issue #65, PRD #56).

Written first (red/green TDD). On top of the stage-1 :class:`~astrolyze.collection.Stack` container
(#64), #65 adds the *physics* of co-addition as **separate, explicit, auditable** steps and the gate
that makes a meaningless average impossible to produce silently (PRD #56 user stories 14-17):

- ``to_common_beam()`` convolves every member to a common (largest) beam, reusing the existing
  convolution machinery and inheriting its **no-super-resolution** guard (a finer-beam request
  raises);
- ``to_velocity_grid()`` resamples every member onto one shared velocity grid;
- ``shift_to_rest(v_sys=...)`` shifts each member to rest velocity, with v_sys resolved per member
  (catalog -> cube metadata -> user value) and a **raise** when none is available (story 17);
- ``coadd(weights=...)`` combines the members into one Cube, **raising** unless they are homogeneous
  (same species/transition/bunit AND a compatible spatial+spectral grid), inverse-variance weighting
  from the noise companions by default (story 15/16).

The fixtures are synthetic tmp-path Zarr stores with **known** data + **known** σ noise companions, so
the coadd arithmetic is checked against hand-computed values (the house style: assert external
behaviour — returned values and raised errors — never implementation details).
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits
from astropy.wcs import WCS

from astrolyze.collection import CoaddError, Stack, StackMember
from astrolyze.core import Cube, NoiseModel
from astrolyze.core.cube import LossyDirectionError

warnings.filterwarnings("ignore", module="spectral_cube")

pytest.importorskip("pyarrow")

REST_CO21 = 230.538e9  # CO(2-1) Hz
PIX_DEG = 2e-4


# --------------------------------------------------------------------------------------
# Synthetic store + companion helpers (known data, known sigma)
# --------------------------------------------------------------------------------------
def _header(*, obj, species, bmaj_arcsec, bunit="K", crval3=0.0, cdelt3=2000.0):
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 170.0
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -PIX_DEG, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", -0.04
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = PIX_DEG, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", crval3
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = cdelt3, 1.0, "m/s"
    h["OBJECT"] = obj
    h["BUNIT"] = bunit
    h["RESTFRQ"] = (REST_CO21, "Hz")
    h["HIERARCH ASTROLYZE SPECIES"] = species
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    beam = radio_beam.Beam(
        major=bmaj_arcsec * u.arcsec, minor=(bmaj_arcsec - 1) * u.arcsec, pa=0 * u.deg
    )
    for k, v in beam.to_header_keywords().items():
        h[k] = v
    return h, beam


def _cube(header, data):
    """A Cube from a header + a known data array (eager FITS path; 3 channels x ny x nx)."""
    from astrolyze.io.access import LoadedData
    from astrolyze.io.schema import Metadata

    loaded = LoadedData(
        data=np.asarray(data, dtype="float32"),
        wcs=WCS(header),
        metadata=Metadata.from_header(header),
        path="<memory>",
        header_string=header.tostring(),
        header=header,
    )
    return Cube.from_loaded(loaded)


def _store_with_noise(tmp_path, name, *, cube, sigma, species="CO", v_sys=None):
    """Write *cube* to a Zarr store with a constant-σ noise companion; return a StackMember.

    The member's cube is reloaded from the written store and stamped with its origin store URI (the
    provenance the noise-weighted coadd reads to find the companion), exactly as a collection-born
    member is. ``sigma`` is the per-pixel constant σ (in the cube's unit) the companion records;
    ``v_sys`` (a velocity Quantity) is stamped onto the reloaded cube's Metadata when given (the
    store-schema systemic velocity the rest-shift resolves)."""
    from dataclasses import replace

    target = tmp_path / name
    target.mkdir(parents=True, exist_ok=True)
    store = cube.to_zarr(target)

    sigma_map = np.full(cube.shape[1:], float(sigma))
    NoiseModel.from_rms_map(cube, sigma_map).to_zarr_companion(store)

    reloaded = Cube.from_zarr(str(store))
    reloaded.metadata = replace(
        reloaded.metadata,
        origin_store_uri=str(store),
        origin_catalog_version="1.0",
        systemic_velocity=v_sys,
    )
    return StackMember(
        cube=reloaded,
        object=cube.metadata.object,
        survey="SURVEY",
        species=species,
        transition="2-1",
        record=None,
    )


def _selection():
    from astrolyze.collection.stack import Selection

    return Selection(kind="position", targets=None, size=None, catalog_version="1.0")


def _stack(members):
    return Stack(members, _selection())


# --------------------------------------------------------------------------------------
# (A) coadd RAISES on heterogeneous members — three separate raise tests
# --------------------------------------------------------------------------------------
def test_coadd_raises_on_mixed_species(tmp_path):
    """Mixed species is a meaningless average — coadd refuses, naming the conflicting axis."""
    h_co, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    h_hi, _ = _header(obj="SRC", species="HI", bmaj_arcsec=5.0)
    co = _store_with_noise(
        tmp_path,
        "co",
        cube=_cube(h_co, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        species="CO",
    )
    hi = _store_with_noise(
        tmp_path,
        "hi",
        cube=_cube(h_hi, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        species="HI",
    )
    with pytest.raises(CoaddError) as exc:
        _stack([co, hi]).coadd()
    message = str(exc.value)
    assert "species" in message
    assert "filter" in message  # names the alignment step that fixes it


def test_coadd_raises_on_mixed_bunit(tmp_path):
    """Mixed brightness units is a meaningless average — coadd refuses, naming bunit."""
    h_k, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0, bunit="K")
    h_jy, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0, bunit="Jy/beam")
    a = _store_with_noise(
        tmp_path, "k", cube=_cube(h_k, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    b = _store_with_noise(
        tmp_path, "jy", cube=_cube(h_jy, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    with pytest.raises(CoaddError) as exc:
        _stack([a, b]).coadd()
    assert "bunit" in str(exc.value)


def test_coadd_raises_on_unaligned_grids(tmp_path):
    """Same species/unit but DIFFERENT spectral grids — coadd refuses until resampled onto one grid."""
    h_a, _ = _header(
        obj="SRC", species="CO", bmaj_arcsec=5.0, crval3=0.0, cdelt3=2000.0
    )
    # A different velocity axis (channels start elsewhere): same shape, incompatible spectral grid.
    h_b, _ = _header(
        obj="SRC", species="CO", bmaj_arcsec=5.0, crval3=20000.0, cdelt3=2000.0
    )
    a = _store_with_noise(
        tmp_path, "a", cube=_cube(h_a, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    b = _store_with_noise(
        tmp_path, "b", cube=_cube(h_b, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    with pytest.raises(CoaddError) as exc:
        _stack([a, b]).coadd()
    message = str(exc.value)
    assert "spectral grid" in message
    assert "to_velocity_grid" in message  # names the alignment step that fixes it


def test_coadd_raises_on_empty_stack():
    """coadd() on an empty stack has nothing to combine — a descriptive raise, not a NaN cube."""
    with pytest.raises(CoaddError):
        _stack([]).coadd()


def test_coadd_raises_on_unknown_weights(tmp_path):
    """An unknown weighting mode is refused, never silently defaulted (ADR-0003)."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m = _store_with_noise(
        tmp_path, "m", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    with pytest.raises(CoaddError) as exc:
        _stack([m]).coadd(weights="bogus")
    assert "bogus" in str(exc.value)


# --------------------------------------------------------------------------------------
# (B) coadd arithmetic vs HAND-COMPUTED values — uniform and inverse-variance weights
# --------------------------------------------------------------------------------------
def test_coadd_inverse_variance_matches_hand_computed(tmp_path):
    """Two members with KNOWN data + KNOWN σ: assert the combined voxel == hand-computed IVW mean.

    Member 1: d=10, σ=1 -> w=1.   Member 2: d=12, σ=2 -> w=0.25.
    IVW value  = (10*1 + 12*0.25) / (1 + 0.25) = 13/1.25 = 10.4
    IVW sigma  = sqrt(1 / (1 + 0.25))           = sqrt(0.8) = 0.8944271909999159
    """
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m1 = _store_with_noise(
        tmp_path, "m1", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    m2 = _store_with_noise(
        tmp_path, "m2", cube=_cube(h, np.full((3, 4, 4), 12.0)), sigma=2.0
    )

    result = _stack([m1, m2]).coadd(weights="noise")
    values = result._data_quantity.to_value(u.K)
    assert np.allclose(values, 10.4)
    # Combined (propagated) noise rides along as a companion product on the result.
    combined_sigma = result._coadd_sigma.to_value(u.K)
    assert np.allclose(combined_sigma, np.sqrt(0.8))


def test_coadd_inverse_variance_three_members(tmp_path):
    """Three members, distinct σ: IVW mean and propagated σ both match the hand computation.

    d = (8, 10, 12), σ = (1, 1, 2) -> w = (1, 1, 0.25).
    value = (8 + 10 + 12*0.25) / 2.25 = 21 / 2.25 = 9.333333333333334
    sigma = sqrt(1 / 2.25)              = 0.6666666666666666
    """
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    members = [
        _store_with_noise(
            tmp_path, "a", cube=_cube(h, np.full((3, 4, 4), 8.0)), sigma=1.0
        ),
        _store_with_noise(
            tmp_path, "b", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
        ),
        _store_with_noise(
            tmp_path, "c", cube=_cube(h, np.full((3, 4, 4), 12.0)), sigma=2.0
        ),
    ]
    result = _stack(members).coadd(weights="noise")
    assert np.allclose(result._data_quantity.to_value(u.K), 21.0 / 2.25)
    assert np.allclose(result._coadd_sigma.to_value(u.K), np.sqrt(1.0 / 2.25))


def test_coadd_uniform_matches_hand_computed(tmp_path):
    """Uniform weights: the combined voxel is the plain mean; σ = sqrt(Σσ²)/N.

    d = (10, 12), σ = (1, 2) -> mean = 11; σ = sqrt(1 + 4)/2 = sqrt(5)/2 = 1.118033988749895
    """
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m1 = _store_with_noise(
        tmp_path, "u1", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    m2 = _store_with_noise(
        tmp_path, "u2", cube=_cube(h, np.full((3, 4, 4), 12.0)), sigma=2.0
    )
    result = _stack([m1, m2]).coadd(weights="uniform")
    assert np.allclose(result._data_quantity.to_value(u.K), 11.0)
    assert np.allclose(result._coadd_sigma.to_value(u.K), np.sqrt(5.0) / 2.0)


def test_coadd_noise_weighting_requires_a_companion(tmp_path):
    """weights='noise' with a member lacking a noise companion raises — IVW cannot be faked."""
    from dataclasses import replace

    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    good = _store_with_noise(
        tmp_path, "g", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    # A member whose origin store carries NO noise companion: write the cube, but no companion.
    target = tmp_path / "nocomp"
    target.mkdir()
    store = _cube(h, np.full((3, 4, 4), 12.0)).to_zarr(target)
    reloaded = Cube.from_zarr(str(store))
    reloaded.metadata = replace(reloaded.metadata, origin_store_uri=str(store))
    bare = StackMember(
        cube=reloaded, object="SRC", survey="S", species="CO", transition="2-1"
    )
    with pytest.raises(CoaddError) as exc:
        _stack([good, bare]).coadd(weights="noise")
    assert "noise companion" in str(exc.value)


# --------------------------------------------------------------------------------------
# (C) shift_to_rest — v_sys resolution order + the missing-v_sys raise (user story 17)
# --------------------------------------------------------------------------------------
def test_shift_to_rest_raises_when_v_sys_missing(tmp_path):
    """No catalog v_sys, none on the cube, none supplied -> a descriptive raise (never guessed)."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m = _store_with_noise(
        tmp_path,
        "novsys",
        cube=_cube(h, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        v_sys=None,
    )
    assert m.cube.metadata.systemic_velocity is None
    with pytest.raises(ValueError) as exc:
        _stack([m]).shift_to_rest()  # no v_sys supplied either
    assert "systemic velocity" in str(exc.value)


def test_shift_to_rest_uses_user_value_when_supplied(tmp_path):
    """With no curated v_sys, the user-supplied value is used: the line at v_sys moves to 0 km/s."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m = _store_with_noise(
        tmp_path, "uv", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0, v_sys=None
    )
    before = m.cube.velocity_axis().to_value(u.km / u.s)  # [0, 2, 4]
    shifted = _stack([m]).shift_to_rest(v_sys=4 * u.km / u.s)
    after = shifted[0].cube.velocity_axis().to_value(u.km / u.s)
    # Each channel relabelled by -v_sys: [0,2,4] - 4 -> [-4,-2,0]. v_sys now sits at rest (0).
    assert np.allclose(after, before - 4.0)


def test_shift_to_rest_prefers_store_v_sys_over_user_value(tmp_path):
    """A curated (store) v_sys wins over a user-supplied scalar (it is the per-source authority)."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m = _store_with_noise(
        tmp_path,
        "sv",
        cube=_cube(h, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        v_sys=2 * u.km / u.s,
    )
    before = m.cube.velocity_axis().to_value(u.km / u.s)
    # Supply a DIFFERENT user value; the curated store v_sys (2 km/s) must take precedence.
    shifted = _stack([m]).shift_to_rest(v_sys=99 * u.km / u.s)
    after = shifted[0].cube.velocity_axis().to_value(u.km / u.s)
    assert np.allclose(after, before - 2.0)


def test_resolve_v_sys_precedence_unit_level(tmp_path):
    """The resolver returns store v_sys when present, else the supplied fallback (precedence)."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    with_store = _store_with_noise(
        tmp_path,
        "ws",
        cube=_cube(h, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        v_sys=7 * u.km / u.s,
    )
    without = _store_with_noise(
        tmp_path, "wo", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0, v_sys=None
    )
    assert with_store.resolve_v_sys(supplied=99 * u.km / u.s) == 7 * u.km / u.s
    assert without.resolve_v_sys(supplied=5 * u.km / u.s) == 5 * u.km / u.s
    assert without.resolve_v_sys(supplied=None) is None


def test_resolve_v_sys_prefers_catalog_value(tmp_path):
    """Tier 1: a curated catalog v_sys_kms (on the Record row) wins over store + user values."""

    class _Row:
        v_sys_kms = 42.0  # the catalog's curated per-source systemic velocity

    class _Record:
        row = _Row()

    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    # The cube ALSO carries a store v_sys (8 km/s); the catalog value must still win.
    member = _store_with_noise(
        tmp_path,
        "cat",
        cube=_cube(h, np.full((3, 4, 4), 10.0)),
        sigma=1.0,
        v_sys=8 * u.km / u.s,
    )
    catalogged = StackMember(
        cube=member.cube,
        object="SRC",
        survey="S",
        species="CO",
        transition="2-1",
        record=_Record(),
    )
    assert catalogged.resolve_v_sys(supplied=99 * u.km / u.s) == 42 * u.km / u.s


# --------------------------------------------------------------------------------------
# (D) to_common_beam — the no-super-resolution guard + the co-addable equal-beam path
# --------------------------------------------------------------------------------------
def test_to_common_beam_raises_on_super_resolution(tmp_path):
    """Asking to_common_beam for a FINER beam than a member's raises (no super-resolution)."""
    h, beam = _header(obj="SRC", species="CO", bmaj_arcsec=10.0)
    m = _store_with_noise(
        tmp_path, "big", cube=_cube(h, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    finer = radio_beam.Beam(major=2 * u.arcsec, minor=1 * u.arcsec, pa=0 * u.deg)
    with pytest.raises(LossyDirectionError):
        _stack([m]).to_common_beam(beam=finer)


def test_to_common_beam_produces_equal_beam_members(tmp_path):
    """to_common_beam brings two differently-beamed members to ONE common beam (>= both)."""
    h_a, beam_a = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    h_b, beam_b = _header(obj="SRC", species="CO", bmaj_arcsec=8.0)
    a = _store_with_noise(
        tmp_path, "a", cube=_cube(h_a, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    b = _store_with_noise(
        tmp_path, "b", cube=_cube(h_b, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    common = _stack([a, b]).to_common_beam()
    beams = [member.beam for member in common]
    # Both members end at the same beam, and it is no finer than either input (no super-resolution).
    assert beams[0] == beams[1]
    assert beams[0].major >= beam_a.major
    assert beams[0].major >= beam_b.major


# --------------------------------------------------------------------------------------
# (E) to_velocity_grid — members land on a common grid that then passes the coadd gate
# --------------------------------------------------------------------------------------
def test_to_velocity_grid_aligns_members_for_coadd(tmp_path):
    """Two members on DIFFERENT velocity axes -> regrid onto one grid -> coadd gate passes."""
    h_a, _ = _header(
        obj="SRC", species="CO", bmaj_arcsec=5.0, crval3=0.0, cdelt3=2000.0
    )
    h_b, _ = _header(
        obj="SRC", species="CO", bmaj_arcsec=5.0, crval3=1000.0, cdelt3=2000.0
    )
    a = _store_with_noise(
        tmp_path, "a", cube=_cube(h_a, np.full((3, 4, 4), 10.0)), sigma=1.0
    )
    b = _store_with_noise(
        tmp_path, "b", cube=_cube(h_b, np.full((3, 4, 4), 12.0)), sigma=1.0
    )
    stack = _stack([a, b])
    # Before alignment the spectral grids differ -> coadd refuses.
    with pytest.raises(CoaddError):
        stack.coadd()
    # Resample both onto one explicit velocity grid; now the gate passes and coadd combines.
    grid = np.array([0.0, 2.0, 4.0]) * u.km / u.s
    aligned = stack.to_velocity_grid(grid)
    for member in aligned:
        assert np.allclose(
            member.cube.velocity_axis().to_value(u.km / u.s), [0.0, 2.0, 4.0]
        )
    result = aligned.coadd(weights="uniform")
    # Per channel: a is [10,10,10] on [0,2,4]; b (native [1,3,5]) regrids to [NaN,12,12] on [0,2,4]
    # (v=0 is off b's coverage -> NaN, zero-weighted). So the uniform combine is [10, 11, 11]:
    # v=0 only a constrains (10); v=2 and v=4 both members (mean(10,12)=11).
    per_channel = np.nanmean(result._data_quantity.to_value(u.K), axis=(1, 2))
    assert np.allclose(per_channel, [10.0, 11.0, 11.0])


# --------------------------------------------------------------------------------------
# (F) END-TO-END: a heterogeneous browse stack -> filter -> align -> coadd
# --------------------------------------------------------------------------------------
def test_end_to_end_filter_align_coadd(tmp_path):
    """The full intended path: a mixed stack is gated, then filter+align make it co-addable.

    Two CO members (different beams + offset velocity grids, σ=1 and σ=2) and one HI member. The
    HI member makes the whole stack heterogeneous; filter(species='CO') narrows to the co-addable
    subset, to_common_beam + to_velocity_grid align it, and coadd yields the hand-computed IVW mean.

    The cubes are 16x16 px with modest beams so the spatial convolution to the common beam leaves
    the uniform interior intact at the centre pixel (edge pixels attenuate, as convolution must);
    the IVW arithmetic is asserted at that beam-invariant centre voxel of the shared channel.
    """
    flat = lambda value: np.full((3, 16, 16), value)  # noqa: E731 - terse fixture data
    h_co_a, _ = _header(obj="SRC", species="CO", bmaj_arcsec=2.0, crval3=0.0)
    h_co_b, _ = _header(obj="SRC", species="CO", bmaj_arcsec=3.0, crval3=1000.0)
    h_hi, _ = _header(obj="SRC", species="HI", bmaj_arcsec=2.0)
    co_a = _store_with_noise(
        tmp_path, "co_a", cube=_cube(h_co_a, flat(10.0)), sigma=1.0, species="CO"
    )
    co_b = _store_with_noise(
        tmp_path, "co_b", cube=_cube(h_co_b, flat(12.0)), sigma=2.0, species="CO"
    )
    hi = _store_with_noise(
        tmp_path, "hi", cube=_cube(h_hi, flat(99.0)), sigma=1.0, species="HI"
    )
    stack = _stack([co_a, co_b, hi])

    # (1) The heterogeneous browse stack refuses to coadd (mixed species).
    assert not stack.is_homogeneous
    with pytest.raises(CoaddError):
        stack.coadd()

    # (2) Filter to the co-addable species, then align beam + velocity grid (explicit, auditable).
    co = stack.filter(species="CO")
    grid = np.array([0.0, 2.0, 4.0]) * u.km / u.s
    aligned = co.to_common_beam().to_velocity_grid(grid)

    # (3) coadd: IVW of d=10 (σ=1, w=1) and d=12 (σ=2, w=0.25) -> 10.4 at the beam-invariant centre
    # voxel of the v=2 km/s channel (covered by both members). σ propagation rides the same weights;
    # the common-beam convolution scales σ but the IVW value is the weighted mean, asserted here.
    result = aligned.coadd(weights="noise")
    assert isinstance(result, Cube)
    centre = result._data_quantity.to_value(u.K)[1, 8, 8]  # v=2 channel, centre pixel
    assert np.isclose(centre, 10.4, atol=1e-3)

    # (4) The result carries combined provenance: which members, alignment steps, catalog version.
    prov = result.metadata.provenance
    assert prov["operation"] == "coadd"
    assert prov["weights"] == "noise"
    assert prov["n_members"] == 2
    assert len(prov["member_origins"]) == 2
    assert any("to_common_beam" in step for step in prov["alignment_steps"])
    assert any("to_velocity_grid" in step for step in prov["alignment_steps"])
    assert prov["catalog_version"] == "1.0"


def test_coadd_homogeneous_aligned_stack_succeeds(tmp_path):
    """The happy path: members already homogeneous + on one grid coadd without any alignment call."""
    h, _ = _header(obj="SRC", species="CO", bmaj_arcsec=5.0)
    m1 = _store_with_noise(
        tmp_path, "h1", cube=_cube(h, np.full((3, 4, 4), 6.0)), sigma=1.0
    )
    m2 = _store_with_noise(
        tmp_path, "h2", cube=_cube(h, np.full((3, 4, 4), 6.0)), sigma=1.0
    )
    result = _stack([m1, m2]).coadd(weights="noise")
    # Identical members (d=6, σ=1 each): IVW mean = 6, combined σ = sqrt(1/2).
    assert np.allclose(result._data_quantity.to_value(u.K), 6.0)
    assert np.allclose(result._coadd_sigma.to_value(u.K), np.sqrt(0.5))
