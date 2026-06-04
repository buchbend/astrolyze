"""Tests for the typed-insufficiency return — the structured "what I'd need" descriptor
(issue #29; ADR-0013 traceable/correctable results, ADR-0003 raise-not-guess).

Written first (red/green TDD). These tests *are* the correctness obligation for making the
existing no-silent-physics guard a first-class, inspectable, returnable value object:

- a small **insufficiency descriptor** type lists each missing-or-ambiguous interpretation
  parameter as a named gap (``name`` + human ``reason`` [+ how to supply it]);
- ``MissingContextError`` *carries* that descriptor, so a raise exposes structured detail and
  not merely a message string;
- a **non-raising probe** (``cube.can_convert_to("K")``) *returns* the descriptor naming the
  gaps, so a caller can branch on what's missing without catching an exception;
- the descriptor can name parameters that do not yet have a schema field (``spectral_frame``,
  ``calibration_scale``, ``eta_mb`` arrive in other issues #24/#25): a gap is just a
  name+reason value object, so the probe can enumerate them from the requested operation
  without requiring the field to exist.

The reference dataset mirrors the tracer-bullet PHANGS cube: NGC 0628, CO(2-1) at
230.538 GHz, a 12"x10" (PA 30 deg) beam, on the radio velocity convention.
"""

import warnings

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.io import Metadata, load
from astrolyze.units import (
    ContextGap,
    Insufficiency,
    MissingContextError,
)

# spectral-cube emits cosmetic warnings (no-beam, stokes) we do not care about here.
warnings.filterwarnings("ignore", module="spectral_cube")

REST = 230.538 * u.GHz  # CO(2-1)
BEAM = radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg)

# The full set of interpretation parameters the descriptor must be able to name — the
# seven listed in the issue. Three of them (spectral_frame / calibration_scale / eta_mb) do
# NOT yet have a Metadata field (they arrive in #24/#25); the descriptor must name them
# anyway, proving a gap is independent of schema shape.
NAMEABLE_GAPS = (
    "rest_frequency",
    "velocity_convention",
    "beam",
    "spectral_frame",
    "calibration_scale",
    "eta_mb",
    "unit",
)


# --------------------------------------------------------------------------------------
# Synthetic cube fixtures (header carries the schema; we knock out fields to make gaps)
# --------------------------------------------------------------------------------------
def _cube_header(*, bunit="Jy/beam"):
    """A 3D cube on the radio velocity convention with the full astrolyze schema."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"] = "NGC0628"
    h["TELESCOP"] = "ALMA"
    h["BUNIT"] = bunit
    h["RESTFRQ"] = (REST.to_value(u.Hz), "Hz")
    for k, v in BEAM.to_header_keywords().items():
        h[k] = v
    h["HIERARCH ASTROLYZE VCONV"] = "radio"
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    return h


@pytest.fixture
def complete_cube(tmp_path):
    """A complete Jy/beam Cube built through the real io.load seam."""
    path = tmp_path / "ngc0628_co21.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), _cube_header())
    return Cube.from_loaded(load(path))


@pytest.fixture
def bare_cube(tmp_path):
    """A Jy/beam Cube with NO rest frequency, convention, or beam (real-archival case).

    The spectral axis is a plain frequency axis (CTYPE3='FREQ', CUNIT3='Hz') so the loader
    reads *no* velocity convention from the WCS either — the convention is genuinely absent,
    not merely missing the dedicated keyword."""
    h = _cube_header()
    del h["RESTFRQ"]
    del h["HIERARCH ASTROLYZE VCONV"]
    h["CTYPE3"], h["CUNIT3"] = "FREQ", "Hz"
    h["CRVAL3"], h["CDELT3"] = REST.to_value(u.Hz), 1.0e6
    for k in BEAM.to_header_keywords():
        del h[k]
    path = tmp_path / "archival.fits"
    fits.writeto(path, np.ones((4, 3, 3), dtype="float32"), h)
    return Cube.from_loaded(load(path))


# --------------------------------------------------------------------------------------
# AC: a first-class insufficiency descriptor type, listing named gaps with name + reason
# --------------------------------------------------------------------------------------
def test_context_gap_carries_name_reason_and_optional_how_to_supply():
    gap = ContextGap(
        name="rest_frequency",
        reason="needed to set the frequency of the Planck/RJ law",
        how_to_supply="set RESTFRQ on the header, or pass rest_frequency=",
    )
    assert gap.name == "rest_frequency"
    assert "Planck" in gap.reason
    assert "rest_frequency=" in gap.how_to_supply
    # how_to_supply is optional (some gaps cannot be supplied generically).
    bare = ContextGap(name="unit", reason="unrecognised unit")
    assert bare.how_to_supply is None


def test_insufficiency_lists_gaps_with_name_and_reason():
    ins = Insufficiency(
        [
            ContextGap("rest_frequency", "needed for the brightness law"),
            ContextGap("calibration_scale", "T_mb vs T_A* changes the value"),
        ]
    )
    assert isinstance(ins, Insufficiency)
    assert ins.names == ["rest_frequency", "calibration_scale"]
    for gap in ins:  # iterable over its gaps
        assert isinstance(gap, ContextGap)
        assert gap.name and gap.reason
    assert ins.get("rest_frequency").reason  # addressable by name
    assert ins.get("not_a_gap") is None


def test_empty_insufficiency_is_falsy_and_a_populated_one_is_truthy():
    # The probe contract: "no gaps" must be a falsy, branchable result.
    assert not Insufficiency([])
    assert Insufficiency([]).is_satisfied is True
    populated = Insufficiency(
        [ContextGap("beam", "Jy/beam needs the beam solid angle")]
    )
    assert populated
    assert populated.is_satisfied is False
    assert "beam" in str(populated)  # renders a legible message (ADR-0013)


@pytest.mark.parametrize("name", NAMEABLE_GAPS)
def test_descriptor_can_name_every_interpretation_parameter(name):
    # Every parameter the issue lists is a buildable gap with a non-empty human reason —
    # including the three that have NO Metadata field yet (#24/#25). A gap is name+reason,
    # decoupled from schema shape, so it can name them today.
    gap = ContextGap.for_parameter(name)
    assert gap.name == name
    assert isinstance(gap.reason, str) and gap.reason.strip()


def test_unrecognised_unit_is_a_first_class_gap():
    # An unrecognised unit is a no-silent-physics case too: name the unit and why, quoting
    # the offending text so the caller sees exactly what was not understood.
    gap = ContextGap.for_parameter("unit", detail="frobnicate/beam")
    assert gap.name == "unit"
    assert "frobnicate/beam" in gap.reason


# --------------------------------------------------------------------------------------
# AC: MissingContextError carries the descriptor (structured detail, not only a string)
# --------------------------------------------------------------------------------------
def test_missing_context_error_carries_the_descriptor():
    ins = Insufficiency([ContextGap("rest_frequency", "needed for the brightness law")])
    err = MissingContextError("missing context", insufficiency=ins)
    assert err.insufficiency is ins
    assert err.insufficiency.names == ["rest_frequency"]
    # Still a plain ValueError-with-message for code that only reads the string.
    assert "missing context" in str(err)


def test_missing_context_error_without_descriptor_still_works():
    # Backward compatible: existing raises that pass only a string keep working.
    err = MissingContextError("just a message")
    assert err.insufficiency is None
    assert str(err) == "just a message"


# --------------------------------------------------------------------------------------
# AC: existing Metadata.missing / ensure_complete preserved, and can produce the descriptor
# --------------------------------------------------------------------------------------
def test_metadata_missing_and_ensure_complete_unchanged():
    m = Metadata(bunit=u.K)  # no rest frequency, no convention
    assert set(m.missing) == {"rest_frequency", "velocity_convention"}
    assert m.is_complete is False
    with pytest.raises(MissingContextError) as exc:
        m.ensure_complete()
    msg = str(exc.value)
    assert "rest_frequency" in msg and "velocity_convention" in msg


def test_metadata_can_produce_the_descriptor():
    m = Metadata(bunit=u.K)
    ins = m.insufficiency()
    assert isinstance(ins, Insufficiency)
    assert set(ins.names) == {"rest_frequency", "velocity_convention"}
    for gap in ins:
        assert gap.reason  # each gap explains itself

    # A complete metadata yields an empty (satisfied) descriptor.
    complete = Metadata(
        rest_frequency=REST,
        velocity_convention="radio",
        bunit=u.K,
    )
    assert complete.insufficiency().is_satisfied is True


def test_ensure_complete_attaches_the_descriptor_to_the_raise():
    # The raise now carries the structured detail too — same gaps as Metadata.missing.
    m = Metadata(bunit=u.K)
    with pytest.raises(MissingContextError) as exc:
        m.ensure_complete()
    ins = exc.value.insufficiency
    assert isinstance(ins, Insufficiency)
    assert set(ins.names) == set(m.missing)


# --------------------------------------------------------------------------------------
# AC: a non-raising probe returns the descriptor naming the gaps instead of raising
# --------------------------------------------------------------------------------------
def test_can_convert_to_returns_empty_descriptor_when_possible(complete_cube):
    # A complete Jy/beam cube with beam + rest frequency: Jy/beam -> MJy/sr is pure
    # geometry and possible, so the probe returns a satisfied (empty) descriptor — and
    # does NOT raise.
    ins = complete_cube.can_convert_to("MJy / sr")
    assert isinstance(ins, Insufficiency)
    assert ins.is_satisfied is True


def test_can_convert_to_K_names_the_gaps_without_raising(bare_cube):
    # Jy/beam -> K needs rest_frequency (brightness law), the beam (Jy/beam geometry), and
    # the genuinely-ambiguous RJ-vs-Planck calibration_scale. The probe RETURNS them; it
    # must NOT raise (that is the whole point of a probe vs .to()).
    ins = bare_cube.can_convert_to("K")
    assert isinstance(ins, Insufficiency)
    assert ins  # truthy: there are gaps
    assert "rest_frequency" in ins.names
    assert "calibration_scale" in ins.names
    assert "beam" in ins.names
    for gap in ins:
        assert gap.reason  # each names *why* it is needed


def test_can_convert_to_K_on_complete_cube_still_flags_the_ambiguous_scale(
    complete_cube,
):
    # Even a fully-complete cube cannot silently pick Rayleigh-Jeans vs Planck: a Jy/beam
    # -> K probe names calibration_scale as the one remaining gap (rest freq + beam are
    # present), so the caller learns the single thing it must decide.
    ins = complete_cube.can_convert_to("K")
    assert "calibration_scale" in ins.names
    assert "rest_frequency" not in ins.names  # the cube carries it
    assert "beam" not in ins.names  # the cube carries it


def test_can_convert_to_flags_unrecognised_target_unit(complete_cube):
    # An unparseable target unit is a no-silent-physics case: name the unit gap, do not
    # raise an opaque astropy error out of a probe.
    ins = complete_cube.can_convert_to("not-a-real-unit")
    assert "unit" in ins.names


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_units_library_code():
    import re
    from pathlib import Path

    import astrolyze.units as units_pkg

    pkg_dir = Path(units_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
