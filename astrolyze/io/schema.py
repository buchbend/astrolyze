"""The astrolyze metadata schema (ADR-0006).

The FITS header is the single source of truth for a dataset's physical context: rest
frequency, velocity convention, beam, distance, species, calibration error, object, and a
free-form name-tag. :class:`Metadata` is the parsed, typed view of that schema.

Two design points worth stating:

- **Lazy enforcement (ADR-0006 ii).** A real archival header missing mandatory context
  (rest frequency / velocity convention) still parses; the object is merely flagged
  *incomplete* (:attr:`Metadata.missing`) and :meth:`Metadata.ensure_complete` raises when an
  operation needs what is absent. We never silently guess (ADR-0003) and never refuse a file.
- **Keyword namespacing.** Standard FITS keywords are used where they exist (``OBJECT``,
  ``TELESCOP``, ``BUNIT``, ``RESTFRQ``, ``BMAJ``/``BMIN``/``BPA``). Everything astrolyze-
  specific lives under ``HIERARCH ASTROLYZE …`` cards — clear, collision-free, round-trips
  cleanly. The velocity convention is also *read* (not guessed) from WCS ``CTYPE3`` when our
  keyword is absent, so real PHANGS/SEDIGISM cubes that state it via ``VRAD``/``VOPT``/``VELO``
  parse as complete.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields

import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.units import (
    MissingContextError,
    VelocityConvention,
    coerce_velocity_convention,
)

# -- header keyword names --------------------------------------------------------------
# Standard FITS keywords (read where the convention already exists).
KEY_OBJECT = "OBJECT"
KEY_TELESCOPE = "TELESCOP"
KEY_BUNIT = "BUNIT"
KEY_RESTFRQ = "RESTFRQ"  # preferred FITS spelling
KEY_RESTFRQ_ALT = "RESTFREQ"  # legacy alias accepted on read

# astrolyze-namespaced HIERARCH cards (accessed by the text after HIERARCH, e.g.
# header["ASTROLYZE VCONV"]).
KEY_VCONV = "ASTROLYZE VCONV"
KEY_SPECIES = "ASTROLYZE SPECIES"
KEY_DISTANCE = "ASTROLYZE DISTANCE"
KEY_DISTUNIT = "ASTROLYZE DISTUNIT"
KEY_CALERR = "ASTROLYZE CALERR"
KEY_NAMETAG = "ASTROLYZE NAMETAG"
KEY_SCHEMA = "ASTROLYZE SCHEMA"

SCHEMA_VERSION = 1
DEFAULT_DISTANCE_UNIT = u.Mpc

# WCS CTYPE codes for a velocity spectral axis -> Doppler convention. Reading what the
# header states is not guessing; this is how real archival cubes declare their convention.
_CTYPE_VELOCITY_CONVENTION = {
    "VRAD": VelocityConvention.RADIO,
    "VOPT": VelocityConvention.OPTICAL,
    "VELO": VelocityConvention.RELATIVISTIC,
}

# Absence of these makes a dataset "incomplete" (ADR-0006 ii): they are the silent-error
# traps (ADR-0003) that gate unit conversion and all velocity work.
REQUIRED_CONTEXT = ("rest_frequency", "velocity_convention")


@dataclass
class Metadata:
    """The typed astrolyze header schema. All fields default to ``None`` so a partial header
    still produces a usable object (lazy enforcement)."""

    object: str | None = None
    telescope: str | None = None
    species: str | None = None
    rest_frequency: u.Quantity | None = None  # frequency
    velocity_convention: VelocityConvention | None = None
    beam: radio_beam.Beam | None = None
    bunit: u.UnitBase | None = None
    distance: u.Quantity | None = None  # length
    calibration_error: float | None = None  # fractional
    name_tag: str | None = None
    #: Provenance seam (ADR-0006 (c)): unused in v1, reserved so the heterogeneity-tracking
    #: layer can attach without breaking the header contract.
    provenance: object | None = field(default=None, repr=False)

    # -- completeness (lazy enforcement) -----------------------------------------------
    @property
    def missing(self) -> list[str]:
        """Mandatory-context fields that are absent (empty when complete)."""
        return [name for name in REQUIRED_CONTEXT if getattr(self, name) is None]

    @property
    def is_complete(self) -> bool:
        return not self.missing

    def ensure_complete(self) -> None:
        """Raise :class:`MissingContextError` naming the missing mandatory context.

        Call this at the door of any operation that needs the physical context (unit
        conversion, velocity work). Never defaults the answer — see ADR-0003/0006."""
        if self.missing:
            raise MissingContextError(
                "dataset is missing mandatory context: "
                + ", ".join(self.missing)
                + " — set it on the header/metadata before this operation (astrolyze never "
                "guesses rest frequency or velocity convention)"
            )

    # -- header (de)serialisation ------------------------------------------------------
    @classmethod
    def from_header(cls, header: fits.Header) -> "Metadata":
        """Parse a :class:`~astropy.io.fits.Header` into the schema. The header is
        authoritative; missing fields simply stay ``None``."""
        rest = header.get(KEY_RESTFRQ, header.get(KEY_RESTFRQ_ALT))
        beam = radio_beam.Beam.from_fits_header(header) if "BMAJ" in header else None

        return cls(
            object=header.get(KEY_OBJECT),
            telescope=header.get(KEY_TELESCOPE),
            species=header.get(KEY_SPECIES),
            rest_frequency=(float(rest) * u.Hz) if rest is not None else None,
            velocity_convention=_read_convention(header),
            beam=beam,
            bunit=_read_unit(header.get(KEY_BUNIT)),
            distance=_read_distance(header),
            calibration_error=_opt_float(header.get(KEY_CALERR)),
            name_tag=header.get(KEY_NAMETAG),
        )

    def to_header(self, base: fits.Header | None = None) -> fits.Header:
        """Write the schema onto a copy of *base* (preserving its WCS), or onto a fresh
        header. Setting keywords replaces rather than appends, so this is idempotent."""
        header = base.copy() if base is not None else fits.Header()
        _put(header, KEY_SCHEMA, (SCHEMA_VERSION, "astrolyze metadata schema version"))

        if self.object is not None:
            _put(header, KEY_OBJECT, self.object)
        if self.telescope is not None:
            _put(header, KEY_TELESCOPE, self.telescope)
        if self.species is not None:
            _put(header, KEY_SPECIES, self.species)
        if self.rest_frequency is not None:
            _put(header, KEY_RESTFRQ, (self.rest_frequency.to_value(u.Hz), "Hz"))
        if self.velocity_convention is not None:
            _put(header, KEY_VCONV, VelocityConvention(self.velocity_convention).value)
        if self.beam is not None:
            for key, value in self.beam.to_header_keywords().items():
                _put(header, key, value)
        if self.bunit is not None:
            _put(header, KEY_BUNIT, str(self.bunit))
        if self.distance is not None:
            _put(header, KEY_DISTANCE, self.distance.value)
            _put(header, KEY_DISTUNIT, str(self.distance.unit))
        if self.calibration_error is not None:
            _put(header, KEY_CALERR, float(self.calibration_error))
        if self.name_tag is not None:
            _put(header, KEY_NAMETAG, self.name_tag)
        return header


# -- header writing --------------------------------------------------------------------
def _put(header: fits.Header, key: str, value) -> None:
    """Set a card, using an explicit ``HIERARCH`` prefix for the long astrolyze keywords so
    astropy writes a HIERARCH card without emitting a VerifyWarning (reading still uses the
    short ``ASTROLYZE …`` form, which astropy resolves to the same card)."""
    header["HIERARCH " + key if key.startswith("ASTROLYZE ") else key] = value


# -- header field parsers --------------------------------------------------------------
def _read_convention(header: fits.Header) -> VelocityConvention | None:
    """The dedicated astrolyze keyword is authoritative; fall back to the WCS spectral
    CTYPE (what real archival cubes use). Returns ``None`` if neither states it."""
    explicit = header.get(KEY_VCONV)
    if explicit is not None:
        return coerce_velocity_convention(explicit)
    for axis in range(1, int(header.get("NAXIS", 0)) + 1):
        code = str(header.get(f"CTYPE{axis}", "")).strip().upper()
        if code in _CTYPE_VELOCITY_CONVENTION:
            return _CTYPE_VELOCITY_CONVENTION[code]
    return None


def _read_unit(value) -> u.UnitBase | None:
    return u.Unit(value) if value else None


def _read_distance(header: fits.Header) -> u.Quantity | None:
    value = header.get(KEY_DISTANCE)
    if value is None:
        return None
    unit = u.Unit(header.get(KEY_DISTUNIT, str(DEFAULT_DISTANCE_UNIT)))
    return float(value) * unit


def _opt_float(value) -> float | None:
    return float(value) if value is not None else None


# Field names in declaration order — handy for callers that diff two Metadata objects.
SCHEMA_FIELDS = tuple(f.name for f in fields(Metadata) if f.name != "provenance")

__all__ = ["Metadata", "REQUIRED_CONTEXT", "SCHEMA_VERSION", "SCHEMA_FIELDS"]
