"""Filename projection (ADR-0006) — the human-browsable mirror of the header.

The FITS header is authoritative; the filename is a *derived* projection of it, regenerated
on save, never hand-maintained. The grammar is harvested from legacy astrolyze-1:

    {source}_{telescope}_{species}_{fluxunit}_{resolution}[.ext]

where ``resolution`` is the beam in arcsec, formatted exactly as the original
``resolutionToString``: circular ``"12.00"``, elliptical ``"12.00x10.00"``, and with a
non-zero position angle ``"…a30.0"`` appended. Any field the (possibly incomplete) metadata
lacks degrades to the token ``"unknown"`` so a name is always producible — the projection
mirrors the header, it never blocks on it.
"""

from __future__ import annotations

import re

import astropy.units as u

from .schema import Metadata

UNKNOWN = "unknown"

# Characters that would break the underscore-delimited grammar.
_UNSAFE = re.compile(r"[\s_/]+")


def _token(value) -> str:
    """A filename-safe token for a string-ish field, or ``unknown`` when absent."""
    if value is None:
        return UNKNOWN
    return _UNSAFE.sub("", str(value)) or UNKNOWN


def _flux_token(unit: u.UnitBase | None) -> str:
    """Brightness unit as a single token: ``Jy / beam`` -> ``Jybeam``, ``MJy / sr`` ->
    ``MJysr``, ``K`` -> ``K``."""
    if unit is None:
        return UNKNOWN
    return str(unit).replace(" ", "").replace("/", "") or UNKNOWN


def _resolution_token(beam) -> str:
    """Beam -> legacy resolution string (arcsec). Mirrors astrolyze-1 ``resolutionToString``."""
    if beam is None:
        return UNKNOWN
    major = beam.major.to_value(u.arcsec)
    minor = beam.minor.to_value(u.arcsec)
    pa = beam.pa.to_value(u.deg)

    circular = f"{major:.2f}" == f"{minor:.2f}"
    shape = f"{major:.2f}" if circular else f"{major:.2f}x{minor:.2f}"
    # A position angle is only meaningful (and only appended) when it is non-zero.
    return shape if f"{pa:.1f}" == "0.0" else f"{shape}a{pa:.1f}"


def project(metadata: Metadata, *, extension: str = "fits") -> str:
    """Return the header-derived filename for *metadata*.

    Parameters
    ----------
    metadata : Metadata
        The parsed schema (the authoritative source the name mirrors).
    extension : str, optional
        File extension without the leading dot (default ``"fits"``).
    """
    stem = "_".join(
        (
            _token(metadata.object),
            _token(metadata.telescope),
            _token(metadata.species),
            _flux_token(metadata.bunit),
            _resolution_token(metadata.beam),
        )
    )
    return f"{stem}.{extension}"


__all__ = ["project"]
