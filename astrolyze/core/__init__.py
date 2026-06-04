"""Core data types: Cube (PPV), Map (2D/moment), Spectrum (1D) — ADR-0004.

Thin wrappers that **compose** (not subclass) spectral-cube / astropy / specutils: each holds
a private upstream handle and the astrolyze :class:`~astrolyze.io.Metadata` that carries the
physical context (beam + rest frequency + velocity convention + provenance). They add the
toolkit's value over the ecosystem:

- **context-carrying type transitions** — ``Cube.moment0() -> Map``, ``cube[:, y, x] ->
  Spectrum`` — so context flows for free;
- the **unit hub** ``.to(unit)``, which routes through :func:`astrolyze.units.convert` with
  the object's own context (ADR-0003c) and never silently guesses the physics;
- a **``.plot()`` seam** onto the viz engine (ADR-0005, issue #6).

Build them from an :class:`~astrolyze.io.LoadedData` via ``Cube.from_loaded(load(path))`` /
``Map.from_loaded(...)``.
"""

from __future__ import annotations

from .cube import Cube, LossyDirectionError
from .map import Map
from .spectrum import Spectrum

__all__ = ["Cube", "LossyDirectionError", "Map", "Spectrum"]
