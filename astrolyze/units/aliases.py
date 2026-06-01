"""Named radio/sub-mm unit aliases (ADR-0003).

Plain ``astropy.units`` objects for the units this field writes by hand a hundred times a
day. They are *not* registered into astropy's global unit registry — importing astrolyze
must have no global side effects (ADR-0005's no-rcParams rule, applied to units), so
``u.Unit("Tmb")`` deliberately still fails. Use these objects directly instead.

``Tmb``/``Ta`` are the Kelvin (Rayleigh-Jeans) brightness-temperature scale; they are
labelled for intent and are identical to ``u.K``.
"""

from __future__ import annotations

import astropy.units as u

#: Main-beam brightness temperature (Rayleigh-Jeans scale, = K).
Tmb = u.K
#: Antenna temperature (Rayleigh-Jeans scale, = K).
Ta = u.K
#: Flux density per beam solid angle.
Jy_beam = u.Jy / u.beam
#: Surface brightness in jansky per steradian.
Jy_sr = u.Jy / u.sr
#: Surface brightness in megajansky per steradian (common for IR/sub-mm maps).
MJy_sr = u.MJy / u.sr
#: Velocity-integrated brightness temperature (integrated intensity).
K_kms = u.K * u.km / u.s

__all__ = ["Tmb", "Ta", "Jy_beam", "Jy_sr", "MJy_sr", "K_kms"]
