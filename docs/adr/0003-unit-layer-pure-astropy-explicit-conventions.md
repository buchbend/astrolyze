# Unit layer: pure-astropy substrate + object-carried context, with explicit convention & rest frequency

Radio/sub-mm unit handling is a top source of silent error (Rayleigh–Jeans vs Planck
brightness temperature, beam-area conversions, and velocity/Doppler convention with an
implicit rest frequency). The toolkit must make these conversions correct, reusable, and
hard to get wrong.

**Decision — two layers, no Quantity subclass:**

- **(a) Pure-astropy units substrate.** A standalone `astrolyze` units module usable with no
  data object: named radio aliases (`Tmb`, `Jy/beam`, `K km/s`, `MJy/sr`, …), a curated
  bundle of equivalencies (brightness temperature, beam angular area, spectral/Doppler),
  and helper converters. Everything is `astropy.units.Quantity`; nothing custom to maintain
  but the domain knowledge.
- **(c) Object-carried context (ergonomic layer).** The data object (`FitsMap`/`Cube`) holds
  beam + rest frequency + velocity convention and exposes `.to(unit)` that supplies the
  right equivalencies automatically. DRY: one conversion implementation (a), two entry
  points (standalone or via the object).
- **NOT (b):** no `Quantity` subclass / context-carrying wrapper — astropy subclassing is
  sharp-edged and the magic is not worth the maintenance.

**Hard rule — no silent physics.** Velocity/Doppler **convention** (radio | optical |
relativistic) and **rest frequency** are **explicit and mandatory** for any conversion that
depends on them; they are never silently defaulted. A missing convention/rest-frequency
raises, it does not guess. This is stricter than astropy's defaults on purpose: a silently
assumed convention or rest frequency is a classic, hard-to-detect bias in exactly this
science.

## Consequences

- Correctness is a **test-suite obligation**, not a code-cleverness one: round-trip tests
  across the unit zoo + at least one hand-checked Planck (non-RJ) value + convention/
  rest-frequency-required tests.
- Slightly more verbose call sites (you must state convention/rest freq) — accepted as the
  price of no silent bias.
- The data object becomes the natural home for provenance (beam, freq, convention),
  reinforcing the I/O-provenance decisions to follow.
