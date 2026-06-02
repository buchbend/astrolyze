# Object model: a thin Cube / Map / Spectrum trio that wraps the ecosystem

The program is about PPV cubes *and* maps *and* spectra, but astrolyze2's `FitsMap` is
2D-only. We need core types that give us the missing PPV `Cube` without re-deriving the
heavy cube machinery (moments, spectral axis, convolution, reproject) that `spectral-cube`
and `astropy` already provide — re-implementing that was the identified "overkill risk".

**Decision:** a small typed trio — **`Cube` (PPV)**, **`Map` (2D image/moment)**,
**`Spectrum` (1D)** — as **thin wrappers that delegate** to `spectral-cube`/`astropy`/
`specutils` for the heavy lifting and add the toolkit's value: object-carried context
(beam + rest frequency + velocity convention, per ADR-0003), the unit-hub `.to()`, a
unified `.plot()`, and provenance. Type transitions flow naturally:
`Cube.moment0() -> Map`, `Cube[:, y, x] -> Spectrum`, carrying context for free.

Rejected: (a) one polymorphic Map branching on dimensionality (2D/3D behaviour diverges
too much — moments/PV/spectral-axis are 3D-only); (c) functions-only over native objects
(maximally DRY but loses the "object knows its context and plots itself" ergonomics that
are the point).

## Riders (decided by recommendation — flag to revisit if wrong)

1. **Compose/wrap, do NOT subclass `SpectralCube`.** Hold a private handle (e.g. `._sc`)
   and delegate. Same reasoning as not subclassing `Quantity` (ADR-0003): less coupling to
   upstream internals, fewer sharp edges.
2. **`Spectrum` is a thin adaptor to `specutils`; ship `Cube` + `Map` first.** Don't
   reinvent spectral fitting/representation — lean on `specutils.Spectrum1D` under the
   `Spectrum` wrapper. `Cube` + `Map` are the v1 priority.

## Consequences

- Cost: wrapper boilerplate, and the discipline to **stay thin** — the temptation to
  reimplement upstream behaviour must be resisted (delegate, don't reinvent).
- The object is the home for context + provenance, reinforcing ADR-0003 and the I/O
  decisions to follow.
- Delegation pattern keeps us resilient to upstream API change (one wrapper to adjust, not
  scattered call sites).
