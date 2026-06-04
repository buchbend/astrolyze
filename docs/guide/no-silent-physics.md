# No silent physics

*Embodies [ADR-0003](../adr/0003-unit-layer-pure-astropy-explicit-conventions.md) (explicit
conventions, no silent physics) and
[ADR-0013](../adr/0013-traceable-correctable-results.md) (traceable, correctable results).*

## What it is

astrolyze refuses to guess the physics that silently biases radio/sub-mm work — the velocity
convention, the rest frequency, the beam, and the Rayleigh-Jeans-vs-Planck brightness scale.
A conversion that lacks the context it needs **raises** rather than defaulting. This guard is
exposed three ways:

- `MissingContextError` — raised when an operation needs context that is absent. It carries an
  optional structured `Insufficiency` descriptor (the field `insufficiency`), so a caller can
  read the gaps without parsing the message string.
- `Insufficiency` / `ContextGap` — the returnable "what I'd need" value objects. An
  `Insufficiency` is an ordered list of `ContextGap`s; each gap has a `name`, a human `reason`,
  and (usually) a `how_to_supply`. An `Insufficiency` is **falsy when empty**, so
  `if insufficiency:` reads as "there are gaps".
- `Cube.can_convert_to(unit)` — a **non-raising probe**. It returns the `Insufficiency`
  describing what is missing to convert this object to `unit`, *without* raising. An empty
  (satisfied) descriptor means the conversion can proceed.

## Why it works this way

- **A silently assumed convention or scale is a hard-to-detect bias (ADR-0003).** The radio,
  optical, and relativistic conventions give different velocities for the same line; the
  Rayleigh-Jeans and Planck scales give different temperatures. Defaulting one of these would
  produce a plausible-looking but wrong number. astrolyze names the choice instead.
- **Refusal as a first-class, inspectable result (ADR-0013).** Making the refusal a returnable
  descriptor — not only an exception — lets a caller ask *"what would I need to do X?"* and
  branch on the answer, without a `try/except`. Results stay traceable and correctable: the
  descriptor renders one legible line per gap.
- **Gaps are decoupled from the schema shape.** A `ContextGap` can name a parameter that has no
  `Metadata` field yet (e.g. the genuinely-ambiguous `calibration_scale` for RJ-vs-Planck),
  because naming a gap is the descriptor's job, not the schema's. The gap list is decided from
  the requested conversion, not from a fixed schema.

## Usage

Probe before acting (no exception to catch):

```python
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))

gaps = cube.can_convert_to("K")     # an Insufficiency
if gaps:                            # falsy when empty -> "there are gaps"
    for gap in gaps:
        print(gap.name, "—", gap.reason, "|", gap.how_to_supply)
else:
    result = cube.to("K")           # safe: the conversion can proceed
```

Or act directly and read the structured detail off the raised error:

```python
from astrolyze.units import MissingContextError

try:
    cube.to("K")
except MissingContextError as exc:
    if exc.insufficiency is not None:
        print(exc.insufficiency.names)     # e.g. ["rest_frequency", "calibration_scale"]
        print(exc.insufficiency.message()) # one legible line per gap
```

The same structured form is available straight off `Metadata` for the mandatory-context gaps:

```python
loaded = load("incomplete_cube.fits")
loaded.metadata.insufficiency()    # an Insufficiency for the missing mandatory context
loaded.metadata.ensure_complete()  # raises MissingContextError carrying that Insufficiency
```

## See also

- [I/O and metadata](io-and-metadata.md) — completeness flags on `Metadata`.
- [Coordinates and validity](coordinates-and-validity.md) — `coordinates.frequency` raises the
  same error when the convention is absent.
