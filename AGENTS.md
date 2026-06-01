# Working with astro data in astrolyze (guide for humans and agents)

This is the house guide for using astrolyze. A coding agent uses astrolyze the **same way a
human does** — by writing astrolyze Python/CLI — producing ordinary, reviewable, reproducible
scripts. There is no AI-only interface (no MCP); the CLI and Python API *are* the interface.

## Principles

1. **Everything goes through astrolyze.** Don't hand-roll matplotlib, unit math, or FITS
   parsing in an analysis. If astrolyze can't do something yet, **extend astrolyze first**
   (with tests, by its conventions), then use it.
2. **Be explicit about physics.** Always supply the velocity convention and rest frequency
   where needed; astrolyze will raise if they're missing — that's intentional, not an error to
   work around.
3. **Use the house display.** Plot via `.plot()` / `astrolyze.viz`; don't restyle globally.
4. **Results are traceable.** Surface how a result was produced (operations + reasoning) so a
   reader can follow and correct it; cite the real artifact behind any claim.
5. **Stay thin.** Prefer delegating to `spectral-cube`/`astropy` over reimplementing.

## The core objects

`Cube` (PPV), `Map` (2D / moment), and `Spectrum` (1D) are thin wrappers that **compose**
spectral-cube / astropy / specutils and carry the physical context (beam + rest frequency +
velocity convention) parsed from the FITS header. Build one from a file through the `io` seam:

```python
from astrolyze.io import load
from astrolyze.core import Cube

cube = Cube.from_loaded(load("ngc0628_co21.fits"))   # <Cube NGC0628 (98, 1600, 1600) [K]>
```

Type transitions carry the context for you:

| operation                | result     |
| ------------------------ | ---------- |
| `cube.moment0()`         | `Map`      |
| `cube.moment(order, axis)` | `Map`    |
| `cube[k]` (channel)      | `Map`      |
| `cube[:, y, x]`          | `Spectrum` |
| `cube[:, y0:y1, x0:x1]`  | `Cube`     |

## The tracer recipe (load → moment0 → convert → plot)

```python
from astrolyze.core import Cube
from astrolyze.io import load

cube = Cube.from_loaded(load("ngc0628_co21.fits"))
integrated = cube.moment0().to("K km/s")   # beam / rest freq / convention come from the object
fig, ax = integrated.plot()                # cividis + WCS axes + beam ellipse + unit colorbar
fig.savefig("ngc0628_co21_moment0.png")
```

`.to(unit)` supplies the object's beam, rest frequency, and velocity convention to the unit
layer so you don't repeat them. It **never** defaults the genuinely-ambiguous
Rayleigh-Jeans-vs-Planck scale — a brightness-temperature conversion needs it stated:

```python
channel = cube[48]                                   # a Jy/beam or K channel Map
channel.to("K", temperature_scale="rayleigh_jeans")  # or "planck"; omitting it raises
```

If the header is missing a rest frequency or velocity convention, the file still **loads**
(so you can inspect archival data), but any operation that needs the missing context raises
`astrolyze.units.MissingContextError` rather than guessing.

## The CLI (same path from the shell)

The CLI exposes the identical spine; an agent reaches for it exactly as a human would.

```bash
astrolyze info  ngc0628_co21.fits                 # metadata schema + completeness (read-only)
astrolyze moment0 ngc0628_co21.fits -u "K km/s" -o ngc0628_mom0.png
astrolyze moment0 ngc0628_co21.fits --temperature-scale planck   # for a brightness conversion
astrolyze --help
```

`info` prints a table of the parsed schema (object, telescope, species, rest frequency,
velocity convention, beam, unit, distance, calibration error) and whether the file is
*complete*. `moment0` runs load → moment0 → `.to(unit)` → plot and writes a house-style PNG
(default `<input>_moment0.png`). A conversion that needs absent context exits non-zero with a
clear message — the CLI surfaces the library's refusal, it doesn't paper over it.

## Try it on real data

A 128×128×50 cutout of the PHANGS-ALMA NGC 628 CO(2-1) cube ships in
`tests/data/ngc0628_co21_cutout.fits.gz` (see `tests/data/PROVENANCE.md`):

```bash
astrolyze info tests/data/ngc0628_co21_cutout.fits.gz
python examples/tracer_ngc628.py tests/data/ngc0628_co21_cutout.fits.gz /tmp/mom0.png
```

Point `$ASTROLYZE_TRACER_CUBE` at a full cube to run the example on real survey data without
an argument.
