# astrolyze

**A toolkit for radio/sub-mm position–position–velocity (PPV) cubes and spectra.**

astrolyze makes data handling, unit conversion, display, and I/O for spectral-line
astronomy *consistent and correct* — so the same cube is always plotted the same way, unit
conversions never silently assume a velocity convention, and every result is traceable. It is
a thin, opinionated layer over the modern Astropy ecosystem (`astropy`, `spectral-cube`,
`radio_beam`, `reproject`, `specutils`), usable by a human or driven identically by a coding
agent.

> **Status: 0.1 pre-alpha — rebuild in progress.** This is the rebuilt astrolyze. The
> original 2012–2016 package lives on at
> [buchbend/astrolyze-legacy](https://github.com/buchbend/astrolyze-legacy) for heritage.

## Design

- **Units are explicit.** Velocity convention (radio/optical/relativistic) and rest frequency
  are required where a conversion needs them — astrolyze raises rather than guessing.
- **Core types carry context.** `Cube` (PPV), `Map` (2D/moment), `Spectrum` (1D) hold their
  beam, rest frequency, and convention; `.to()` and `.plot()` use that context.
- **One display style.** WCS axes, a drawn beam, a colorbar labelled in the data's units;
  `cividis` by default. Plotting never mutates your global matplotlib settings.
- **Self-describing data.** The FITS header is authoritative; astrolyze-written files get a
  metadata-derived, browsable filename.

## Install

```bash
pip install -e ".[dev]"   # development install
pytest                    # run the test suite
```

## License

BSD-3-Clause. See [LICENSE](LICENSE).
