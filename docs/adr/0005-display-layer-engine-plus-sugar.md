# Display layer: free-function engine + thin object-method sugar, style applied locally

"Never re-implement how spectra are plotted or maps are displayed" is a core toolkit goal.
astrolyze2's `plot_image` (WCSAxes + beam ellipse) is the seed. We need one plotting
implementation that is both ergonomic and composable, with a consistent house style.

**Decision — (c) engine + sugar:**

- **Free plotting functions are the engine:** `plot_map(map, ax=...)`,
  `plot_spectrum(spec, ax=...)`, `plot_cube(...)`, `plot_pv(...)` — take an object + an
  `ax`, composable into any figure/grid (needed for multi-panel / Stack overviews).
- **Thin object methods are sugar:** `map.plot()` / `cube.plot()` / `spectrum.plot()` just
  call the corresponding free function. One implementation, two ergonomics (same DRY
  pattern as the unit layer, ADR-0003). Object context means colorbar units, beam, and WCS
  are auto-correct.

**House style — shipped, applied locally, NEVER global on import:**

- Ship a matplotlib style sheet (`astrolyze.mplstyle`) + conventions (always draw the beam,
  WCS ticks, colorbar with units from the object, perceptually-uniform default cmap).
- Apply it **locally** via a context manager (`with astrolyze.style(): ...`) inside the plot
  functions. **The library must never mutate global `rcParams` at import** — hijacking the
  user's matplotlib is a cardinal library sin. (This is the surprising/deliberate bit: a
  reader expecting `import astrolyze` to "just style everything" will find it does not, on
  purpose.)

## Riders (decided)

1. **matplotlib-only for v1, publication-first.** Leave a *possible* backend seam (so
   interactive plotly/glue could slot in later) but **do not build it** (YAGNI).
2. **Default colormap: `cividis`** — perceptually-uniform and prints well on white (radio
   astronomers don't always want dark themes). **Always overridable** by the caller.

## Consequences

- Plotting logic lives in functions, not entangled in the data classes; the classes stay
  about data + context.
- Consistent house look without global side effects; users keep full control of their
  matplotlib state.
- The unbuilt backend seam is a design constraint on the function signatures, not code to
  write now.
