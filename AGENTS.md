# Working with astro data in astrolyze (guide for humans and agents)

This is the house guide for using astrolyze. A coding agent uses astrolyze the **same way a
human does** — by writing astrolyze Python/CLI — producing ordinary, reviewable, reproducible
scripts. There is no AI-only interface.

> Placeholder — filled out in issue #0006 once the API exists. The principles below are fixed.

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
