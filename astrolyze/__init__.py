"""astrolyze — a toolkit for radio/sub-mm PPV cubes and spectra.

Importing astrolyze has no side effects: it does not configure logging handlers and
never mutates matplotlib's global rcParams (house style is applied locally by the viz
layer). See docs/adr for the design decisions behind the package.
"""

__version__ = "0.1.0.dev0"

__all__ = ["__version__"]
