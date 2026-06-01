"""astrolyze — a toolkit for radio/sub-mm PPV cubes and spectra.

Importing astrolyze has no side effects: it does not configure logging handlers and
never mutates matplotlib's global rcParams (house style is applied locally by the viz
layer). See docs/adr for the design decisions behind the package.
"""

__version__ = "0.1.0.dev0"

__all__ = ["__version__", "style"]


def __getattr__(name):
    """Resolve ``astrolyze.style`` lazily so plain ``import astrolyze`` stays matplotlib-free.

    The house-style context manager lives in :mod:`astrolyze.viz`, which imports matplotlib.
    Exposing it via PEP 562 module ``__getattr__`` keeps that cost off the bare package import
    while still letting users write ``with astrolyze.style(): ...`` (ADR-0005)."""
    if name == "style":
        from astrolyze.viz import style

        return style
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
