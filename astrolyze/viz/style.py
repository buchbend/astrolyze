"""The house style, applied **locally** (ADR-0005).

The shipped sheet ``astrolyze.mplstyle`` is never installed into matplotlib's global
``rcParams`` on import — hijacking the user's matplotlib is a cardinal library sin. Instead
the plot functions wrap their drawing in :func:`style`, a context manager over
:func:`matplotlib.pyplot.style.context`, which saves and restores ``rcParams`` exactly. The
user keeps full control of their matplotlib state; the toolkit gets a consistent look only
where it draws.
"""

from __future__ import annotations

from contextlib import contextmanager
from importlib.resources import files

import matplotlib.pyplot as plt

#: Path to the shipped style sheet (package data, see pyproject ``package-data``).
STYLE_PATH = files("astrolyze.viz") / "astrolyze.mplstyle"


@contextmanager
def style(extra=None):
    """Apply the astrolyze house style for the duration of the ``with`` block.

    Usage: ``with astrolyze.style(): ...``. On exit matplotlib's global ``rcParams`` are
    restored to exactly what they were — astrolyze never leaves the house style behind.

    Parameters
    ----------
    extra :
        Optional additional style(s) stacked on top of the house sheet — a style name, a
        path, an rcParams ``dict``, or an iterable of those (same vocabulary as
        :func:`matplotlib.style.use`). Lets a caller tweak the look without losing the base.
    """
    styles = [str(STYLE_PATH)]
    if extra is not None:
        if isinstance(extra, (str, dict)):
            styles.append(extra)
        else:
            styles.extend(extra)
    with plt.style.context(styles):
        yield
