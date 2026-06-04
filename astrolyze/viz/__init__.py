"""Display: plotting engine + house style (ADR-0005).

Free-function engine (``plot_map(map, ax=...)``, ``plot_spectrum``, ``plot_cube``) with thin
object-method sugar (``obj.plot()`` in :mod:`astrolyze.core`). The house style is shipped
(``astrolyze.mplstyle``) and applied **locally** by :func:`style` — astrolyze never mutates
matplotlib's global rcParams on import. Default colormap is cividis, always overridable.

Importing this module pulls in matplotlib (it is the plotting layer); importing the top-level
``astrolyze`` package does not — ``astrolyze.style`` resolves here lazily on first access.
"""

from .engine import DEFAULT_CMAP, plot_cube, plot_map, plot_noise, plot_spectrum
from .style import style

__all__ = [
    "plot_map",
    "plot_spectrum",
    "plot_cube",
    "plot_noise",
    "style",
    "DEFAULT_CMAP",
]
