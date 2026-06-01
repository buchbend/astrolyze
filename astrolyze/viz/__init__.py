"""Display: plotting engine + house style.

Free-function engine (`plot_map(map, ax=...)`) with thin object-method sugar. The house
style sheet is applied locally (context manager); astrolyze never mutates global
matplotlib rcParams on import. Default colormap is cividis, always overridable.
Implemented in issue #0005.
"""
