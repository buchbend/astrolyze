"""Scaffold smoke tests — prove the package imports cleanly with no side effects.

These guard issue #0001's acceptance: install works, import succeeds, and importing
astrolyze does not hijack the user's matplotlib (ADR-0005).
"""


def test_import_and_version():
    import astrolyze

    assert isinstance(astrolyze.__version__, str)


def test_subpackages_import():
    import astrolyze.core  # noqa: F401
    import astrolyze.io  # noqa: F401
    import astrolyze.units  # noqa: F401
    import astrolyze.viz  # noqa: F401


def test_import_does_not_mutate_global_rcparams():
    import matplotlib

    before = matplotlib.rcParams.copy()
    import importlib

    import astrolyze

    importlib.reload(astrolyze)
    after = matplotlib.rcParams.copy()
    assert before == after, "importing astrolyze must not mutate global rcParams"
