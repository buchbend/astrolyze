"""Tests for astrolyze.experiment.layout — the fixed experiment skeleton (issue #10, ADR-0009).

Written first (red/green TDD). These are behavioral tests over the public surface — what a
caller observes on disk and from the value object — never private attributes:

- ``Experiment.init`` creates the exact ADR-0009 skeleton and a default ``config.toml``;
- ``init`` is idempotent (a second call is a no-op, raises nothing, and never clobbers a
  ``config.toml`` the user has already edited);
- ``Experiment(root)`` exposes the resolved path of every skeleton location;
- ``role_of`` classifies a path as raw / interim / processed / output (and ``None`` when the
  path is outside those four data/output trees);
- settings load from ``config.toml`` through dynaconf (the manifest DB URL is readable),
  proving the config seam later slices read from is wired.

``raw/`` is sacred: nothing in this subpackage writes to or renames it. That guarantee is
enforced in the ingest slice; here we only build the skeleton and never touch ``raw/`` after
creating it.
"""

import re
from pathlib import Path

from astrolyze.experiment import Experiment, Role, role_of

# The exact ADR-0009 skeleton, as paths relative to the experiment root.
SKELETON_DIRS = (
    "data/raw",
    "data/interim",
    "data/processed",
    "outputs/figures",
    "outputs/tables",
    "logs",
)


# --------------------------------------------------------------------------------------
# init: creates the fixed skeleton + a default config.toml
# --------------------------------------------------------------------------------------
def test_init_creates_full_skeleton(tmp_path):
    root = tmp_path / "study"
    Experiment.init(root)
    for relative in SKELETON_DIRS:
        assert (root / relative).is_dir(), f"missing skeleton dir: {relative}"


def test_init_writes_default_config(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    assert experiment.config.is_file()
    assert experiment.config.name == "config.toml"


def test_init_returns_experiment_rooted_at_directory(tmp_path):
    root = tmp_path / "study"
    experiment = Experiment.init(root)
    assert isinstance(experiment, Experiment)
    assert experiment.root == root


# --------------------------------------------------------------------------------------
# init: idempotent — second call is a no-op and never clobbers an edited config
# --------------------------------------------------------------------------------------
def test_init_is_idempotent(tmp_path):
    root = tmp_path / "study"
    Experiment.init(root)
    # A second init must not raise and must leave the skeleton intact.
    Experiment.init(root)
    for relative in SKELETON_DIRS:
        assert (root / relative).is_dir()


def test_init_preserves_an_existing_config(tmp_path):
    """Re-running init must not overwrite a config.toml the user has already edited —
    idempotent means the user's settings survive a second scaffold."""
    experiment = Experiment.init(tmp_path / "study")
    experiment.config.write_text('# edited by hand\n[manifest]\ndb_url = "sqlite:///mine.db"\n')
    Experiment.init(experiment.root)
    assert "mine.db" in experiment.config.read_text()


# --------------------------------------------------------------------------------------
# Experiment(root): resolves the path of every skeleton location
# --------------------------------------------------------------------------------------
def test_experiment_exposes_resolved_paths(tmp_path):
    root = tmp_path / "study"
    experiment = Experiment(root)
    assert experiment.raw == root / "data" / "raw"
    assert experiment.interim == root / "data" / "interim"
    assert experiment.processed == root / "data" / "processed"
    assert experiment.figures == root / "outputs" / "figures"
    assert experiment.tables == root / "outputs" / "tables"
    assert experiment.logs == root / "logs"
    assert experiment.config == root / "config.toml"


def test_experiment_accepts_string_root(tmp_path):
    experiment = Experiment(str(tmp_path / "study"))
    assert isinstance(experiment.root, Path)


# --------------------------------------------------------------------------------------
# role_of: classify a path as raw / interim / processed / output
# --------------------------------------------------------------------------------------
def test_role_of_classifies_each_tree(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    assert role_of(experiment, experiment.raw / "archival.fits") is Role.RAW
    assert role_of(experiment, experiment.interim / "harmonised.fits") is Role.INTERIM
    assert role_of(experiment, experiment.processed / "final.fits") is Role.PROCESSED
    assert role_of(experiment, experiment.figures / "mom0.png") is Role.OUTPUT
    assert role_of(experiment, experiment.tables / "fluxes.ecsv") is Role.OUTPUT


def test_role_of_returns_none_outside_the_data_and_output_trees(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    assert role_of(experiment, experiment.logs / "run.jsonl") is None
    assert role_of(experiment, experiment.config) is None
    assert role_of(experiment, tmp_path / "elsewhere" / "file.fits") is None


def test_role_of_is_available_as_a_method(tmp_path):
    experiment = Experiment(tmp_path / "study")
    assert experiment.role_of(experiment.raw / "x.fits") is Role.RAW


# --------------------------------------------------------------------------------------
# config: settings load from config.toml via dynaconf (the seam later slices read)
# --------------------------------------------------------------------------------------
def test_settings_load_from_config(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    settings = experiment.settings
    # The default config records a sqlite manifest DB URL; later slices read it from here.
    assert "sqlite" in str(settings.get("manifest.db_url", ""))


def test_settings_reflect_edited_config(tmp_path):
    experiment = Experiment.init(tmp_path / "study")
    experiment.config.write_text('[manifest]\ndb_url = "sqlite:///custom.db"\n')
    assert "custom.db" in str(Experiment(experiment.root).settings.get("manifest.db_url", ""))


# --------------------------------------------------------------------------------------
# House rule: no SystemExit / print in library code (ADR-0005 / the legacy sin)
# --------------------------------------------------------------------------------------
def test_no_systemexit_or_print_in_library_code():
    import astrolyze.experiment as experiment_pkg

    pkg_dir = Path(experiment_pkg.__file__).parent
    for src_file in pkg_dir.glob("*.py"):
        src = src_file.read_text()
        assert "SystemExit" not in src, f"{src_file.name} uses SystemExit"
        assert not re.search(r"\bprint\s*\(", src), f"{src_file.name} calls print()"
