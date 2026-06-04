"""Tests for astrolyze.experiment.narrative — the narrative *offer* (issue #14, ADR-0010).

Written first (red/green TDD). Behavioural tests over the public surface — what a caller
observes from :func:`narrate` and the note on disk — never internals:

- ``narrate(experiment)`` scaffolds a human-readable **markdown** note in the experiment and
  returns its path;
- the note **references real run-log artifacts** — the run-log file of the run it documents and
  the operations that actually ran (ADR-0010/0013), so a narrative is checkable against what ran;
- it is **offered, never enforced**: idempotent and it never overwrites a note you have started
  editing (re-running ``narrate`` opens the existing note, it does not clobber your prose);
- a specific run can be targeted by id; with **no** run recorded yet the offer still produces a
  note (it never crashes — narrative is never a required artifact).

A run is produced the ordinary way: open a :class:`RunLog` and run the spine on a tiny
synthetic cube, exactly as a human/agent would (no private hooks).
"""

import warnings

import matplotlib

matplotlib.use("Agg")  # headless: the spine calls .plot()

import numpy as np
import pytest
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.experiment import Experiment, RunLog, narrate
from astrolyze.io import load

warnings.filterwarnings("ignore", module="spectral_cube")

REST_HZ = 230.538e9  # CO(2-1)


def _synthetic_k_cube(path):
    """A small schema-complete brightness-temperature cube (K, radio convention)."""
    h = fits.Header()
    h["CTYPE1"], h["CRVAL1"] = "RA---SIN", 24.174
    h["CDELT1"], h["CRPIX1"], h["CUNIT1"] = -2e-4, 1.0, "deg"
    h["CTYPE2"], h["CRVAL2"] = "DEC--SIN", 15.784
    h["CDELT2"], h["CRPIX2"], h["CUNIT2"] = 2e-4, 1.0, "deg"
    h["CTYPE3"], h["CRVAL3"] = "VRAD", 0.0
    h["CDELT3"], h["CRPIX3"], h["CUNIT3"] = 2000.0, 1.0, "m/s"
    h["OBJECT"], h["TELESCOP"], h["BUNIT"] = "NGC0628", "ALMA", "K"
    h["RESTFRQ"] = (REST_HZ, "Hz")
    h["BMAJ"], h["BMIN"], h["BPA"] = 3.33e-4, 2.78e-4, 30.0
    data = np.arange(5 * 6 * 6, dtype="float32").reshape((5, 6, 6))
    fits.writeto(path, data, h)
    return path


@pytest.fixture
def experiment(tmp_path):
    return Experiment.init(tmp_path / "study")


def _run_the_spine(experiment, tmp_path, name="ngc0628_co21.fits") -> RunLog:
    """Run load -> moment0 -> .to -> plot with a run open; return the (closed) RunLog."""
    path = _synthetic_k_cube(tmp_path / name)
    with RunLog.open(experiment) as run:
        Cube.from_loaded(load(path)).moment0().to("K km/s").plot()
    return run


# --------------------------------------------------------------------------------------
# narrate scaffolds a markdown note in the experiment and returns its path
# --------------------------------------------------------------------------------------
def test_narrate_scaffolds_a_markdown_note(experiment, tmp_path):
    _run_the_spine(experiment, tmp_path)
    note = narrate(experiment)
    assert note.suffix == ".md"
    assert note.is_file()
    # the note lives inside the experiment (legible, browsable next to the run logs)
    assert note.resolve().is_relative_to(experiment.root.resolve())


# --------------------------------------------------------------------------------------
# the note references the REAL run-log artifacts (the run file + the ops that ran)
# --------------------------------------------------------------------------------------
def test_note_references_the_run_log_file(experiment, tmp_path):
    run = _run_the_spine(experiment, tmp_path)
    text = narrate(experiment).read_text()
    # The narrative points at the machine record of what actually ran (ADR-0010/0013).
    assert run.path.name in text


def test_note_references_the_recorded_operations(experiment, tmp_path):
    _run_the_spine(experiment, tmp_path)
    text = narrate(experiment).read_text()
    # The ops the run actually performed are surfaced so the account is checkable.
    for op in ("load", "moment", "to", "plot"):
        assert op in text


# --------------------------------------------------------------------------------------
# offered, never enforced: idempotent, never clobbers a note you have started editing
# --------------------------------------------------------------------------------------
def test_narrate_is_idempotent_and_preserves_edits(experiment, tmp_path):
    _run_the_spine(experiment, tmp_path)
    note = narrate(experiment)
    note.write_text("# My findings\n\nThe inner ring lights up in CO(2-1).\n")
    again = narrate(experiment)
    assert again == note
    # Re-running the offer must not overwrite the human's prose.
    assert "inner ring lights up" in again.read_text()


# --------------------------------------------------------------------------------------
# a specific run can be targeted by id
# --------------------------------------------------------------------------------------
def test_narrate_targets_a_specific_run(experiment, tmp_path):
    first = _run_the_spine(experiment, tmp_path, name="first.fits")
    second = _run_the_spine(experiment, tmp_path, name="second.fits")
    assert first.run_id != second.run_id

    note_first = narrate(experiment, run_id=first.run_id)
    assert first.path.name in note_first.read_text()
    # default (no run_id) documents the latest run, not the first
    note_latest = narrate(experiment)
    assert second.path.name in note_latest.read_text()
    assert note_first != note_latest


# --------------------------------------------------------------------------------------
# the offer never crashes when nothing has been recorded yet
# --------------------------------------------------------------------------------------
def test_narrate_with_no_runs_still_offers_a_note(experiment):
    note = narrate(experiment)  # no run has been opened in this fresh experiment
    assert note.suffix == ".md"
    assert note.is_file()
