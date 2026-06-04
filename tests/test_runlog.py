"""Tests for astrolyze.experiment.runlog — the always-on run log (issue #13, ADR-0010).

Written first (red/green TDD). Behavioural tests over the public surface — what a caller
observes: the per-run JSONL file on disk, the parsed records and their required keys, in-order
accumulation, append-only growth, and the inactive no-op — never private attributes:

- ``RunLog.open(experiment)`` starts a run and creates a per-run ``.jsonl`` file in ``logs/``;
- ``emit`` appends one parseable JSON record carrying op / params / inputs / outputs / software
  version / data version / timestamp;
- multiple ops accumulate, in order, in the same run file (append-only — never overwritten);
- separate runs get their own file/id so sessions don't tangle;
- with **no** active run ``emit`` is a no-op, so the library stays usable outside an experiment
  (the tracer spine runs unchanged);
- the existing operations (``load`` / ``save`` / ``moment*`` / ``.to`` / ``.plot``) each emit
  through the seam while a run is active.

Synthetic FITS are built header-first (mirroring tests/test_io.py and tests/test_tracer.py);
matplotlib runs head-less (Agg) because ``.plot()`` is exercised.
"""

import json
import warnings

import matplotlib

matplotlib.use("Agg")  # headless: .plot() is exercised below

import numpy as np
import pytest
import astropy.units as u
import radio_beam
from astropy.io import fits

from astrolyze.core import Cube
from astrolyze.experiment import Experiment
from astrolyze.experiment.runlog import RunLog, emit
from astrolyze.io import Metadata, load, save
from astrolyze.units import VelocityConvention

warnings.filterwarnings("ignore", module="spectral_cube")

REST_HZ = 230.538e9  # CO(2-1)


# --------------------------------------------------------------------------------------
# fixtures: an experiment + a schema-complete synthetic cube
# --------------------------------------------------------------------------------------
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


def _lines(run):
    text = run.path.read_text()
    return [line for line in text.splitlines() if line.strip()]


def _records(run):
    return [json.loads(line) for line in _lines(run)]


# --------------------------------------------------------------------------------------
# RunLog.open: starts a run, creating a per-run JSONL file in logs/
# --------------------------------------------------------------------------------------
def test_open_creates_per_run_jsonl_file_in_logs(experiment):
    with RunLog.open(experiment) as run:
        assert run.path.parent == experiment.logs
        assert run.path.suffix == ".jsonl"
        assert run.path.is_file()
        assert run.run_id  # a non-empty run id


def test_open_run_id_is_in_the_filename(experiment):
    with RunLog.open(experiment) as run:
        assert run.run_id in run.path.name


# --------------------------------------------------------------------------------------
# emit: one parseable JSON record carrying every required key
# --------------------------------------------------------------------------------------
def test_emit_appends_a_record_with_required_keys(experiment):
    with RunLog.open(experiment) as run:
        emit("load", inputs=["data/raw/x.fits"], outputs=[])

    records = _records(run)
    assert len(records) == 1
    record = records[0]
    for key in ("op", "params", "inputs", "outputs", "timestamp", "software", "data"):
        assert key in record, f"record is missing the {key!r} key"

    assert record["op"] == "load"
    assert record["inputs"] == ["data/raw/x.fits"]
    # software version: at least astrolyze's own version is captured.
    assert "astrolyze" in record["software"]
    # data version: the io schema version rides along.
    assert "schema_version" in record["data"]
    # timestamp is ISO-8601 (parseable) and run id is carried.
    from datetime import datetime

    assert datetime.fromisoformat(record["timestamp"])
    assert record["run_id"] == run.run_id


def test_emit_returns_the_record_when_a_run_is_active(experiment):
    with RunLog.open(experiment):
        record = emit("save", outputs=["data/interim/x.fits"])
    assert record is not None
    assert record["op"] == "save"


# --------------------------------------------------------------------------------------
# multiple ops accumulate, in order, append-only (never overwritten)
# --------------------------------------------------------------------------------------
def test_multiple_ops_accumulate_in_order(experiment):
    with RunLog.open(experiment) as run:
        emit("load", inputs=["a.fits"])
        emit("moment", params={"order": 0})
        emit("to", params={"unit": "K km/s"})

    assert [r["op"] for r in _records(run)] == ["load", "moment", "to"]


def test_records_are_append_only(experiment):
    with RunLog.open(experiment) as run:
        emit("load")
        first = run.path.read_text()
        emit("save")
        after = run.path.read_text()

    # The earlier record is preserved verbatim; the new one is appended after it.
    assert after.startswith(first)
    assert len(_lines(run)) == 2


# --------------------------------------------------------------------------------------
# each run gets its own file/id so separate sessions don't tangle
# --------------------------------------------------------------------------------------
def test_separate_runs_use_separate_files(experiment):
    with RunLog.open(experiment) as run1:
        emit("load")
    with RunLog.open(experiment) as run2:
        emit("save")

    assert run1.path != run2.path
    assert run1.run_id != run2.run_id
    assert [r["op"] for r in _records(run1)] == ["load"]
    assert [r["op"] for r in _records(run2)] == ["save"]


# --------------------------------------------------------------------------------------
# no active run: emit is a no-op (library usable outside an experiment)
# --------------------------------------------------------------------------------------
def test_emit_is_a_noop_without_an_active_run():
    # Outside any RunLog.open block, emit must write nothing and return None.
    assert emit("load", inputs=["a.fits"]) is None


def test_run_is_inactive_again_after_the_block(experiment):
    with RunLog.open(experiment):
        pass
    # The context closed: emit is a no-op once more.
    assert emit("load") is None


def test_spine_outside_experiment_emits_nothing(tmp_path):
    """The tracer spine run with no active run must not raise and must log nothing."""
    path = _synthetic_k_cube(tmp_path / "ngc0628_co21.fits")
    cube = Cube.from_loaded(load(path))
    cube.moment0().to("K km/s")  # no active run -> every emit is a silent no-op
    assert emit("noop") is None


# --------------------------------------------------------------------------------------
# the existing operations emit through the seam while a run is active
# --------------------------------------------------------------------------------------
def test_operations_emit_records_when_a_run_is_active(experiment, tmp_path):
    path = _synthetic_k_cube(tmp_path / "ngc0628_co21.fits")
    with RunLog.open(experiment) as run:
        cube = Cube.from_loaded(load(path))
        mom0 = cube.moment0().to("K km/s")
        mom0.plot()

    ops = [r["op"] for r in _records(run)]
    assert "load" in ops
    assert "moment" in ops
    assert "to" in ops
    assert "plot" in ops


def test_load_records_the_input_path(experiment, tmp_path):
    path = _synthetic_k_cube(tmp_path / "ngc0628_co21.fits")
    with RunLog.open(experiment) as run:
        load(path)
    load_records = [r for r in _records(run) if r["op"] == "load"]
    assert len(load_records) == 1
    assert str(path) in load_records[0]["inputs"]


def test_save_records_the_written_path(experiment, tmp_path):
    metadata = Metadata(
        object="NGC0628",
        rest_frequency=230.538 * u.GHz,
        velocity_convention=VelocityConvention.RADIO,
        beam=radio_beam.Beam(major=12 * u.arcsec, minor=10 * u.arcsec, pa=30 * u.deg),
        bunit=u.K,
        species="CO21",
    )
    with RunLog.open(experiment) as run:
        out_path = save(np.zeros((3, 3), dtype="float32"), metadata, experiment.interim)

    save_records = [r for r in _records(run) if r["op"] == "save"]
    assert len(save_records) == 1
    assert str(out_path) in save_records[0]["outputs"]


def test_moment_records_its_order(experiment, tmp_path):
    path = _synthetic_k_cube(tmp_path / "ngc0628_co21.fits")
    with RunLog.open(experiment) as run:
        Cube.from_loaded(load(path)).moment0()
    moment_records = [r for r in _records(run) if r["op"] == "moment"]
    assert len(moment_records) == 1
    assert moment_records[0]["params"]["order"] == 0
