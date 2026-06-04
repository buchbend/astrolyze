"""End-to-end acceptance for the experiment layer (issue #14, PRD 0002).

Drives the committed example script's ``main()`` — the whole path a human/agent runs:

    init -> place a cube in raw/ -> ingest -> analyse (run log open) -> manifest list -> narrate

and asserts the observable outcome of the full slice working together:

- the complete cube in ``raw/`` is **ingested and registered** in the manifest;
- the analysis ran **with a run log active**, so ``logs/`` holds the JSONL records for the ops
  that ran (load / moment / to / save / plot), in order;
- the derived product landed in ``processed/`` and the figure in ``figures/``;
- the **narrative offer** produced a note referencing the real run-log artifacts.

This is the capstone test: it exercises layout + ingest + manifest + runlog + narrative through
one ordinary script, not their internals.
"""

import importlib.util
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: the example saves a figure

import numpy as np
import pytest
from astropy.io import fits

from astrolyze.experiment import Manifest

warnings.filterwarnings("ignore", module="spectral_cube")

REST_HZ = 230.538e9  # CO(2-1)
EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "experiment_ngc628.py"


def _load_example_main():
    """Import ``main`` from the example script (not a package) by file path."""
    spec = importlib.util.spec_from_file_location("experiment_ngc628_example", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.main


def _complete_cube(path):
    """A schema-complete brightness-temperature cube (K, radio convention, beam, species)."""
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
    h["HIERARCH ASTROLYZE SPECIES"] = "CO21"
    fits.writeto(path, np.arange(5 * 6 * 6, dtype="float32").reshape((5, 6, 6)), h)
    return path


@pytest.fixture
def ran_experiment(tmp_path):
    """Run the example end-to-end and hand back the resulting Experiment."""
    main = _load_example_main()
    cube = _complete_cube(tmp_path / "ngc0628_co21.fits")
    return main(tmp_path / "study", cube)


# --------------------------------------------------------------------------------------
# ingest -> manifest: the complete cube is registered
# --------------------------------------------------------------------------------------
def test_ingested_dataset_is_in_the_manifest(ran_experiment):
    records = Manifest.for_experiment(ran_experiment).all()
    assert [r.metadata.object for r in records] == ["NGC0628"]
    assert records[0].source_path == "data/raw/ngc0628_co21.fits"
    assert records[0].metadata.species == "CO21"


# --------------------------------------------------------------------------------------
# run log: logs/ holds the JSONL records for the computed ops, in order
# --------------------------------------------------------------------------------------
def test_logs_hold_the_runlog_records_for_the_computation(ran_experiment):
    run_files = list(ran_experiment.logs.glob("run-*.jsonl"))
    assert len(run_files) == 1

    records = [
        json.loads(line)
        for line in run_files[0].read_text().splitlines()
        if line.strip()
    ]
    ops = [r["op"] for r in records]
    # the whole analysed spine was captured automatically, in execution order
    assert ops == ["load", "moment", "to", "save", "plot"]
    # every record carries the required envelope
    for record in records:
        for key in (
            "op",
            "params",
            "inputs",
            "outputs",
            "software",
            "data",
            "timestamp",
        ):
            assert key in record


# --------------------------------------------------------------------------------------
# derived product + figure landed in the right trees
# --------------------------------------------------------------------------------------
def test_derived_product_and_figure_are_written(ran_experiment):
    processed = list(ran_experiment.processed.glob("*.fits"))
    figures = list(ran_experiment.figures.glob("*.png"))
    assert len(processed) == 1 and processed[0].stat().st_size > 0
    assert len(figures) == 1 and figures[0].stat().st_size > 0
    # raw/ is sacred: the one file we put in is still there, unmodified in name
    assert [p.name for p in ran_experiment.raw.glob("*.fits")] == ["ngc0628_co21.fits"]


# --------------------------------------------------------------------------------------
# narrative offer: a note referencing the real run-log artifacts
# --------------------------------------------------------------------------------------
def test_narrative_note_references_real_artifacts(ran_experiment):
    notes = list(ran_experiment.logs.glob("*.md"))
    assert len(notes) == 1
    text = notes[0].read_text()
    run_file = next(ran_experiment.logs.glob("run-*.jsonl"))
    assert run_file.name in text  # points at the machine record of what ran
    for op in ("load", "moment", "to", "save", "plot"):
        assert op in text
