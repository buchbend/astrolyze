#!/usr/bin/env python
"""The whole experiment path, end-to-end (PRD 0002 acceptance, issue #14).

One analysis from an empty directory to a documented result, written exactly as a human (or a
coding agent, identically — ADR-0011/0012) would write it:

    init  ->  place a cube in raw/  ->  ingest  ->  analyse (with a run log open)  ->
    manifest list  ->  narrate

Every step is an ordinary astrolyze call; the corresponding shell one-liners are
``astrolyze init / ingest / manifest list / narrate`` plus the ``moment0`` tracer. Two things
ride along for free and are the point of the experiment layer:

- **the merciless gate** — ``ingest`` validates ``raw/`` against the metadata schema and only
  registers what carries the mandatory physical context (rest frequency, velocity convention);
- **the always-on run log** — with a run open, ``load`` / ``moment`` / ``.to`` / ``save`` /
  ``.plot`` each append a record to ``logs/``, so *what ran* is captured without anyone
  remembering to log it. ``narrate`` then offers (never forces) a human note over that record.

Run it self-contained on the committed NGC 628 CO(2-1) cutout::

    python examples/experiment_ngc628.py                       # uses tests/data + a temp dir
    python examples/experiment_ngc628.py ./ngc628-study        # keep the experiment here
    python examples/experiment_ngc628.py ./ngc628-study /path/to/cube.fits
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

from astrolyze.core import Cube
from astrolyze.experiment import Experiment, Manifest, ingest, narrate
from astrolyze.experiment.runlog import RunLog
from astrolyze.io import load, save

# The committed real-data cutout (see tests/data/PROVENANCE.md) — the default input, so the
# example runs with no arguments at all.
DEFAULT_CUBE = (
    Path(__file__).resolve().parent.parent
    / "tests"
    / "data"
    / "ngc0628_co21_cutout.fits.gz"
)


def main(experiment_dir: str | Path, cube_path: str | Path) -> Experiment:
    cube_path = Path(cube_path)

    # 1. Scaffold the fixed experiment skeleton (idempotent). data/{raw,interim,processed},
    #    outputs/{figures,tables}, logs/, config.toml — the same shape every study uses.
    exp = Experiment.init(experiment_dir)
    print(f"init      {exp.root}")

    # 2. Place the cube in raw/. raw/ is sacred: we (the user) put data in; astrolyze only ever
    #    reads it, never renames or rewrites it. Real inputs keep their upstream names.
    raw_cube = exp.raw / cube_path.name
    if not raw_cube.exists():
        shutil.copy2(cube_path, raw_cube)
    print(
        f"placed    {raw_cube.relative_to(exp.root)}  (raw/ is sacred — read-only to astrolyze)"
    )

    # 3. Merciless ingest: validate every raw file's header and register what is complete. A
    #    file missing rest frequency / velocity convention would be *rejected* here, naming the
    #    gap — nothing incomplete reaches the manifest.
    report = ingest(exp)
    print(f"ingest    {report.n_accepted} accepted, {report.n_rejected} rejected")
    for rejected in report.rejected:
        print(
            f"            rejected {rejected.source_path}: missing {rejected.missing}"
        )

    # 4. Analyse with a run log open. From here every operation appends a record to logs/ — the
    #    machine record of what actually ran (load -> moment0 -> .to -> save -> plot).
    with RunLog.open(exp) as run:
        cube = Cube.from_loaded(load(raw_cube))
        mom0 = cube.moment0().to(
            "K km/s"
        )  # beam / rest freq / convention come from the cube

        # a derived product in processed/, named by the header projection (ADR-0006)
        derived = save(np.asarray(mom0.data.value), mom0.metadata, exp.processed)
        # a house-style figure in figures/ (cividis on WCS axes, beam ellipse, unit colorbar)
        fig, _ = mom0.plot()
        figure = exp.figures / f"{raw_cube.stem.split('.')[0]}_moment0.png"
        fig.savefig(figure, bbox_inches="tight", dpi=150)
    print(f"analysed  run {run.run_id}")
    print(
        f"            run log -> {run.path.relative_to(exp.root)} ({len(run.entries())} records)"
    )
    print(f"            derived -> {derived.relative_to(exp.root)}")
    print(f"            figure  -> {figure.relative_to(exp.root)}")

    # 5. The manifest: see what is in the experiment at a glance, queryable — no grepping
    #    filenames. (Read-only; ingest generates and keeps it in sync.)
    manifest = Manifest.for_experiment(exp)
    print("manifest  registered datasets:")
    for record in manifest.all():
        m = record.metadata
        print(f"            #{record.id} {record.source_path}  {m.object} {m.species}")

    # 6. The narrative offer: a markdown note scaffolded over the run log, referencing the real
    #    artifacts above. Offered, never enforced — fill in the "why" if it is worth writing.
    note = narrate(exp)
    print(
        f"narrate   {note.relative_to(exp.root)}  (optional — open it to add the why)"
    )

    return exp


if __name__ == "__main__":
    where = sys.argv[1] if len(sys.argv) > 1 else None
    cube = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_CUBE
    if not Path(cube).exists():
        raise SystemExit(
            f"cube not found: {cube}\n"
            "usage: python examples/experiment_ngc628.py [experiment_dir] [cube.fits]\n"
            "   (defaults to the committed tests/data cutout)"
        )
    # Default to a throwaway temp experiment so the example is runnable with zero arguments;
    # pass a directory to keep the scaffolded experiment around.
    experiment_dir = where or tempfile.mkdtemp(prefix="ngc628-study-")
    main(experiment_dir, cube)
    if where is None:
        print(
            f"\n(experiment scaffolded under {experiment_dir} — pass a path to keep one in place)"
        )
