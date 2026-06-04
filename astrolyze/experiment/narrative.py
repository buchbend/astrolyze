"""The narrative *offer* (ADR-0010).

The run log captures *what* ran, automatically. The narrative is the other half — the human
"why I did this / what I found / what it means" — and it is deliberately **offered, never
enforced**: a mandated narrative rots into box-ticking, so astrolyze only ever *invites* one.

:func:`narrate` scaffolds (or re-opens) a small markdown note for a run, pre-filled with a
reference to that run's machine record and the operations it actually performed, plus empty
prose sections to fill in. Two properties make the offer safe to ignore and safe to repeat:

- **It references real run-log artifacts.** The note points at the run's ``run-<id>.jsonl`` and
  lists the ops/outputs that run recorded, so any account written here is checkable against what
  actually ran (ADR-0010/0013 — verify-before-claiming). The note is *not* the source of truth;
  the run log is.
- **It never clobbers your prose.** Re-running :func:`narrate` returns the existing note
  untouched (idempotent), so the offer is repeatable and the human's writing is never lost. The
  note lives beside the run logs in ``logs/`` — browsable, and obviously paired with its run.

Minimal by design: composing a rich narrative is a future-UI concern (ADR-0010), so this is a
thin scaffold, not an editor.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from . import runlog

if TYPE_CHECKING:
    from .layout import Experiment

# The general note used when an experiment has no run recorded yet — the offer still stands.
_UNANCHORED_NOTE = "narrative.md"

# The invitation header, written once into a fresh note and meant to be deleted by the author.
_OFFER_COMMENT = """\
<!--
This note is an *offer*, never a requirement (ADR-0010). Write the scientific "why / what I
found / what it means" below when it is genuinely worth writing — astrolyze never mandates it.
Keep your account checkable: point at the real run-log / figure artifacts referenced above (the
run log is the machine record of what actually ran). Delete this comment when you start.
-->"""

_PROSE_SECTIONS = "## Why\n\n\n## What I found\n\n\n## What it means\n"


def narrate(experiment: "Experiment", *, run_id: str | None = None) -> Path:
    """Scaffold (or re-open) the markdown narrative note for a run, and return its path.

    By default the note documents the **latest** run; pass *run_id* to target a specific one.
    The note is created beside the run logs in ``logs/`` and pre-filled with a reference to the
    run's ``run-<id>.jsonl`` and the operations it recorded. It is an **offer**: idempotent, and
    an existing note is returned untouched (never overwriting prose you have started). When the
    experiment has no run recorded yet, a general note is still produced — the offer never
    fails, because narrative is never a required artifact.
    """
    logs = experiment.logs
    logs.mkdir(parents=True, exist_ok=True)

    run_path = _resolve_run(logs, run_id)
    note_path = (
        run_path.with_suffix(".md") if run_path is not None else logs / _UNANCHORED_NOTE
    )
    if note_path.exists():
        return note_path  # never clobber an existing (possibly edited) note — open, don't overwrite
    note_path.write_text(_scaffold(experiment, run_path), encoding="utf-8")
    return note_path


# -- helpers ---------------------------------------------------------------------------
def _resolve_run(logs: Path, run_id: str | None) -> Path | None:
    """The run-log file to document: the one named by *run_id*, else the most recent. Run-log
    filenames lead with a UTC timestamp, so the lexicographic max is the latest run. Returns
    ``None`` when no run has been recorded yet."""
    if run_id is not None:
        return logs / f"run-{run_id}.jsonl"
    runs = sorted(logs.glob("run-*.jsonl"))
    return runs[-1] if runs else None


def _scaffold(experiment: "Experiment", run_path: Path | None) -> str:
    """The starter markdown for a fresh note: a header, the run-log reference block (real
    artifacts), the offer comment, and empty prose sections."""
    if run_path is None:
        header = f"# Narrative — {experiment.root.name or experiment.root}"
        reference = (
            "- **Experiment:** `"
            + str(experiment.root)
            + "`\n- **Run log:** none recorded yet — run an analysis with a run open "
            "(`RunLog.open(experiment)`), then re-run `astrolyze narrate` to document it."
        )
        return f"{header}\n\n{reference}\n\n{_OFFER_COMMENT}\n\n{_PROSE_SECTIONS}"

    label = run_path.stem  # "run-<id>"
    entries = runlog.read(run_path)
    operations = _ordered_ops(entries)
    artifacts = _output_artifacts(entries, experiment.root)

    header = f"# Narrative — {label}"
    reference = "\n".join(
        (
            f"- **Experiment:** `{experiment.root}`",
            # The note sits in logs/ beside the run log, so a bare filename is the right link.
            f"- **Run log:** [`{run_path.name}`]({run_path.name}) — the machine record of what "
            "actually ran",
            "- **Operations recorded:** "
            + (" → ".join(operations) if operations else "— (none recorded)"),
            "- **Output artifacts:** "
            + (
                ", ".join(f"`{a}`" for a in artifacts)
                if artifacts
                else "— (none recorded)"
            ),
        )
    )
    return f"{header}\n\n{reference}\n\n{_OFFER_COMMENT}\n\n{_PROSE_SECTIONS}"


def _ordered_ops(entries: list[dict]) -> list[str]:
    """The distinct operations the run recorded, in first-seen order (load, moment, to, …)."""
    seen: list[str] = []
    for entry in entries:
        op = entry.get("op")
        if op is not None and op not in seen:
            seen.append(op)
    return seen


def _output_artifacts(entries: list[dict], root: Path) -> list[str]:
    """Every output path the run recorded, de-duplicated and made relative to the experiment
    where possible (so the note reads in experiment-local terms, like the manifest)."""
    artifacts: list[str] = []
    for entry in entries:
        for output in entry.get("outputs", []):
            ref = _relativize(output, root)
            if ref not in artifacts:
                artifacts.append(ref)
    return artifacts


def _relativize(output, root: Path) -> str:
    """``output`` as a path relative to the experiment *root* when it lies inside it, else
    unchanged (a manifest id or an external path passes through)."""
    try:
        return Path(str(output)).resolve().relative_to(Path(root).resolve()).as_posix()
    except (ValueError, OSError):
        return str(output)


__all__ = ["narrate"]
