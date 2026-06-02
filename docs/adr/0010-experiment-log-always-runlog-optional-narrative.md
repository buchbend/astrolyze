# Experiment log: always auto-capture what ran; narrative offered, never enforced

ADR-0008 requires experiments to be legible — reconstructable after the fact. Two failure
modes are opposite: a pure machine log records *what ran* but not *why*; a mandated
narrative captures *why* but rots into box-ticking and isn't reliably maintained.

**Decision:**

- **"What ran" is always captured, automatically.** A machine-written, append-only
  structured run log (per run, in the experiment's `logs/`) records every astrolyze
  operation: params, inputs/outputs, software/data versions, timestamps. This is
  non-negotiable and does **not** depend on an agent remembering to write it — astrolyze
  operations emit structured records (via the `logging` already wired in astrolyze2).
- **Narrative is offered, not enforced.** The human-readable "why / what was found / what it
  means" account is an **offer**, never a required artifact, and never an inline obligation
  the agent must satisfy each session (that is what rots). The natural place to *compose*
  the narrative is the **future web UI** (the DB-backed manifest seam, ADR-0009) — it can
  surface the run log and invite a narrative on top.

This operationalizes verify-before-claiming: because the run log always records the real
artifacts, any narrative later written about results can be checked against what actually
ran — claims must point at real run-log/figure artifacts.

## Considered / rejected

- **Enforce a maintained narrative `LOG.md` inline** — rejected: enforcement breeds
  box-ticking and rot; the value of narrative is real only when voluntarily written.
- **Run log only** — rejected as the *whole* answer: it loses the scientific "why"; we keep
  the door open via the offer, not by mandate.

## Consequences

- Reproducibility ("what ran") is guaranteed and automatic; understanding ("why") is
  encouraged and tooled-for, not forced.
- Requires the run-log record format (likely per-run JSONL in `logs/`) — a build-time
  detail; keep it UI-readable so the future frontend can render and annotate it.
- The narrative-offer is itself a future-UI feature; not built now (YAGNI).
