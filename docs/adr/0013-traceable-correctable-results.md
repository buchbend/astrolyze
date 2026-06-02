# Traceable, correctable results — legibility over enforcement

The problem to solve is being **outpaced**: a result appears and the user cannot easily see
how it came to be, so they cannot follow along or correct it. The fix is **transparency of
derivation**, NOT gating. Enforcement that slows the user down is the wrong tool — it is the
"stones in the way" anti-pattern (ADR-0011: astrolyze enables, never obstructs). The burden
sits on the assistant/tooling to *surface* how a result was produced, never on the user to
wade through checkpoints.

## Decision

**Every result carries its derivation, legibly, so the user can follow and correct.**

- **Traceability:** when a result/figure/number is produced, the *how* is visible alongside
  it — the operations and inputs that produced it (from the always-on Run Log, ADR-0010) and
  the reasoning behind the choices. The user can answer "how did this come to be?" without
  asking.
- **Correctability:** because the derivation is visible, the user can intervene and correct —
  catch a wrong assumption, a wrong approach, a misread — at the point it matters.
- **Verify-before-claiming (kept):** a stated result must cite the real artifact it came from
  (figure/table/run-log entry), so claims are checkable, not asserted. (This session's
  lesson; in memory.)
- **No user-facing friction:** this does **not** gate the user's tool calls or force
  reason-before-act checkpoints on them. Legibility is produced *as a byproduct* of doing the
  work through astrolyze (which logs itself, ADR-0010/0011), not by adding stop-gates.

## What changed from the first draft of this ADR

The first version centred on **enforcement/pacing** (a hard hook + reason-before-act at
"forks"). Rejected by the user: that slows the user down, which is the opposite of the goal.
The goal is to *follow easily and correct*, achieved by **surfacing derivation**, not by
gating. Pacing/gating is out; traceability + correctability is in.

## How (lightweight, build-time detail)

- Lean on the always-on Run Log (ADR-0010) for "what ran" — it already exists with no user
  burden.
- The assistant surfaces the reasoning/derivation with results (and in the experiment's
  `logs/`), in plain language, citing artifacts — so the trail is followable and correctable.
- Keep it plain markdown + the existing log; standalone, no lore dependency (lore is passive
  cross-session capture, alpha — different concern).
- If any automation is added later, it must reduce the assistant's opacity, never add user
  friction.

## Consequences

- Legibility is the assistant's job, surfaced by default; the user is never the one slowed.
- Reconciles with ADR-0008 (legibility pillar) and ADR-0011 (enabling, not obstructing) — and
  drops the friction the first draft wrongly introduced.
- Exact surfacing format is a build-time detail, deliberately light.
