# All data work goes through astrolyze; formalize-before-use; astrolyze is AI-independent

This is the operating principle of the whole collaboration — the generalization of the
hard wall (ADR-0001) from *code promotion* to *every data operation*.

**Decision:**

1. **astrolyze is primary and AI-independent.** It must be fully usable by a human with no
   agent in the loop — a clean Python API + a CLI (typer/rich). The AI layer is an *optional*
   accelerator on top, never a requirement and never a parallel path.
2. **No ad-hoc data twiddling.** All data handling — human or AI — goes **through astrolyze**
   and is **logged by astrolyze** (ADR-0010). There is no off-to-the-side numpy hack on a
   cube that escapes the toolkit and the log.
3. **Formalize before use.** To do something astrolyze can't yet do, you **extend astrolyze
   first** (with tests, by its conventions), then use it. New capability is added to the
   tool, not improvised in an experiment. Proven, reused extensions graduate per the rule of
   three (ADR-0001).
4. **AI does it as a human would.** When present, an agent uses astrolyze the same way a
   human does — by writing astrolyze Python/CLI pipelines, runbooks, notebooks — producing
   the same legible, logged, reproducible artifacts. One interface, used identically by human
   and AI.
5. **Conventions live in the tool, not in policing prompts.** Correctness is enforced *in
   astrolyze* (e.g. explicit convention/rest-frequency or it raises, ADR-0003), so no narrow
   "cop" skills (e.g. `unit-check`) are needed. The "how we work with astro data" knowledge
   lives in guidelines; the guardrails live in code.
6. **astrolyze enables, it does not obstruct.** The house way must be the path of *least*
   resistance — ergonomic enough that reaching around it is never tempting. (The merciless
   ingest gate, ADR-0009, is not a contradiction: it enforces data *quality* so downstream
   work is trustworthy; it is a gate at the door, not friction in daily use.)

## Corrections this forces

- **Ingest is astrolyze functionality, not an AI skill** (revises the Q14 framing). An agent
  may *drive* ingest as a human would, but the ingest logic lives in the toolkit/CLI.
- The AI-layer artifact list collapses from "many narrow skills" toward "a usage guide + the
  conventions enforced in code."

## Consequences

- Because astrolyze logs itself and raises on bad physics, no AI-specific logging path or
  validation cop is needed — a major simplification.
- Discipline cost: genuinely new analysis requires a tool extension first (slower than a
  one-off hack), bought back as durability, reuse, and trust.
- The agent-integration *mechanism* (usage guide vs. MCP) is decided separately; current
  lean is guide + CLI, no MCP (an MCP would be an AI-only second interface, against
  principle #1).
