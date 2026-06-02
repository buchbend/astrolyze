# Agent integration: usage guide + CLI/API, no MCP

Given "one interface, AI uses astrolyze as a human would" (ADR-0011), agents integrate by
**writing astrolyze Python/CLI exactly as a human does** — producing the same legible,
logged, reproducible artifacts (scripts, runbooks, notebooks). There is **no MCP / no
AI-only interface**.

astrolyze therefore ships, as first-class deliverables:
- a clean Python API + a typer/rich **CLI**,
- honest, accurate **docstrings**,
- an agent-facing **usage guide** (a `AGENTS.md` / CLAUDE.md-style "how we work with and
  display astro data" doc) that teaches agents the house way and to reach for astrolyze by
  default (ADR-0008 agent-native).

## Why not an MCP

An MCP would be an **AI-only second interface** — humans wouldn't use it — which violates
ADR-0011 principle #1 (astrolyze AI-independent, one interface) and reintroduces the very
human/AI divergence ADR-0008 exists to prevent. An agent writing an astrolyze script *is*
the artifact we want; an MCP makes the agent's actions less like a human's, not more.

## Consequences

- The guide + docstrings are part of astrolyze's deliverable surface, not afterthoughts;
  they must stay honest (verify-before-claiming applies to docs too).
- **No MCP now and no seam obligation**; if a concrete need ever earns one (rule of three),
  it must be a thin transport wrapping the *same* CLI/API — never a parallel logic path.
- Reuse existing skills (`gildas-class`, `heterodyne-calibration`, `radio-astronomer`,
  `tdd`) by wiring astrolyze conventions into them; do not build narrow policing skills
  (correctness is enforced in astrolyze itself — ADR-0011 #5).
