# Design: alpha-release tester onboarding (a clean **research template repo** alpha users fork)

**Date:** 2026-06-06 (decisions finalized 2026-06-08)
**Status:** Ready for review. All core decisions locked (§3); remaining items in §10 are build-time verifications, not design forks.
**Related:**
- `multiomics_explorer/docs/superpowers/specs/2026-06-01-explorer-package-release-design.md` — the explorer becomes a git-installable package; the template consumes it transitively via `uv sync`.
- `multiomics_biocypher_kg/plans/alpha_release.md` — Track A (lab box `132.75.249.47`); owns KG hosting + the canonical **connection** guide (`docs/kg_mcp_guide.md`, §2.7 / Workstream C) + `/release-kg`.

## 1. Overview & goal

Get lab alpha testers from "empty laptop" to "running their first KG-backed analysis in Claude Code" with the least friction. Delivery model (decided 2026-06-08):

**A new, clean `multiomics_research_template` repo that testers fork. One clone of their fork carries everything Claude-facing; `uv sync` carries everything Python.** Forking (rather than bare-cloning) gives per-tester attribution and a place to commit usage logs + analyses back — see §9. The current `multiomics_research` repo is frozen as a **historical artifact** (retains the maintainer's analyses + full history) and is no longer the tester target.

The key insight that drove this: the three pieces a tester needs are *different kinds of thing* and arrive by different mechanisms — but a **single clone + `uv sync`** delivers all of them:

| Piece | What it is | Delivered how | Lands where |
|---|---|---|---|
| **KG** | Neo4j DB on the lab box | URL + user + password | a gitignored `.env` in the repo (or shell/OS env) |
| **explorer** | Python code (library + MCP server) | template's `pyproject.toml` → `uv sync` | the local **venv** (Python imports it; `multiomics-kg-mcp` runs from it) |
| **research skills** | Markdown + hooks + `.mcp.json` (instructions for Claude) | **committed in the template repo** | read by Claude Code **from the cloned folder** — no venv, no install |
| **analyses** | the tester's research output | they create them | `analyses/` in their clone |

A marketplace-distributed plugin was evaluated and **rejected for the alpha** (§3) — skills-in-the-cloned-repo is simpler and matches the existing `uv sync` mental model. Revisit only if cross-tester skill-update management becomes painful.

This spec covers **only** the template/onboarding side. It does not re-specify the explorer package (its own spec) or KG hosting / the connection guide (KG plan).

## 2. Scope

**In scope:**
- A new **clean `multiomics_research_template`** repo: toolkit only (skills, recipes, hooks, scripts, `.mcp.json`, `pyproject.toml`, README, CLAUDE.md) + an **empty `analyses/` scaffold**. No analyses, no maintainer-specific paths.
- The skill-load mechanism that auto-loads on clone in the VS Code extension (§4 — skills move to `.claude/skills/`; one-time smoke test).
- `pyproject.toml`: depend on the explorer via a **git tag**, not a local editable path.
- `.mcp.json`: `uv run multiomics-kg-mcp` (from the synced venv), **no env block** — creds reach the explorer via a gitignored `.env` at the repo root.
- `.claude/settings.json`: portable hook paths (`${CLAUDE_PROJECT_DIR}`); author-specific allows stripped.
- A **preflight script** (`scripts/preflight.sh`): the DOA gate — `kg_release_info` (version contract) + a 2-line Python API smoke call.
- A **tester onboarding README** that owns clone→sync→creds→preflight→start-analysis and **links** to the KG connection guide for network/credential specifics.
- **Per-platform credential-setup instructions in the README** — *how* to make the MCP see `NEO4J_*` (not the values): a gitignored `.env` (primary) vs shell-profile / OS env, covering **Linux / Win11 / Win11+Remote-SSH**. (Values + network specifics stay in the KG guide.)
- Freezing the current repo as a historical artifact.

**Out of scope:**
- The explorer package itself (its spec) — the template only pins and consumes it.
- KG hosting, firewall, shared-credential distribution, the canonical connection guide — KG plan owns these.
- A marketplace/plugin distribution path — rejected for the alpha (§3).
- Migrating the maintainer's analyses anywhere — they stay in the frozen current repo.
- Per-analysis scaffolding automation — the existing step-protocol scaffolding stands.

## 3. Locked decisions

| Decision | Choice | Rationale |
|---|---|---|
| Delivery model | **Skills ride in the repo as files (not a marketplace plugin)** | Matches the tester's `uv sync` mental model. Skills are markdown read by Claude Code from the workspace — no install step, no venv. One clone of the fork + `uv sync` delivers everything. |
| Tester repo | **Fork** the new clean `multiomics_research_template` | A per-tester remote gives attribution for free — *who* = the fork owner (no `ALPHA_TESTER_ID`). The hook commits usage logs into the fork and the tester's analyses are captured too — both feed this repo's usage-logging/eval purpose. Template updates via standard `git pull upstream` (analyses + logs live in tester-only paths → clean merges). **Public template + public forks; analyses *and* logs both public** (decided 2026-06-08, open-science posture — §9). (Earlier "clone, not fork" is superseded: forking now pays for itself via attribution + log/analysis collection.) |
| Current repo | **Frozen as a historical artifact** | Retains the maintainer's analyses + full history. No longer the tester target; no migration needed. (Resolves the prior open question of where the maintainer's analyses live: they stay put.) |
| Marketplace/plugin | **Rejected for the alpha** | A git-marketplace plugin works (verified in the VS Code extension) and bundles skills+hooks+MCP, but it adds an install step + version-bump discipline the alpha doesn't need. Skills-in-repo is simpler. Recorded as a future option if multi-tester update management gets painful. |
| Explorer dependency | **`git = "…/multiomics_explorer.git", tag = "v0.1.0-alpha.1"`** in the template `pyproject.toml` | Reproducible `uv sync`; no sibling clone needed. The current editable `path = "../multiomics_explorer"` is a maintainer-only convenience that breaks any clone. **One** pin (pyproject); the MCP just runs whatever's in the venv. |
| MCP registration | Committed `.mcp.json`: `command: uv`, `args: ["run","multiomics-kg-mcp"]`, **no `env` block** | `uv run` from the workspace runs with CWD = repo root, so the explorer reads creds from `.env` / process env directly. The `${VAR}` env block is *omitted on purpose*: an unset `${VAR}` makes Claude Code fail to parse `.mcp.json` (explorer spec §6) — exactly the GUI/remote launch failure mode. Dropping it lets `.env` work uniformly (§6-credentials row, §10 #3). **(Supersedes the earlier `${VAR}`-env-block decision.)** |
| Credentials | **Primary: a gitignored `.env` in the repo root** (`NEO4J_URI` / `NEO4J_USERNAME` / `NEO4J_PASSWORD`, optional `NEO4J_DATABASE`). Alternative: shell-profile / OS env. | The explorer (pydantic-settings, `env_file: ".env"` — a relative path) reads `.env` from the process **CWD = repo root**, *not* from the venv/site-packages where the explorer code is installed — so the `.env` lives in the tester's cloned repo root. **Uniform across Linux / Win11 / Win11+Remote-SSH**, sidestepping VS Code's GUI/remote env-inheritance problem (the root cause of §10 #3). Mirrors KG plan §2.5 (shared `explorer` read login, out-of-band); `NEO4J_USER` accepted as alias. Refines the explorer-spec "env-vars = install model / `.env` = dev": testers **clone** the template → in-clone → `.env` is the right tool. |
| Toolkit packaging | **Not a buildable package** — drop the `[tool.hatch...wheel]` include block | The template is a *cloned workspace*, not a pip artifact. `pyproject.toml` exists only to declare deps for `uv sync`. |
| Portable settings | `.claude/settings.json` hooks use `${CLAUDE_PROJECT_DIR}`; author-specific allows pruned | Current absolute `/home/osnat/…` hook paths and analysis-specific `Bash(...)` allows break any clone / are meaningless to testers. |
| DOA gate | `scripts/preflight.sh`: `kg_release_info` (version contract) + a 2-line Python API smoke call | One green/red answer before opening a research chat (KG plan acceptance #6). |
| Doc ownership | KG `docs/kg_mcp_guide.md` is the **canonical connection guide**; the template README links to it | One source of truth per concern; avoids 3-repo drift on "how to connect." |

## 4. Skill-load mechanism (verify before locking)

The template must make the research skills **auto-load when a tester clones it and opens the folder in the VS Code extension — with no `/plugin install` step.**

**Decided 2026-06-08: move skills into `.claude/skills/`** (project-scoped, auto-load on workspace trust — the most reliable "clone and it just works" path) and drop the `.claude-plugin/plugin.json` manifest + root `skills/` layout. Recipes likewise move under `.claude/skills/`.

**Still smoke-test once** (the "does the skill actually load?" check): scaffold the template, open in the extension, trust the workspace, confirm the research-methodology skill appears in a new chat. This needs no KG/explorer — skill loading is independent — so it can be verified immediately, before the gated end-to-end dogfood.

## 5. Clean template contents

**Ships:** `.claude/skills/` (research-methodology + recipes), `hooks/`, `scripts/` (incl. `preflight.sh`), `evals/` (optional for testers — decide), `.claude/settings.json` (portable), `.mcp.json` (tester form), `pyproject.toml` (git-tag explorer, no wheel build), `CLAUDE.md`, tester `README.md`, empty `analyses/` (`.gitkeep` + `analyses/README.md` describing the per-analysis structure the step-protocol expects).

**Concrete cleanups vs the current repo (verified 2026-06-06):**
- [pyproject.toml](../../../pyproject.toml) — editable `../multiomics_explorer` source → git tag; remove the `[tool.hatch.build.targets.wheel]` block.
- [.mcp.json](../../../.mcp.json) — replace `${MULTIOMICS_EXPLORER_DIR}` + `--directory` with `uv run multiomics-kg-mcp` (**no env block** — creds via `.env`).
- Add a committed **`.env.example`** (the three `NEO4J_*` keys, no values) and ensure `.env` is gitignored — the tester copies `.env.example` → `.env` and fills in operator-provided values.
- [.claude/settings.json](../../../.claude/settings.json) — hook `command` → `${CLAUDE_PROJECT_DIR}/hooks/log-mcp-usage.sh`; drop the three analysis-script `Bash(...)` allows; keep the MCP allows + usage-logging hooks.
- **Org/owner** — canonical owner is **`wosnat`** (decided 2026-06-08); the new repo is `github.com/wosnat/multiomics_research_template`. (The current [plugin.json](../../../.claude-plugin/plugin.json) `github.com/osnat/…` was wrong; the manifest is being dropped anyway per §4.)
- [README.md](../../../README.md) — currently "18 tools / 15 organisms / 76 experiments"; live KG is 40 tools / 37 organisms / 197 experiments. Write a fresh tester README (§7); don't restate counts that drift — link the KG guide.

## 6. Preflight ("DOA") gate

`scripts/preflight.sh` runs, in order:
1. **Version contract** — `kg_release_info`: assert the synced explorer satisfies the live KG's `Schema_info.mcp_min_version`; surface KG identity + ok/warn/unknown.
2. **Python API smoke** — a 2-line `GraphConnection()` + `gene_overview([...])` proving import + Neo4j auth + a real round-trip.

Output: one green/red summary with remediation hints for the three common failures — env var unset, off-subnet, version mismatch. Red = don't proceed to an analysis chat; report it.

## 7. Tester onboarding README (outline)

Owns clone→sync→creds→preflight→start-analysis; links to the KG guide for connection specifics.

```
# multiomics_research_template (alpha)

[One paragraph: clone this, point it at the lab KG, run analyses in Claude Code.
Carries the research skills + MCP wiring; uv sync pulls the explorer tools.]

## Prerequisites
- git, Python 3.11+, uv, Claude Code, VS Code
- Lab-subnet access to the KG box → see <KG connection guide link>

## 1. Fork, clone & install
   fork <template url> → git clone <your fork> && cd multiomics_research_template
   git remote add upstream <template url>     # for pulling updates later
   uv sync                       # pulls the explorer (pinned tag) + analysis deps

## 2. Set KG credentials (values from the operator → <KG guide link>; never commit them)
   Recommended — all platforms: a gitignored .env in the repo root
       cp .env.example .env   # then fill in:
       NEO4J_URI=bolt://132.75.249.47:17687
       NEO4J_USERNAME=explorer
       NEO4J_PASSWORD=…
   The MCP reads .env automatically (it runs from the repo root). Identical on Linux, Win11, Remote-SSH.

   Alternative — shell / OS environment (if you'd rather not keep a creds file):
     • Linux:       export NEO4J_* in ~/.bashrc, then launch VS Code from that terminal (`code .`)
                    — a GUI-launched VS Code won't inherit them.
     • Win11:       set them as User Environment Variables (PowerShell:
                    [Environment]::SetEnvironmentVariable("NEO4J_URI","…","User")), then restart VS Code.
     • Win11+Remote-SSH: set them on the REMOTE host (the MCP runs there), not on Windows
                    — or just use the .env method, which is simplest here.

## 3. Open in Claude Code (VS Code) & trust the workspace
   skills auto-load; the MCP registers (.mcp.json). Confirm in /mcp.

## 4. Preflight (DOA)
   ./scripts/preflight.sh         # version contract + API smoke → green/red

## 5. Start an analysis
   new chat; the research-methodology skill loads; work lives in analyses/ in YOUR clone.

## Updating
   git pull upstream main && uv sync     (preflight tells you when you're behind)
```

**Do NOT duplicate** from the KG guide: lab Bolt URI rationale, firewall/subnet checks, `Test-NetConnection`/`nc`, shared-credential handling. Link instead.

## 8. Verification & acceptance

1. **Clean-clone dogfood.** From a fresh dir (or a second machine), fork→clone→`uv sync`→set env→open in the extension→trust→`./scripts/preflight.sh` green, **with no sibling explorer clone present**. Proves the git-tag dependency path.
2. **Skill auto-loads on clone** in the VS Code extension without `/plugin install` (§4 mechanism confirmed).
3. **No maintainer leakage.** `grep -rn '/home/osnat'` empty; no `analyses/` content; no `.env`. `git ls-files analyses/` shows only scaffold.
4. **MCP reads creds from `.env`** at the repo root and connects — verified on all 3 platforms (MCP CWD = repo root; on Remote-SSH the `.env` is in the remote workspace). The shell/OS-env alternative also works when the launching process carries the vars.
5. **Preflight catches the three common failures** — env unset, off-subnet, version mismatch — each with a distinct message.
6. **First analysis runs** — skill loads, one MCP tool returns real data, an analysis dir scaffolds under the clone.
7. **Doc non-duplication** — README defers connection specifics to the KG guide via links.

## 9. Versioning, tracking & updates

The clone model still needs explicit versioning — "just `git pull`" tells you neither who is on what, nor *when* a tester should pull. Methodology:

**Template versioning.** Semver git tags on the template repo (`v0.1.0-alpha.N`), mirroring `/release-kg` discipline. A `CHANGELOG.md` (Keep a Changelog) accumulates notable changes under `[Unreleased]`, cut to `[version]` at release; a GitHub Release per tag. The tag (via `git describe`) or a `VERSION` file is what preflight reads.

**The version triple.** Each tester runs three independently-versioned things; surface all three:
- **template** `vX.Y.Z` — skills + wiring (git tag)
- **explorer** `vA.B.C` — pinned in `pyproject.toml`, in the venv
- **KG** — `Schema_info.version`, live from the box (`kg_release_info`)

`scripts/preflight.sh` prints the triple in one line so any tester reports their exact stack. Template→explorer is deterministic (read the pin); KG comes from `kg_release_info`.

**When to pull — preflight is the trip-wire.** Preflight does a lightweight `git fetch --tags` and compares the local tag to the latest release; if behind, it prints the exact command + CHANGELOG pointer:
```
⚠ template v0.1.0-alpha.1 — latest is v0.1.0-alpha.3
   update:  git pull upstream main && uv sync   (changes: CHANGELOG.md)
```
This fires at the point of use, not in a channel testers might miss. Backed by the CHANGELOG, watchable GitHub Releases, and an out-of-band announce for big/breaking cuts.

**How to pull — always two steps as one:** `git pull upstream main && uv sync`. The pull brings new skills + the bumped explorer pin; `uv sync` installs the new explorer. Preflight verifies the synced explorer matches the pin (catches "pulled but forgot to sync").

**Who is using what — attribution via git (no id needed).** All testers share the one `explorer` DB login (KG plan §2.5), so the DB cannot distinguish them. Instead of a client-side id, use the fork: the usage-logging hook writes JSONL **into the fork** (e.g. `usage/*.jsonl`, un-ignored) with the **version triple stamped per line**; testers commit it alongside their per-step analysis commits, so logs ride along on every push. *Who* = the fork owner; *what version* = stamped in the line; *their analyses* are captured too. The maintainer aggregates by pulling forks (publicly listable via the GitHub forks API). Requires: change the hook to write into the workspace (today it writes `~/.claude/logs/multiomics-kg-usage.jsonl`, gitignored), un-ignore the in-repo log path, and document the commit/push cadence (piggyback on step commits).

**Visibility — public, all of it (decided 2026-06-08).** Public template, public forks, and **both analyses and usage logs are pushed publicly** — open-science posture, chosen with the trade-off explicit. Consequence acknowledged: testers' in-progress, unpublished analyses (including findings targeted for publication) are world-visible and indexable pre-publication; this is a one-way door. Onboarding should state this plainly so testers know their forks are public before they start. If a specific analysis ever needs to stay private, that tester gitignores its dir locally — the default is public.

**Release flow.** A short `/release-template` checklist mirroring `/release-kg`: bump explorer pin if needed → bump template VERSION → update CHANGELOG → commit → tag `v…` → push `--follow-tags` → GitHub Release → announce.

## 10. Open questions / dependencies

1. ~~**Skill-load mechanism**~~ — **RESOLVED 2026-06-08: move skills to `.claude/skills/`** (§4). Remaining: the one-time smoke test that it auto-loads on a bare clone in the VS Code extension (no deps — do first).
2. ~~**GitHub org/owner**~~ — **RESOLVED 2026-06-08: `wosnat`** (§5). New repo: `github.com/wosnat/multiomics_research_template`.
3. **Credential delivery — verify on 3 target platforms** (Linux / Win11 / Win11+Remote-SSH, the supported environments). Largely *resolved by the `.env`-primary decision* (§3, §6): a `.env` at the repo root is read by the MCP uniformly, dodging the VS Code env-inheritance differences. **Remaining verification** = confirm the MCP process actually launches with **CWD = repo root** in each environment so `.env` is found — especially **Remote-SSH** (the MCP runs on the remote host; the `.env` must be in the *remote* workspace). The shell/OS-env *alternative* still carries the per-platform inheritance caveats documented in the README (§7 step 2); the `.env` path is recommended precisely to avoid them. Fold the CWD check into the §4 smoke test.
4. **`evals/` in the template?** — testers probably don't need the eval suite; decide whether it ships or stays maintainer-only.
5. **Cross-repo sequencing** — end-to-end dogfood (§8.1) gates on (a) the explorer `v0.1.0-alpha.1` tag existing (explorer spec not yet implemented) and (b) the KG deployed on the lab box (KG plan items #7–#14). The template + config + preflight can be built now.
6. **Win11 ergonomics** — confirm `uv sync`, env-var setting, and the workspace-trust + `/mcp` flow on a clean Win11 box before inviting testers (mirrors KG plan §6.4(4)).
7. **Where the maintainer's *future* analyses go** — presumably a fresh clone of the template (maintainer becomes a normal tester). Not a blocker; the frozen repo holds the past.
8. **Staleness-check mechanism** (§9) — offline `git fetch --tags` + tag compare (no auth, works on lab subnet) vs a GitHub API call (needs network/rate-limit) vs none (manual/announce-only). Recommend the offline tag compare.
9. ~~**Fork visibility / privacy**~~ — **RESOLVED 2026-06-08: fully public** (template, forks, analyses, logs). Open-science posture; one-way-door consequence acknowledged (§9). Onboarding must state forks are public up front. Remaining: tester consent/awareness wording in the README.
