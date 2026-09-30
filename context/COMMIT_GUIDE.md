# CellMetPro UI — Commit Message Guide

Format: `<type>(<scope>): <short description>`

- Description is lowercase, imperative mood ("add" not "added" / "adds")
- No period at the end
- Max 72 characters on the first line
- Breaking changes: append `!` after type/scope, or add `BREAKING CHANGE:` in footer

---

## Types

| Type       | When to use                                              | Version impact |
|------------|----------------------------------------------------------|----------------|
| `feat`     | A new feature visible to the user or API consumer       | minor bump     |
| `fix`      | A bug fix                                                | patch bump     |
| `perf`     | Performance improvement with no behavior change          | patch bump     |
| `refactor` | Code restructure — no new feature, no bug fix            | none           |
| `test`     | Adding or correcting tests                               | none           |
| `docs`     | Documentation only (README, CONTEXT, comments)           | none           |
| `chore`    | Maintenance — dependency updates, config, tooling        | none           |
| `build`    | Changes to the build system (pyproject.toml, vite, etc.) | none           |
| `ci`       | CI/CD workflow changes (.github/workflows)               | none           |
| `style`    | Formatting only — no logic change (whitespace, quotes)   | none           |
| `revert`   | Reverts a previous commit                                | none           |

A `!` after the type (or `BREAKING CHANGE:` in the footer) triggers a **major bump**:
```
feat(server)!: change job status response shape
```

---

## Scopes

| Scope       | What it covers                                              |
|-------------|-------------------------------------------------------------|
| `server`    | FastAPI server code (`packages/server/cellmetpro_server/`)  |
| `desktop`   | React frontend source (`packages/desktop/src/`)             |
| `electron`  | Electron shell (`packages/desktop/electron/`)               |
| `tooling`   | Root-level config: ESLint, ruff, mypy, pre-commit, pnpm     |
| `ci`        | GitHub Actions workflows (`.github/workflows/`)             |
| `docker`    | Dockerfile and docker-compose (`docker/`)                   |
| `deps`      | Dependency updates (pyproject.toml, package.json)           |
| `tests`     | Test files or test infrastructure                           |
| `context`   | Project documentation in `context/` (ROADMAP, CONTEXT, etc.)|
| `release`   | Version bumps, changelogs, tagging                          |

Scope is optional but strongly encouraged. Omit it only for truly cross-cutting changes.

---

## Examples

```
feat(server): add health check endpoint
fix(server): correct prefix missing leading slash in system router
chore(tooling): configure pre-commit hooks
chore(deps): add pydantic-settings to server dependencies
test(server): add tests for health and version endpoints
refactor(server): extract settings into config module
docs(context): add commit guide
ci: add lint and test jobs to ci.yml
feat(desktop): add file upload drag-and-drop zone
fix(electron): prevent zombie process on unexpected server crash
feat(server)!: change job progress WebSocket message shape
```

---

## Multi-line commits (when the why needs explaining)

```
fix(server): handle port-in-use error on server startup

Previously the server crashed with an unhandled OSError when the
default port 8000 was already occupied. Now it logs a clear message
and exits with code 1 so Electron can surface a user-facing error.
```

The body answers **why**, not what. The diff already shows what changed.
