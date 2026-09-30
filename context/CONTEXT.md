# BioSeqFlow — Project Context

## What this is

A single, self-contained desktop application that wraps the industry-standard NGS toolchain
(FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts, pyDESeq2, SnpEff) for RNA-seq
differential expression and DNA variant calling behind a clean GUI. One repo, one product: there
is no separate Nextflow pipeline or standalone Python library to build and keep in sync first —
the app directly orchestrates each tool as a managed subprocess.

This project makes those capabilities accessible to scientists who are not comfortable with the
command line or managing bioinformatics tool installations. The goal is a
**production-quality, open-source, cross-platform desktop application** — not a demo or internal
tool. Every architectural and code decision should reflect that standard.

---

## The developer

**Oumar Ndiaye** — bioinformatics engineer, author of CellMetPro.

- Comfortable with Python: built CellMetPro end-to-end (COMPASS algorithm, FBA, scRNA-seq pipeline,
  CLI, test suite, CI/CD, PyPI)
- Strong JavaScript/TypeScript: React, Vue, Node.js; has shipped production JS projects
- Already built CellMetPro UI with the same architecture — this is the second OmicsPilot desktop
  application, applying established patterns to a new domain
- Not a fan of Java; avoid unless there is a genuinely compelling reason
- This project deepens expertise in:
  - The NGS toolchain itself (FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts,
    pyDESeq2, SnpEff) — learned hands-on while building each subprocess wrapper, not as a
    separate prerequisite project
  - Managing long-running external processes from a web server (hours-scale, not seconds)
  - WebSocket progress streaming from tool-specific log formats
  - Reference genome management (large file downloads, indexing, persistence)
  - Cross-platform tool availability and environment validation
  - Docker images bundling multiple bioinformatics tools
- Quality bar: same as CellMetPro UI — no shortcuts, no "clean up later." Code is typed, tested,
  documented, and CI-gated.

---

## Architecture decision

### Two-piece design

```
┌─────────────────────────────────────────────────────────────┐
│  Piece 1: bioseqflow-server  (Python, pip-installable)      │
│  FastAPI application that wraps NGS tool subprocesses        │
│  Runs wherever the user wants: local machine, HPC, cloud VM  │
└────────────────────────┬────────────────────────────────────┘
                         │ REST + WebSocket
┌────────────────────────▼────────────────────────────────────┐
│  Piece 2: BioSeqFlow Desktop  (Electron + React)            │
│  The React frontend works standalone in any browser too      │
│  Electron is a thin shell: process management + native UX   │
└─────────────────────────────────────────────────────────────┘
```

### Two operating modes (same codebase, config-switchable)

**Mode 1 — Local (default for most users)**
- User downloads and double-clicks the desktop app
- Electron checks for required NGS tools at startup; shows a setup guide if any are missing
- Electron spawns the Python server as a child process and manages its lifecycle
- User never sees a terminal or a Nextflow command

**Mode 2 — Remote (for power users and labs)**
- User (or their sysadmin) runs `pip install bioseqflow-server && bioseqflow-server start` on a
  remote machine (HPC cluster, cloud VM, lab workstation) where all NGS tools are already installed
- User opens the desktop app, goes to Settings, enters the server URL and an auth token
- The React frontend connects to the remote server; all computation happens there
- The desktop app is purely a UI client in this mode — particularly useful for HPC where GATK and
  STAR benefit from many cores and large memory

### Why this architecture

- Open source and zero hosting costs — users bring their own compute
- Non-technical users get a double-click install experience
- Power users get full HPC/cloud flexibility
- The React frontend is backend-agnostic: it only talks to a URL
- The API layer makes BioSeqFlow accessible to future R clients or third-party tools
- Consistent with CellMetPro UI: same monorepo layout, same tech stack, same quality standards

### Key difference from CellMetPro UI

CellMetPro UI wraps a Python library — computation happens inside the server process. BioSeqFlow
wraps external binaries (STAR, GATK, etc.) as managed subprocesses. This changes three things:

| Aspect | CellMetPro UI | BioSeqFlow |
|---|---|---|
| Computation | Python library call | External subprocess calls |
| Progress | Task events | Parsed stdout/stderr streams |
| Runtime | Minutes | Hours (STAR ~20min, GATK ~2hr for human) |
| Tool deps | Python packages only | External binaries + reference genomes |
| File sizes | h5ad files (<2GB typical) | FASTQ pairs (can be >50GB each) |

---

## Monorepo structure

```
bioseqflow/                            # GitHub repo: omicspilot/bioseqflow
├── packages/
│   ├── server/                        # pip-installable Python package
│   │   ├── bioseqflow_server/
│   │   │   ├── main.py                # FastAPI app entrypoint
│   │   │   ├── routers/               # one router per domain
│   │   │   │   ├── projects.py        # project + sample management
│   │   │   │   ├── references.py      # reference genome download + indexing
│   │   │   │   ├── qc.py              # FastQC + fastp
│   │   │   │   ├── rnaseq.py          # STAR + featureCounts + pyDESeq2
│   │   │   │   ├── variants.py        # BWA-MEM2 + samtools + GATK + SnpEff
│   │   │   │   ├── results.py         # serve result files + export
│   │   │   │   └── io.py              # file upload + validation
│   │   │   ├── process_manager.py     # subprocess wrapper + log streaming
│   │   │   ├── tool_check.py          # verify external tool availability
│   │   │   ├── jobs.py                # background job management
│   │   │   ├── auth.py                # token auth for remote mode
│   │   │   ├── ws.py                  # WebSocket progress streaming
│   │   │   └── parsers/               # tool-specific log parsers
│   │   │       ├── fastqc.py          # parse FastQC JSON output
│   │   │       ├── fastp.py           # parse fastp JSON output
│   │   │       ├── star.py            # parse STAR Log.final.out
│   │   │       ├── samtools.py        # parse flagstat / markdup stats
│   │   │       ├── gatk.py            # parse HaplotypeCaller progress
│   │   │       └── snpeff.py          # parse SnpEff summary stats
│   │   ├── pyproject.toml
│   │   └── tests/
│   │
│   └── desktop/                       # Electron + React
│       ├── src/                       # React app (also works in browser)
│       │   ├── app/                   # routing, layout
│       │   ├── features/              # one folder per feature domain
│       │   │   ├── setup/             # tool availability check + setup guide
│       │   │   ├── projects/          # project + sample management
│       │   │   ├── references/        # reference genome manager
│       │   │   ├── pipeline/          # pipeline config + launch
│       │   │   ├── progress/          # live run monitoring
│       │   │   ├── qc/                # QC results dashboard
│       │   │   ├── rnaseq/            # DE results, volcano, PCA
│       │   │   ├── variants/          # variant table, impact charts
│       │   │   ├── export/            # download results
│       │   │   └── settings/          # server URL, auth, resources
│       │   ├── api/                   # typed API client (auto-generated from OpenAPI)
│       │   ├── store/                 # Zustand global state
│       │   └── components/            # shared UI components
│       ├── electron/                  # Electron shell
│       │   ├── main.ts                # main process: window, server lifecycle
│       │   ├── preload.ts             # context bridge
│       │   └── server-manager.ts      # spawn / kill Python server
│       ├── package.json
│       └── vite.config.ts
│
├── docker/
│   └── Dockerfile                     # bioconda base + all NGS tools pre-installed
│   └── docker-compose.yml             # server + reference genome volume
├── .github/
│   └── workflows/
│       ├── ci.yml                     # lint + test on every PR
│       └── release.yml                # build installers + publish on tag
└── README.md
```

---

## Technology stack

### Backend (`packages/server`)

| Concern | Choice | Reason |
|---|---|---|
| Framework | **FastAPI** | Async, WebSocket native, auto OpenAPI docs |
| Database | **SQLAlchemy 2.x async** + **Alembic** + **aiosqlite** | Typed models, tracked migrations, portable to Postgres |
| Task queue | **FastAPI BackgroundTasks** (start), migrate to **ARQ** if needed | Avoids Celery overhead early on |
| Process management | **asyncio.create_subprocess_exec** | Non-blocking subprocess with streamed stdout/stderr |
| Progress streaming | **WebSocket** per job ID | Real-time progress bars from long-running tools |
| Auth (remote mode) | **Bearer token** (static, user-generated) | Simple, no OAuth overhead for v1 |
| Testing | **pytest** + **httpx** (async test client) | Consistent with CellMetPro |
| Linting | **ruff** + **mypy** | ruff handles both linting and formatting |

### Frontend (`packages/desktop`)

| Concern | Choice | Reason |
|---|---|---|
| Framework | **React 18 + TypeScript** | Same as CellMetPro UI, known stack |
| Build | **Vite** | Fast HMR, small bundles |
| State | **Zustand** | Minimal, TypeScript-first, no boilerplate |
| Server state | **TanStack Query (React Query)** | Caching, background refetch, loading states |
| UI components | **shadcn/ui** (Radix + Tailwind) | Same as CellMetPro UI, shared design tokens |
| Routing | **React Router v6** | File-based routes inside Vite |
| Charts | **Plotly.js** | Scientific plots: volcano, PCA, heatmap, variant distributions |
| API client | **Auto-generated from OpenAPI spec** (openapi-ts) | Types always in sync with backend |
| Forms | **React Hook Form + Zod** | Type-safe validation for pipeline config |
| Testing | **Vitest** + **Testing Library** | Fast, Vite-native |

### Electron shell

| Concern | Choice |
|---|---|
| Electron version | Latest stable |
| Bundler | **electron-builder** |
| IPC | Context bridge (no `nodeIntegration`) |
| Server management | `child_process.spawn` with stdout/stderr piped to log panel |
| Auto-updater | **electron-updater** (GitHub Releases) |
| Distribution | `.dmg` (macOS), `.exe` NSIS installer (Windows), `.AppImage` (Linux) |

---

## Key product decisions

### Project model
- Every user session is organised into **projects** — named workspaces that group samples, pipeline
  runs, and results
- Projects persist across server restarts via SQLAlchemy + SQLite
- A user can have multiple concurrent projects and switch between them
- Projects are the top-level API resource: all samples, runs, and results are project-scoped
- Data directory: `~/.bioseqflow/projects/{project_id}/` holds files on disk; metadata in SQLite

### Sample model
- Samples are FASTQ file pairs (R1 + R2 for paired-end, single file for single-end)
- Metadata per sample: name, condition label (used for DE contrasts), sequencing type (RNA-seq / WGS / WES)
- Samples belong to a project; a project can have multiple samples processed in parallel
- Large file uploads are streamed; files stored under the project directory

### Reference genome model
- Reference genomes are shared across projects (stored in `~/.bioseqflow/references/`)
- Each reference has: FASTA, GTF annotation, STAR index, BWA-MEM2 index, SnpEff database
- Indexing is a tracked background job with WebSocket progress (STAR index for human: ~30 min)
- Common genomes (human GRCh38, mouse GRCm39) available for one-click download
- Custom genomes can be registered by providing a FASTA + GTF path

### Pipeline run model
- Every analysis is a **pipeline run** with a UUID, scoped to a project
- Runs carry a `pipeline_type` field: `rnaseq` or `variantcalling`
- A run consists of ordered **steps**: QC → Trimming → Alignment → Quantification/Calling → Analysis
- Run states: `queued → running → complete | failed | cancelled`
- Step states tracked independently: a failed GATK step does not erase completed alignment
- Progress streamed via WebSocket: `{ run_id, step, tool, progress, message }`
- Users can resume a failed run from the last completed step

### Job model
- Each pipeline step is a **job** wrapping one external tool process
- Jobs are managed by the process manager: spawn, stream stdout/stderr, parse progress, handle exit codes
- Each tool has a dedicated log parser that extracts structured progress from raw output
- Tool exit codes are mapped to human-readable error messages in the UI
- Logs are persisted in SQLite and viewable in the UI even after the run completes

### Normalised progress event schema
All tools emit progress through a single WebSocket event schema:
```json
{
  "run_id": "uuid",
  "step": "alignment",
  "tool": "STAR",
  "progress": 0.42,
  "message": "Mapping reads to reference...",
  "timestamp": "2026-07-17T10:23:11Z"
}
```
Each log parser (STAR, GATK, fastp, etc.) is responsible for translating raw tool output into
this schema. The frontend is tool-agnostic — it only consumes normalised events.

### Version management
- The server reports tool versions at startup (`GET /tools/versions`)
- The desktop app displays tool versions in the status bar
- Remote mode: user connects to a server with pre-installed tools; tool versions shown in UI
- This enables reproducibility: a run record stores the tool versions used

### Auth model (remote mode only)
- Server started with `bioseqflow-server start --token <token>`
- Client sends `Authorization: Bearer <token>` on every request
- Local mode: no auth (localhost only, token ignored)

---

## Quality standards (non-negotiable)

- **TypeScript strict mode** throughout the frontend — no `any`
- **mypy** on the backend with `--strict` or `--ignore-missing-imports` at minimum
- **100% of API routes have tests** (pytest + httpx)
- **CI blocks merges** on failing tests, lint, or type errors
- **No `console.log` or `print` in production code** — structured logging only
- **Semantic versioning** — server and desktop versions are kept in sync
- **Changelog** maintained for every release (Keep a Changelog format)
- **Accessibility**: all interactive UI elements keyboard-navigable, ARIA labels on custom components

---

## Relationship to OmicsPilot

- `bioseqflow-server` calls NGS tools as subprocesses — it does not reimplement the analysis logic
  of industry-standard tools (FastQC, STAR, GATK, etc. remain the source of truth for their own
  algorithms)
- The desktop app is version-agnostic as long as the server API contract is maintained
- Breaking API changes follow semver (major bump)
- BioSeqFlow and CellMetPro share the same OmicsPilot brand, design tokens, and monorepo
  conventions — they are sibling products, not competing tools

---

## Open source

- License: **MIT** (consistent with CellMetPro)
- GitHub: `omicspilot/bioseqflow`
- Issues and PRs welcome from day one
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and issue templates from the start
