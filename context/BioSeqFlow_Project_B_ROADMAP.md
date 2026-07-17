# BioSeqFlow — Implementation Roadmap

This file is the single source of truth for the full implementation plan.
Each step is tackled linearly. Check off items as they are completed.

> **Mission**: Make NGS analysis accessible to any researcher, regardless of command-line experience. BioSeqFlow wraps the full NGS toolchain (FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts, pyDESeq2, SnpEff) behind a clean desktop GUI — the same architecture as CellMetPro.
>
> **Prerequisite**: `bioseqflow-core` (Project A) must be complete. You need to understand what every tool does, what it outputs, and what can go wrong before wrapping it in a GUI.

---

## Phase 0 — Monorepo & Tooling Foundation

> Mirror the CellMetPro monorepo structure exactly for consistency across OmicsPilot.

- [ ] **0.1** — Monorepo structure with `pnpm` workspaces (`packages/server`, `packages/desktop`)
- [ ] **0.2** — Python environment: `uv`, `pyproject.toml`, `ruff` / `mypy`
- [ ] **0.3** — Node/TS environment: `pnpm`, `tsconfig`, `eslint` / `prettier`
- [ ] **0.4** — Git hygiene: `.gitignore`, `commitlint`, conventional commits
- [ ] **0.5** — Pre-commit hooks (Python + JS in one repo)
- [ ] **0.6** — Shared OmicsPilot design tokens: import colour palette, typography, component styles from CellMetPro

---

## Phase 1 — FastAPI Server (`packages/server`)

> The server wraps NGS tools as managed subprocesses, streams their output over WebSocket, and exposes results via REST. This is the core difference from CellMetPro: instead of calling a Python library, you are orchestrating long-running external processes.

- [ ] **1.1** — FastAPI app skeleton, project layout, Uvicorn
- [ ] **1.2** — Health check + `GET /version` endpoint
- [ ] **1.3** — Tool availability check: verify FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts, SnpEff are on PATH; return JSON status per tool
- [ ] **1.4** — File I/O router: multipart FASTQ upload, reference genome registration, temp directory management
- [ ] **1.5** — Process manager: subprocess wrapper with stdout/stderr streaming, PID tracking, graceful kill
- [ ] **1.6** — Background jobs system: job UUID → state model (QUEUED / RUNNING / COMPLETE / FAILED / CANCELLED), job queue with concurrency limit
- [ ] **1.7** — WebSocket progress streaming: parse tool-specific log formats → normalised progress events → `ws.py`

  > Each tool emits progress differently. Define a normalised event schema:
  > ```json
  > { "job_id": "...", "tool": "STAR", "stage": "alignment", "progress": 0.42, "message": "Mapping reads..." }
  > ```

- [ ] **1.8** — Database layer: SQLAlchemy async engine, Alembic setup
  - Models: `Project`, `Sample`, `ReferenceGenome`, `PipelineRun`, `RunStep`, `ResultFile`
  - First migration
- [ ] **1.9** — Project management router: CRUD for projects and samples, project-scoped pipeline runs
- [ ] **1.10** — QC router: wrap FastQC + fastp, parse JSON outputs, return structured QC metrics
- [ ] **1.11** — RNA-seq router:
  - STAR alignment (stream progress via WebSocket)
  - featureCounts quantification
  - pyDESeq2 differential expression (condition metadata from project)
  - Return: count matrix, DE results table, PCA coordinates, volcano plot data
- [ ] **1.12** — Variant calling router:
  - BWA-MEM2 alignment (stream progress)
  - samtools sort + markdup + index
  - GATK HaplotypeCaller (stream progress)
  - SnpEff annotation
  - Return: annotated VCF summary, variant table, impact distribution
- [ ] **1.13** — Reference genome router: download, index (STAR + BWA-MEM2), register in DB
- [ ] **1.14** — Results router: serve result files, pre-computed plot data, export to CSV/VCF
- [ ] **1.15** — Auth middleware: Bearer token (remote mode only)
- [ ] **1.16** — Error handling: tool-specific exit codes → human-readable messages; structured logging; CORS
- [ ] **1.17** — pytest + httpx test suite (100% route coverage target)
- [ ] **1.18** — OpenAPI spec export (feeds the frontend type generator)

---

## Phase 2 — React Frontend (`packages/desktop/src`)

> Same stack as CellMetPro: React + TypeScript + Vite + Tailwind + shadcn/ui + Zustand + TanStack Query + Plotly.js.

### 2.0 — Scaffold & Shared Setup
- [ ] **2.1** — Vite + React + TypeScript scaffold
- [ ] **2.2** — Tailwind CSS + shadcn/ui setup, OmicsPilot design tokens applied
- [ ] **2.3** — React Router v6: root layout, route definitions
- [ ] **2.4** — Zustand stores: connection state, active project state, active run state, tool availability state
- [ ] **2.5** — Auto-generated API client from OpenAPI spec (`openapi-ts`)
- [ ] **2.6** — TanStack Query: data fetching layer, query/mutation patterns

### 2.1 — Tool Status & Setup
- [ ] **2.7** — Tool availability panel: show green/red status per tool on first launch
- [ ] **2.8** — Setup guide: instructions for installing missing tools via conda/mamba
- [ ] **2.9** — Reference genome manager: browse registered genomes, download new ones, track indexing progress

### 2.2 — Project & Sample Management
- [ ] **2.10** — Project list view: create, rename, archive projects
- [ ] **2.11** — Sample manager: upload FASTQ files (drag-and-drop, paired-end support), sample metadata form
- [ ] **2.12** — Sample sheet preview: validate pairs, show file sizes, FASTQ format check

### 2.3 — Pipeline Configuration & Launch
- [ ] **2.13** — Pipeline selector: choose RNA-seq or Variant Calling track (or both)
- [ ] **2.14** — RNA-seq config form: reference genome, annotation GTF, strandness, condition assignment per sample, DE contrast definition
- [ ] **2.15** — Variant calling config form: reference genome, calling mode (germline/somatic), target BED (optional), SnpEff database
- [ ] **2.16** — Resource config: threads per tool, memory limits, output directory
- [ ] **2.17** — Run launch: validate config → submit to server → redirect to progress view

### 2.4 — Progress Monitoring
- [ ] **2.18** — Run progress view: WebSocket consumer, per-tool progress bars, stage labels, elapsed time
- [ ] **2.19** — Live log panel: stream raw tool output, collapsible per tool
- [ ] **2.20** — Run controls: pause queue, cancel run, resume from last completed step
- [ ] **2.21** — Run history: list past runs per project, status badges, link to results

### 2.5 — QC Results
- [ ] **2.22** — QC dashboard: per-sample quality cards (pass/warn/fail per FastQC module)
- [ ] **2.23** — Per-base quality plot: interactive Plotly line chart, before/after trimming overlay
- [ ] **2.24** — GC content distribution: per-sample, flag outliers
- [ ] **2.25** — Adapter content panel: show detected adapters, trimming summary
- [ ] **2.26** — MultiQC-style aggregate table: sortable, filterable, flag failing samples

### 2.6 — RNA-seq Results
- [ ] **2.27** — Count matrix preview: gene × sample heatmap (top 50 by variance), log-normalised
- [ ] **2.28** — PCA plot: samples coloured by condition, interactive tooltips
- [ ] **2.29** — Volcano plot: log2FC vs -log10(padj), interactive hover (gene name, values), significance thresholds configurable
- [ ] **2.30** — DE results table: sortable by padj/log2FC, filterable, gene name links to Ensembl
- [ ] **2.31** — Alignment statistics panel: % uniquely mapped, % multi-mapped, % unmapped per sample

### 2.7 — Variant Calling Results
- [ ] **2.32** — Alignment statistics: % duplicates, mean depth, % mapped per sample
- [ ] **2.33** — Variant summary: SNP count, indel count, Ti/Tv ratio, % PASS
- [ ] **2.34** — Impact distribution chart: HIGH / MODERATE / LOW / MODIFIER counts (stacked bar)
- [ ] **2.35** — Variant table: sortable/filterable by gene, impact, consequence, AF, DP; ClinVar badge if known pathogenic
- [ ] **2.36** — Variant detail panel: click row → show full annotation, gene context, external links (ClinVar, gnomAD, OMIM)

### 2.8 — Export & Settings
- [ ] **2.37** — Export panel: download DE results CSV, annotated VCF, QC report HTML, plots as SVG/PNG
- [ ] **2.38** — Settings page: server URL, auth token, mode switch (local/remote), default resource config
- [ ] **2.39** — Accessibility audit: keyboard navigation, ARIA labels, reduced motion
- [ ] **2.40** — Vitest + Testing Library test suite

---

## Phase 3 — Electron Shell (`packages/desktop/electron`)

> Identical pattern to CellMetPro. The server manager is more complex here because it must also verify external tool availability at startup.

- [ ] **3.1** — Electron main process: `BrowserWindow`, security hardening, CSP headers
- [ ] **3.2** — Context bridge / preload: IPC surface design
- [ ] **3.3** — Server manager: spawn Python FastAPI server, pipe logs, kill on quit
- [ ] **3.4** — Tool availability check on startup: run `fastqc --version`, `STAR --version` etc., show setup dialog if any missing
- [ ] **3.5** — Local mode vs remote mode switching
- [ ] **3.6** — Native dialogs: file open (FASTQ, BED, GTF), folder selection, system notifications on run complete
- [ ] **3.7** — Auto-updater (`electron-updater` + GitHub Releases)
- [ ] **3.8** — electron-builder config: `.dmg` (macOS), `.exe` NSIS (Windows), `.AppImage` (Linux)

---

## Phase 4 — Integration & End-to-End

> Wire all three pieces together and harden the seams. Test with real public datasets.

- [ ] **4.1** — Local mode full flow: Electron spawns server → React connects → run RNA-seq pipeline on test data end to end
- [ ] **4.2** — Local mode full flow: variant calling track end to end
- [ ] **4.3** — Remote mode full flow: React → auth token → remote server running on a separate machine
- [ ] **4.4** — Error states: tool not found, reference genome missing, FASTQ file corrupted, GATK OOM, server unreachable, version mismatch
- [ ] **4.5** — Large file handling: test with full-size FASTQ files (>10GB), verify no memory issues in upload and processing
- [ ] **4.6** — E2E testing strategy: Playwright + Electron driver, test happy path for both tracks

---

## Phase 5 — CI/CD, Docker & Release

> Ship it like a real open-source product.

- [ ] **5.1** — GitHub Actions `ci.yml`: lint + test on every PR (Python + TypeScript)
- [ ] **5.2** — GitHub Actions `release.yml`: build installers + publish on tag push
- [ ] **5.3** — Docker: `omicspilot/bioseqflow-server` image with all NGS tools pre-installed (bioconda base)
- [ ] **5.4** — Docker Compose: server + reference genome volume for self-hosted deployment
- [ ] **5.5** — Community files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue templates, bug report form
- [ ] **5.6** — Changelog (`CHANGELOG.md`, Keep a Changelog format) + semantic versioning strategy
- [ ] **5.7** — OmicsPilot website: add BioSeqFlow to the tools listing alongside CellMetPro

---

## Architecture Notes

### How tools are wrapped

Every NGS tool runs as a managed subprocess. The server:
1. Validates input files and parameters before launching
2. Constructs the command with safe argument escaping
3. Launches as a subprocess with stdout/stderr piped
4. Parses tool-specific output line by line for progress
5. Emits normalised WebSocket events to the frontend
6. Captures exit code and maps to success/failure/partial states
7. Stores result file paths in the database

Example process flow for STAR alignment:
```
Client POST /api/runs/{id}/align
  → server validates BAM exists, reference indexed
  → creates RunStep record (status=RUNNING)
  → subprocess.Popen(["STAR", "--runThreadN", "4", ...])
  → read stdout line by line
  → parse "... % mapped ..." lines
  → emit WebSocket: { tool: "STAR", progress: 0.6, message: "Mapping..." }
  → on exit code 0: update RunStep to COMPLETE, register output BAM
  → on exit code != 0: update to FAILED, store stderr in DB
```

### Key difference from CellMetPro

| Aspect | CellMetPro | BioSeqFlow |
|---|---|---|
| Core computation | Python library call (cellmetpro) | External subprocess calls |
| Progress | Task completion events | Parsed stdout/stderr streams |
| Input files | h5ad (single file, <2GB) | FASTQ pairs (can be >50GB each) |
| Runtime | Minutes | Hours (STAR ~20min, GATK ~2hr) |
| Tool dependencies | Python packages only | External binaries (STAR, GATK, etc.) |
| Reference data | None | Reference genomes (several GB) |

### Reference genome management

Reference genomes are large (human GRCh38 = 3.2GB FASTA, STAR index = 27GB). The GUI must:
- Provide a download manager for common genomes (human, mouse, rat)
- Support user-provided custom FASTA + GTF
- Track indexing progress (STAR index takes ~30min for human)
- Store indexed references in a persistent location (not temp)
- Allow reuse across projects

---

## Legend

| Symbol | Meaning |
|---|---|
| `[ ]` | Not started |
| `[~]` | In progress |
| `[x]` | Done |
