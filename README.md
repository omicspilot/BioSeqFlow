<div align="center">
  <h1>BioSeqFlow</h1>
  <h2>Desktop GUI for NGS Analysis</h2>
  <p>
    <b>RNA-seq differential expression and DNA variant calling — wrapping the industry-standard
    NGS toolchain (FastQC, fastp, STAR, BWA-MEM2, samtools, GATK, featureCounts, pyDESeq2, SnpEff)
    behind a clean, cross-platform desktop app. No command line required.</b>
  </p>
  <p>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
    <a href="https://codecov.io/gh/omicspilot/BioSeqFlow"><img src="https://codecov.io/gh/omicspilot/BioSeqFlow/graph/badge.svg" alt="Codecov"></a>
    <br/>
    <a href="https://github.com/omicspilot/BioSeqFlow/actions/workflows/ci.yml"><img src="https://github.com/omicspilot/BioSeqFlow/actions/workflows/ci.yml/badge.svg" alt="CI workflow"></a>
    <a href="https://github.com/omicspilot/BioSeqFlow/commits/main"><img src="https://img.shields.io/github/last-commit/omicspilot/BioSeqFlow/main" alt="GitHub last commit"></a>
  </p>
</div>

---

## Overview

BioSeqFlow is a single, self-contained desktop application that makes NGS analysis accessible to
researchers who aren't comfortable with the command line or managing bioinformatics tool
installations. One repo, one product: there's no separate pipeline or standalone library to
install — the app directly orchestrates each NGS tool as a managed subprocess.

**Two tracks:**
- **RNA-seq** — QC → trimming → STAR alignment → featureCounts quantification → pyDESeq2
  differential expression (count matrix, DE results, volcano plot)
- **Variant calling** — QC → trimming → BWA-MEM2 alignment → samtools cleanup → GATK
  HaplotypeCaller → SnpEff annotation (annotated VCF, variant table)

See [context/CONTEXT.md](context/CONTEXT.md) for the full architecture and
[context/ROADMAP.md](context/ROADMAP.md) for the live implementation plan and current status.

## Project Status

Early development — Phase 0 (monorepo & tooling foundation) is complete; Phase 1 (the FastAPI
server) is in progress. Nothing here is runnable as a product yet. Check `context/ROADMAP.md` for
what's actually done versus planned.

## Development Setup

```bash
git clone https://github.com/omicspilot/BioSeqFlow.git
cd BioSeqFlow

pnpm install                    # desktop tooling (packages/desktop)

cd packages/server
make dev                        # python env (uv) + git hooks
make check                      # lint + type check
make test                       # test suite
```

## Project Structure

```
bioseqflow/
├── packages/
│   ├── server/                        # bioseqflow-server: pip-installable FastAPI package
│   │   ├── bioseqflow_server/
│   │   │   ├── toolkit/               # internal QC/preprocessing implementation (not a public API)
│   │   │   ├── routers/               # REST endpoints, one file per domain
│   │   │   ├── parsers/               # tool-specific log/output parsers
│   │   │   └── main.py                # FastAPI app entrypoint
│   │   ├── tests/
│   │   ├── pyproject.toml
│   │   └── Makefile
│   └── desktop/                       # Electron + React GUI
│       ├── src/
│       └── electron/
└── context/                           # architecture docs + implementation roadmap
```

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass and coverage is maintained
5. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## Support

- **Discussions**: [GitHub Discussions](https://github.com/omicspilot/BioSeqFlow/discussions)
- **Issues**: [GitHub Issues](https://github.com/omicspilot/BioSeqFlow/issues)
