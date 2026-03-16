# Ramanujan Discovery Lab

[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![Lean](https://img.shields.io/badge/Lean-4.28.0-0B5CAD?style=for-the-badge)](https://lean-lang.org/)
[![Pages](https://img.shields.io/website?style=for-the-badge&url=https%3A%2F%2Fpeter941221.github.io%2Framanujan-discovery-lab%2F&label=GitHub%20Pages)](https://peter941221.github.io/ramanujan-discovery-lab/)
[![Workflow](https://img.shields.io/github/actions/workflow/status/peter941221/ramanujan-discovery-lab/research-pages.yml?style=for-the-badge&label=Research%20Pages)](https://github.com/peter941221/ramanujan-discovery-lab/actions/workflows/research-pages.yml)
[![License](https://img.shields.io/github/license/peter941221/ramanujan-discovery-lab?style=for-the-badge)](LICENSE)

```text
+----------------------------------------------------------------------------+
|                         Ramanujan Discovery Lab                            |
|                                                                            |
|                 1                                                          |
| F(q) = -----------------------     discover -> verify -> report -> site    |
|               a1(q)                                                        |
|        b0(q) + --------------     analyze -> research -> identify         |
|                     a2(q)                                                  |
|               b1(q) + -------- ... formalize -> Lean 4 -> public notes    |
|                          ...                                               |
|                                                                            |
| A CLI-first research workspace for Ramanujan-style q-continued fractions.  |
+----------------------------------------------------------------------------+
```

*Mathematical beauty is treated as a clue, not as a proof.*


- Repository: `https://github.com/peter941221/ramanujan-discovery-lab`
- GitHub Pages: `https://peter941221.github.io/ramanujan-discovery-lab/`
- Stack: `Python 3.11+`, `mpmath`, `sympy`, `Lean 4`
- Automation: GitHub Actions reruns discovery, verification, reporting, and Pages deployment on `main`
- Posture: conservative, reproducible, benchmark-backed, proof-aware

```text
+==========================================================================+
|                                Contents                                  |
+==========================================================================+
```

- [Why This Repository Exists](#why-this-repository-exists)
- [Research Posture](#research-posture)
- [Project Shape](#project-shape)
- [Pipeline](#pipeline)
- [Quickstart](#quickstart)
- [Current Mathematical Focus](#current-mathematical-focus)
- [Lean Proof Layer](#lean-proof-layer)
- [Benchmark Catalog](#benchmark-catalog)
- [Repository Map](#repository-map)
- [Status Semantics](#status-semantics)
- [Current Limits](#current-limits)
- [Public Artifacts](#public-artifacts)

```text
+==========================================================================+
|                           Why This Repository Exists                     |
+==========================================================================+
```

## Why This Repository Exists

```text
Goal
|-- rediscover classical Ramanujan-style q-continued fractions
|-- surface stable unexplained candidates from a small interpretable search box
|-- audit those candidates with symbolic, numerical, and formal tools
`-- publish artifacts that stay honest about what is known and what is not
```

This repository is a local-first research workspace for discovering and
stress-testing Ramanujan-style `q`-continued-fraction templates. It is designed
to do two things well at the same time:

1. recover classical benchmark identities from an independently checkable catalog
2. investigate nearby unexplained candidates without prematurely calling them
   “new identities”

The current award-track direction is centered on a single hero case,
`cb60fd71d1d7`, with the long-term goal of reaching:

```text
unexplained candidate
      |
      +-- exact recognition / closed form
      +-- theorem-grade proof
      +-- Lean formalization
      `-- literature closure
  |
             `--> only then: publishable novelty posture
```

```text
+==========================================================================+
|                              Research Posture                            |
+==========================================================================+
```

## Research Posture

```text
candidate
  |
  +-- exact benchmark match ----------> known
  |
  +-- strong noncanonical benchmark --> known_variant
  |
  +-- internal regression control ----> fixture
  |
  `-- stable but unmatched -----------> review
                                         |
                                         `-- review != novelty claim
```

This repository is intentionally conservative:

- `review` means “numerically stable and not matched by the current built-in
  catalog”; it does **not** mean “new formula”
- closed-world benchmark novelty is not the same as literature novelty
- bounded transform scans are evidence, not final proofs
- public-facing notes stay cautious until the theorem gate and literature gate
  are both materially satisfied

```text
+==========================================================================+
|                               Project Shape                              |
+==========================================================================+
```

## Project Shape

```text
Ramanujan Discovery Lab
|-- src/ramanujan_discovery/
|   |-- discovery.py        search bounded template families
|   |-- verification.py     high-precision rescoring and classification
|   |-- reporting.py        markdown + GitHub Pages artifacts
|   |-- analysis.py         focused candidate notes
|   |-- research.py         heavier symbolic and transform probes
|   |-- identification.py   bounded algebraic relation guesses
|   `-- formalization.py    theorem-prep notes + Lean module generation
|
|-- tests/                  regression checks
|-- proofs/                 Lean 4 workspace
|-- results/                generated candidate snapshots
|-- notes/
|   |-- hero/               hero-case research and formalization notes
|   `-- review/             review-set audits and comparison notes
|-- docs/                   static public site for GitHub Pages
`-- README.md               public-facing project overview
```

```text
+==========================================================================+
|                                  Pipeline                                |
+==========================================================================+
```

## Pipeline

```text
discover ---> verify ---> report ---> site
    |            |
    |            +----> analyze ---> research ---> identify
    |                                   |
    |                                   `----------> formalize ---> Lean 4
    |
    `---- bounded, interpretable search over q-continued-fraction templates
```

### Core Stages

- `discover`
  Search a constrained family of `q`-continued-fraction templates.
- `verify`
  Recompute promising templates at higher precision and classify them.
- `report`
  Emit a Markdown research snapshot from verified candidates.
- `site`
  Render a static gallery for GitHub Pages.

### Deep-Dive Stages

- `analyze`
  Produce a focused note for one candidate with symbolic `q`-series deltas.
- `research`
  Run heavier probes: ratio-series fits, Euler-product exponent extraction,
  structured fit attempts, and transform audits.
- `identify`
  Attempt small algebraic relation guesses against the nearest benchmark.
- `formalize`
  Produce formalization-prep notes and optionally auto-generate Lean modules.

```text
+==========================================================================+
|                                 Quickstart                               |
+==========================================================================+
```

## Quickstart

```powershell
python -m pip install -e .[dev]
$env:PYTHONPATH='src'
python -m ramanujan_discovery discover --depth 36 --precision 80 --budget-hours 0.1 --out results/candidates.jsonl
python -m ramanujan_discovery verify --in results/candidates.jsonl --precision 160 --out results/verified.jsonl
python -m ramanujan_discovery report --in results/verified.jsonl --out results/report.md
python -m ramanujan_discovery site --in results/verified.jsonl --out-dir docs
```

Then open `docs/index.html` locally or publish `docs/` with GitHub Pages.

### Hero-Case Deep Dive

```powershell
$env:PYTHONPATH='src'
python -m ramanujan_discovery analyze --in results/verified.jsonl --candidate-id cb60fd71d1d7 --stdout-format unicode --out notes/hero/HERO_CASE_CB60FD71D1D7.md
python -m ramanujan_discovery research --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 151 --out notes/hero/CB60FD71D1D7_RESEARCH_NOTE.md
python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md --lean-out proofs/Proofs/Generated/Cb60fd71d1d7.lean
python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md
```

### Fast Smoke Validation

```powershell
$env:PYTHONPATH='src'
python -m ramanujan_discovery research --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 12 --series-order 61 --smoke --out tmp/research-smoke.md
python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --smoke --out tmp/formalize-smoke.md --lean-out tmp/formalize-smoke.lean
python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --smoke --out tmp/identify-smoke.md
```

`--smoke` keeps the same output sections while reducing bounded symbolic search
coverage for quicker local validation.

```text
+==========================================================================+
|                         Current Mathematical Focus                       |
+==========================================================================+
```

## Current Mathematical Focus

```text
RR(q^3) neighborhood
|
|-- e2cc74240b6f   plain step branch
`-- cb60fd71d1d7  hybrid branch / current hero case
```

The current hero case is `cb60fd71d1d7`.

Current conservative read:

- it tracks `RR(q^3)` through `q^9`
- it first diverges at `q^12`
- it appears structurally richer than the plain `RR(q^3)` step-branch candidate
  `e2cc74240b6f`
- it remains an **unexplained candidate**, not a settled novelty claim

Key notes:

- [Hero case summary](notes/hero/CB60FD71D1D7_PUBLIC_SUMMARY.md)
- [Hero case study](notes/hero/CB60FD71D1D7_CASE_STUDY.md)
- [Hero transform audit](notes/hero/CB60FD71D1D7_TRANSFORM_AUDIT.md)
- [Hero exact subsequence obstruction](notes/hero/CB60FD71D1D7_EXACT_SUBSEQUENCE_OBSTRUCTION.md)
- [Hero Heine `cor2cf` obstruction](notes/hero/CB60FD71D1D7_HEINE_COR2CF_OBSTRUCTION.md)
- [Hero identification note](notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md)
- [Hero formalization note](notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md)
- [Hero novelty gate](notes/hero/CB60FD71D1D7_NOVELTY_GATE.md)
- [Hero bibliography matrix](notes/hero/CB60FD71D1D7_BIBLIOGRAPHY_MATRIX.md)
- [Hero literature log](notes/hero/CB60FD71D1D7_LITERATURE_LOG.md)
- [Hero public article draft](notes/hero/CB60FD71D1D7_PUBLIC_ARTICLE.md)

```text
+==========================================================================+
|                              Lean Proof Layer                            |
+==========================================================================+
```

## Lean Proof Layer

```text
proofs/
|-- Proofs/GeneralizedCF.lean
|-- Proofs/HeroCaseObjects.lean
|-- Proofs/HeroCaseLocal.lean
|-- Proofs/HeroCasePage43.lean
|-- Proofs/HeroCaseHeineCor2cf.lean
|-- Proofs/HeroCaseSubsequence.lean
|-- Proofs/HeroCaseSubsequenceExact.lean
|-- Proofs/HeroCaseBauerMuir.lean
|-- Proofs/RationalEquivalence.lean
`-- Proofs/Generated/Cb60fd71d1d7.lean
```

The Lean workspace is still obstruction-first rather than final-identity-first,
but it is already useful:

- `Proofs/GeneralizedCF.lean`
  formalizes finite-truncation continuants and convergent recurrences
- `Proofs/HeroCaseLocal.lean`
  proves exact local mismatch lemmas and convergent-factor reduction facts
- `Proofs/HeroCasePage43.lean`
  formalizes bounded page-43 family exclusions
- `Proofs/HeroCaseSubsequence.lean` and `Proofs/HeroCaseSubsequenceExact.lean`
  cover bounded and stronger exact subsequence obstructions
- `Proofs/HeroCaseBauerMuir.lean`
  mirrors the bounded Bauer-Muir search layer
- `Proofs/RationalEquivalence.lean`
  proves the reverse equivalence witness over `RatFunc Rat`
- `formalize --lean-out ...`
  can generate candidate-specific Lean modules

Build the proof layer with:

```powershell
Set-Location proofs
lake build
lake env lean Proofs/Generated/Cb60fd71d1d7.lean
```

```text
+==========================================================================+
|                              Benchmark Catalog                           |
+==========================================================================+
```

## Benchmark Catalog

```text
benchmark catalog
|-- Rogers-Ramanujan: q, q^2, q^3, q^4
|-- Ramanujan cubic: q, q^2, q^3
`-- internal fixtures: shifted / denominator-perturbed
```

Current built-in benchmark families:

- Rogers-Ramanujan normalized continued fraction
- Rogers-Ramanujan family benchmarks at `q^2`, `q^3`, and `q^4`
- Ramanujan cubic continued fraction
- Ramanujan cubic family benchmarks at `q^2` and `q^3`
- two internal fixtures for regression control

```text
+==========================================================================+
|                               Repository Map                             |
+==========================================================================+
```

## Repository Map

```text
.
|-- README.md
|-- AGENTS.md
|-- MEMORY.md
|-- RUNBOOK.md
|-- TECHNICAL_DESIGN.md
|-- notes/
|   |-- hero/
|   `-- review/
|-- src/
|   `-- ramanujan_discovery/
|       |-- __main__.py
|       |-- cli.py
|       |-- discovery.py
|       |-- verification.py
|       |-- reporting.py
|       |-- analysis.py
|       |-- research.py
|       |-- identification.py
|       `-- formalization.py
|-- tests/
|-- proofs/
|-- results/
`-- docs/
```

Important entrypoints:

- `src/ramanujan_discovery/cli.py`
- `src/ramanujan_discovery/__main__.py`
- `notes/hero/`
- `notes/review/`
- `proofs/`
- `tests/`
- `docs/`

```text
+==========================================================================+
|                              Status Semantics                            |
+==========================================================================+
```

## Status Semantics

- `known`
  exact canonical match to a classical benchmark
- `known_variant`
  strong numerical match to a classical benchmark but not canonical
- `fixture`
  internal regression target
- `review`
  stable unexplained candidate that does not strongly match the current
  benchmark catalog

Review audit:

- [Review audit note](notes/review/REVIEW_AUDIT.md)
- Pages view: `review-audit.html`

```text
+==========================================================================+
|                               Current Limits                             |
+==========================================================================+
```

## Current Limits

```text
what this repo does not yet claim
|-- full literature closure
|-- final theorem for the hero identity
|-- final Lean proof of the target source identity
`-- broad open-ended search over all Ramanujan-style families
```

Current constraints:

- novelty is only checked against the built-in benchmark catalog
- the search family is intentionally small and biased toward interpretable
  templates
- `formalize` currently emits theorem-prep assets rather than a completed final
  source proof
- bounded scans and symbolic fits are valuable evidence, but they are not the
  same thing as a theorem

```text
+==========================================================================+
|                              Public Artifacts                            |
+==========================================================================+
```

## Public Artifacts

- GitHub Pages snapshot: `docs/` or
  `https://peter941221.github.io/ramanujan-discovery-lab/`
- candidate report: `results/report.md`
- hero-case long-form notes:
  `notes/hero/HERO_CASE_CB60FD71D1D7.md`,
  `notes/hero/CB60FD71D1D7_RESEARCH_NOTE.md`,
  `notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md`

---

```text
Beauty -> pattern
pattern -> conjecture
conjecture -> proof
proof -> publication
```

This project tries to keep all four stages in one place.
