# Ramanujan Discovery Lab

`Ramanujan Discovery Lab` is a local-first research tool for rediscovering and stress-testing Ramanujan-style `q`-continued fractions.

- Repository: `https://github.com/peter941221/ramanujan-discovery-lab`
- Pages: `https://peter941221.github.io/ramanujan-discovery-lab/`

The current release is intentionally narrow:

- CLI-first workflow
- high-precision numerical evaluation with `mpmath`
- benchmark verification against independent infinite-product formulas
- GitHub Pages-friendly static site generation
- GitHub Actions workflow that reruns discovery, verification, reporting, and Pages deployment on `main`

## What It Does

The pipeline has four stages:

1. `discover`
   Search a constrained family of `q`-continued-fraction templates.
2. `verify`
   Recompute promising templates at higher precision and classify them.
3. `report`
   Emit a Markdown snapshot for research notes or release artifacts.
4. `site`
   Generate a static gallery suitable for GitHub Pages.

Optional research stage:

5. `analyze`
   Write a focused Markdown note for one candidate, including structural deltas, symbolic `q`-series comparisons, and terminal-friendly formula output.
6. `research`
   Run heavier candidate probes (higher-order ratio series, Euler-product exponent extraction, and targeted literature coefficient checks).
7. `formalize`
   Write a formalization-prep note for one candidate, separating exact local lemmas from bounded search evidence.
8. `identify`
   Attempt small algebraic relation guesses between a candidate and its nearest benchmark (step-reduced if possible).

Current benchmark catalog:

- Rogers-Ramanujan normalized continued fraction
- Rogers-Ramanujan family benchmarks at `q^2`, `q^3`, and `q^4`
- Ramanujan cubic continued fraction
- Ramanujan cubic family benchmarks at `q^2` and `q^3`
- Conservative review audit for numerically stable unexplained candidates

## Quickstart

```powershell
python -m pip install -e .[dev]
$env:PYTHONPATH='src'
python -m ramanujan_discovery discover --depth 36 --precision 80 --budget-hours 0.1 --out results/candidates.jsonl
python -m ramanujan_discovery verify --in results/candidates.jsonl --precision 160 --out results/verified.jsonl
python -m ramanujan_discovery report --in results/verified.jsonl --out results/report.md
python -m ramanujan_discovery analyze --in results/verified.jsonl --candidate-id cb60fd71d1d7 --stdout-format unicode --out HERO_CASE_CB60FD71D1D7.md
python -m ramanujan_discovery research --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 151 --out CB60FD71D1D7_RESEARCH_NOTE.md
python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out CB60FD71D1D7_FORMALIZATION_NOTE.md --lean-out proofs/Proofs/Generated/Cb60fd71d1d7.lean
python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out CB60FD71D1D7_IDENTIFICATION_NOTE.md
python -m ramanujan_discovery site --in results/verified.jsonl --out-dir docs
```

`analyze` stdout formats:

- `unicode`: terminal-friendly pretty math via SymPy
- `latex`: copyable LaTeX string for Markdown or papers
- `plain`: ASCII fallback
- `none`: suppress terminal summary

Fast local validation:

```powershell
python -m ramanujan_discovery research --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 12 --series-order 61 --smoke --out tmp/research-smoke.md
python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --smoke --out tmp/formalize-smoke.md --lean-out tmp/formalize-smoke.lean
python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --smoke --out tmp/identify-smoke.md
```

`--smoke` keeps the same output sections but trims bounded symbolic searches so
the command is suitable for quick regression checks.

Then open `docs/index.html` locally or publish the `docs/` folder with GitHub Pages.

## Lean Proof Layer

The repository now includes an initial Lean 4 proof workspace under `proofs/`.
Its current scope is intentionally narrow but no longer purely ad hoc:

- `Proofs/GeneralizedCF.lean` formalizes finite-truncation continuants and the
  convergent recurrence
- `Proofs/HeroCaseObjects.lean` centralizes the canonical hero-case objects,
  nearby source objects, and generic reverse-scale definitions
- `Proofs/HeroCaseLocal.lean` machine-checks exact local obstruction lemmas for
  `cb60fd71d1d7`, including polynomial-level no-match theorems and an exact
  convergent-factor reduction theorem
- `Proofs/HeroCasePage43.lean` formalizes bounded page-43 family exclusions for
  the current hero case in Laurent polynomials
- `Proofs/HeroCaseSubsequence.lean` and `Proofs/HeroCaseBauerMuir.lean` mirror
  the bounded subsequence and Bauer-Muir scans as computation-checked theorems
- the current hero case also admits an explicit reverse equivalence
  transformation from the reduced fraction back to the original reciprocal
  object
- `Proofs/RationalEquivalence.lean` now proves that reverse equivalence witness
  over `RatFunc Rat` for all stages of the shared hero-case object
- `formalize` and `research` now emit that exact reduction/equivalence witness
  explicitly; the reverse transform scales are rational functions, so the next
  formal step likely needs a fraction-field coefficient layer
- `formalize --lean-out ...` can now auto-generate a candidate-specific Lean
  proof module such as `proofs/Proofs/Generated/Cb60fd71d1d7.lean`

```powershell
Set-Location proofs
lake build
lake env lean Proofs/Generated/Cb60fd71d1d7.lean
```

## Status Semantics

- `known`: exact canonical match to a classical benchmark
- `known_variant`: strong numerical match to a classical benchmark but not canonical
- `fixture`: internal regression target
- `review`: stable unexplained candidate that does not strongly match the benchmark catalog

Review audit:

- [REVIEW_AUDIT.md](REVIEW_AUDIT.md)
- Pages view: `review-audit.html`

## Current Limits

- novelty is only checked against the built-in benchmark catalog
- there is not yet a full formal proof of a final continued-fraction identity
- `formalize` now emits a theorem-prep note plus an optional Lean proof module, but the formal layer is still obstruction-first rather than a final source proof
- the search family is still intentionally small and biased toward interpretable templates
