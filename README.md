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
python -m ramanujan_discovery site --in results/verified.jsonl --out-dir docs
```

Then open `docs/index.html` locally or publish the `docs/` folder with GitHub Pages.

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
- there is no formal proof layer yet
- the search family is still intentionally small and biased toward interpretable templates
