# AI-Friendly Publishing Plan

This plan assumes the goal is:

```text
maximize discoverability
+ keep claims conservative
+ stay easy for humans and machines to read, cite, and summarize
```

## Best Stack

Use all four together, in this order:

1. `GitHub Pages`
2. `arXiv`
3. `Zenodo`
4. `Hugging Face Hub`

That stack gives you:

- human-readable HTML
- stable academic citation
- DOI-based indexing
- strong AI/ML audience reach

## Recommended Order

### 1. GitHub Pages

Publish first on the existing project site:

- fastest to update
- static HTML is easy for search engines and AI tools to parse
- lets you link every supporting artifact directly

Use it for:

- the public article
- the hero summary
- the transform audit
- the identification note
- the repo + proof links

Recommended page to publish:

- `notes/hero/CB60FD71D1D7_PUBLIC_ARTICLE.md`

### 2. arXiv

Use arXiv for the "serious research object" version.

Recommended framing:

```text
case study / methods note
```

not

```text
new identity announcement
```

Suggested title direction:

- `An Audit-First Case Study in Experimental Mathematics for q-Continued Fractions`

Suggested positioning:

- discovery pipeline
- conservative benchmark-relative classification
- one hero case with exact local Lean results

### 3. Zenodo

Archive the repo release on Zenodo and mint a DOI.

Use it for:

- repository snapshot
- release tarball
- exact artifact version cited by the article

This is the cleanest way to make the software/artifact citable even before the
mathematical identity is resolved.

### 4. Hugging Face Hub

Use a `dataset` or `space` style repo as a distribution mirror for:

- `results/verified.jsonl`
- the public article
- audit notes
- optional small demo or browser for candidate records

This is valuable because the AI/ML crowd already uses Hugging Face as a
discovery surface, and its metadata/tagging system is machine-friendly.

## Platform Ranking

### Tier A

1. `arXiv`
   - strongest for academic legitimacy and agent-friendly scholarly discovery
2. `GitHub Pages`
   - strongest for easy HTML reading, linking, and fast iteration
3. `Zenodo`
   - strongest for DOI, metadata, and archival citation
4. `Hugging Face Hub`
   - strongest for AI-native audience and metadata-driven discoverability

### Tier B

5. `OpenAlex` via DOI/arXiv ingestion
   - not a publishing platform, but important as a machine-readable downstream index
6. `DataCite` via Zenodo DOI metadata
   - also not a direct publishing platform, but critical for downstream discovery

### Tier C

7. `Hashnode` or `DEV`
   - fine as blog mirrors, but weaker than the four above for citation-grade research artifacts
8. `X / Reddit / Hacker News`
   - useful for traffic bursts, not for durable scholarly presence

## Where I Would Actually Publish This

If we want a practical, AI-friendly rollout, I would do:

1. publish the HTML version on GitHub Pages
2. submit a short methods/case-study version to arXiv
3. create a Zenodo release with DOI for the exact repo snapshot
4. mirror the public article and structured outputs on Hugging Face Hub
5. only then amplify on X / Reddit / Hacker News / LessWrong if desired

## Copy Blocks

### Short Blurb

We built a CLI-first search pipeline for Ramanujan-style `q`-continued fractions
and found a strong Rogers-Ramanujan-adjacent candidate. We are not claiming a
new identity yet. The interesting part is the audit workflow: symbolic
elimination, conservative benchmark-relative claims, and Lean-checked exact
local structure around one hard example.

### One-Line Tagline

Experimental mathematics, but with stronger audit discipline and theorem-grade
local checkpoints.

### Best Audience Fit

- experimental mathematics
- formalized mathematics / Lean
- AI for mathematics
- open science / reproducible research tooling

