# `cb60fd71d1d7` Novelty Gate

Status date: `2026-03-16`

## Decision

`cb60fd71d1d7` should **not** yet be labeled a new discovery.

The strongest honest current label is:

```text
high-value unexplained Ramanujan-style candidate
with exact local structure
and multiple nearby-family eliminations
```

## Why The Current Evidence Is Strong

The project now has all of the following in place:

1. Stable numerical behavior
   - the candidate survives higher-precision verification
   - it remains benchmark-adjacent to `RR(q^3)` rather than collapsing under
     re-checks

2. Clean symbolic object
   - in `t = q^3`, the reciprocal object is

```text
C(t) = 1 + K_(n>=1) (t^n + t^(2n)) / (1 + t^n)
```

3. Exact structural proof layer
   - `Proofs/HeroCaseObjects.lean`
   - `Proofs/HeroCaseLocal.lean`
   - `Proofs/RationalEquivalence.lean`

4. Exact or computation-checked elimination layers
   - bounded page-43 monomial-substitution exclusions
   - bounded arithmetic-subsequence contraction exclusions
   - bounded 1/2/3-step Bauer-Muir exclusions
   - small-degree algebraic relation search against the nearest benchmark

5. Public audit posture
   - the repo already separates:
     - exact theorem-grade claims
     - bounded-search evidence
     - open novelty risk

## Why This Is Still Not Enough

Three gaps remain.

### 1. Literature closure is still incomplete

Nearby primary sources have been checked, including:

- Bowman, Mc Laughlin, Wyshinski, 2019 arXiv version of the page-43 / `H`
  families:
  `https://arxiv.org/abs/1901.00584`
- Lee, Mc Laughlin, Sohn, 2019 arXiv version on Bauer-Muir and generalized
  Rogers-Ramanujan-type fractions:
  `https://arxiv.org/abs/1906.11991`
- Berndt, Rebaka, 2025 survey:
  `https://arxiv.org/abs/2512.19952`
- McLaughlin, Monks, Reid, 2025:
  `https://doi.org/10.5281/zenodo.14580993`

This is enough to say "the obvious nearby explanations did not match."

It is **not** enough to say "nobody has published this identity."

### 2. Several exclusion layers are still bounded

The current subsequence and Bauer-Muir exclusions are strong and reproducible,
but they are still bounded by:

- finite transform depth
- finite modifier family
- exact rational sample points
- finite truncation depth

That is excellent audit evidence, but it is not yet a theorem excluding the
entire natural transform class.

### 3. There is no final source theorem

The project does not yet have either of:

- a positive theorem identifying `cb60fd71d1d7` with a known source object
- a negative theorem excluding a sufficiently broad, explicitly defined source
  class that specialists would accept as the relevant neighborhood

Without one of those, "new discovery" would still be ahead of the evidence.

## Current Safe Claim

As of `2026-03-16`, the safe research claim is:

```text
Within the currently audited nearby Ramanujan-style families and bounded
transform classes, cb60fd71d1d7 remains an unexplained and unusually structured
candidate. We do not yet have enough evidence to promote it to a novelty claim.
```

## What Would Be Enough

At least one of the following two routes should be completed.

### Route A: Positive identification

1. find a final source family
2. prove the exact source identity
3. downgrade the object from "mystery" to "recognized variant"

### Route B: Strong novelty support

1. define a broad natural source class precisely
2. prove or computation-check that class at a level substantially beyond the
   current bounded scans
3. run a wider bibliography pass across classical and modern RR / cubic /
   Gollnitz-Gordon / transform papers
4. get external mathematical review before public novelty wording

## Practical Next Actions

1. Upgrade one bounded exclusion family into an all-parameter theorem.
2. Expand the literature log into a dated bibliography matrix by source family.
3. Keep public wording at "unexplained candidate" until steps 1 and 2 are both
   materially stronger.
