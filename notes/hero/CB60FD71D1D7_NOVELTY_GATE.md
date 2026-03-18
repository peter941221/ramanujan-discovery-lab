# `cb60fd71d1d7` Novelty Gate

Status date: `2026-03-17`

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
   - `Proofs/HeroCaseFinalIdentity.lean`, which now packages the current
     exact waypoint from finite-convergent rational equivalence plus the
     page-43 nearest-shift-cube exclusion layer, even though it still does not
     contain a positive final identity theorem

4. Exact or computation-checked elimination layers
     - bounded page-43 monomial-substitution exclusions
     - an exact zero-shift `f2/gcf3` `n`-dependent equivalence obstruction,
       where the necessary residual polynomial forces `a = 0`, `b = 0`, then
       `lambda = 1`, but still leaves a surviving `m^2` coefficient `t`
     - an exact zero-shift `f4/gcf2` `n`-dependent equivalence obstruction,
       where the necessary residual polynomial again forces `a = 0`, `b = 0`,
       then `lambda = 1`, but still leaves a surviving `m^2` coefficient `t`
     - an exact unit-`a`-shift extension of those same `f2/gcf3` and `f4/gcf2`
       lanes, where the forced parameter specializations still leave a
       surviving nonzero `m^2` coefficient `t`
     - an exact unit-`b`-shift extension of those same `f2/gcf3` and `f4/gcf2`
       lanes, where the forced parameter specializations again still leave a
       surviving nonzero `m^2` coefficient `t`
     - an exact mixed unit-`a` / unit-`b` extension of those same
       `f2/gcf3` and `f4/gcf2` lanes, where the forced parameter
       specializations still end with the same surviving nonzero `m^2`
       coefficient `t`
     - an exact mixed unit-`a` / unit-`lambda` extension of those same
       `f2/gcf3` and `f4/gcf2` lanes, where the surviving `m^1`
       coefficient is again `lambda*t^2 - t`, so no constant `lambda` can
       make either lane vanish identically
     - an exact mixed unit-`b` / unit-`lambda` extension of those same
       `f2/gcf3` and `f4/gcf2` lanes, where the surviving `m^1`
       coefficient is again `lambda*t^2 - t`, so no constant `lambda` can
       make either lane vanish identically
     - an exact mixed unit-`a` / unit-`b` / unit-`lambda` extension of those
       same `f2/gcf3` and `f4/gcf2` lanes, where the surviving `m^1`
       coefficient is again `lambda*t^2 - t`, so no constant `lambda` can
       make either lane vanish identically
     - an exact unit-`lambda`-shift extension of those same `f2/gcf3` and
       `f4/gcf2` lanes, where the surviving `m^1` coefficient is
       `lambda*t^2 - t`, so no constant `lambda` can make either nearest shift
       lane vanish identically
     - collectively, those eight nearby cases now close the full nearest-shift
       cube at the current exact page-43 theorem layer
     - bounded arithmetic-subsequence contraction exclusions, now accompanied by a
       stronger exact obstruction for the RR and cubic source lanes, plus a
       machine-checked Lean proof layer:
       - `Proofs/HeroCaseSubsequenceExact.lean`
     - bounded 1/2/3-step Bauer-Muir exclusions
     - small-degree algebraic relation search against the nearest benchmark,
       including:
       - direct polynomial boxes in `(C, B1)` and low-degree RR-tower prefixes
       - ratio-object polynomial prefix boxes for `F = candidate / RR(q^3)`
       - ratio-object multiplicative RR-tower boxes of the form `F = prod_i B_i^(e_i)`
       - ratio-object one-layer fractional-linear RR-tower boxes
       - ratio-object two-layer single-basis fractional-linear RR-tower boxes
       - named source-family multiplicative boxes for
         `F = candidate / RR(q^3)` built from `RR`, `cubic`, `GG`, and `S`
       - named source-family fractional-linear boxes for
         `F = candidate / RR(q^3)` built from the same `RR`, `cubic`, `GG`,
         and `S` basis order
       - named source-family two-layer fractional-linear boxes for
         `F = candidate / RR(q^3)` built from the same source-family prefix
         order, still with `0` hits through prefixes ending at `cubic`, `GG`,
         and `S`
       - per-family powered source-family ladders, still with `0` hits in low-degree
         polynomial, multiplicative, fractional-linear, family-preserving
         two-layer fractional-linear, within-family quotient-ladder, and
         quotient-ladder two-layer, mixed quotient-basis, and mixed
         quotient-basis two-layer boxes for
         `RR/RR2/RR3/RR4`, `cubic/cubic2/cubic3/cubic4`,
         `GG/GG2/GG3/GG4/GG5/GG7/GG11`, and `S/S2/S3/S4`
       - an explicit `GG` / `S` transform-template pass over direct objects,
         reciprocals, and pairwise quotients in those same ladders, still with
         `0` exact hits

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
- classical Gordon / Hirschhorn citations are now pinned more explicitly in the
  project notes:
  - Basil Gordon, *Some continued fractions of the Rogers-Ramanujan type*,
    Duke Mathematical Journal 32 (1965), 741-748,
    https://doi.org/10.1215/S0012-7094-65-03278-3
  - M. D. Hirschhorn, *A continued fraction*,
    Duke Mathematical Journal 41 (1974), 27-33,
    https://doi.org/10.1215/S0012-7094-74-04104-0
  - M. D. Hirschhorn, *A continued fraction of Ramanujan*,
    Journal of the Australian Mathematical Society Series A 29 (1980), 80-86,
    https://doi.org/10.1017/S1446788700020954
- a modern structural Hirschhorn paper:
  - Bhatnagar, Ismail, *Orthogonal polynomials associated with a continued
    fraction of Hirschhorn*,
    https://arxiv.org/abs/1901.09985
- a direct GG power-relation paper:
  - Vasuki, Srivatsa Kumar, *Certain identities for
    Ramanujan-Göllnitz-Gordon continued fraction*,
    https://doi.org/10.1016/j.cam.2005.03.038
- a notebook-wide survey source:
  - Bhatnagar, *Ramanujan's `q`-continued fractions*,
    https://arxiv.org/abs/2208.12656
- Berndt, Rebaka, 2025 survey:
  `https://arxiv.org/abs/2512.19952`

This is enough to say "the obvious nearby explanations did not match."

It is **not** enough to say "nobody has published this identity."

### 2. Several exclusion layers are still bounded

The current subsequence and Bauer-Muir exclusions are strong and reproducible,
but they are still bounded by:

- finite transform depth
- finite modifier family
- finite generic relation-template families around the RR tower
- exact rational sample points
- finite truncation depth

That is excellent audit evidence, but it is not yet a theorem excluding the
entire natural transform class.

Note: for the two closest source lanes (reduced RR reciprocal and reduced cubic
reciprocal), the arithmetic-subsequence obstruction is now exact and
machine-checked in Lean. The remaining boundedness mostly concerns other source
classes, the Bauer-Muir / transform neighborhoods, and source families not
captured by the current low-complexity RR-tower ansatz boxes.

### 3. There is no final source theorem

The project does not yet have either of:

- a positive theorem identifying `cb60fd71d1d7` with a known source object
- a negative theorem excluding a sufficiently broad, explicitly defined source
  class that specialists would accept as the relevant neighborhood

Without one of those, "new discovery" would still be ahead of the evidence.

## Current Safe Claim

As of `2026-03-17`, the safe research claim is:

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

1. Move beyond the now-formalized nearest-shift page-43 cube into either a
   broader all-parameter theorem for that neighborhood or the next natural
   exact transform template that still preserves source-family meaning.
2. Move beyond powered family ladders into explicit reciprocation / quotient /
   transform templates for the `GG` and `S` orbits, keeping the family meaning
   visible instead of flattening everything into anonymous basis elements.
3. Continue literature closure across the Gordon / Hirschhorn / Weber orbit,
   especially sources that could hide transformed presentations of `GG(q)` or
   `S(q)` rather than just their base product formulas.
4. Keep public wording at "unexplained candidate" until steps 1-3 are all
   materially stronger.

## Companion Notes

- `CB60FD71D1D7_EXACT_SUBSEQUENCE_OBSTRUCTION.md`
- `CB60FD71D1D7_HEINE_COR2CF_OBSTRUCTION.md`
- `CB60FD71D1D7_BIBLIOGRAPHY_MATRIX.md`
- `CB60FD71D1D7_LITERATURE_LOG.md`
