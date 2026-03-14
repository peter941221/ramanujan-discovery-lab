# `cb60fd71d1d7` Literature Log

Audit date: `2026-03-13`

## Goal

Check whether the current hero-case candidate

```text
1 / (1 + (q^3 + q^6)/(1 + q^3 + (q^6 + q^12)/(1 + q^6 + ...)))
```

is already explained by a nearby primary-source family.

## Primary Sources Checked

1. Bowman, Mc Laughlin, Wyshinski, 2006
   - arXiv source: https://arxiv.org/abs/1901.00584
   - Useful hits:
     - four-parameter `q`-continued fraction families `H(a,b,c,d,q)` and `H_1(a,b,c,d,q)`
     - Ramanujan page-43 ratio family
       `G(a q, \lambda q; b; q) / G(a, \lambda; b; q)`
     - nearby corollaries for:
       - Rogers-Ramanujan-type denominator perturbations
       - Ramanujan cubic continued fractions

2. Lee, Mc Laughlin, Sohn, 2020
   - arXiv source: https://arxiv.org/abs/1906.11991
   - Useful hits:
     - Bauer-Muir transformation links among generalized Rogers-Ramanujan-type fractions
     - Ramanujan page-43 continued fractions `gcf2`, `gcf3`
     - Hirschhorn-type transformed families `gcf4`

## What Survived

- The candidate is close to the `RR(q^3)` benchmark family.
- The candidate also reuses the cubic-style extra numerator pattern `q^(3n) + q^(6n)`.
- The nearest primary-source neighborhood is therefore:
  - generalized Rogers-Ramanujan page-43 ratio families
  - denominator-perturbed Rogers-Ramanujan families
  - Ramanujan cubic continued fractions
  - Bauer-Muir transformed relatives of those families

## What Was Ruled Out

- No direct constant-parameter specialization of Ramanujan Entry 6.4.4 / `gcf3` was found.
- For `t = q^3`, matching

```text
lambda * t^n - a*b * t^(2n) = t^n + t^(2n)
1 + t^n * (a*t + b) = 1 + t^n
```

with fixed complex constants `a`, `b`, `lambda` has no solution.

## Search Outcome

- No exact primary-source hit for the mixed pattern

```text
(q^(3n) + q^(6n)) / (1 + q^(3n))
```

was found in the sources checked.
- The closest primary-source matches remain family-level, not exact-pattern matches.

## Transform Audit Update

- A dedicated audit note now records the first explicit transform eliminations:
  - `CB60FD71D1D7_TRANSFORM_AUDIT.md`
- Current eliminations now include:
  - no direct constant-parameter match to Lee-Mc Laughlin-Sohn `f2`, `f3`, `f4`, or `cor2cf`
  - no direct constant-parameter match to Bowman-Mc Laughlin-Wyshinski `H_1`
  - no simple odd-part origin from Bowman-Mc Laughlin-Wyshinski `H_2` or `H_3`
  - no direct 1-step Bauer-Muir origin from either the RR reciprocal or the cubic reciprocal in reduced `t = q^3` form
  - no hit in the current tiny low-complexity 1-step / 2-step / 3-step Bauer-Muir scan built from those RR/cubic reciprocals
- Additional signals recorded in the transform audit:
  - `cb60fd71d1d7 / RR(q^3)` appears to be a pure `t=q^3` series in the visible symbolic window
  - the only viable two-step contraction branch around the `cor2cf` specialization that preserves `b0=1` fails at the first numerator shape
  - simple odd/even-part contractions of the cubic reciprocal were also reconstructed and do not match the candidate's partial denominators

## Current Interpretation

- This is not enough to claim novelty.
- It is enough to keep `cb60fd71d1d7` as the best current hero case.
- The next literature pass should be narrower:
  - contraction / even-odd part identities
  - longer or less rigid Bauer-Muir chains beyond the direct and tiny 3-step RR/cubic paths already ruled out
  - page-43 Ramanujan ratio families after nontrivial substitutions
