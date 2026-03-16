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

## 2026-03-15 Update: Extra Sources Scanned

The goal of this pass was to widen the net toward older Rogers-Ramanujan-type
continued fractions with `1 + q^n` partial denominators and classical modular
invariant presentations.

1. Bowman, Mc Laughlin, Wyshinski, 2005 (Hirschhorn special issue)
   - PDF mirror: https://www.wcupa.edu/sciences-mathematics/mathematics/documents/hirschhornnov08.pdf
   - Relevant hits:
     - reprints Ramanujan notebook page-43 ratio family (Entry 6.4.4) in a
       form consistent with the 2006 arXiv version
     - includes two classical `1+q^n`-denominator fractions:
       - `S(q)` with product `(q^2;q^3)_∞/(q;q^3)_∞`
       - the Ramanujan--Gollnitz--Gordon continued fraction `GG(q)`
   - Outcome: no exact stage pattern match to
     `(t^n + t^(2n)) / (1 + t^n)` was found in the explicit formulas scanned.

2. Adiga, Kim, et al., 2017 (Ramanujan--Weber class invariants)
   - PDF: https://www.filomat.org/index.php/filomat/article/download/4026/4026/31-13-2-4026.pdf
   - Relevant hits:
     - continued fractions tied to Weber class invariants and the
       Ramanujan--Gollnitz--Gordon neighborhood
     - functional equations among those invariants
   - Outcome: noted as a good next target for transform-style recognition, but
     no direct coefficient-level template match was extracted in this scan.

3. Ismail--Stanton, 2003 (orthogonal polynomial viewpoint)
   - PDF: http://www-users.math.umn.edu/~isman/papers/ramanujan.pdf
   - Reason: many Rogers-Ramanujan-type continued fractions can be recognized
     as J-fractions / Stieltjes transforms with structured recurrence data.
   - Outcome: no direct match was found in a first text search for the exact
     `1+q^n` pattern, but the framework is relevant for deeper classification.

## Updated Interpretation

- The candidate still has no identified source identity.
- The `1 + q^n` partial denominator feature is common in older literature
  (e.g., `S(q)`, `GG(q)`), but the exact mixed numerator
  `q^(3n) + q^(6n)` structure remains unlocated.

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

## 2026-03-16 Update: Novelty Gate Reminder

This pass was not aimed at finding one more nearby family. It was aimed at
checking whether the project is already in a position to upgrade the hero case
from "best unexplained object" to "new discovery".

Additional external sources checked:

1. Berndt, Rebaka, 2025 survey
   - arXiv: https://arxiv.org/abs/2512.19952
   - Relevance:
     - recent survey-level snapshot of Ramanujan-style continued fractions and
       related identities
     - strong reminder that the surrounding literature is still active and not
       frozen at the older notebook / page-43 neighborhood
   - Outcome:
     - no direct exact-pattern hit for
       `(q^(3n) + q^(6n)) / (1 + q^(3n))`
     - but this source widens, rather than closes, the novelty burden

2. McLaughlin, Monks, Reid, 2025
   - INTEGERS article: https://doi.org/10.5281/zenodo.14580993
   - Relevance:
     - recent primary-source work deriving further Ramanujan-style continued
       fraction identities in the Rogers-Ramanujan / Gollnitz-Gordon orbit
     - confirms that "nearby known families" remain an actively developing
       literature area rather than a closed historical list
   - Outcome:
     - no direct exact-pattern hit was recovered in this scan
     - but the existence of recent identities means the project still needs a
       broader bibliography pass before using novelty language

Updated interpretation as of `2026-03-16`:

- `cb60fd71d1d7` remains the strongest current unexplained candidate.
- The exact pattern still has no source hit in the scanned nearby literature.
- That is still **not** enough to call it a new discovery.
- The honest threshold is now better stated as:
  - either identify a final source identity
  - or build a much stronger named-class exclusion plus a substantially wider
    literature closure pass

See also:

- `CB60FD71D1D7_NOVELTY_GATE.md`
