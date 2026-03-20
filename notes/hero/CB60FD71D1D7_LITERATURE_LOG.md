# `cb60fd71d1d7` Literature Log

Audit date: `2026-03-19`

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

## 2026-03-16 Update: Encoded `1+q^n` Literature Bases

This pass turned the older `1 + q^n` denominator references into exact project
benchmarks rather than keeping them only as prose reminders.

1. `S(q)` benchmark
   - Encoded as `hirschhorn_s_normalized`.
   - Continued fraction model:

```text
1 / (1 - q/(1 + q + (-q^3)/(1 + q^2 + (-q^5)/(1 + q^3 + ...))))
```

   - Verified against the product
     `(q^2; q^3)_∞ / (q; q^3)_∞`.

2. `GG(q)` benchmark
   - Encoded as `gollnitz_gordon_normalized`.
   - Normalized object:

```text
1 / (1 + q + q^2/(1 + q^3 + q^4/(1 + q^5 + ...)))
```

   - Verified against the product
     `(q, q^7; q^8)_∞ / ((q^3; q^8)_∞ (q^5; q^8)_∞)`.

3. Outcome for the hero case
   - The identification note now scans the multiplicative correction
     `candidate / RR(q^3)` against the named source-family basis
     `RR`, `cubic`, `GG`, `S`.
   - Result: `0` hits in the current exact multiplicative prefix boxes.

Interpretation:

- This narrows the "simple named source-family product correction" explanation.
- It still does **not** rule out transformed, shifted, reciprocated, or more
  nonlinear Gordon/Hirschhorn presentations.

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

## 2026-03-16 Update: Verified Post-2020 RR Slice

This pass was aimed at replacing vague "recent literature exists" language with
exactly titled, link-verified primary sources. The point was not to inflate the
novelty story; it was to make the novelty burden more honest.

Verified external sources checked this pass:

1. Aricheta, Guadalupe, 2024
   - Title: "A Remark on Modularity of Certain Products of the Rogers-Ramanujan Continued Fraction"
   - arXiv: https://arxiv.org/abs/2405.06678
   - Relevance:
     - recent modular-product work directly in the RR continued-fraction orbit
     - increases the burden on any claim that a nearby product-style identity is new
   - Outcome:
     - no direct exact-pattern hit for
       `(q^(3n) + q^(6n)) / (1 + q^(3n))`
     - but this is exactly the kind of modern source that must be cleared before
       prize-level novelty language is credible

2. Aricheta, Guadalupe, 2024
   - Title: "On Certain Order 10 Modular Functions Involving the Rogers-Ramanujan Continued Fraction"
   - arXiv: https://arxiv.org/abs/2404.05756
   - Relevance:
     - recent modular-function identities and structure around the RR orbit
   - Outcome:
     - no exact coefficient-level source match extracted in this pass

3. Akkarapakam, Morton, 2024
   - Title: "A remark on modular equations involving the Rogers-Ramanujan continued fraction via 5-dissections"
   - arXiv: https://arxiv.org/abs/2410.14149
   - Relevance:
     - modern modular-equation lane very close to the RR identity ecosystem
   - Outcome:
     - no direct hero-case source identity was recovered from the scanned formulas

4. Yamamoto, 2024
   - Title: "Proof and Generalization of Conjectures of Ramanujan Machine"
   - arXiv: https://arxiv.org/abs/2403.09729
   - Relevance:
     - explicit machine-discovery / proof workflow on RR-adjacent identities
     - reminder that discovery-style pipelines can and do land publishable identities in this neighborhood
   - Outcome:
     - no exact hit for the current hero pattern
     - but it raises the standard for how strong the project's own novelty proof must be

5. Baruah, Sarma, 2025
   - Title: "Sign Patterns and Congruences for Certain Infinite Products involving the Rogers-Ramanujan Continued Fraction"
   - arXiv: https://arxiv.org/abs/2503.08517
   - Relevance:
     - recent RR-adjacent product behavior and arithmetic structure
   - Outcome:
     - no direct source match
     - confirms the RR-product side is actively developing

6. Ghoshal, Jana, 2025
   - Title: "Affirmation of Certain Conjectures on the Rogers-Ramanujan Continued Fraction"
   - arXiv: https://arxiv.org/abs/2503.17950
   - Relevance:
     - recent RR continued-fraction conjecture follow-up
   - Outcome:
     - no direct source match
     - again confirms the surrounding literature is current rather than closed

7. Berndt, Rebaka, 2025 survey
   - Title: "The Rogers-Ramanujan continued fraction"
   - arXiv: https://arxiv.org/abs/2512.19952
   - Relevance:
     - recent survey-level snapshot of Ramanujan-style continued fractions and
       related identities
   - Outcome:
     - no direct exact-pattern hit for the hero case
     - but this source widens, rather than closes, the novelty burden

Bibliographic hygiene note:

- A previously mentioned `2025` McLaughlin / Monks / Reid placeholder was not
  kept in the verified set for this pass because its exact title normalization
  was not confirmed from a primary source.
- That lead may still matter later, but it should not be counted as a closed
  citation until the exact paper metadata is pinned down.

Updated interpretation as of `2026-03-16`:

- `cb60fd71d1d7` remains the strongest current unexplained candidate.
- The exact pattern still has no source hit in the scanned nearby literature.
- That is still **not** enough to call it a new discovery.
- The honest threshold is now better stated as:
  - either identify a final source identity
  - or build a much stronger named-class exclusion plus a substantially wider
    literature closure pass

## 2026-03-17 Update: Classical Citation Hygiene + Parameterized Family Pass

This pass did three things at once:

1. It pinned the classical Gordon / Hirschhorn citation spine more explicitly in
   the project notes:
   - Basil Gordon, *Some continued fractions of the Rogers-Ramanujan type*,
     Duke Mathematical Journal 32 (1965), 741-748
   - The Gordon paper now also has a verified DOI-backed landing page:
     https://doi.org/10.1215/S0012-7094-65-03278-3
   - M. D. Hirschhorn, *A continued fraction*,
     Duke Mathematical Journal 41 (1974), 27-33
   - M. D. Hirschhorn, *A continued fraction of Ramanujan*,
     Journal of the Australian Mathematical Society Series A 29 (1980), 80-86
   - The 1974 paper now has a verified DOI-backed landing page:
     https://doi.org/10.1215/S0012-7094-74-04104-0
   - The 1980 paper now also has a verified Cambridge/DOI landing page:
     https://doi.org/10.1017/S1446788700020954

2. It upgraded the source-family recognition tooling from only mixed-family
   boxes to explicit per-family powered ladders:
   - `RR`, `RR2`, `RR3`, `RR4`
   - `cubic`, `cubic2`, `cubic3`, `cubic4`
   - `GG`, `GG2`, `GG3`, `GG4`
   - `S`, `S2`, `S3`, `S4`

3. It paired that recognition pass with a wider theorem-grade page-43
   exclusion layer:
   - the first mixed unit-`a` / unit-`b` exact equivalence lanes for
     `f2/gcf3` and `f4/gcf2` are now ruled out, and both end with the same
     surviving `m^2` coefficient `t` after the forced parameter
     specializations
   - the first mixed unit-`a` / unit-`lambda`, mixed unit-`b` /
     unit-`lambda`, and full mixed unit-`a` / unit-`b` / unit-`lambda`
     exact-equivalence lanes are also now ruled out; all three fail earlier at
     the same surviving `m^1` coefficient `lambda*t^2 - t`

Outcome of this pass:

- The older Gordon / Hirschhorn citation spine is now less vague.
- The project now has DOI-backed landing pages for Gordon 1965, Hirschhorn
  1974, and Hirschhorn 1980, so the classical citation spine is less fragile
  than before.
- The new per-family powered ladders still give `0` hits in both the
  multiplicative and fractional-linear boxes for the `GG` and `S` orbits, so
  the Gordon/Hirschhorn line remains unexplained rather than identified.
- A new explicit direct / reciprocal / quotient template pass inside the same
  `GG/GG2/GG3/GG4` and `S/S2/S3/S4` ladders also gives `0` exact hits, so the
  Gordon/Hirschhorn line still does not collapse to the most obvious quotient
  or reciprocation interpretations either.
- The Weber / class-invariant side remains open as a transform-recognition
  target rather than a recognized source identity.
- The exact page-43 nearest-shift cube is now wider, but it still does **not**
  amount to a full all-parameter theorem for that neighborhood.

## 2026-03-17 Update: Hirschhorn / GG / Notebook Deepening

This pass widened the literature map on the specific orbit that still feels
closest to the hero case: the `GG` / `S` / Ramanujan-notebook neighborhood.

1. Bhatnagar, Ismail, 2019
   - Title: "Orthogonal polynomials associated with a continued fraction of Hirschhorn"
   - arXiv: https://arxiv.org/abs/1901.09985
   - Relevance:
     - treats Hirschhorn's continued fraction as a structural object, not only
       as a product identity
     - includes convergents, generating functions, orthogonality, and the
       Stieltjes transform
     - explicitly notes RR and Ramanujan-generalization special cases
   - Outcome:
     - no exact hero-pattern source identity extracted
     - but this confirms the Hirschhorn orbit has richer analytic
       presentations than the current benchmark/product scans alone

2. Vasuki, Srivatsa Kumar, 2006
   - Title: "Certain identities for Ramanujan-Göllnitz-Gordon continued fraction"
   - DOI: https://doi.org/10.1016/j.cam.2005.03.038
   - Relevance:
     - gives explicit identities relating `H(q)` to `H(q^3)`, `H(q^5)`,
       `H(q^7)`, and `H(q^11)`
     - directly widens the transform/power neighborhood around the encoded
       `GG(q)` benchmark
   - Outcome:
     - no exact coefficient-level hero match extracted
     - but it shows the GG orbit contains nontrivial power/substitution
       identities beyond the simple direct / reciprocal / quotient templates
       already scanned locally

3. Bhatnagar, 2022
   - Title: "Ramanujan's `q`-continued fractions"
   - arXiv: https://arxiv.org/abs/2208.12656
   - Relevance:
     - reorganizes the notebook `q`-continued-fraction landscape rather than
       viewing entries only in notebook order
     - useful as a literature-closure source, not merely as a single-family
       identity paper
   - Outcome:
     - no exact hero-pattern source identity extracted in this pass
     - but it widens the burden for any future novelty claim because the
       notebook-level `q`-continued-fraction literature is broader than the
       currently encoded benchmark families

Updated interpretation after this pass:

- The project's literature map is materially better on the Hirschhorn / GG /
  Ramanujan-notebook side.
- The hero case still does **not** collapse to an obvious `GG` / `S` /
  Hirschhorn power relation or notebook entry.
- The identification note now also includes a per-family low-degree polynomial
  scan over the same powered `RR`, `cubic`, `GG`, and `S` ladders; for
  `cb60fd71d1d7`, those family-preserving polynomial boxes still show `0` hits
  through the checked prefixes.
- The same identification note now also includes a within-family quotient
  ladder `Qk = Tk / T1` over those same powered `RR`, `cubic`, `GG`, and `S`
  ladders; for `cb60fd71d1d7`, those quotient-ladder polynomial /
  multiplicative / fractional-linear boxes still show `0` hits through the
  checked prefixes.
- The same identification note now also includes quotient-ladder two-layer
  fractional-linear boxes over those same powered `RR`, `cubic`, `GG`, and `S`
  ladders; for `cb60fd71d1d7`, those quotient-coordinate nonlinear correction
  boxes also still show `0` hits through the checked prefixes.
- The same identification note now also includes mixed quotient-basis boxes
  that combine the family base object with the quotient ladder over those same
  powered `RR`, `cubic`, `GG`, and `S` ladders; for `cb60fd71d1d7`, those
  mixed quotient-basis polynomial / multiplicative / fractional-linear boxes
  also still show `0` hits through the checked prefixes.
- The same identification note now also includes mixed quotient-basis two-layer
  fractional-linear boxes over those same powered `RR`, `cubic`, `GG`, and `S`
  ladders; for `cb60fd71d1d7`, those stronger base-plus-quotient nonlinear
  correction boxes also still show `0` hits through the checked prefixes.
- On the `GG` side, the family-preserving ladder now explicitly includes the
  literature-motivated powers `GG5`, `GG7`, and `GG11` suggested by the
  Vasuki--Srivatsa Kumar identity paper, and those still give `0` hits in the
  checked polynomial / multiplicative / fractional-linear boxes.
- The Gordon / Hirschhorn citation spine is now slightly less fragile because
  the Gordon 1965, Hirschhorn 1974, and Hirschhorn 1980 papers all now have
  verified DOI-backed landing pages in the notes.
- That is still not literature closure, only a stronger reason to keep the
  public label at "unexplained candidate."

See also:

- `CB60FD71D1D7_EXACT_SUBSEQUENCE_OBSTRUCTION.md`
- `CB60FD71D1D7_BIBLIOGRAPHY_MATRIX.md`
- `CB60FD71D1D7_NOVELTY_GATE.md`

## 2026-03-18 Update: GG Modular-Equation Deep Research

This pass was narrower than the previous bibliography broadening.

The goal here was not to add more rows for their own sake, but to answer one
practical question:

```text
What should the next local recognition lane be?
```

### Sources checked in this pass

1. Chan, Huang, 1997
   - Title: "On the Ramanujan-Göllnitz-Gordon continued fraction"
   - Link: https://mrc.sdu.edu.cn/ziliao/8.pdf
   - Relevance:
     - treats `H(q)` as more than a product identity
     - derives exact identities tying `H(q)` to `H(-q)` and `H(q^2)`
     - explains how modular equations generate relations between `H(q)` and
       `H(q^n)`, with explicit attention to the `q^3` and `q^4` lanes
   - Outcome:
     - no final hero identity extracted directly
     - but this is the clearest source-level evidence so far that the current
       local `GG` scans are still missing an important **structured** lane

2. Cho, Koo, Park, 2009
   - Title: "Arithmetic of the Ramanujan-Göllnitz-Gordon continued fraction"
   - DOI: https://doi.org/10.1016/j.jnt.2008.09.018
   - Relevance:
     - extends the modular-equation story to all odd primes
     - treats the GG object as a modular unit with arithmetic structure
     - suggests that low-complexity source-faithful modular relations may exist
       far beyond the current direct / reciprocal / quotient boxes
   - Outcome:
     - no exact hero-pattern source identity extracted in this pass
     - but this pushes the next implementation target toward a
       modular-equation-aware `GG` lane rather than another anonymous basis
       expansion

3. Yuttanan, 2012
   - Title: "New properties for the Ramanujan-Göllnitz-Gordon continued fraction"
   - Link: https://eudml.org/doc/278973
   - Relevance:
     - adds later exact structural properties around the same GG object
     - confirms that the orbit still carries nontrivial coefficient and
       evaluation structure beyond the earlier identity papers
   - Outcome:
     - no direct hero identification extracted in this pass
     - but it weakens any temptation to treat the GG orbit as "already
       exhausted" by the current product / quotient / ladder scans

### Updated interpretation after this pass

- The strongest next local move is now clearer:
  build a literature-driven modular-equation recognition lane for the `GG`
  orbit.
- That lane should preserve source meaning explicitly:
  `GG(q)`, `GG(-q)`, `GG(q^2)`, `GG(q^3)`, `GG(q^4)`, and then odd-prime
  descendants only if needed.
- The current negative `GG` results are still valuable, but they now look more
  like a statement about the **limitations of the current boxes** than about
  the full exhaustion of the Gordon / Hirschhorn orbit.
- The Weber / class-invariant side remains open, but the first implementation
  priority is now the `GG` modular-equation lane because the literature there is
  more explicit and more immediately codable.

### Practical next action

The next implementation pass should **not** start by widening the generic basis
catalog again.

It should start by encoding a small source-faithful modular template family
whose variables are drawn directly from the Chan--Huang / Cho--Koo--Park
objects.

## 2026-03-19 Update: Tail-Family-First `GG` Quotient Coordinates

This pass asked a more local question:

```text
if the reduced coefficients force a stationary tail family T(x),
do the first concrete tail objects U(t^2), U(t^3), U(t^4)
or their gap-normalized residuals fall into the exact GG quotient coordinates
highlighted by Chan--Huang?
```

### Local objects checked

- exact tail-family ladder:
  - `U_t2 = T(t^2)/(1+t^2)`
  - `U_t3 = T(t^3)/(1+t^3)`
  - `U_t4 = T(t^4)/(1+t^4)`
- gap-normalized residual ladder through depth `3`
- explicit `GG` quotient coordinates:
  - `Q_neg = GG(-t) / GG(t)`
  - `Q_2 = GG(t^2) / GG(t)`
  - `Q_3 = GG(t^3) / GG(t)`
  - `Q_4 = GG(t^4) / GG(t)`

### Outcome

- all `12` sampled tail-family objects still give `0` hits in the current
  `GG/Weber modular-equation` lane
- the narrowest exact quotient-coordinate sublane also gives `0` hits:
  - no exact Chan--Huang hit on `Q_3`
  - no exact Chan--Huang hit on `Q_4`
- the note now records obstruction witnesses rather than only a yes/no result:
  each sampled tail object reports the first visible residual coefficient where
  the exact `Q_3` / `Q_4` modular-equation template fails

### Updated interpretation

- this does **not** weaken the `GG` direction
- it does weaken the idea that the next gain comes from a broader `GG` prefix
  box
- the current best local move is now:

```text
keep the GG / Weber orbit
but narrow the next lane toward exact quotient-coordinate elimination
```

In other words:

- stay in the same mathematical neighborhood
- switch from a wider hallway scan
- to a smaller number of more theorem-shaped locks

### Weighted modular-curve follow-up

The new weighted local diagnostic is:

```text
W_34 = Q_3^3 / Q_4^2
```

This still does **not** identify the tail-family ladder, but it sharpens the
literature reading:

- Chan--Huang gives exact quotient coordinates `Q_3`, `Q_4`
- Cho--Koo--Park treats the GG object through arithmetic / modular-curve style
  structure at odd-prime levels
- computational modular-polynomial notes around the same orbit suggest that the
  right coordinates can require a level shift before they become clean modular
  functions, rather than appearing directly in the raw `Q_3`, `Q_4` pair

So the best literature-driven next step is no longer “search wider over GG
prefixes,” but “look for a deeper Weber / modular-curve Hauptmodul-like
parameter compatible with the observed `3:2` weighting.”

## 2026-03-19 Update: Latest Literature Cross-Check On Direction

This pass asked a simpler strategic question:

```text
after the new three-lane implementation
(direct modular-unit / eta, periodic-point, tail-operator),
does the latest primary literature support the current direction,
or should the project pivot?
```

### Latest primary-source signals

1. RR modular-unit / eta signal
   - Russelle Guadalupe, 2024 / 2025:
     *Modularity of certain products of the Rogers-Ramanujan continued fraction*
     (`https://arxiv.org/abs/2405.06678`)
   - Main signal:
     products built from `r(t)` and `r(2t)` can generate full modular-function
     fields and are proved with `eta`-quotient / generalized `eta`-quotient
     methods.
   - Direction implication:
     the direct modular-unit / eta lane is not an arbitrary box; it matches the
     way recent RR literature actually compresses continued fractions into
     theorem-grade modular objects.

2. RR modular-function / modular-equation signal
   - Russelle Guadalupe, 2024 / 2026:
     *Ramanujan's continued fractions of order 10 as modular functions*
     (`https://arxiv.org/abs/2404.05756`)
   - Nayandeep Deka Baruah, Pranjal Talukdar, 2024:
     *Identities for the Rogers-Ramanujan Continued Fraction*
     (`https://arxiv.org/abs/2410.17110`)
   - Russelle Guadalupe, 2024:
     *A remark on modular equations involving the Rogers-Ramanujan continued
     fraction via 5-dissections* (`https://arxiv.org/abs/2410.14149`)
   - Main signal:
     recent RR papers keep pushing toward exact modular equations among
     explicitly named coordinates such as `R(q)`, `R(q^2)`, `R(q^3)`, `R(q^4)`,
     `R(q^6)`, `R(q^8)`, `R(q^16)`, `R(q^20)`, rather than toward anonymous
     low-degree fit boxes.
   - Direction implication:
     keep explicit RR / GG quotient-coordinate and modular-curve lanes active;
     do not let broad generic prefix growth become the main bet again.

3. Ramanujan periodic-point / algebraic-function signal
   - Akkarapakam, Morton, 2024:
     *Periodic points of algebraic functions related to a continued fraction of
     Ramanujan* (`https://nyjm.albany.edu/j/2024/30-36.html`)
   - Main signal:
     a Ramanujan continued fraction can be characterized through periodic points
     of a fixed algebraic function.
   - Direction implication:
     the new Morton-inspired periodic-point lane is literature-faithful and
     worth keeping, but its current `0` hits on the sampled tail-family ladder
     mean it is a diagnostic obstruction lane for now, not the main positive
     recognition trunk.

4. Machine-discovery proof-method signal
   - Yamamoto, 2024:
     *Proof and generalization of conjectures of Ramanujan Machine*
     (`https://arxiv.org/abs/2403.09729`)
   - Chao Wang, 2026:
     *A Formal Proof of a Continued Fraction Conjecture for π Originating from
     the Ramanujan Machine* (`https://arxiv.org/abs/2601.08461`)
   - Chao Wang, 2026:
     *Analytic Proof of a Quartic Continued Fraction Identity for 8/pi^2 via
     Operator Decoupling* (`https://arxiv.org/abs/2602.03027`)
   - Main signal:
     recent proof workflows do use recurrences, differential equations,
     equivalence transformations, and operator decompositions, but usually after
     the right analytic kernel or structural source object is already in hand.
   - Direction implication:
     the tail-operator lane is still worth keeping, but mainly as a theorem /
     proof scaffold; it should not replace source-recognition lanes as the main
     local priority after the current `0`-hit affine q-difference pass.

### Updated direction verdict

The latest literature does **not** support a pivot away from the current
award-track road.

It supports a more precise ordering:

1. keep direct modular-unit / eta recognition as the main positive-recognition
   trunk
2. keep exact RR / GG quotient-coordinate and modular-curve lanes as the main
   named-family refinement
3. keep periodic-point / algebraic-function lanes as literature-driven
   diagnostic obstructions
4. keep operator lanes as theorem / proof scaffolding, not as the main
   replacement for source recognition

### 2026-03-19 Addendum: Published-status cross-check after `G2_W34`

The latest local result on the deeper weighted `GG` lane is:

```text
F / W_34 - 1 first fails at t^1 with coefficient -1
G_W34 = (1 - F / W_34) / t
G2_W34 = (G_W34 - 1) / (-4*t^2)
```

and both `G_W34` and `G2_W34` still give `0` hits in the checked small
eta-quotient, modular-unit / eta, and one-core `RR/GG` source-family
eta-correction boxes.

That result is still consistent with the newest primary-source picture:

- Guadalupe's `2405.06678` was revised on `2025-09-16` and published in
  `Ramanujan J. 68 (2025)`; it still points toward eta-quotient /
  generalized-eta modular recognition rather than anonymous fit boxes
- Aricheta--Guadalupe's `2404.05756` was revised on `2025-03-13` and published
  in `J. Number Theory 278 (2026)`; it still points toward named modular
  coordinates and modular equations at many levels
- Deka Baruah--Talukdar's `2410.17110` (submitted `2024-10-22`) still points
  toward explicit RR coordinate identities among named powers and levels

So the literature-driven reading stays the same:

- do **not** pivot away from modular-unit / eta plus named coordinate lanes
- do **not** promote the current operator lane above source recognition
- do move the next positive-recognition attempt deeper than the first weighted
  coordinate `W_34` and its first two normalized corrections `G_W34`,
  `G2_W34`

## 2026-03-20 Update: Phase-1 Literature Exit And The Weber-Schlafli Hand-Off

This pass asked the exact Phase-1 question from the execution plan:

```text
after freezing the local tail-family conclusion,
which primary-source spine most sharply determines the next coding target?
```

### Verified source spine

1. Chan, Huang, 1997
   - link:
     `https://mrc.sdu.edu.cn/ziliao/8.pdf`
   - verified role:
     the exact `GG` quotient-coordinate source spine around the `q^3` and
     `q^4` lanes.

2. Cho, Koo, Park, 2009
   - link:
     `https://doi.org/10.1016/j.jnt.2008.09.018`
   - verified role:
     extends the same orbit to all odd primes by computing affine models of
     modular curves `X(Γ)` with `Γ = Γ_1(8) ∩ Γ_0(16p)`, and treats `v(τ)` as a
     modular unit over `Z`.
   - direction consequence:
     the next lane should look like a modular-curve coordinate lane, not an
     anonymous prefix fit.

3. Akkarapakam, Morton, 2024
   - link:
     `https://nyjm.albany.edu/j/2024/30-36.html`
   - verified role:
     gives a deeper Weber-Schlafli coordinate change for the relevant
     continued-fraction object.
   - the key exact coordinate relation recorded there is:

```text
2*p(8τ) = 1/v(τ) - v(τ)
```

   - the same paper then gives low-degree exact relations in the deeper
     coordinates, including:

```text
p(τ)^2 p(2τ)^2 + p(τ)^2 - 2 p(2τ) = 0
b(τ)^4 = p(τ)^8 + 16 p(τ)^4
b(τ/2)^2 = b(τ) + 4
```

   - direction consequence:
     after the misses on the raw `Q_3`, `Q_4`, `W_34`, `G_W34`, and `G2_W34`
     lane, the first literature-backed next coordinate should be `p(8τ)`, not
     another raw quotient combination.

4. Yui, Zagier, 1997
   - link:
     `https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1090/S0025-5718-97-00854-5/fulltext.pdf`
   - verified role:
     classical Weber modular-function singular-value anchor behind the later
     modular-curve and class-invariant story.
   - direction consequence:
     moving from raw `GG` quotient coordinates to Weber-style modular-function
     coordinates is a source-backed shift, not a random coordinate trick.

5. Berndt, Chan, Zhang, 1997
   - link:
     `https://mrc.sdu.edu.cn/ziliao/10.pdf`
   - verified role:
     explicit Ramanujan-Weber / classical-Weber dictionary:
     Ramanujan-Weber `G_n` / `g_n` are identified there with the classical
     Weber `f` / `f1` normalization.
   - direction consequence:
     once paired with Yui--Zagier's classical Weber `f`, `f1`, `f2` trio,
     the current `g12_ws` / `p12_ws` / `G_f2_ws` shell is no longer just a
     suggestive algebra gadget; it becomes a source-backed classical Weber
     trio in project normalization.

### Phase-1 exit conclusion

Phase 1 now exits with a sharper local target:

```text
keep the GG / Weber neighborhood
stop widening anonymous boxes
move next to the Weber-Schlafli coordinate lane
```

Concretely, the next implementation target should start with:

```text
P = p(8τ), where 2P = 1/v - v
```

and only pull in the companion coordinate

```text
B = b(4τ)
```

if the one-coordinate `P` lane still misses.

### Why this beats another quotient-box widening step

Analogy:

- the raw `Q_3 / Q_4` lane was like checking the front-facing symmetry of a
  machine
- the Weber-Schlafli lane is like opening the casing and reading the gearbox
  coordinates directly

The current literature now says the gearbox view is the more faithful next
representation.

## 2026-03-20 Update: Classical Weber Trio Naming Closure

This pass asked a naming-closure question rather than a new scan question:

```text
after the Ramanujan-Weber class-invariant compression,
can the current `g12_ws / p12_ws / G_f2_ws` shell be named from sources
rather than only treated as an internally useful coordinate package?
```

Primary source spine used for this closure:

1. Berndt, Chan, Zhang, 1997
   - role:
     identifies Ramanujan-Weber `G_n` / `g_n` with the classical Weber
     `f` / `f1` normalization.

2. Yui, Zagier, 1997
   - role:
     supplies the classical Weber `f`, `f1`, `f2` trio as the ambient
     modular-function source language.

### Outcome

- the current `g12_ws` / `p12_ws` / `G_f2_ws` shell should now be read as a
  named classical Weber trio in project normalization, not as an anonymous
  product gadget
- correspondingly, the branch
  `Q_gp_ws -> X_g_ws -> G_X_ws`
  is now best described as the source-faithful classical Weber
  quotient/template lane
- this is a naming and closure gain, not a positive final identity:
  the product follow-up branch
  `G_f2_ws -> H_f2_ws`
  still stays negative in the first theorem-shaped boxes, and the quotient
  branch remains the cleaner constructive trunk

### Direction consequence

The next move inside the same orbit should still prioritize:

- `X_g_ws`
- `G_X_ws`
- `H_X_ws`
- the later quotient-return bridge

rather than pivoting the whole search toward the classical product lane alone.

## 2026-03-20 Update: Companion `B = b(4τ)` Pass Closes Too

The first `P = p(8τ)` pass was already enough to justify staying in the named
`GG / Weber` orbit.

This follow-up asks the next natural literature-backed question:

```text
if P_ws = (1/Y - Y) / 2 still misses,
does the direct companion coordinate B = b(4τ)
stabilize the same tail-family ladder?
```

The same Akkarapakam--Morton formulas give the coding bridge:

```text
b(τ)^4 = p(τ)^8 + 16 p(τ)^4
b(τ/2)^2 = b(τ) + 4
```

so the local companion coordinate is:

```text
B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)
```

and the first exact companion templates are:

```text
B_ws^2 - B_ws,2 - 4 = 0
B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0
```

Current exact outcome:

- all `12` sampled tail-family objects still give `0` hits
- unlike the `P_ws` lane, the `B_ws` lane does not spread into several
  low-order classes
- instead, it collapses immediately to one universal constant-term
  obstruction:

```text
B_ws^2 - B_ws,2 - 4 -> (t^0, -2)
B_ws,2^4 - P_ws^8 - 16*P_ws^4 -> (t^0, 16)
```

Direction consequence:

- this keeps the same named `GG / Weber` orbit active
- but it says the direct `P -> B` companion closure is already exact-blocked
- so the next source-faithful step should be another named Weber coordinate or
  a theorem-grade obstruction package, not a wider anonymous scan

## 2026-03-20 Update: Return-Bridge Literature Closure On `Q_XK_ws` / `L_XK_ws`

This pass asked a narrower closure question inside the same named orbit:

```text
after the Weber class-invariant compression and the direct return bridge
K_XR_ws -> H_X_ws,
do the derived return-bridge objects Q_XK_ws or L_XK_ws
fall back into the named GG modular-equation coordinates?
```

The exact local objects checked were:

```text
D_XK_ws = K_XR_ws - H_X_ws
Q_XK_ws = K_XR_ws / H_X_ws
L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)
```

and the source-faithful comparison box was the same literature-driven `GG`
modular-equation lane built from the named Chan--Huang coordinates:

```text
GG(t), GG(-t), GG(t^2), GG(t^3), GG(t^4)
```

together with their direct, quotient, and mixed quotient-coordinate templates.

### Outcome

- `Q_XK_ws` still gives `0` hits in the checked named `GG` modular-equation
  box:
  - exact templates: `0`
  - exact polynomial templates: `0`
  - quotient exact polynomial templates: `0`
  - direct polynomial / multiplicative / fractional-linear / two-layer boxes:
    all `0`
  - quotient polynomial / multiplicative / fractional-linear / two-layer
    boxes: all `0`
  - mixed quotient polynomial / multiplicative / fractional-linear / two-layer
    boxes: all `0`
- `L_XK_ws` gives the same verdict: every checked named `GG` modular-equation
  box remains at `0` hits

### Updated interpretation

- this does **not** say the `GG / Weber` orbit is exhausted
- it does say the current return-bridge objects are still best read as
  **derived bridge coordinates**, not as recovered named literature
  coordinates
- in particular, the honest label for `Q_XK_ws` and `L_XK_ws` is now:

```text
derived return bridge
```

rather than:

```text
identified GG modular-equation coordinate
```

### Practical consequence

The next source-faithful move should stay in the same named orbit but push one
step deeper:

- either another named Weber coordinate
- or a theorem-grade obstruction package that explains why the return bridge
  fails to close back inside the current named `GG` coordinates
