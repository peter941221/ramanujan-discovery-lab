# `cb60fd71d1d7` Shortlist (Identity Hunting)

This is a deliberately narrow checklist for the next literature + recognition pass.
It is **not** a novelty claim.

## Working Object

In step-reduced variable `t = q^3`, the reciprocal continued fraction is

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n).
```

The closest built-in benchmark reciprocal is the Rogers-Ramanujan ladder

```text
B1(t) = 1 + K_{n>=1} t^n / 1.
```

## What We Already Know (Internal Audit Summary)

- `cb60fd71d1d7` matches `RR(q^3)` through `q^9` and diverges first at `q^12`.
- No direct constant-parameter specialization was found in the nearest page-43
  ratio families or nearby Bauer--Muir chain families (see `CB60FD71D1D7_TRANSFORM_AUDIT.md`).
- A small-degree algebraic relation guess between the reciprocals `C(t)` and
  `B1(t)` did not find any `P(C,B1)=0` up to total degree `4` modulo `t^90`
  and the current low-degree RR-tower prefix scan also found no
  candidate-dependent relation using
  `B1(t^2), B1(t^3), B1(t^4), B1(t^5), B1(t^6), B1(t^12), B1(t^20)`
  (see `CB60FD71D1D7_IDENTIFICATION_NOTE.md`).
- The same note now also records that the multiplicative correction object
  `candidate / RR(q^3)` has no candidate-dependent RR-tower prefix relation in
  the current low-degree (`1` / `2`) search boxes, and no exact small-integer
  multiplicative RR-tower product relation
  `F = prod_i B_i^(e_i)` was found in the current benchmark-prefix scan either.
- The named source-family multiplicative scan for
  `F = candidate / RR(q^3)` also reports `0` hits in the current prefix boxes
  built from `RR`, `cubic`, `GG`, and `S`, so wiring the Gordon/Hirschhorn line
  into the benchmark family did not collapse the hero case to an obvious
  literature basis product.
- The richer named source-family fractional-linear scan
  `F = (1 + sum a_i (S_i - 1)) / (1 + sum b_i (S_i - 1))` also reports `0`
  hits through the same `RR/cubic/GG/S` prefix order.
- The second-ring named source-family two-layer scan also reports `0` hits for
  prefixes ending at `cubic`, `GG`, and `S`.
- A new per-family powered-ladder pass also reports `0` hits in low-degree
  polynomial, multiplicative, fractional-linear, quotient-ladder, and
  quotient-ladder two-layer, mixed quotient-basis, and mixed quotient-basis
  two-layer boxes for
  `RR/RR2/RR3/RR4`, `cubic/cubic2/cubic3/cubic4`,
  `GG/GG2/GG3/GG4/GG5/GG7/GG11`, and `S/S2/S3/S4`.
- A new explicit `GG` / `S` transform-template pass also reports `0` exact hits
  for direct objects, reciprocals, and pairwise quotients inside the same
  `GG/GG2/GG3/GG4` and `S/S2/S3/S4` ladders.
- The exact page-43 equivalence layer now rules out the zero-shift, the first
  unit-`a`, the first unit-`b`, and the first unit-`lambda` nearest
  `f2/gcf3` and `f4/gcf2` lanes; the unit-`a`, unit-`b`, and first mixed
  unit-`a` / unit-`b` lanes still end with a surviving `m^2` term `t`, while
  the unit-`lambda`, mixed unit-`a` / unit-`lambda`, and mixed unit-`b` /
  unit-`lambda`, and mixed unit-`a` / unit-`b` / unit-`lambda` lanes already
  fail at the surviving `m^1` term
  `lambda*t^2 - t`.
- Taken together, those eight nearby exact lanes now close the full
  nearest-shift cube for the current page-43 equivalence audit.

## Shortlist: Next Source Lanes To Deepen

The point of this list is not breadth, but "most likely to contain a close cousin"
of a `1 + q^n`-denominator continued fraction with a two-term numerator pattern.

1. Basil Gordon, 1965
   - Exact citation:
     *Some continued fractions of the Rogers-Ramanujan type*,
     Duke Mathematical Journal 32 (1965), 741-748.
   - DOI: https://doi.org/10.1215/S0012-7094-65-03278-3
   - Current state: the nearby Gordon/Hirschhorn orbit is now partially encoded
     in the benchmark layer through `GG(q)` and `S(q)`.
   - Goal: revisit Gordon-style variants, reciprocal forms, and transform
     presentations around those base objects rather than just re-encoding the
     base objects themselves.

2. Hirschhorn, 1974 and 1980
   - Exact citations:
     - M. D. Hirschhorn, *A continued fraction*,
       Duke Mathematical Journal 41 (1974), 27-33.
     - M. D. Hirschhorn, *A continued fraction of Ramanujan*,
       Journal of the Australian Mathematical Society Series A 29 (1980), 80-86.
     - Bhatnagar, Ismail, *Orthogonal polynomials associated with a continued
       fraction of Hirschhorn*, arXiv: https://arxiv.org/abs/1901.09985
   - Current state: `S(q)` is now a verified normalized project benchmark, and
     the hero ratio object still has no simple multiplicative scan hit against
     the `RR/cubic/GG/S` basis.
   - Goal: look for shifted, reciprocated, or equivalence-transformed
     Hirschhorn-family presentations that could still hide the hero pattern,
     especially analytic / orthogonal-polynomial presentations that are not
     obvious from the product formulas alone.

3. Berndt (Ramanujan notebooks) and Slater-style catalogs
   - Goal: search for continued fraction entries indexed by modulus `3` or
     nearby (since `t=q^3` purity is strong in this case).

4. Ramanujan--Weber / class-invariant papers
   - Goal: check whether `C(t)` can be expressed as a rational function of a
     known Weber/Ramanujan continued fraction (even if not an eta-quotient itself).

5. GG power-relation papers
   - Exact citation:
     - Vasuki, Srivatsa Kumar, *Certain identities for
       Ramanujan-Göllnitz-Gordon continued fraction*,
       DOI: https://doi.org/10.1016/j.cam.2005.03.038
     - Chan, Huang, *On the Ramanujan-Göllnitz-Gordon continued fraction*,
       https://mrc.sdu.edu.cn/ziliao/8.pdf
     - Cho, Koo, Park, *Arithmetic of the Ramanujan-Göllnitz-Gordon continued fraction*,
       DOI: https://doi.org/10.1016/j.jnt.2008.09.018
   - Goal: test whether the hero case can hide in a nontrivial powered /
     substituted `GG` presentation before moving to even larger transform boxes.
   - Current scan status: the explicit literature-motivated powers
     `q^5`, `q^7`, and `q^11` are now in the family-preserving `GG` ladder, and
     still give `0` hits in the checked polynomial / multiplicative /
     fractional-linear / quotient-ladder / quotient-ladder two-layer /
     mixed quotient-basis / mixed quotient-basis two-layer boxes.
     A newer source-faithful modular-equation box now also keeps
     `GG(t)`, `GG(-t)`, `GG(t^2)`, `GG(t^3)`, `GG(t^4)`, `GG(t^5)`,
     `GG(t^7)`, `GG(t^11)` and the quotient coordinates against `GG(t)`
     explicit, and that lane still reports `0` hits in its exact-template,
     polynomial, multiplicative, fractional-linear, two-layer, and
     quotient-coordinate scans.

6. Page-43 substitutions with explicit source-family meaning
   - Goal: move beyond the currently formalized zero-shift, single-parameter,
     and mixed two-parameter nearest lanes and test the next low-complexity
     parameter shifts or transform templates that still preserve interpretable
     `f2/gcf3` or `f4/gcf2` provenance.

## Automation Targets (If A Source Family Looks Promising)

These are concrete expansions to the current tooling that would make a new
"candidate origin" test cheap:

1. Broaden the RR-tower relation search
   - Keep the current reciprocal, ratio-object polynomial, ratio-object
     multiplicative, ratio-object fractional-linear, two-layer single-basis
     nonlinear, named source-family scans, and the new per-family powered
     source-family scans.
   - The direct / reciprocal / quotient `GG` / `S` template box is now covered;
     the low-degree per-family polynomial box, the within-family quotient
     ladder, the mixed quotient-basis box, and its two-layer follow-up are now
     covered too;
     the next jump should be richer transform templates for those same
     literature families before trying larger raw degree boxes.
   - On the `GG` side specifically, the next best jump is no longer another
     anonymous mixed-family box; it is an odd-prime modular-equation lane that
     preserves `GG(t^p)` / `GG(t)` quotient meaning for the primes already
     emphasized by the literature.

2. Extend the source-family-specific exact lane
   - The zero-shift, first unit-`a`, first unit-`b`, first mixed unit-`a` /
     unit-`b`, first mixed unit-`a` / unit-`lambda`, first mixed unit-`b` /
     unit-`lambda`, first mixed unit-`a` / unit-`b` / unit-`lambda`, and first
     unit-`lambda` nearest `f2/gcf3` and `f4/gcf2` `n`-dependent equivalence
     obstructions are now executable and formalized in Lean.
   - The next proof-facing step is to widen that exact page-43 layer beyond the
     currently formalized nearest lanes, ideally by proving an all-parameter
     theorem for a natural page-43 neighborhood instead of just one more
     bounded nearby specialization.

3. Add a second-ring Bauer--Muir chain search
   - Current scan uses a tiny fixed modifier family; if a literature family is
     found, implement its exact `w_n` templates and search within that.

4. Keep the `1+q^n`-denominator benchmark family literature-backed
   - `GG(q)` and `S(q)` are now wired as normalized benchmarks with verified
     product formulas.
   - Only add more Gordon/Hirschhorn-style bases when the exact continued
     fraction and product formula are both sourced confidently.
