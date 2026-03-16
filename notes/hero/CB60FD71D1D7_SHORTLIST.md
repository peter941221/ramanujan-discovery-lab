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

## Shortlist: Next Sources To Check

The point of this list is not breadth, but "most likely to contain a close cousin"
of a `1 + q^n`-denominator continued fraction with a two-term numerator pattern.

1. Basil Gordon, 1965
   - Exact citation:
     *Some continued fractions of the Rogers-Ramanujan type*,
     Duke Mathematical Journal 32 (1965), 741-748.
   - Goal: find any RR-like continued fractions whose partial denominators are
     `1 + q^n` (or can be equivalence-transformed into that form).

2. Hirschhorn, 1974 and 1980
   - Exact citations:
     - M. D. Hirschhorn, *A continued fraction*,
       Duke Mathematical Journal 41 (1974), 27-33.
     - M. D. Hirschhorn, *A continued fraction of Ramanujan*,
       Journal of the Australian Mathematical Society Series A 29 (1980), 80-86.
   - Goal: locate fractions where numerators involve both `q^n` and `q^(2n)`
     (or a fixed-shift variant that could become our numerator after an index shift).

3. Berndt (Ramanujan notebooks) and Slater-style catalogs
   - Goal: search for continued fraction entries indexed by modulus `3` or
     nearby (since `t=q^3` purity is strong in this case).

4. Ramanujan--Weber / class-invariant papers
   - Goal: check whether `C(t)` can be expressed as a rational function of a
     known Weber/Ramanujan continued fraction (even if not an eta-quotient itself).

## Automation Targets (If A Source Family Looks Promising)

These are concrete expansions to the current tooling that would make a new
"candidate origin" test cheap:

1. Broaden the RR-tower relation search
   - Keep the current reciprocal, ratio-object polynomial, ratio-object
     multiplicative, ratio-object fractional-linear, and two-layer single-basis
     nonlinear scans, then move to source-family-specific relation templates
     before trying larger raw degree boxes.

2. Extend the source-family-specific exact lane
   - The zero-shift `f2/gcf3` and `f4/gcf2` `n`-dependent equivalence
     obstructions are now executable and formalized in Lean.
   - The next proof-facing step is to widen that exact page-43 layer beyond the
     current zero-shift nearest lanes rather than re-formalizing the same case.

3. Add a second-ring Bauer--Muir chain search
   - Current scan uses a tiny fixed modifier family; if a literature family is
     found, implement its exact `w_n` templates and search within that.

4. Add a dedicated `1+q^n`-denominator benchmark family (only with a reliable source)
   - This would let the discovery pipeline classify "Hirschhorn/Gordon-style"
     neighbors early, instead of treating them as unexplained review cases.
