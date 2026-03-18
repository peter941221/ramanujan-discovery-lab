# Case Study: `cb60fd71d1d7` (Audit-First Experimental Mathematics)

This note is written as a "public-facing but conservative" narrative of what the
project can already do on a single hard example, without claiming a new identity.

## The Candidate

The current hero-case candidate is

```text
1 / (1 + (q^3 + q^6)/(1 + q^3 + (q^6 + q^12)/(1 + q^6 + ...))).
```

It was produced by the discovery pipeline as a numerically stable template whose
nearest built-in benchmark is the Rogers-Ramanujan family at `q^3`.

## Why It Is Interesting (But Not Yet "A New Formula")

- Symbolically, it matches `RR(q^3)` through `q^9`.
- The first divergence from `RR(q^3)` is at `q^12`.
- It survives several rings of "obvious" explanations:
  constant-parameter specializations of the nearest published families,
  direct odd/even contractions, and short constrained Bauer--Muir chains.

Those eliminations do not prove novelty. They just justify why the candidate is
worth treating as a high-value audit target.

## The Clean Step-Reduced Model

In reduced variable `t = q^3`, it is natural to study the **reciprocal**
continued fraction (the `1 + ...` object):

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n).
```

The nearest benchmark reciprocal is the Rogers-Ramanujan ladder:

```text
B1(t) = 1 + K_{n>=1} t^n / 1.
```

## What Is Already Exact (Machine-Checkable)

The proof layer currently focuses on theorem-grade local statements:

1. Finite-truncation convergent recurrences for generalized continued fractions.
2. Exact convergent-factor reduction:
   the `n`th convergent numerator and denominator of `C(t)` both carry factor
   `1 + t^n`, inducing a second "reduced-by-factor" continued fraction.
3. Equivalence transformations:
   the reverse scales needed to reconstruct the original coefficients from the
   reduced-by-factor model can be expressed in a fraction-field coefficient layer.
4. Local obstruction lemmas:
   several nearby transform/contraction origins are ruled out at the coefficient level.
5. Award-track waypoint packaging:
   the Lean scaffold `proofs/Proofs/HeroCaseFinalIdentity.lean` now packages
   the current exact waypoint as named theorems and a reusable
   `currentExactWaypointCertificate`, instead of acting as a purely
   comment-only placeholder.

The repository's Lean files reflect this scope (finite convergents and exact
coefficient equalities), not a final analytic identity.

## What Has Been Tried For Identification (Negative Results)

Three recognition layers have now been run on the step-reduced reciprocals or
their ratio objects:

1. Euler-product fingerprints:
   the ratio `candidate / RR(q^3)` and the reciprocal `C(t)` both yield dense,
   rapidly growing Euler-product exponents, making a small eta/Pochhammer tweak
   explanation unlikely in the searched boxes.
2. Generic RR-tower relation boxes:
   no polynomial relation `P(C,B1)=0` was found up to total degree `4`, checked
   modulo `t^90`; the broader RR-tower pass still found no candidate-dependent
   low-degree relation in the scanned prefixes built from
   `B1(t^2), B1(t^3), B1(t^4), B1(t^5), B1(t^6), B1(t^12), B1(t^20)`.
   The same low-degree prefix scan also found no candidate-dependent relation
   for the multiplicative correction object `candidate / RR(q^3)` against the
   same RR tower, and the exact multiplicative ansatz
   `F = prod_i B_i^(e_i)` also has no small-integer hit in the current prefix
   scan. Even the richer structured box
   `F = (1 + sum a_i (B_i - 1)) / (1 + sum b_i (B_i - 1))` still gives no hit
   through the same RR-tower prefixes. A second-ring nonlinear pass with two
   single-basis fractional-linear factors still gives no hit there either.
3. Named source-family and exact nearby-lane checks:
   the older `1 + q^n` denominator literature is now represented directly by
   the normalized source-family bases `RR`, `cubic`, `GG`, and `S`, but the
   exact multiplicative correction box
   `F = prod_i S_i^(e_i)` still has no hit for the ratio
   `F = candidate / RR(q^3)`. The richer source-family fractional-linear box
   `F = (1 + sum a_i (S_i - 1)) / (1 + sum b_i (S_i - 1))` also has no hit in
   the current prefix scan through `RR`, `cubic`, `GG`, and `S`. Even the
   second-ring source-family two-layer box built from two single-basis
   fractional-linear factors still has no hit through the prefixes ending at
   `cubic`, `GG`, and `S`. A new per-family powered-ladder pass also reports
   `0` hits in low-degree polynomial, multiplicative, fractional-linear,
   within-family quotient-ladder, quotient-ladder two-layer, mixed
   quotient-basis, and mixed quotient-basis two-layer boxes for
   `RR/RR2/RR3/RR4`, `cubic/cubic2/cubic3/cubic4`,
   `GG/GG2/GG3/GG4/GG5/GG7/GG11`, and `S/S2/S3/S4`, so the
   Gordon/Hirschhorn line still does not explain the hero ratio object even
   before families are mixed together.
   On the page-43 side, the zero-shift `f2/gcf3` and `f4/gcf2`
   `n`-dependent equivalence lanes are ruled out exactly: in both cases the
   necessary residual polynomial forces `a = b = 0`, then `lambda = 1`, but
   still leaves a surviving `m^2` coefficient `t`. The first nearby
   unit-`a`-shift lanes also fail exactly and end with the same surviving
   `m^2` coefficient `t` after the forced parameter specializations. The first nearby
   unit-`b`-shift lanes now fail in exactly the same way. The first mixed
   unit-`a` / unit-`b` nearby lanes also fail with that same final surviving
   `m^2` coefficient `t`. The first mixed unit-`a` / unit-`lambda` and
   unit-`b` / unit-`lambda` nearby lanes fail earlier, with the same surviving
   `m^1` coefficient `lambda*t^2 - t`. The first full mixed
   unit-`a` / unit-`b` / unit-`lambda` nearby lanes also fail at that same
   surviving `m^1` coefficient. Taken together, these eight nearby exact cases
   now close the full nearest-shift cube for the current page-43 equivalence
   audit. The first nearby
   unit-`lambda`-shift lanes fail even earlier: after the same forced
   specialization, the surviving `m^1` coefficient becomes `lambda*t^2 - t`,
   so no constant `lambda` can make those lanes vanish identically.
   A newer literature-driven `GG` modular-equation lane also now keeps the
   sign and substitution structure explicit instead of flattening everything
   into anonymous basis elements: it scans the basis
   `GG(t)`, `GG(-t)`, `GG(t^2)`, `GG(t^3)`, `GG(t^4)`, `GG(t^5)`,
   `GG(t^7)`, `GG(t^11)` together with the quotient coordinates
   `GG(-t)/GG(t)`, `GG(t^2)/GG(t)`, `GG(t^3)/GG(t)`, `GG(t^4)/GG(t)`,
   `GG(t^5)/GG(t)`, `GG(t^7)/GG(t)`, `GG(t^11)/GG(t)`, and the mixed
   quotient-coordinate prefixes that keep `GG(t)` explicit while allowing
   corrections in those quotient variables.
   That same lane now also checks the first exact Chan--Huang
   modular-equation polynomials in the `q^3` and `q^4` lanes, both in direct
   form against `GG(t^3)`, `GG(t^4)` and in quotient-coordinate form against
   `GG(t^3)/GG(t)`, `GG(t^4)/GG(t)`.
   That more source-faithful box still gives `0` hits in its exact template,
   low-degree polynomial, multiplicative, fractional-linear, two-layer
   fractional-linear, quotient-coordinate, mixed quotient-coordinate, and
   exact Chan--Huang direct / quotient polynomial scans.

These are not impossibility theorems. They only remove a class of "too simple"
closed forms.

## What Would Make This Publishable As A New Identity

To cross the bar for a genuine new-formula paper, we would need one of:

1. A sourced identity statement: `C(t) = known_object(t)` (or `candidate = ...`)
   with a reliable reference chain, plus proof.
2. A uniquely characterizing theorem that implies a new identity:
   for example, a functional equation + normalization that pins down the function.

At the moment, we have strong audit evidence and machine-checked local structure,
but no final source object.

## What Is Publishable Today (More Realistically)

Even without a final identity, the current state is strong enough for a "methods"
write-up if framed honestly:

- automated search + verification pipeline for q-continued fractions
- conservative benchmark-relative audit semantics
- a worked "hero case" showing how computational eliminations and Lean-checked
  local lemmas can interact

This is closer to an experimental-math + formalization case study than a
classical "new continued fraction identity" announcement.
