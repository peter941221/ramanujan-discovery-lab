# `cb60fd71d1d7` Direction Review

## Question

After the new RHS uniqueness pass, is the current main direction still correct?

## Short Answer

Yes on the **main road**, no on the **current box shape**.

- Yes on the main road:
  pushing toward a theorem-facing right-hand side for the hero ratio object is
  still the correct award-track direction.
- No on the current box shape:
  the new bounded RHS scans now show that several of the simplest “small
  equation” doors are closed.

## What Changed In This Round

The identification layer now checks three theorem-facing lanes on the ratio
object side:

1. direct RHS uniqueness on
   `F(t) = candidate / RR(q^3)`
2. one-source-core corrected uniqueness on
   `G(t) = F(t) / S(t)` with `S in {RR, cubic, GG, S}`
3. two-source-core corrected uniqueness on
   `H(t) = F(t) / (S1(t) * S2(t))`

The bounded boxes now checked are:

- self-polynomial:
  `P(t, X(t), X(t^m)) = 0`
  with `m in {2,3,4,5,6}`, `deg_(X, X_m) <= 3`, `deg_t <= 3`
- self-fractional-linear:
  `X(t) = (A(t) + B(t)*(X(t^m) - 1)) / (C(t) + D(t)*(X(t^m) - 1))`
  with `deg_t <= 3` on the direct ratio object
- corrected-object scans:
  the same low-complexity polynomial / fractional-linear lanes after removing
  one or two nearby source cores

Current result:

```text
Direct RHS boxes: 0 hits
One-core corrected RHS boxes: 0 hits
Two-core corrected RHS boxes: 0 hits
```

## Why The Direction Is Still Correct

Analogy:

- We know the left side of the equation very well now.
- The new pass tried several “small locks” on the right side.
- None of those locks opened.
- That means the problem is probably not “wrong corridor”.
- It means the true lock is more structured than the first key ring.

In award-track terms, this is still the correct trunk because:

1. the prize bar is a theorem bar, so a compact defining equation or a source
   object with a uniqueness theorem still has the highest value
2. the direct family-expansion roads have already accumulated many exact and
   bounded no-hit results
3. the new one-core / two-core correction pass says the hero case is not
   immediately explained by “nearby known source object + tiny self-equation”

## What The New `0` Hits Actually Mean

They do **not** mean the RHS-first direction failed.

They mean something more specific:

- the hero object is probably not captured by a very low-degree self-polynomial
  over `t` and `X(t^m)`
- it is also probably not captured by the first low-degree
  fractional-linear-recursive box
- dividing out one or two nearby named source cores still does not drop the
  candidate into those same tiny boxes

So the current situation is closer to:

```text
left side known
right side unknown
small-box RHS guesses ruled out
structured RHS generation still needed
```

## Recommended Direction Now

### 1) Keep The RHS-First Direction

Do not go back to broad family expansion as the main bet.

The direct theorem-facing road is still the best road because it is the only
road that can naturally end in:

- a uniqueness statement
- a modularity statement
- a final Lean-formalizable theorem

### 2) Stop Expecting Tiny Generic Boxes To Be Enough

The new evidence suggests the next gain probably will **not** come from merely
turning:

- `deg <= 2` into `deg <= 3`
- or `m <= 4` into `m <= 6`

Those were worth testing, and now they count as useful eliminations.
But the next box should be more structured, not just larger.

### 3) Prefer Structured Intermediate Objects

The best next RHS work should probably come from one of these lanes:

1. derive an intermediate object from the exact rational-equivalence /
   convergent-reduction story and search equations for that object instead of
   for raw `F`
2. extract recurrence or transfer equations forced by the continued-fraction
   coefficients themselves and convert those into candidate-defining functional
   equations
3. use modular-unit / eta-quotient architecture to propose the RHS shape first,
   then verify it, instead of asking a generic low-degree box to guess it

## Current Verdict

```text
Direction verdict:
RHS-first is still correct.

Box verdict:
the first low-complexity direct / one-core / two-core uniqueness boxes are now
better viewed as ruled-out neighborhoods, not as the final destination.
```
