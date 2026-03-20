# `cb60fd71d1d7` Direction Review

## Question

After the new RHS uniqueness pass, is the current main direction still correct?

## Short Answer

Yes on the **main road**, with one **local steering adjustment** on the box shape.

- Yes on the main road:
  pushing toward a theorem-facing right-hand side for the hero ratio object is
  still the correct award-track direction.
- The steering adjustment:
  the next `GG` pass should now favor exact or near-exact quotient-coordinate
  lanes over another broad prefix-box widening.

## What Changed In This Round

The identification layer now checks three theorem-facing lanes on the ratio
object side:

1. direct RHS uniqueness on
   `F(t) = candidate / RR(q^3)`
2. one-source-core corrected uniqueness on
   `G(t) = F(t) / S(t)` with `S in {RR, cubic, GG, S}`
3. two-source-core corrected uniqueness on
   `H(t) = F(t) / (S1(t) * S2(t))`
4. rational-equivalence reduced-object uniqueness on
   the reduced reciprocal bridge object `R(t)` and the reduced ratio
   `F_red(t) = B1(t) / R(t)`
5. exact reduced-coefficient transfer extraction on
   the stationary tail family forced by the reduced coefficients from stage `3`
   onward
6. anchored tail-object scans on
   `T_tail = T(t^2)` and its first visible normalization
   `U_tail = T_tail / (1 + t^2)`
7. reciprocal-normalized next-tail scans on
   `R_tail = (1 + t^3) / T(t^3)`
8. first-gap tail-residual scans on
   `G_tail = (U_tail - 1) / t^3` and
   `H_tail = (1 - R_tail) / t^4`
9. reduced-ratio reverse-scale transfer scans on
   `F_red = F_red(t^m)^a * prod_r (-t^r; t^m)_inf^{e_r}`
10. reduced-ratio reverse-scale transfer scans with a small eta tail on
    `F_red = F_red(t^m)^a * prod_r (-t^r; t^m)_inf^{e_r} * eta_tail`
11. second-gap tail-residual scans on
    `G2_tail = (G_tail - 1) / t^a` and
    `H2_tail = (H_tail - 1) / t^b`
    followed by one-core source-family eta-correction checks
12. one-core reduced-ratio reverse-scale mixed scans on
    `G_red = F_red / T`
    and
    `G_red = G_red(t^m)^a * prod_r (-t^r; t^m)_inf^{e_r} * eta_tail`

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
Reduced-tail transfer law: exact hit
Anchored tail boxes: 0 hits
Normalized anchored tail boxes: 0 hits
Reciprocal-normalized next-tail boxes: 0 hits
Gap-normalized tail boxes: 0 hits
Gap-normalized tail source-core eta boxes: 0 hits
Reduced-object bridge boxes: 0 hits
Reduced-ratio modular-unit boxes: 0 hits
Reduced-object Mahler/transfer boxes: 0 hits
Reduced-ratio plus-product boxes: 0 hits
Reduced-ratio signed-product boxes: 0 hits
Reduced-ratio plus-Pochhammer boxes: 0 hits
Reduced-ratio plus-Pochhammer + eta boxes: 0 hits
Second-gap tail source-core eta boxes: 0 hits
One-core reduced-ratio plus-Pochhammer + eta boxes: 0 hits
Reduced-ratio signed-eta transfer boxes: 0 hits
```

The new positive structural fact is:

```text
For n >= 3, let x_n = b_n_red - 1.
Then b_n_red = 1 + x_n,
     a_n_red = x_n*(t + x_n),
     x_{n+1} = t*x_n.

So the stationary tail family satisfies
T(x) = 1 + x + x*(t + x)/T(t*x),
equivalently
T(x)*T(t*x) - (1 + x)*T(t*x) - x*(t + x) = 0.
```

And even after anchoring that family at the first stationary state and trying
the first visible normalization,

```text
T_tail = T(t^2)
U_tail = T_tail / (1 + t^2)
```

the same bounded self-polynomial / self-fractional-linear / eta /
product-style boxes still give `0` hits.

Pushing one exact step deeper along the same transfer chain,

```text
R_tail = (1 + t^3) / T(t^3)
```

still gives `0` hits in that same bounded self / eta / product family.

After stripping the first visible nonzero gaps from those two exact tail
residuals,

```text
G_tail = (U_tail - 1) / t^3
H_tail = (1 - R_tail) / t^4
```

the same bounded self / eta / product family still gives `0` hits, and even
the one-core source-flavored box

```text
X_tail = T * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

still gives `0` hits for both `G_tail` and `H_tail` across nearby
`RR/cubic/GG/S` basis ladders through powers `2,3,4`.

And even after stripping one more visible gap from those first-gap residuals,

```text
G2_tail = (G_tail - 1) / t^a
H2_tail = (H_tail - 1) / t^b
```

the same one-core source-family eta-correction question still gives `0` hits
for both deeper residuals in the same nearby basis ladders.

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
4. the new reduced-object bridge pass says even the exact convergent-factor
   reduction does not immediately drop the problem into a tiny self-polynomial,
   tiny self-fractional-linear, finite-product, or small eta-quotient box
5. the same reduced-object bridge now gives a positive exact gain:
   from stage `3` onward, the coefficients collapse into a stationary tail
   family with transfer law
   `T(x) = 1 + x + x*(t + x)/T(t*x)`,
   so the search is no longer “blind RHS guessing” on the reduced lane
6. even after anchoring that transfer law at `x = t^2` and dividing by the
   first visible factor `1 + t^2`, the concrete tail objects still avoid the
   same small functional / modular-unit boxes
7. even after moving one exact step further to the reciprocal-normalized
   next-tail object `R_tail = (1 + t^3) / T(t^3)`, the same bounded boxes
   still stay empty
8. even after stripping the first visible nonzero gaps from those exact tail
   residuals, the same bounded self / eta / product boxes still stay empty
9. and even after asking whether those first-gap residuals look like one
   nearby named source-family ladder times a small eta tail, both `G_tail` and
   `H_tail` still give `0` hits
10. the next deeper source-informed passes still stay empty:
   bounded multi-level Mahler/transfer boxes and bounded `(1+t^r)` self-quotient
   product boxes both still give `0` hits on the reduced lanes
11. even after mixing the two simplest modular-unit building blocks into one
   bounded signed-product lane
   `prod (1-t^r)^{a_r} (1+t^r)^{b_r}`,
   the reduced ratio still gives `0` hits
12. even after matching the reverse-scale arithmetic-progression shape more
    directly with periodic plus-Pochhammer blocks `(-t^r; t^m)_inf`, the
    reduced ratio still gives `0` hits
13. even after combining those periodic plus-Pochhammer blocks with the
    simplest self-copy lane and a small eta tail, the reduced ratio still
    gives `0` hits
14. even after stripping a second visible gap from the first-gap tail residuals,
    the one-core source-family eta-correction lane still stays empty
15. even after removing one nearby source-family core from the reduced ratio
    and then asking for that same reverse-scale-aware plus-Pochhammer + eta
    self-transfer box, the bounded lane still gives `0` hits
16. even after allowing the reduced ratio to carry its own shifted copy together
    with signed modular-unit pieces and eta pieces in one bounded
    signed-eta transfer lane, the reduced ratio still gives `0` hits
17. even after moving from the reduced-ratio side to the tail-family-first lane
    and scanning the sampled `U(t^2)`, `U(t^3)`, `U(t^4)` objects plus their
    gap-normalized residuals through depth `3`, the one-core source-family eta
    lane still gives `0` hits on all `12` samples
18. and even after adding a literature-driven `GG/Weber` lane on those same
    `12` tail samples, the direct / quotient / mixed quotient modular-equation
    boxes still give `0` hits; in particular, the exact Chan--Huang
    quotient-coordinate lane on `Q_3 = GG(t^3)/GG(t)` and
    `Q_4 = GG(t^4)/GG(t)` still gives `0` hits

## What The New `0` Hits Actually Mean

They do **not** mean the RHS-first direction failed.

They mean something more specific:

- the hero object is probably not captured by a very low-degree self-polynomial
  over `t` and `X(t^m)`
- it is also probably not captured by the first low-degree
  fractional-linear-recursive box
- dividing out one or two nearby named source cores still does not drop the
  candidate into those same tiny boxes
- even after moving to the exact reduced reciprocal object suggested by the
  rational-equivalence proof layer, the reduced object and the reduced ratio
  still avoid the current small functional / modular-unit boxes
- but the reduced coefficients are no longer shapeless:
  they now point to a canonical stationary tail object `T(x)` with an exact
  transfer equation, so we have a more natural “right-hand-side scaffold”
  than raw `F_red` alone
- and even the first concrete incarnations of that scaffold,
  `T_tail = T(t^2)` and `U_tail = T_tail / (1 + t^2)`,
  still do not fall into the current low-complexity self / eta / product boxes
- and even one exact step deeper,
  `R_tail = (1 + t^3) / T(t^3)`,
  still does not fall into those same bounded boxes
- and even after first-gap normalization,
  `G_tail = (U_tail - 1) / t^3` and
  `H_tail = (1 - R_tail) / t^4`,
  those residuals still do not fall into the same bounded self / eta /
  product boxes
- and even after asking whether those first-gap residuals look like one nearby
  source-family core times a small eta tail, both still return `0` hits across
  the scanned `RR/cubic/GG/S` basis ladders
- and even after stripping one more visible gap from those residuals, the same
  one-core source-family eta-correction lane still stays empty
- even after adding a more source-informed recursive lane (`Mahler/transfer`)
  and a more source-informed modular-unit lane (products built from `1+t^r`,
  matching the reverse-scale shape), the reduced lanes still avoid the current
  bounded boxes
- even after letting the reduced ratio mix the `1-t^r` and `1+t^r` building
  blocks in one signed modular-unit box, the bounded lane still stays empty
- even after replacing those finite `1+t^r` pieces by the more global
  arithmetic-progression blocks `(-t^r; t^m)_inf`, the bounded reduced-ratio
  lane still stays empty
- even after combining that plus-Pochhammer lane with the simplest self-copy
  and a small eta tail, the bounded reduced-ratio transfer lane still stays
  empty
- even after removing one nearby source-family core from `F_red`, that same
  reverse-scale-aware plus-Pochhammer + eta transfer lane still stays empty
- even after combining that signed modular-unit box with the simplest self-copy
  lane `F_red(t^m)` and small eta pieces, the bounded transfer lane still stays
  empty

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

### 2.5) Keep The `GG` Orbit, But Narrow The Coordinates

The new tail-family-first `GG/Weber` pass suggests a small tactical change:

- keep the `GG/Weber` orbit as the main nearby named orbit
- stop treating the next step as “more `GG` prefixes”
- instead treat it as “fewer, more exact quotient coordinates”

Analogy:

- broad prefix scans are like shaking every door in the hallway
- quotient-coordinate exact lanes are like taking the blueprint for the one
  lock family that still matters

So the right adjustment is not a pivot away from `GG`.
It is a pivot from:

```text
broader GG box
```

to:

```text
exact / near-exact GG quotient-coordinate elimination
```

### 2.6) The Obstruction Classes Compress, But Only At Leading Order

The new `Q_3 / Q_4` witness table has a useful extra pattern:

```text
R_Q3(Y) = 2*rho(Y)*t^m(Y) + ...
R_Q4(Y) = 3*rho(Y)*t^m(Y) + ...
```

for every currently sampled tail-family object `Y`.

So the classes do unify one layer deeper than the raw table suggests.
The natural tail-derived quotient coordinate is:

```text
Theta_tail(Y) = (m(Y), rho(Y))
```

Analogy:

- `Q_3` and `Q_4` are like two cameras pointed at the same obstruction
- `Theta_tail` is the shared object casting the shadow
- but only the first silhouette is stable so far; the finer contour still
  changes from class to class

That means:

- yes, the local `Q_3 / Q_4` lane is more structured than a plain no-hit table
- no, the current renormalization is still not strong enough to replace the
  source search by itself
- the first weighted follow-up `W_34 = Q_3^3 / Q_4^2` and its logarithmic form
  `3*log(Q_3) - 2*log(Q_4)` also still miss on the sampled tail-family ladder
  and do not produce degree-`<= 2` polynomial or one-coordinate
  fractional-linear closure
- the first weighted correction `F / W_34` and the normalized residual
  `G_W34` now also still miss in the checked small eta-quotient,
  modular-unit / eta, and one-core `RR/GG` source-family eta-correction boxes
- the refreshed hero note now also shows the deeper residual
  `G2_W34 = (G_W34 - 1) / (-4*t^2)`, and that second-normalized layer still
  misses in the same checked small eta-quotient, modular-unit / eta, and
  one-core `RR/GG` source-family eta-correction boxes
- so the next move should not be another anonymous box widening
- it should be a deeper Weber / modular-curve parameter built to respect the
  same `3:2` weighting, not the raw pair `Q_3`, `Q_4` alone

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

This round partially advanced all three:

- the exact reduced bridge object is now inside `identify`
- the reduced coefficient profile / reverse scales are now surfaced directly in
  the identification note
- the reduced coefficient profile now also collapses into an exact stationary
  tail-transfer law, giving a more canonical intermediate object than the raw
  reduced ratio
- the first anchored tail objects `T_tail` and `U_tail` now also get their own
  bounded self / eta / product scan layer, and those boxes still return `0`
- the next reciprocal-normalized tail object `R_tail = (1 + t^3) / T(t^3)` now
  also gets the same bounded scan layer, and that box family still returns `0`
- the first-gap residual objects
  `G_tail = (U_tail - 1) / t^3` and `H_tail = (1 - R_tail) / t^4`
  now also get both the same bounded self / eta / product scan layer and a
  one-core source-family eta-correction lane, and both still return `0`
- the second-gap residual objects derived from `G_tail` and `H_tail` now also
  get that one-core source-family eta-correction lane, and that deeper layer
  still returns `0`
- the reduced ratio now also gets a bounded modular-unit lane
- the reduced ratio now also gets a bounded plus-product lane based on
  `1+t^r`, and both the reduced object and reduced ratio now get a bounded
  multi-level Mahler/transfer lane
- the reduced ratio now also gets a bounded periodic plus-Pochhammer lane based
  on the reverse-scale arithmetic progression factors, and that lane still
  returns `0`
- the reduced ratio now also gets a bounded periodic plus-Pochhammer + eta
  transfer lane, and that lane still returns `0`
- the reduced ratio now also gets a one-core source-family version of that same
  plus-Pochhammer + eta transfer lane, and that lane still returns `0`
- the dedicated tail-family note now also checks a source-faithful
  `GG/Weber modular-equation` lane directly on the sampled `U(x)` ladder and
  its gap-normalized residuals, and that lane still returns `0`
- within that same tail-family note, the narrowest exact quotient-coordinate
  lane on `Q_3` and `Q_4` also still returns `0`

So the next gain probably needs a more bespoke object than `F`, `F/S`,
`F/(S1*S2)`, `R`, or `F_red` alone; the most natural next candidate is the
tail family `T(x)` itself or a deeper renormalization of `T(x)` beyond the
first visible factors `1 + t^2` and `1 + t^3`, matched against a smaller
`GG` quotient-coordinate orbit rather than a broader anonymous prefix box.

## Current Verdict

```text
Direction verdict:
RHS-first is still correct.

Box verdict:
the first low-complexity direct / one-core / two-core uniqueness boxes are now
better viewed as ruled-out neighborhoods, not as the final destination.

Adjustment verdict:
keep the `GG/Weber` neighborhood, but narrow the next move toward exact
quotient-coordinate / modular-curve style elimination on `Q_3`, `Q_4`, and
their tail-derived analogues; in practice this now means using the shared
leading obstruction scalar behind the `Q_3/Q_4` pair as a diagnostic waypoint,
then moving toward a deeper weighted modular coordinate rather than broader
anonymous prefix scans. The new `G2_W34` miss means the next coordinate should
be deeper than the first two weighted correction layers, not just a slightly
larger anonymous box around `W_34`.

New structural verdict:
the reduced lane now has a canonical exact tail-transfer equation, but that
equation has not yet been converted into the final closed-form RHS.
```

## Latest Literature Check

The newest primary-source read strengthens the same conclusion.

- recent RR papers keep compressing continued fractions into modular-function,
  modular-equation, and eta-quotient coordinates, not into anonymous low-degree
  fit boxes
- the new periodic-point paper of Akkarapakam--Morton supports keeping an
  algebraic-function lane, but our current Morton-inspired box is still a
  diagnostic no-hit lane rather than a positive recognition lane
- recent Ramanujan-Machine proof papers support keeping an operator lane, but
  mostly as a proof device once the right kernel or source object is known

So the tactical order should now be read as:

1. primary:
   direct modular-unit / eta recognition on `F`, `T(x)`, or deeper tail
   renormalizations
2. secondary:
   exact RR / GG quotient-coordinate or modular-curve lanes
3. tertiary:
   tail-operator / q-difference lanes as theorem scaffolding

Analogy:

- lane 1 is trying to identify the metal of the key
- lane 2 is checking the lock family stamped on the door
- lane 3 is reverse-engineering the hinge mechanics after the right door is
  already in view

The project should keep all three lanes alive, but not give them equal rank.

### 2026-03-19 Addendum: After The `G2_W34` Layer

The deeper weighted-correction result does not change that ordering.

- local result:
  `F / W_34 - 1` first fails at `t^1` with coefficient `-1`, so the current
  hero normalization is `G_W34 = (1 - F / W_34) / t`, and the refreshed hero
  note now also shows `G2_W34 = (G_W34 - 1) / (-4*t^2)`
- scan result:
  both `G_W34` and `G2_W34` still give `0` hits in the checked small
  eta-quotient, modular-unit / eta, and one-core `RR/GG` source-family
  eta-correction boxes
- strategy result:
  this argues against both an operator-first pivot and another anonymous box
  growth step; it keeps the project on the same modular-function /
  quotient-coordinate road, but pushes the next positive-recognition attempt
  deeper than `W_34`, `G_W34`, and `G2_W34`
