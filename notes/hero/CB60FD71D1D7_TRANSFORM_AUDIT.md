# `cb60fd71d1d7` Transform Audit

Audit date: `2026-03-13`

## Working Model

Set `t = q^3` and study the reciprocal continued fraction

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n).
```

This is the cleanest form for comparing against published generalized Rogers-Ramanujan type families, because many primary-source formulas are written as `1 + ...` continued fractions.

## Primary Sources Used

1. Bowman, Mc Laughlin, Wyshinski, 2006
   - `A generalization of Ramanujan's q-continued fraction formula`
   - arXiv: https://arxiv.org/abs/1901.00584
   - Relevant objects:
     - the four-parameter families `H(a,b,c,d,q)` and `H_1(a,b,c,d,q)`
     - canonical odd-part / contraction machinery
     - the `H_2` and `H_3` families obtained from odd parts

2. Lee, Mc Laughlin, Sohn, 2020
   - `Applications of the Heine and Bauer-Muir transformations to Rogers-Ramanujan type continued fractions`
   - arXiv: https://arxiv.org/abs/1906.11991
   - Relevant objects:
     - the Bauer-Muir chain `f1 -> f2 -> f3 -> f4`
     - Ramanujan page-43 forms `gcf2`, `gcf3`
     - the Heine-derived alternating family in Corollary `cor2cf`

## Direct Constant-Parameter Matches Ruled Out

### 1) `f2` / `gcf3`

With base variable `t`, the family has stage-`n` term

```text
(lambda*t^n - a*b*t^(2n)) / (1 + b*t^n + a*t^(n+1)).
```

Matching this against

```text
(t^n + t^(2n)) / (1 + t^n)
```

forces

```text
lambda = 1
b = 1
a = 0
a*b = -1
```

which is impossible.

### 2) `f3`

The `f3` numerators are of the form

```text
lambda*t^n + b,
```

so every stage carries a constant term `b`. The candidate numerators have no constant term at all. First-stage coefficient matching already fails.

### 3) `f4` / `gcf2`

With base variable `t`, the family has stage-`n` term

```text
(a*t + lambda*t^n) / (1 - a*t + b*t^n).
```

The first numerator can only be linear in `t`, but the candidate already needs `t + t^2` at stage 1. Constant-parameter matching fails immediately.

### 4) `H_1(a,b,c,d,t)`

The Bowman-Mc Laughlin-Wyshinski family has stage-`n` term

```text
(c*t^(n-1) - a*b*t^(2n-1)) / (d + (a+b)*t^n).
```

To match the candidate, one would need the numerator to behave like `t^n + t^(2n)`, which shifts both powers one step too high. Coefficient matching on the first two stages gives no solution.

### 5) Heine-derived `cor2cf`

The alternating family from Lee-Mc Laughlin-Sohn has the first two nontrivial stages

```text
(lambda*t - a*b*t^2)/(1 + a*t^2)
(b*t + lambda*t^2)/(1 + a*t^3).
```

The candidate uses the same structural rule at every stage, while `cor2cf` alternates two different numerator shapes. Matching the first two stages simultaneously has no constant-parameter solution.

## Simple Contraction / Odd-Part Origins Ruled Out

### 6) Odd part of `H_2`

From the odd-part formula quoted and used in Bowman-Mc Laughlin-Wyshinski, the `H_2` odd part has stage-`n` term

```text
(-a*t^n*(b*t^n + e)) / (1 + e + b*t^n + a*t^(n+1)).
```

Matching with the candidate forces incompatible conditions:

```text
a*b = -1
a*e = -1
b = 1
a = 0
e = 0
```

So the candidate is not a direct odd part of `H_2`.

### 7) Odd part of `H_3`

The analogous odd part has stage-`n` term

```text
(-b*t^n*(a*t^n + e)) / (1 + e + b*t^n + a*t^(n+1)).
```

Matching forces another inconsistent system:

```text
a*b = -1
b*e = -1
b = 1
a = 0
e = 0
```

So the candidate is not a direct odd part of `H_3` either.

## Interim Conclusion

This audit rules out the most obvious transformed origins:

- no direct constant-parameter specialization of `f2`, `f3`, `f4`, or `H_1`
- no direct match to the Heine-derived `cor2cf` alternating family
- no simple odd-part contraction coming from `H_2` or `H_3`
- no low-complexity RR-tower correction of the form

```text
F = candidate / RR(q^3)
```

  inside the current one-layer or two-layer fractional-linear prefix boxes built
  from `B1(t^k) - 1`
- no zero-shift `f2` / `gcf3` origin even after allowing arbitrary
  `n`-dependent equivalence factors in the exact necessary-identity lane

That does not prove novelty. It only means the candidate is not sitting in the first ring of standard published transforms around the Rogers-Ramanujan and Ramanujan-cubic neighborhoods already checked.

## What Still Remains Open

- an `n`-dependent equivalence transformation rather than a fixed-parameter specialization
- a multi-step contraction of a more complicated alternating continued fraction
- a nontrivial infinite Bauer-Muir chain not visible from first-stage coefficient matching
- a source outside the two audited papers

## `n`-Dependent Equivalence Transformation Audit

Let

```text
C(t) = 1 + K_{n>=1} a_n / b_n
```

with

```text
a_n = t^n + t^(2n) = t^n(1+t^n)
b_n = 1 + t^n.
```

Suppose an arbitrary equivalence transformation with nonzero factors `r_n`
maps `C(t)` to a target family

```text
1 + K_{n>=1} alpha_n / beta_n.
```

Then

```text
alpha_n = r_(n-1) r_n a_n
beta_n = r_n b_n.
```

Eliminating `r_n` gives the necessary identity

```text
alpha_n (1 + t^(n-1)) = t^n beta_(n-1) beta_n.
```

This is stronger than direct coefficient matching, because it already allows the scaling factors `r_n` to depend on `n`.

### RR Family

For the pure Rogers-Ramanujan ladder,

```text
alpha_n = t^n
beta_n = 1.
```

The necessary identity becomes

```text
t^n (1+t^(n-1)) = t^n,
```

which is impossible for nonzero `t`.

### `f2` / `gcf3`

Write `m = t^(n-1)`. Then

```text
alpha_n = lambda*m*t - a*b*m^2*t^2
beta_(n-1) = 1 + b*m + a*m*t
beta_n = 1 + b*m*t + a*m*t^2.
```

The necessary identity expands in powers of `m` with leading coefficient

```text
-t^2 (a^2 t^2 + 2ab t + ab + b^2).
```

For this to vanish identically, one is forced to `a = 0` and `b = 0`. But then the `m^1` coefficient forces `lambda = 1`, while the `m^2` coefficient becomes `t`, still nonzero. So even arbitrary `r_n` cannot move the candidate into the `f2` / `gcf3` family.

### `f3`

The expanded necessary identity has constant term

```text
b.
```

So `b = 0`. Then the `m^1` coefficient forces `lambda = 1`, the `m^3`
coefficient forces `a = 0`, and the remaining `m^2` coefficient is again `t`.
Thus no `n`-dependent equivalence transformation sends the candidate into `f3`.

### `f4` / `gcf2`

The expanded necessary identity has constant term

```text
a t,
```

so `a = 0`. The `m^3` coefficient then becomes

```text
-b^2 t^2,
```

so `b = 0`. The `m^1` coefficient forces `lambda = 1`, and the `m^2`
coefficient again survives as `t`. Therefore the candidate is not equivalent to `f4` / `gcf2` under arbitrary `r_n`.

### `H_1`

For the four-parameter `H_1` family, the `m^1` coefficient of the necessary
identity is

```text
c - d^2 t.
```

Vanishing identically forces `c = 0` and `d = 0`, which already degenerates the denominator side. So the candidate is not even `n`-dependently equivalent to `H_1`.

## Two-Step Contraction Audit Around `cor2cf`

The Heine-derived alternating family `cor2cf` is the natural place to look for a hidden two-step contraction.

To compare with the candidate reciprocal, note first that `cor2cf` starts with

```text
1 + a t.
```

The candidate starts with `1`. Since equivalence transformations and canonical
contractions do not alter that initial term, any contraction path landing on the
candidate must satisfy

```text
a = 0.
```

So only the specialized two-parameter alternating family remains:

```text
1 + (lambda t)/1 + (b t + lambda t^2)/1 + (lambda t^3)/1 + (b t^2 + lambda t^4)/1 + ...
```

### First-Step Odd and Even Parts

For `a = 0`, the first odd part begins

```text
1 + lambda t - lambda t^2 (b + lambda t) / (1 + b t + lambda t^2 + lambda t^3) + ...
```

Its initial term is `1 + lambda t`, not `1`, unless `lambda = 0`. But `lambda = 0`
kills the nontrivial Rogers-Ramanujan-adjacent branch. So the odd-first branch is out.

The first even part begins

```text
1 + (lambda t) / (1 + b t + lambda t^2) - lambda t^4 (b + lambda t) / (1 + b t^2 + lambda t^3 + lambda t^4) + ...
```

This keeps initial term `1`, so it is the only viable first contraction branch.
But already its first numerator is only linear in `t`, whereas the candidate
needs `t + t^2` at the first nontrivial stage.

### Second-Step Contractions of the Even Part

The odd part of that even part starts with

```text
1 + b t + lambda t + lambda t^2 + ...
```

so its initial term is no longer `1`. That branch is out immediately.

The even part of the even part is the only surviving two-step branch, and its
first nontrivial numerator is

```text
lambda t (1 + b t^2 + lambda t^3 + lambda t^4).
```

This numerator has no `t^2` term at all, so it cannot equal the candidate's
first numerator

```text
t + t^2
```

for any constants `b` and `lambda`.

## Direct 1-Step Bauer-Muir Audit

Another natural hypothesis is that the reciprocal candidate

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n)
```

comes from a direct 1-step Bauer-Muir transform of a nearby source fraction.

For a source reciprocal

```text
S(t) = 1 + K_{n>=1} a_n / b_n
```

and a 1-step Bauer-Muir modification by `w_n`, the transformed data begin with

```text
b0' = b0 + w0
b1' = b1 + w1
a1' = a1 - w0 b1'.
```

So whenever both source and target have `b0 = 1`, matching the initial term
forces `w0 = 0`, and then the first transformed numerator is locked to

```text
a1' = a1.
```

### RR Reciprocal

For the Rogers-Ramanujan reciprocal,

```text
S_RR(t) = 1 + K_{n>=1} t^n / 1,
```

the first source numerator is

```text
a1 = t.
```

But the candidate needs

```text
a1' = t + t^2.
```

Since `w0 = 0` is forced by `b0' = 1`, the first transformed numerator stays
`t`. Therefore the candidate is not a direct 1-step Bauer-Muir transform of the
RR reciprocal.

### Cubic Reciprocal

For the Ramanujan cubic reciprocal,

```text
S_cubic(t) = 1 + K_{n>=1} (t^n + t^(2n)) / 1,
```

the first numerator already matches the candidate. But matching the first two
partial denominators

```text
1 + t, 1 + t^2
```

forces

```text
w0 = 0
w1 = t
w2 = t^2.
```

Under those forced choices, the second transformed numerator becomes

```text
a2' = t^2 + t^4 - t(1 + t^2) = -t + t^2 - t^3 + t^4,
```

while the candidate needs

```text
t^2 + t^4.
```

So the candidate is not a direct 1-step Bauer-Muir transform of the cubic
reciprocal either.

## Constrained Two-Step Bauer-Muir Search

As a small computational extension of the direct obstruction above, the current
`research` tooling now scans a fixed low-complexity Bauer-Muir pattern family
over one-step and two-step chains from the nearby RR and cubic reciprocals.

Search shape per step:

```text
0
±(t^n - 1)
±(t^(2n) - 1)
±(t^n - t^(2n))
```

with scales `1` and `2`, giving `13` modifiers per step.

Current result for the reduced candidate reciprocal:

- RR source:
  - `13` one-step chains checked
  - `169` two-step chains checked
  - `2197` three-step chains checked
  - no hits
- cubic source:
  - `13` one-step chains checked
  - `169` two-step chains checked
  - `2197` three-step chains checked
  - no hits

These checks are still heuristic rather than exhaustive. They are verified at
exact rational sample points `t = 1/10, 1/7, 1/5` for reciprocal coefficients
through depth `4`. Still, they make the "tiny hidden low-complexity Bauer-Muir
chain" explanation less plausible.

## Updated Conclusion

The candidate survives a noticeably stronger audit than before:

- not a direct constant-parameter specialization of the main nearby families
- not an arbitrary `n`-dependent equivalence transformation of `RR`, `f2`, `f3`, `f4`, or `H_1`
- not a direct 1-step Bauer-Muir transform of either the RR reciprocal or the cubic reciprocal
- not any 1-step, 2-step, or 3-step chain inside the current tiny low-complexity Bauer-Muir pattern family built from those RR/cubic reciprocals
- not a one-step odd/even contraction of the `cor2cf` specialization relevant to the Rogers-Ramanujan neighborhood
- not the only remaining viable two-step contraction branch of that `cor2cf` specialization

The remaining open space is therefore narrower than a generic "maybe a simple transform did it" explanation.

## Extra Signal: The Ratio Is a Pure `t = q^3` Series

Using symbolic `q`-series expansion (depth `40`, order `60`), the multiplicative ratio

```text
cb60fd71d1d7 / RR(q^3)
```

has support only on exponents `q^(3k)` (no `q^(3k+1)` or `q^(3k+2)` terms appear in the same truncation window).

The first visible terms are:

```text
1 - q^12 + 2 q^15 - q^18 - q^21 + 4 q^24 - 9 q^27 + 10 q^30 - 2 q^33 - 14 q^36 + ...
```

This is consistent with the template-level expectation that the deviation from `RR(q^3)` should be expressible as a series in `t = q^3` alone.

As an additional sanity check, expressing this `t`-series as the unique Euler product

```text
R(t) = ∏_{n>=1} (1 - t^n)^{c_n},
```

the inferred exponents `c_n` up to the visible order are dense and grow quickly in magnitude (for example, among `n=1..31`, only 3 of the `c_n` vanish). This makes a "small finite eta/product tweak" explanation unlikely for this ratio.

## Cubic Odd/Even Parts Ruled Out

Because the reciprocal of `cb60fd71d1d7` uses the same numerator pattern as the
Ramanujan cubic continued fraction family (but with an added `1+t^n` in the partial
denominators), a natural hypothesis is:

```text
C(t) might be an odd/even part (canonical contraction) of the cubic reciprocal.
```

Direct reconstruction of the odd and even parts of the cubic reciprocal (treating it as
`1 + K_{n>=1} (t^n + t^(2n)) / 1`) shows:

- the odd part has initial term `d0 = 1 + t + t^2` (already incompatible with `C(t)` which has `d0 = 1`)
- the even part does have `d0 = 1` and the first numerator `c1 = t + t^2`, but the first denominator is
  `d1 = 1 + t^2 + t^4`, not `1 + t`

So `cb60fd71d1d7` is not a simple odd/even contraction of the cubic continued fraction reciprocal.
