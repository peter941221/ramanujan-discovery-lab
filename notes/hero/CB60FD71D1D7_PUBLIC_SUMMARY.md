# `cb60fd71d1d7` Public Summary

## What It Is

`cb60fd71d1d7` is the current strongest `review` candidate in the Ramanujan Discovery Lab search pipeline.

Its template is

```text
1 / (1 + (q^3 + q^6)/(1 + q^3 + (q^6 + q^12)/(1 + q^6 + ...)))
```

and its nearest built-in benchmark is `rogers_ramanujan_q3_normalized`.

## Why It Matters

- It is not random search noise.
- Symbolically, it matches `RR(q^3)` through `q^9`.
- The first visible divergence from `RR(q^3)` occurs at `q^12`.
- It is materially closer to `RR(q^3)` than to `cubic(q^3)` in the first visible symbolic orders.

## Current Best Description

The candidate looks like a structured Rogers-Ramanujan-adjacent hybrid perturbation:

- it keeps the `RR(q^3)` main numerator progression
- it adds a cubic-style extra numerator term
- it also adds a Rogers-Ramanujan-type denominator perturbation

## What Has Been Ruled Out

- No direct constant-parameter specialization of Ramanujan Entry 6.4.4 / the nearest page-43 ratio family has been found.
- No exact match for the mixed pattern

```text
(q^(3n) + q^(6n)) / (1 + q^(3n))
```

has been found in the primary sources checked so far.
- No source-faithful `GG` modular-equation hit has been found in the first
  explicit box built from `GG(q^3)`, `GG(-q^3)`, `GG(q^6)`, `GG(q^9)`,
  `GG(q^12)`, `GG(q^15)`, `GG(q^21)`, `GG(q^33)` and the corresponding
  quotient coordinates against `GG(q^3)`, even after adding the mixed
  quotient-coordinate pass that keeps `GG(q^3)` explicit while scanning
  low-degree corrections in `GG(-q^3)/GG(q^3)` and `GG(q^{3p})/GG(q^3)`.

## What Has Not Been Claimed

- This is not a proof of novelty.
- This is not yet a publishable new-formula claim.
- It is still best treated as a high-value audit target.

## Current Working Formula

To low order,

```text
candidate = RR(q^3) * (1 - q^12 + 2 q^15 - q^18) + O(q^31)
```

## Sources

- Bowman, Mc Laughlin, Wyshinski, 2006: https://arxiv.org/abs/1901.00584
- Lee, Mc Laughlin, Sohn, 2020: https://arxiv.org/abs/1906.11991
