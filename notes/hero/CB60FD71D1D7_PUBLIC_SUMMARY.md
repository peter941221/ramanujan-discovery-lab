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
- The first exact Chan--Huang modular-equation polynomials in the `q^3` and
  `q^4` lanes also still give `0` hits, both in the direct
  `F` vs. `GG(q^{3m})` form and in the quotient-coordinate
  `F` vs. `GG(q^{3m})/GG(q^3)` form.
- The next source-faithful Weber-Schlafli `P/B` branch also still gives `0`
  hits:
  `P_ws = (1/F - F) / 2` misses in repeated low-order classes, while the
  direct companion
  `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)` already fails in one uniform
  constant-term class across all `12` sampled tail-family objects.
- A deeper Ramanujan-Weber class-invariant compression now also exists on the
  same normalized `GG` variable:
  `g12_ws = 4*t*(Z_g - 1/Z_g)` and
  `p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`, with
  `Z_g = ((1 - t*F^2)^2) / (4*t*F^2)`.
  On the true `GG` source these recover the explicit named objects
  `(t^2; t^4)_inf^12` and `(-t^2; t^4)_inf^12`, but on the current hero
  tail-family ladder they still give `0` hits after the first direct-template,
  eta / modular-unit, and plus-Pochhammer correction passes.
- That same Weber pair is now tied together by an exact bridge
  `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
  This lets the project keep `G_g12_ws` as the current primary eta-side
  residual, treat `G_p12_ws` as the constrained companion, and insert the
  tighter quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2` together with the
  exact bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0`,
  `Q_gp_ws = p12_ws / g12_ws`.
  The same compression now also exposes the template-normalized coordinate
  `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`, which is
  exact `1` on true normalized `GG`; on the hero sample it first differs from
  `1` at `t^1` with coefficient `4`, and its first follow-up
  `H_X_ws = (G_X_ws - 1) / (4*t^1)` still stays outside the same first
  theorem-shaped closure boxes.
  The same lane now also compares the two hero-side normalized follow-ups
  directly:
  `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws` both first
  differ at `t^2` with coefficient `-24`, and no degree-`<= 3` polynomial or
  one-coordinate fractional-linear bridge was found in that first comparison
  box.
  The quotient side of that same bridge now also carries a stripped follow-up
  `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, which on the hero sample first
  differs from `1` at `t^1` with coefficient `2` and still gives `0` hits in
  the first self-polynomial, self-fractional-linear, self-quotient
  finite-product, eta / modular-unit, and plus-Pochhammer boxes.
  That stripped quotient lane is now also compared back to the
  template-normalized branch:
  `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws` both first
  differ at `t^1` with coefficient `-5/2`, and the next quotient follow-up
  `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)` still gives `0` hits in the same
  first narrow theorem-shaped closure boxes.
  After that compression, the next positive-recognition step still stays
  focused on the residual quotient `R_gp_ws = G_p12_ws / G_g12_ws`, which
  still gives `0` hits in the first narrow self-polynomial,
  self-fractional-linear, self-quotient finite-product, eta / modular-unit,
  and plus-Pochhammer boxes.
  The same lane now also keeps a normalized follow-up
  `H_gp_ws = (R_gp_ws - 1) / (96*t^3)` under watch, and that follow-up also
  remains outside the first theorem-shaped closure boxes.
- That direct Weber companion obstruction is now also packaged in Lean under
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`, so this is no longer
  only a notebook-level negative result.
- The new focused Weber residual bridge is likewise packaged under
  `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean` and exposed
  from `Proofs/HeroCaseFinalIdentity.lean` as the current theorem-shaped
  hand-off after the `g/p` compression.

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
- Chan, Huang, 1997: https://mrc.sdu.edu.cn/ziliao/8.pdf
