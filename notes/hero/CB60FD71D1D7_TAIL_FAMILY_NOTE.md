# Tail-Family Note: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized`
- Variable view: `t = q^3`
- Tail stages checked: `3, 4, 5`
- Max gap depth checked: `3`

## Exact Tail Family

- From stage `3` onward, the reduced coefficients collapse into one stationary tail family.

```text
For n >= 3, let x_n = b_n_red - 1.
Then b_n_red = 1 + x_n, a_n_red = x_n*(t + x_n), x_{n+1} = t*x_n.
T(x) = 1 + x + (x*(t + x))/T(t*x)
T(x)*T(t*x) - (1 + x)*T(t*x) - x*(t + x) = 0
```

## Variable-Level Source-Core Recognition Lane

We now treat the normalized tail family itself as the main intermediate object:

```text
U(x) = T(x) / (1 + x)
```

- For each sampled state `x = t^k`, we scan `U(x)` and its repeated gap-normalized residuals.
- Unlike the broader `identify` note, this lane keeps only the source-driven one-core eta-correction question:

```text
Y = S * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- We also keep two more source-faithful direct recognition lanes on the same sampled objects:

```text
Y = prod_r (1 - t^r)^{a_r} * prod_r (1 + t^r)^{b_r} * prod_{d|N} (t^d; t^d)_inf^{e_d}
f(Y, Y_2) = 0,   g(Y^2, Y_2^2) = 0,   f(Y, (1-Y_2)/(1+Y_2)) = 0
```

- The first is a direct modular-unit / eta lane.
- The second is a Morton-2024-inspired periodic-point / algebraic-function lane built from the exact low-degree polynomials attached to the GG/Weber orbit.
- Phase 2 now also opens the first deeper Weber-Schlafli coordinate lane on the same sampled objects:

```text
P_ws = (1/Y - Y) / 2
B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)
P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0
B_ws^2 - B_ws,2 - 4 = 0
B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0
```

- Phase 3 now adds the first Ramanujan-Weber class-invariant compression on the same normalized `GG` variable:

```text
Z_g = ((1 - t*Y^2)^2) / (4*t*Y^2)
g12_ws = 4*t*(Z_g - 1/Z_g)
p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)
g12_ws ?= (t^2; t^4)_inf^12
p12_ws ?= (-t^2; t^4)_inf^12
R_gp_ws = G_p12_ws / G_g12_ws
g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0
```

- Phase 4 now treats the eta-side residual `G_g12_ws` as the primary Weber hand-off, uses the exact algebraic bridge between `g12_ws` and `p12_ws` to keep `G_p12_ws` as a constrained companion, and then probes the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws` before reopening any broader search box.

- We also check a `GG/Weber modular-equation` lane on the same sampled `U(x)` objects and their gap residuals:

```text
Y = T
Y = 1 / T
Y = T_i / T_j
P(Y, T_i) = 0
Y = prod_i T_i^e_i
Y = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))
```

- GG base benchmark: `gollnitz_gordon_normalized`
- GG basis ladder: `GG = GG(t)`, `GGneg = GG(-t)`, `GG2 = GG(t^2)`, `GG3 = GG(t^3)`, `GG4 = GG(t^4)`
- Preferred quotient coordinates: `Q_neg = GG(-t) / GG(t)`, `Q_2 = GG(t^2) / GG(t)`, `Q_3 = GG(t^3) / GG(t)`, `Q_4 = GG(t^4) / GG(t)`
- The narrowest exact quotient-coordinate lane keeps special attention on `Q_3 = GG(t^3)/GG(t)` and `Q_4 = GG(t^4)/GG(t)`, because those are the Chan--Huang exact modular-equation coordinates.
- Each sample below reports direct, quotient-coordinate, and mixed quotient-coordinate prefix summaries for this literature-driven lane.

### `U_t2`

- Start stage: `3`
- Gap depth: `0`
- State: `t^2`

```text
U_t2 = T(t^2) / (1 + t^2)
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^6` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (4*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `2`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-1`.
- GG normalized weighted correction `G_W34`: G_W34 = (1 - F / W_34) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-5*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t2_g1`

- Start stage: `3`
- Gap depth: `1`
- State: `t^2`

```text
U_t2_g1 = (U_t2 - 1) / t^3
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `8`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (8*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `11/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `192`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (192*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `11/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `4`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `6`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^2` with coefficient `-2`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^2` with coefficient `-2`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^2` with coefficient `-2`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / (-2*t^2).
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-3/2*t^1).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t2_g2`

- Start stage: `3`
- Gap depth: `2`
- State: `t^2`

```text
U_t2_g2 = (U_t2_g1 - 1) / t^1
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `-1`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^2` with coefficient `4`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^2` with coefficient `4`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^2` with coefficient `-8`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (-8*t^2)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^2` with coefficient `-163/6`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^4` with coefficient `-192`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (-192*t^4)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^2` with coefficient `-4`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^2` with coefficient `-6`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-2`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-2`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-2`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / (-2*t^1).
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (1/2*t^1).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t2_g3`

- Start stage: `3`
- Gap depth: `3`
- State: `t^2`

```text
U_t2_g3 = (1 - U_t2_g2) / t^1
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `8`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-6`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-6`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `12`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (12*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `47/6`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `6`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `288`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (288*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `47/6`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `6`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `9`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `1`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (1 - G_W34) / t^1.
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t3`

- Start stage: `4`
- Gap depth: `0`
- State: `t^3`

```text
U_t3 = T(t^3) / (1 + t^3)
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^8` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (4*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `2`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-1`.
- GG normalized weighted correction `G_W34`: G_W34 = (1 - F / W_34) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-4*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t3_g1`

- Start stage: `4`
- Gap depth: `1`
- State: `t^3`

```text
U_t3_g1 = (U_t3 - 1) / t^4
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^4` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (4*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `11/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `2`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `11/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-1`.
- GG normalized weighted correction `G_W34`: G_W34 = (1 - F / W_34) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (1 - G_W34) / t^1.
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t3_g2`

- Start stage: `4`
- Gap depth: `2`
- State: `t^3`

```text
U_t3_g2 = (U_t3_g1 - 1) / t^2
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `-1`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^2` with coefficient `2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^2` with coefficient `2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^2` with coefficient `-4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (-4*t^2)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `6`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^2` with coefficient `-151/6`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^4` with coefficient `-96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (-96*t^4)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `6`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^2` with coefficient `-2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^2` with coefficient `-3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-2`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-2`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-2`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / (-2*t^1).
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-3/2*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t3_g3`

- Start stage: `4`
- Gap depth: `3`
- State: `t^3`

```text
U_t3_g3 = (1 - U_t3_g2) / t^1
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `8`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (8*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `7`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `192`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (192*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `7`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `4`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `6`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^2` with coefficient `1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^2` with coefficient `1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^2` with coefficient `1`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / t^2.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (2*t^1).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t4`

- Start stage: `5`
- Gap depth: `0`
- State: `t^4`

```text
U_t4 = T(t^4) / (1 + t^4)
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^10` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (4*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `2`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-1`.
- GG normalized weighted correction `G_W34`: G_W34 = (1 - F / W_34) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-4*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t4_g1`

- Start stage: `5`
- Gap depth: `1`
- State: `t^4`

```text
U_t4_g1 = (U_t4 - 1) / t^5
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^6` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (4*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `2`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `9/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-1`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-1`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-1`.
- GG normalized weighted correction `G_W34`: G_W34 = (1 - F / W_34) / t^1.
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-5*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t4_g2`

- Start stage: `5`
- Gap depth: `2`
- State: `t^4`

```text
U_t4_g2 = (U_t4_g1 - 1) / t^3
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `-1`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^2` with coefficient `2`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^2` with coefficient `2`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^2` with coefficient `-4`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (-4*t^2)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `5`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^2` with coefficient `-151/6`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^4` with coefficient `-96`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (-96*t^4)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `5`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^2` with coefficient `-2`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^2` with coefficient `-3`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^1` with coefficient `-2`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^1` with coefficient `-2`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^1` with coefficient `-2`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / (-2*t^1).
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-2*t^2).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

### `U_t4_g3`

- Start stage: `5`
- Gap depth: `3`
- State: `t^4`

```text
U_t4_g3 = (1 - U_t4_g2) / t^1
```

- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`.
- Direct eta-quotient templates: `0` / `4` hit boxes.
- Direct modular-unit / eta templates: `0` / `12` hit boxes.
- Morton periodic-point / algebraic-function templates: `0` / `4` exact hits.
- Morton obstruction witnesses: `Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`` first fails at `t^0` with coefficient `2`; `Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`` first fails at `t^0` with coefficient `4`; `Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`` first fails at `t^0` with coefficient `1`; `Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`` first fails at `t^0` with coefficient `1`.
- Morton Weber-Schlafli coordinate `P_ws`: `P_ws = (1/F - F) / 2`.
- Morton Weber-Schlafli templates on `P_ws`: `0` / `1` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`` first fails at `t^2` with coefficient `3`.
- Morton Weber-Schlafli coordinate `B_ws`: `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`.
- Morton Weber-Schlafli templates on `B_ws`: `0` / `2` exact hits.
- Morton Weber-Schlafli obstruction witnesses: `Morton Weber companion template `B^2 - B_2 - 4`` first fails at `t^0` with coefficient `-2`; `Morton Weber companion template `B_2^4 - P^8 - 16*P^4`` first fails at `t^0` with coefficient `16`.
- Weber class-invariant coordinate `g12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), g12_ws = 4*t*(Z_g - 1/Z_g)`.
- Weber class-invariant template on `g12_ws`: `(t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber g-coordinate template` first differs from `(t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_g12_ws`: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber class-invariant coordinate `p12_ws`: `Z_g = ((1 - t*F^2)^2) / (4*t*F^2), p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`.
- Weber class-invariant template on `p12_ws`: `(-t^2; t^4)_inf^12`.
- Weber class-invariant obstruction witness: `Chan--Huang Weber G-coordinate template` first differs from `(-t^2; t^4)_inf^12` at `t^1` with coefficient `-4`.
- Weber class-invariant correction `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber class-invariant correction eta templates: `0` / `4` hit boxes.
- Weber class-invariant correction modular-unit / eta templates: `0` / `12` hit boxes.
- Weber class-invariant correction plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber class-invariant correction plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual bridge keeps `G_g12_ws` as the current primary residual: `G_g12_ws = g12_ws / (t^2; t^4)_inf^12`.
- Weber residual bridge reason: The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion.
- Weber residual companion `G_p12_ws`: `G_p12_ws = p12_ws / (-t^2; t^4)_inf^12`.
- Weber residual exact coordinate bridge: `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
- Weber residual exact coordinate bridge verdict: matches through the checked truncation.
- Weber residual exact residual bridge: `(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0`.
- Weber residual quotient-coordinate `X_g_ws`: `X_g_ws = 16*t^2 / g12_ws^2`.
- Weber residual exact quotient-coordinate bridge: `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws`.
- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`.
- Weber quotient-coordinate template bridge: `G_X_ws*G_g12_ws^2 - 1 = 0`.
- Weber quotient-coordinate template normalization `G_X_ws`: `G_X_ws - 1` first fails at `t^1` with coefficient `8`.
- Weber quotient-coordinate template self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate template self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate template eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate template modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate template plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws = (G_X_ws - 1) / (8*t^1)`.
- Weber quotient-coordinate normalized follow-up `H_X_ws`: `H_X_ws - 1` first fails at `t^1` with coefficient `13/2`.
- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber quotient-coordinate normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized eta templates: `0` / `4` hit boxes.
- Weber quotient-coordinate normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber normalized follow-up bridge difference `D_XR_ws`: `D_XR_ws = H_gp_ws - H_X_ws`.
- Weber normalized follow-up bridge difference `D_XR_ws`: first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws = H_gp_ws / H_X_ws`.
- Weber normalized follow-up bridge quotient `Q_XR_ws`: `Q_XR_ws - 1` first fails at `t^2` with coefficient `-24`.
- Weber normalized follow-up bridge polynomial boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge fractional-linear box: `0` / `1` hit boxes.
- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`.
- Weber normalized follow-up bridge quotient normalized follow-up `K_XR_ws`: `K_XR_ws - 1` first fails at `t^1` with coefficient `4`.
- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized eta templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: `0` / `9` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: `0` / `9` hit boxes.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws`.
- Weber residual quotient diagnostic `R_gp_ws`: `R_gp_ws = G_p12_ws / G_g12_ws - 1` first fails at `t^3` with coefficient `192`.
- Weber residual quotient self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual quotient self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual quotient self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual quotient eta templates: `0` / `4` hit boxes.
- Weber residual quotient modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual quotient plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual quotient plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws = (R_gp_ws - 1) / (192*t^3)`.
- Weber residual normalized follow-up `H_gp_ws`: `H_gp_ws - 1` first fails at `t^1` with coefficient `13/2`.
- Weber residual normalized self-polynomial uniqueness boxes: `0` / `18` hit boxes.
- Weber residual normalized self-fractional-linear uniqueness boxes: `0` / `9` hit boxes.
- Weber residual normalized self-quotient finite-product boxes: `0` / `3` hit boxes.
- Weber residual normalized eta templates: `0` / `4` hit boxes.
- Weber residual normalized modular-unit / eta templates: `0` / `12` hit boxes.
- Weber residual normalized plus-Pochhammer templates: `0` / `3` hit boxes.
- Weber residual normalized plus-Pochhammer + eta templates: `0` / `12` hit boxes.
- GG direct / reciprocal / quotient templates: `0` / `30` exact template hits.
- GG direct exact modular-equation templates: `0` / `2` exact Chan--Huang direct hits.
- GG quotient exact modular-equation templates: `0` / `2` exact Chan--Huang quotient-coordinate hits.
- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.
- GG exact quotient-coordinate obstruction witnesses: `Chan--Huang Cor. 3.2(i) on (F, Q_3)` first fails at `t^1` with coefficient `4`; `Chan--Huang Cor. 3.2(ii) on (F, Q_4)` first fails at `t^1` with coefficient `6`.
- GG weighted quotient-coordinate diagnostic `W_34 = Q_3^3 / Q_4^2`: `F - W_34` first fails at `t^3` with coefficient `4`; `log(F) - (3*log(Q_3) - 2*log(Q_4))` first fails at `t^3` with coefficient `4`.
- GG weighted correction `F / W_34`: `F / W_34 - 1` first fails at `t^3` with coefficient `4`.
- GG normalized weighted correction `G_W34`: G_W34 = (F / W_34 - 1) / (4*t^3).
- GG second normalized weighted correction `G2_W34`: G2_W34 = (G_W34 - 1) / (-3/4*t^1).
- GG weighted quotient-coordinate `W_34`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found.
- GG weighted correction `F / W_34`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes.
- GG normalized weighted correction `G_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes.
- GG normalized weighted correction `G_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes.
- GG second normalized weighted correction `G2_W34`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases.
- GG second normalized weighted correction `G2_W34`: no explicit GG transform-template eta-correction hit was found in the checked small boxes.
- GG second normalized weighted correction `G2_W34` quotient-coordinate prefixes: polynomial `0` / `6` hit prefixes; skipped `2`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG second normalized weighted correction `G2_W34` mixed quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `3`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG direct prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.
- GG quotient-coordinate prefixes: polynomial `0` / `7` hit prefixes; skipped `1`; multiplicative `0` / `4` hit prefixes; fractional-linear `0` / `4` hit prefixes; two-layer fractional-linear `0` / `3` hit prefixes.
- GG mixed quotient-coordinate prefixes: polynomial `0` / `8` hit prefixes; skipped `2`; multiplicative `0` / `5` hit prefixes; fractional-linear `0` / `5` hit prefixes; two-layer fractional-linear `0` / `4` hit prefixes.

## Tail Verdict

- Samples checked: `12`
- Source-core eta hits found: `0`
- Direct eta-quotient sample hits found: `0`
- Direct modular-unit / eta sample hits found: `0`
- GG/Weber modular-equation sample hits found: `0`
- GG exact quotient-coordinate sample hits found: `0`
- Morton periodic-point / algebraic-function sample hits found: `0`
- Morton Weber-Schlafli sample hits found: `0`
- Weber g-class-invariant sample hits found: `0`
- Weber G-class-invariant sample hits found: `0`
- Weber residual-quotient sample hits found: `0`
- Weber residual-follow-up sample hits found: `0`
- Current reading: the tail-family ladder remains structurally informative, but the sampled `U(x)` objects and their deeper gap residuals still do not collapse into the first direct eta / modular-unit boxes, the first nearby one-core eta-correction boxes, the direct Morton algebraic-function templates, the first Weber-Schlafli coordinate / companion templates, the first Ramanujan-Weber class-invariant compression boxes, the focused Weber residual-quotient box, the normalized Weber residual-follow-up box, or the first literature-driven GG/Weber modular-equation boxes.

- Build elapsed seconds: `565.03`
