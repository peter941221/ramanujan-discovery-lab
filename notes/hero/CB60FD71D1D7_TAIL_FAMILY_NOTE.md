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
- Current reading: the tail-family ladder remains structurally informative, but the sampled `U(x)` objects and their deeper gap residuals still do not collapse into the first direct eta / modular-unit boxes, the first nearby one-core eta-correction boxes, the direct Morton algebraic-function templates, the first Weber-Schlafli coordinate / companion templates, or the first literature-driven GG/Weber modular-equation boxes.

- Build elapsed seconds: `344.15`
