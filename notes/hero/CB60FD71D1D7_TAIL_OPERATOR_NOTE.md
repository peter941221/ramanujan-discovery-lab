# Tail-Operator Note: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized`
- Variable view: `t = q^3`
- Tail stages checked: `3, 4, 5`
- Max gap depth checked: `3`

## Exact Tail Family

- The operator lane starts from the same exact stationary tail law rather than from a fresh anonymous ansatz.

```text
For n >= 3, let x_n = b_n_red - 1.
Then b_n_red = 1 + x_n, a_n_red = x_n*(t + x_n), x_{n+1} = t*x_n.
T(x) = 1 + x + (x*(t + x))/T(t*x)
T(x)*T(t*x) - (1 + x)*T(t*x) - x*(t + x) = 0
```

## Operator Lane

We now ask a more recurrence-first question on the sampled tail ladder:

```text
A_0(t) + A_1(t)*Y(t) + A_2(t)*Y(t^m) + A_3(t)*Y(t^(m^2)) = 0
```

- This is the current affine q-difference / Mahler-style operator box.
- The goal is not another family label; it is a compact operator statement that could later support an operator-factorization or uniqueness proof.

### `U_t2`

- Start stage: `3`
- Gap depth: `0`
- State: `t^2`

```text
U_t2 = T(t^2) / (1 + t^2)
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t2_g1`

- Start stage: `3`
- Gap depth: `1`
- State: `t^2`

```text
U_t2_g1 = (U_t2 - 1) / t^3
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t2_g2`

- Start stage: `3`
- Gap depth: `2`
- State: `t^2`

```text
U_t2_g2 = (U_t2_g1 - 1) / t^1
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t2_g3`

- Start stage: `3`
- Gap depth: `3`
- State: `t^2`

```text
U_t2_g3 = (1 - U_t2_g2) / t^1
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t3`

- Start stage: `4`
- Gap depth: `0`
- State: `t^3`

```text
U_t3 = T(t^3) / (1 + t^3)
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t3_g1`

- Start stage: `4`
- Gap depth: `1`
- State: `t^3`

```text
U_t3_g1 = (U_t3 - 1) / t^4
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t3_g2`

- Start stage: `4`
- Gap depth: `2`
- State: `t^3`

```text
U_t3_g2 = (U_t3_g1 - 1) / t^2
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t3_g3`

- Start stage: `4`
- Gap depth: `3`
- State: `t^3`

```text
U_t3_g3 = (1 - U_t3_g2) / t^1
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t4`

- Start stage: `5`
- Gap depth: `0`
- State: `t^4`

```text
U_t4 = T(t^4) / (1 + t^4)
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t4_g1`

- Start stage: `5`
- Gap depth: `1`
- State: `t^4`

```text
U_t4_g1 = (U_t4 - 1) / t^5
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t4_g2`

- Start stage: `5`
- Gap depth: `2`
- State: `t^4`

```text
U_t4_g2 = (U_t4_g1 - 1) / t^3
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

### `U_t4_g3`

- Start stage: `5`
- Gap depth: `3`
- State: `t^4`

```text
U_t4_g3 = (1 - U_t4_g2) / t^1
```

- Moduli checked: `m=2`, `m=3`
- Recurrence depths checked: `levels=2`, `levels=3`
- Polynomial t-degrees checked: `deg_t=1`, `deg_t=2`, `deg_t=3`
- No affine q-difference / Mahler operator hit was found in the scanned box.

## Operator Verdict

- Samples checked: `12`
- Total affine q-difference / Mahler hits found: `0`
- Current reading: the exact tail law is strong enough to justify an operator-first endgame, but the first low-degree affine q-difference box is still mostly a diagnostic lane rather than a final theorem.
- Build elapsed seconds: `32.70`
