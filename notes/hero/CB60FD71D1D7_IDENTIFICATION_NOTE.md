# Identification Note: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized`
- Variable view: `t = q^3`
- Depth: `40`
- Series order: `90`
- Polynomial relation search: total degree `<= 4`
- Build elapsed seconds before final render: `1749.75`

## Build Timing

- `series-and-benchmark-setup`: `0.00`
- `rhs-uniqueness-search`: `69.33`
- `source-family-scans`: `414.93`
- `cross-family-functional-scans`: `816.92`
- `explicit-gg-family-scans`: `375.27`
- `benchmark-tower-scans`: `73.30`
- `final-render`: `13.43`

## Objects

- Candidate template:
  - `qcf:top=1;base=1;ns=1;nr=1;nk=1;nes=1;ner=2;nek=2;dc=1;ds=1;dr=1;dk=1`
- Benchmark template:
  - `qcf:top=1;base=1;ns=1;nr=1;nk=1;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`

We run the relation search on the **reciprocal** continued fractions (the `1 + ...` objects):

- `C = 1 / candidate`
- `B1 = 1 / rogers_ramanujan_q3_normalized`

## Result

No nontrivial polynomial relation

```text
P(C, B1) = 0
```

was found in the search box `total degree <= 4` when checked modulo `t^90`.

## Extra Multivariate Search

We also tried a small multivariate search that includes benchmark power substitutions:

- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

Skipped multivariate relation search:

```text
underdetermined polynomial relation search: 715 monomials > 90 constraints (increase order, lower the search box, or reduce variables)
```

## RHS Uniqueness Search

We also ran a theorem-facing search directly on the ratio object, looking for a compact right-hand-side defining equation for:

- `F = candidate / rogers_ramanujan_q3_normalized`
- Moduli checked: `m=2`, `m=3`, `m=4`, `m=5`, `m=6`
- Self-polynomial boxes: `deg_(F,G) <= 1`, `deg_(F,G) <= 2`, `deg_(F,G) <= 3`
- `t`-degree boxes: `deg_t <= 1`, `deg_t <= 2`, `deg_t <= 3`
- Self-fractional-linear `t`-degree boxes: `deg_t <= 1`, `deg_t <= 2`, `deg_t <= 3`

### Polynomial Functional Box

```text
P(t, F(t), F(t^m)) = 0
```

No candidate-dependent self-polynomial uniqueness relation was found in the scanned box.

### Fractional-Linear Functional Box

```text
F(t) = (A(t) + B(t)*(F(t^m) - 1)) / (C(t) + D(t)*(F(t^m) - 1))
```

No candidate-dependent self-fractional-linear uniqueness relation was found in the scanned box.

### One-Source-Core Correction Objects

We then stripped a single nearby source core and repeated the uniqueness search on the residual correction object:

```text
G = F / S
P(t, G(t), G(t^m)) = 0
G(t) = (A(t) + B(t)*(G(t^m) - 1)) / (C(t) + D(t)*(G(t^m) - 1))
```

- Source cores checked: `RR`, `cubic`, `GG`, `S`
- Correction objects checked: `4`

No one-core self-polynomial correction hit was found in the scanned box.

No one-core self-fractional-linear correction hit was found in the scanned box.

### Two-Source-Core Correction Objects

We also stripped products of two nearby source cores and repeated the same bounded uniqueness search on the residual object:

```text
H = F / (S1 * S2)
P(t, H(t), H(t^m)) = 0
H(t) = (A(t) + B(t)*(H(t^m) - 1)) / (C(t) + D(t)*(H(t^m) - 1))
```

- Two-core correction objects checked: `6`

No two-core self-polynomial correction hit was found in the scanned box.

No two-core self-fractional-linear correction hit was found in the scanned box.

### Rational-Equivalence Reduced Object

We also switched from the raw hero reciprocal to the reduced-by-factor object coming from the exact convergent-factor / rational-equivalence bridge:

```text
R = reduced hero reciprocal object
F_red = B1 / R    where    B1 = 1 / rogers_ramanujan_q3_normalized
a1_red = t
b1_red = 1
a2_red = t^2
b2_red = t + 1
a3_red = t^4 + t^3
b3_red = t^2 + 1
r1 = t + 1
r2 = t^2/(t + 1) + 1/(t + 1)
r3 = t^3/(t^2 + 1) + 1/(t^2 + 1)
```

- This is the first source-informed RHS lane: it uses the exact reduction/equivalence bridge rather than a generic low-degree guess on the raw ratio object.
- To keep this bridge tractable inside `identify`, the reduced-object lane is currently bounded at depth `12` and order `36`.

- From stage `3` onward, the reduced coefficients collapse into one stationary tail family.

```text
For n >= 3, let x_n = b_n_red - 1.
Then b_n_red = 1 + x_n, a_n_red = x_n*(t + x_n), x_{n+1} = t*x_n.
T(x) = 1 + x + (x*(t + x))/T(t*x)
T(x)*T(t*x) - (1 + x)*T(t*x) - x*(t + x) = 0
```

- We also anchored that tail law at the first stationary stage and scanned the concrete tail object itself.

```text
T_tail = T(t^2)
U_tail = T_tail / (1 + t^2)
```

No anchored-tail hit was found in the scanned self-polynomial, self-fractional-linear, eta-quotient, finite-product, plus-product, signed-product, or signed-eta boxes.

No normalized anchored-tail hit was found in the same scanned box family after dividing by the visible factor `1 + x`.

No reduced-object self-polynomial hit was found in the scanned box.

No reduced-object self-fractional-linear hit was found in the scanned box.

No reduced-object Mahler/transfer hit was found in the scanned box.

No reduced-ratio self-polynomial hit was found in the scanned box.

No reduced-ratio self-fractional-linear hit was found in the scanned box.

No reduced-ratio Mahler/transfer hit was found in the scanned box.

No reduced-ratio self-quotient finite-product hit was found in the scanned box.

No reduced-ratio self-quotient plus-product hit was found in the scanned box.

No reduced-ratio self-quotient signed-product hit was found in the scanned box.

No reduced-ratio self signed-eta transfer hit was found in the scanned box.

No reduced-ratio eta-quotient hit was found in the scanned box.


## Benchmark Power-Tower Prefix Scan

We also ran a structured low-degree scan against prefixes of the benchmark power tower:

- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

- Prefixes checked: `(C, B1, B2)`, then `(C, B1, B2, B3)`, and so on through the final listed power.
- Degrees scanned here are intentionally low (`1` and `2`) to stay in candidate-dependent, well-determined boxes.

No candidate-dependent relation was found in any scanned prefix box.

- `total degree <= 1`: no hit for prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.
- `total degree <= 2`: no hit for prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.

## Ratio-Object Source-Family Multiplicative Scan

We also searched for exact multiplicative corrections built from nearby named source families:

```text
F = prod_i S_i^e_i
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- The source-family bases are evaluated in the same variable view used above.
- `RR = rogers_ramanujan_normalized`
- `cubic = ramanujan_cubic_normalized`
- `GG = gollnitz_gordon_normalized`
- `S = hirschhorn_s_normalized`

- Prefixes are scanned in that order, solving exact integer exponents from the log-series constraints and then verifying by exact series re-expansion.

No source-family multiplicative relation was found in any scanned prefix box.

- No hit for source-family prefixes ending at `RR`, `cubic`, `GG`, `S`.

## Ratio-Object Source-Family Fractional-Linear Scan

We also searched for low-complexity fractional-linear corrections built from nearby named source families:

```text
F = (1 + sum a_i*(S_i - 1)) / (1 + sum b_i*(S_i - 1))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- The source-family bases are evaluated in the same variable view used above.
- `RR = rogers_ramanujan_normalized`
- `cubic = ramanujan_cubic_normalized`
- `GG = gollnitz_gordon_normalized`
- `S = hirschhorn_s_normalized`

- Prefixes are scanned in that order, solving an exact linear system for the numerator and denominator correction coefficients in each source-family box.

No source-family fractional-linear relation was found in any scanned prefix box.

- No hit for source-family fractional-linear prefixes ending at `RR`, `cubic`, `GG`, `S`.

## Ratio-Object Source-Family Two-Layer Fractional-Linear Scan

We then expanded to a second-ring nonlinear box built from two single-basis factors drawn from the named source-family prefixes:

```text
F = ((1 + a0*(X1 - 1)) / (1 + b0*(Y1 - 1))) * ((1 + a1*(X2 - 1)) / (1 + b1*(Y2 - 1)))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- The source-family bases are evaluated in the same variable view used above.
- `RR = rogers_ramanujan_normalized`
- `cubic = ramanujan_cubic_normalized`
- `GG = gollnitz_gordon_normalized`
- `S = hirschhorn_s_normalized`

- Prefixes checked: `(RR, cubic)`, then `(RR, cubic, GG)`, and so on through the final listed source-family basis.
- Each prefix scans low-complexity two-factor templates and verifies any candidate hit by exact series re-expansion.

No source-family two-layer fractional-linear relation was found in any scanned prefix box.

- No hit for source-family two-layer prefixes ending at `cubic`, `GG`, `S`.

## Ratio-Object Parameterized Source-Family Power Scan

We also scanned short power ladders inside each named source family so the family meaning stays explicit:

```text
P(F, T_i) = 0
F = prod_i T_i^e_i
F = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))
F = prod_j (1 + a_j*(T_r(j) - 1)) / (1 + b_j*(T_s(j) - 1))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Each family is scanned separately, using the base object together with powered substitutions in the same variable view.
- This keeps the Gordon/Hirschhorn family labels explicit instead of collapsing them into one anonymous mixed basis.
- The low-degree polynomial box is motivated by literature where `GG` / Hirschhorn-type objects can satisfy nontrivial power-substitution identities without reducing to a pure product or a simple quotient.
- We now also include a family-preserving two-layer fractional-linear box, so simple nonlinear corrections stay inside one literature family instead of mixing labels.
- We also scan a within-family quotient ladder `Qk = Tk / T1`, which is often a more natural coordinate for power-substitution identities than the raw powered objects themselves.
- That quotient ladder now also gets its own two-layer fractional-linear pass, so the first nonlinear corrections can stay in quotient coordinates without crossing families.
- We also scan a mixed quotient basis built from the family base object together with the quotient ladder, so relations of the form `T1 * correction(Q2, Q3, ...)` can surface without mixing literature families.
- That mixed quotient basis now also gets its own two-layer fractional-linear pass, so the first nonlinear corrections can use both the base object and quotient coordinates while still staying in one family.
- The exact powered labels are listed separately inside each family subsection, because the literature-motivated ladders are now family-specific.

### `RR` Family

- Base benchmark: `rogers_ramanujan_normalized`
- Basis ladder: `RR = RR(t)`, `RR2 = RR(t^2)`, `RR3 = RR(t^3)`, `RR4 = RR(t^4)`
- Polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Polynomial `total degree <= 1`: no hit for prefixes ending at `RR`, `RR2`, `RR3`, `RR4`.
- Polynomial `total degree <= 2`: no hit for prefixes ending at `RR`, `RR2`, `RR3`, `RR4`.
- Multiplicative scan: no hit for prefixes ending at `RR`, `RR2`, `RR3`, `RR4`.
- Fractional-linear scan: no hit for prefixes ending at `RR`, `RR2`, `RR3`, `RR4`.
- Two-layer fractional-linear scan: no hit for prefixes ending at `RR2`, `RR3`, `RR4`.
- Quotient ladder: `Q2 = RR2 / RR`, `Q3 = RR3 / RR`, `Q4 = RR4 / RR`
- Quotient-ladder polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Quotient-ladder polynomial `total degree <= 1`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder polynomial `total degree <= 2`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder multiplicative scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder two-layer fractional-linear scan: no hit for prefixes ending at `Q3`, `Q4`.
- Mixed quotient basis: `RR = RR(t)`, `Q2 = RR2 / RR`, `Q3 = RR3 / RR`, `Q4 = RR4 / RR`
- Mixed-quotient polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Mixed-quotient polynomial `total degree <= 1`: no hit for prefixes ending at `RR`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient polynomial `total degree <= 2`: no hit for prefixes ending at `RR`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient multiplicative scan: no hit for prefixes ending at `RR`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient fractional-linear scan: no hit for prefixes ending at `RR`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient two-layer fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
### `cubic` Family

- Base benchmark: `ramanujan_cubic_normalized`
- Basis ladder: `cubic = cubic(t)`, `cubic2 = cubic(t^2)`, `cubic3 = cubic(t^3)`, `cubic4 = cubic(t^4)`
- Polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Polynomial `total degree <= 1`: no hit for prefixes ending at `cubic`, `cubic2`, `cubic3`, `cubic4`.
- Polynomial `total degree <= 2`: no hit for prefixes ending at `cubic`, `cubic2`, `cubic3`, `cubic4`.
- Multiplicative scan: no hit for prefixes ending at `cubic`, `cubic2`, `cubic3`, `cubic4`.
- Fractional-linear scan: no hit for prefixes ending at `cubic`, `cubic2`, `cubic3`, `cubic4`.
- Two-layer fractional-linear scan: no hit for prefixes ending at `cubic2`, `cubic3`, `cubic4`.
- Quotient ladder: `Q2 = cubic2 / cubic`, `Q3 = cubic3 / cubic`, `Q4 = cubic4 / cubic`
- Quotient-ladder polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Quotient-ladder polynomial `total degree <= 1`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder polynomial `total degree <= 2`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder multiplicative scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder two-layer fractional-linear scan: no hit for prefixes ending at `Q3`, `Q4`.
- Mixed quotient basis: `cubic = cubic(t)`, `Q2 = cubic2 / cubic`, `Q3 = cubic3 / cubic`, `Q4 = cubic4 / cubic`
- Mixed-quotient polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Mixed-quotient polynomial `total degree <= 1`: no hit for prefixes ending at `cubic`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient polynomial `total degree <= 2`: no hit for prefixes ending at `cubic`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient multiplicative scan: no hit for prefixes ending at `cubic`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient fractional-linear scan: no hit for prefixes ending at `cubic`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient two-layer fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
### `GG` Family

- Base benchmark: `gollnitz_gordon_normalized`
- Basis ladder: `GG = GG(t)`, `GG2 = GG(t^2)`, `GG3 = GG(t^3)`, `GG4 = GG(t^4)`, `GG5 = GG(t^5)`, `GG7 = GG(t^7)`, `GG11 = GG(t^11)`
- Polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Polynomial `total degree <= 1`: no hit for prefixes ending at `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Polynomial `total degree <= 2`: no hit for prefixes ending at `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Multiplicative scan: no hit for prefixes ending at `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Fractional-linear scan: no hit for prefixes ending at `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Two-layer fractional-linear scan: no hit for prefixes ending at `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Quotient ladder: `Q2 = GG2 / GG`, `Q3 = GG3 / GG`, `Q4 = GG4 / GG`, `Q5 = GG5 / GG`, `Q7 = GG7 / GG`, `Q11 = GG11 / GG`
- Quotient-ladder polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Quotient-ladder polynomial `total degree <= 1`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Quotient-ladder polynomial `total degree <= 2`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Quotient-ladder multiplicative scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Quotient-ladder fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Quotient-ladder two-layer fractional-linear scan: no hit for prefixes ending at `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Mixed quotient basis: `GG = GG(t)`, `Q2 = GG2 / GG`, `Q3 = GG3 / GG`, `Q4 = GG4 / GG`, `Q5 = GG5 / GG`, `Q7 = GG7 / GG`, `Q11 = GG11 / GG`
- Mixed-quotient polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Mixed-quotient polynomial `total degree <= 1`: no hit for prefixes ending at `GG`, `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Mixed-quotient polynomial `total degree <= 2`: no hit for prefixes ending at `GG`, `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Mixed-quotient multiplicative scan: no hit for prefixes ending at `GG`, `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Mixed-quotient fractional-linear scan: no hit for prefixes ending at `GG`, `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
- Mixed-quotient two-layer fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
### `S` Family

- Base benchmark: `hirschhorn_s_normalized`
- Basis ladder: `S = S(t)`, `S2 = S(t^2)`, `S3 = S(t^3)`, `S4 = S(t^4)`
- Polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Polynomial `total degree <= 1`: no hit for prefixes ending at `S`, `S2`, `S3`, `S4`.
- Polynomial `total degree <= 2`: no hit for prefixes ending at `S`, `S2`, `S3`, `S4`.
- Multiplicative scan: no hit for prefixes ending at `S`, `S2`, `S3`, `S4`.
- Fractional-linear scan: no hit for prefixes ending at `S`, `S2`, `S3`, `S4`.
- Two-layer fractional-linear scan: no hit for prefixes ending at `S2`, `S3`, `S4`.
- Quotient ladder: `Q2 = S2 / S`, `Q3 = S3 / S`, `Q4 = S4 / S`
- Quotient-ladder polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Quotient-ladder polynomial `total degree <= 1`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder polynomial `total degree <= 2`: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder multiplicative scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.
- Quotient-ladder two-layer fractional-linear scan: no hit for prefixes ending at `Q3`, `Q4`.
- Mixed quotient basis: `S = S(t)`, `Q2 = S2 / S`, `Q3 = S3 / S`, `Q4 = S4 / S`
- Mixed-quotient polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.
- Mixed-quotient polynomial `total degree <= 1`: no hit for prefixes ending at `S`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient polynomial `total degree <= 2`: no hit for prefixes ending at `S`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient multiplicative scan: no hit for prefixes ending at `S`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient fractional-linear scan: no hit for prefixes ending at `S`, `Q2`, `Q3`, `Q4`.
- Mixed-quotient two-layer fractional-linear scan: no hit for prefixes ending at `Q2`, `Q3`, `Q4`.

## Ratio-Object Source-Family Eta-Correction Scan

We also checked whether the ratio object can be written as one nearby source-family basis object times a small eta-quotient correction:

```text
F = T * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a more direct closed-form recognition lane than the polynomial / fractional-linear boxes above.
- Eta levels checked: `N=1`, `N=2`, `N=3`, `N=4`, `N=5`, `N=6`, `N=12`, `N=20`

### `RR` Eta-Correction Box

- Base benchmark: `rogers_ramanujan_normalized`
- Raw basis choices: `RR`, `RR2`, `RR3`, `RR4`
- Raw-basis eta-correction scan: no hit for basis choices `RR`, `RR2`, `RR3`, `RR4`.
- Quotient basis choices: `Q2 = RR2 / RR`, `Q3 = RR3 / RR`, `Q4 = RR4 / RR`
- Quotient-basis eta-correction scan: no hit for basis choices `Q2`, `Q3`, `Q4`.
### `cubic` Eta-Correction Box

- Base benchmark: `ramanujan_cubic_normalized`
- Raw basis choices: `cubic`, `cubic2`, `cubic3`, `cubic4`
- Raw-basis eta-correction scan: no hit for basis choices `cubic`, `cubic2`, `cubic3`, `cubic4`.
- Quotient basis choices: `Q2 = cubic2 / cubic`, `Q3 = cubic3 / cubic`, `Q4 = cubic4 / cubic`
- Quotient-basis eta-correction scan: no hit for basis choices `Q2`, `Q3`, `Q4`.
### `GG` Eta-Correction Box

- Base benchmark: `gollnitz_gordon_normalized`
- Raw basis choices: `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`
- Raw-basis eta-correction scan: no hit for basis choices `GG`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Quotient basis choices: `Q2 = GG2 / GG`, `Q3 = GG3 / GG`, `Q4 = GG4 / GG`, `Q5 = GG5 / GG`, `Q7 = GG7 / GG`, `Q11 = GG11 / GG`
- Quotient-basis eta-correction scan: no hit for basis choices `Q2`, `Q3`, `Q4`, `Q5`, `Q7`, `Q11`.
### `S` Eta-Correction Box

- Base benchmark: `hirschhorn_s_normalized`
- Raw basis choices: `S`, `S2`, `S3`, `S4`
- Raw-basis eta-correction scan: no hit for basis choices `S`, `S2`, `S3`, `S4`.
- Quotient basis choices: `Q2 = S2 / S`, `Q3 = S3 / S`, `Q4 = S4 / S`
- Quotient-basis eta-correction scan: no hit for basis choices `Q2`, `Q3`, `Q4`.

## Ratio-Object Two-Core Source-Family Eta Scan

We also checked a low-complexity hybrid source box built from two raw basis objects from different nearby families together with a small eta tail:

```text
F = T1 * T2 * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Here `T1` and `T2` come from distinct named-family raw ladders, and each source-core exponent is restricted to `±1`.
- Basis pairs checked: `132`
- Total pair-level boxes checked: `1056`
- Family-pair split: `GG×S` -> `28` pair(s), `RR×GG` -> `28` pair(s), `RR×S` -> `16` pair(s), `RR×cubic` -> `16` pair(s), `cubic×GG` -> `28` pair(s), `cubic×S` -> `16` pair(s)

No cross-family two-core eta-correction hit was found in the scanned box.


## Ratio-Object Quotient-Core Source-Family Eta Scan

We also checked a hybrid source box where one nearby family contributes a quotient core and a second family contributes one raw basis object, again with a small eta tail:

```text
F = Q * T * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Here `Q = T_k / T_1` comes from one named family, `T` comes from a different family's raw ladder, and both source-core exponents are restricted to `±1`.
- Quotient/raw basis pairs checked: `207`
- Total pair-level boxes checked: `1656`
- Quotient/raw family split: `GG->RR` -> `24` pair(s), `GG->S` -> `24` pair(s), `GG->cubic` -> `24` pair(s), `RR->GG` -> `21` pair(s), `RR->S` -> `12` pair(s), `RR->cubic` -> `12` pair(s), `S->GG` -> `21` pair(s), `S->RR` -> `12` pair(s), `S->cubic` -> `12` pair(s), `cubic->GG` -> `21` pair(s), `cubic->RR` -> `12` pair(s), `cubic->S` -> `12` pair(s)

No cross-family quotient-core eta-correction hit was found in the scanned box.


## Ratio-Object Two-Quotient-Core Source-Family Eta Scan

We also checked a quotient-only hybrid source box where two distinct nearby families each contribute one quotient core, again with a small eta tail:

```text
F = Q1 * Q2 * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Here `Q1 = T_i / T_1` and `Q2 = U_j / U_1` come from distinct named families, and both quotient-core exponents are restricted to `±1`.
- Quotient-pair basis pairs checked: `81`
- Total pair-level boxes checked: `648`
- Quotient-pair family split: `GG×S` -> `18` pair(s), `RR×GG` -> `18` pair(s), `RR×S` -> `9` pair(s), `RR×cubic` -> `9` pair(s), `cubic×GG` -> `18` pair(s), `cubic×S` -> `9` pair(s)

No cross-family two-quotient-core eta-correction hit was found in the scanned box.


## Ratio-Object Two-Quotient-Core Self-Quotient Finite-Product Scan

We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object with a compact finite-product self-quotient equation:

```text
G = F / (Q1 * Q2)
G(t) / G(t^m) = prod_{r=1}^{m-1} (1 - t^r)^{e_r}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a source-aware Mahler-style lane: a hit would point to a recursive product correction after factoring out two nearby quotient cores.
- Moduli checked: `m=2`, `m=3`, `m=4`
- Quotient-pair basis pairs checked: `81`
- Total pair-level boxes checked: `243`
- Quotient-pair family split: `GG×S` -> `18` pair(s), `RR×GG` -> `18` pair(s), `RR×S` -> `9` pair(s), `RR×cubic` -> `9` pair(s), `cubic×GG` -> `18` pair(s), `cubic×S` -> `9` pair(s)

No cross-family two-quotient-core finite-product self-quotient hit was found in the scanned box.


## Ratio-Object Two-Quotient-Core Self-Eta Functional Scan

We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-complexity self-eta functional equation:

```text
G = F / (Q1 * Q2)
G(t) = G(t^m)^a * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a product-theoretic theorem-facing lane: a hit would identify the residual correction through a recursive self-eta equation after factoring out two nearby quotient cores.
- Moduli checked: `m=2`, `m=3`, `m=4`
- Eta levels checked: `N=2`, `N=3`, `N=4`
- Quotient-pair basis pairs checked: `81`
- Total pair-level boxes checked: `729`
- Quotient-pair family split: `GG×S` -> `18` pair(s), `RR×GG` -> `18` pair(s), `RR×S` -> `9` pair(s), `RR×cubic` -> `9` pair(s), `cubic×GG` -> `18` pair(s), `cubic×S` -> `9` pair(s)

No cross-family two-quotient-core self-eta functional-equation hit was found in the scanned box.


## Ratio-Object Two-Quotient-Core Self-Fractional-Linear Scan

We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-complexity self-fractional-linear equation with a small eta tail:

```text
G = F / (Q1 * Q2)
G(t) = (1 + a*(G(t^m) - 1) + ... ) / (1 + b*(G(t^m) - 1) + ... )
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a theorem-facing nonlinear lane: a hit would give a compact recursive rational equation for the residual correction after factoring out two nearby quotient cores.
- Moduli checked: `m=2`, `m=3`, `m=4`
- Eta levels checked: `N=2`, `N=3`, `N=4`
- Quotient-pair basis pairs checked: `81`
- Total pair-level boxes checked: `729`
- Quotient-pair family split: `GG×S` -> `18` pair(s), `RR×GG` -> `18` pair(s), `RR×S` -> `9` pair(s), `RR×cubic` -> `9` pair(s), `cubic×GG` -> `18` pair(s), `cubic×S` -> `9` pair(s)

No cross-family two-quotient-core self-fractional-linear hit was found in the scanned box.


## Ratio-Object Two-Quotient-Core Self-Polynomial Functional Scan

We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-degree algebraic self-functional equation:

```text
G = F / (Q1 * Q2)
P(G(t), G(t^m)) = 0
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a theorem-facing lane: a hit would suggest a compact defining functional equation for the residual correction after factoring out two nearby quotient cores.
- Moduli checked: `m=2`, `m=3`, `m=4`
- Degrees checked: `total degree <= 1`, `total degree <= 2`
- Quotient-pair basis pairs checked: `81`
- Total pair-level boxes checked: `486`
- Quotient-pair family split: `GG×S` -> `18` pair(s), `RR×GG` -> `18` pair(s), `RR×S` -> `9` pair(s), `RR×cubic` -> `9` pair(s), `cubic×GG` -> `18` pair(s), `cubic×S` -> `9` pair(s)

No cross-family two-quotient-core self-polynomial functional-equation hit was found in the scanned box.


## Ratio-Object Explicit GG/S Transform Template Scan

We also checked a smaller family-meaning-preserving box tailored to the Gordon/Hirschhorn orbit:

```text
F = T
F = 1 / T
F = T_i / T_j
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Families checked here are the literature-family ladders `GG` and `S`.
- This does not enlarge the algebraic search box much; it makes the reciprocal / quotient interpretations explicit in the note.

### `GG` Explicit Transform Box

- Base benchmark: `gollnitz_gordon_normalized`
- Basis ladder: `GG = GG(t)`, `GG2 = GG(t^2)`, `GG3 = GG(t^3)`, `GG4 = GG(t^4)`, `GG5 = GG(t^5)`, `GG7 = GG(t^7)`, `GG11 = GG(t^11)`
- Templates checked: `56` exact direct / reciprocal / quotient templates.
- No exact direct / reciprocal / quotient template hit was found in this family.

### `S` Explicit Transform Box

- Base benchmark: `hirschhorn_s_normalized`
- Basis ladder: `S = S(t)`, `S2 = S(t^2)`, `S3 = S(t^3)`, `S4 = S(t^4)`
- Templates checked: `20` exact direct / reciprocal / quotient templates.
- No exact direct / reciprocal / quotient template hit was found in this family.


## Ratio-Object GG Modular-Equation Template Scan

We also checked a narrower literature-driven `GG` box motivated by the modular-equation papers of Chan--Huang and Cho--Koo--Park:

```text
F = T
F = 1 / T
F = T_i / T_j
P(F, T_i) = 0
F = prod_i T_i^e_i
F = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Base benchmark: `gollnitz_gordon_normalized`
- This lane keeps the sign and substitution objects explicit instead of flattening them into a larger anonymous basis box.
- The literature-motivated basis here starts with `GG(t)`, `GG(-t)`, `GG(t^2)`, `GG(t^3)`, and `GG(t^4)`, and in the full profile it also includes the odd-prime descendants suggested by the GG modular-equation papers.
- We also run a mixed quotient-coordinate pass that keeps `GG(t)` explicit while letting the correction move in quotient coordinates such as `GG(-t)/GG(t)` and `GG(t^p)/GG(t)`.
- Basis ladder: `GG = GG(t)`, `GGneg = GG(-t)`, `GG2 = GG(t^2)`, `GG3 = GG(t^3)`, `GG4 = GG(t^4)`, `GG5 = GG(t^5)`, `GG7 = GG(t^7)`, `GG11 = GG(t^11)`
- Exact direct / reciprocal / quotient templates checked: `72`.
- No exact direct / reciprocal / quotient template hit was found in this modular-equation box.
- Exact literature polynomial templates checked: `Chan--Huang Cor. 3.2(i) on (F, GG3)`, `Chan--Huang Cor. 3.2(ii) on (F, GG4)`.
- No exact Chan--Huang direct modular-equation polynomial hit was found.

- Polynomial scan: no candidate-dependent hit was found in the checked modular-equation prefixes.
- Polynomial `total degree <= 1`: no hit for prefixes ending at `GG`, `GGneg`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Polynomial `total degree <= 2`: no hit for prefixes ending at `GG`, `GGneg`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Multiplicative scan: no hit for modular-equation prefixes ending at `GG`, `GGneg`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Fractional-linear scan: no hit for modular-equation prefixes ending at `GG`, `GGneg`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Two-layer fractional-linear scan: no hit for modular-equation prefixes ending at `GGneg`, `GG2`, `GG3`, `GG4`, `GG5`, `GG7`, `GG11`.
- Quotient basis: `Q_neg = GG(-t) / GG(t)`, `Q_2 = GG(t^2) / GG(t)`, `Q_3 = GG(t^3) / GG(t)`, `Q_4 = GG(t^4) / GG(t)`, `Q_5 = GG(t^5) / GG(t)`, `Q_7 = GG(t^7) / GG(t)`, `Q_11 = GG(t^11) / GG(t)`
- Exact quotient-coordinate literature templates checked: `Chan--Huang Cor. 3.2(i) on (F, Q_3)`, `Chan--Huang Cor. 3.2(ii) on (F, Q_4)`.
- No exact Chan--Huang quotient-coordinate modular-equation polynomial hit was found.
- Quotient-coordinate polynomial scan: no candidate-dependent hit was found in the checked quotient prefixes.
- Quotient-coordinate polynomial `total degree <= 1`: no hit for prefixes ending at `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Quotient-coordinate polynomial `total degree <= 2`: no hit for prefixes ending at `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Quotient-coordinate multiplicative scan: no hit for prefixes ending at `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Quotient-coordinate fractional-linear scan: no hit for prefixes ending at `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Quotient-coordinate two-layer fractional-linear scan: no hit for prefixes ending at `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Mixed quotient basis: `GG = GG(t)`, `Q_neg = GG(-t) / GG(t)`, `Q_2 = GG(t^2) / GG(t)`, `Q_3 = GG(t^3) / GG(t)`, `Q_4 = GG(t^4) / GG(t)`, `Q_5 = GG(t^5) / GG(t)`, `Q_7 = GG(t^7) / GG(t)`, `Q_11 = GG(t^11) / GG(t)`
- Mixed quotient-coordinate polynomial scan: no candidate-dependent hit was found in the checked prefixes.
- Mixed quotient-coordinate polynomial `total degree <= 1`: no hit for prefixes ending at `GG`, `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Mixed quotient-coordinate polynomial `total degree <= 2`: no hit for prefixes ending at `GG`, `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Mixed quotient-coordinate multiplicative scan: no hit for prefixes ending at `GG`, `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Mixed quotient-coordinate fractional-linear scan: no hit for prefixes ending at `GG`, `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.
- Mixed quotient-coordinate two-layer fractional-linear scan: no hit for prefixes ending at `Q_neg`, `Q_2`, `Q_3`, `Q_4`, `Q_5`, `Q_7`, `Q_11`.


## Ratio-Object Explicit GG/S Template Eta-Correction Scan

We also checked whether one explicit Gordon/Hirschhorn-orbit template times a small eta tail explains the ratio object:

```text
F = T * prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- Here `T` ranges over the exact direct / reciprocal / quotient templates from the preceding GG/S transform box.
- Eta levels checked: `N=1`, `N=2`, `N=3`, `N=4`, `N=5`, `N=6`, `N=12`, `N=20`

### `GG` Explicit Eta-Correction Box

- Base benchmark: `gollnitz_gordon_normalized`
- Templates checked: `56` explicit direct / reciprocal / quotient templates.
- No explicit-template eta-correction hit was found in this family.

### `S` Explicit Eta-Correction Box

- Base benchmark: `hirschhorn_s_normalized`
- Templates checked: `20` explicit direct / reciprocal / quotient templates.
- No explicit-template eta-correction hit was found in this family.


## Ratio-Object RR-Tower Prefix Scan

We also scanned the multiplicative correction object against prefixes of the benchmark tower:

- `F = candidate / rogers_ramanujan_q3_normalized`
- `B1 = rogers_ramanujan_q3_normalized`
- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

- Prefixes checked: `(F, B1, B2)`, then `(F, B1, B2, B3)`, and so on through the final listed power.
- Degrees scanned here are intentionally low (`1` and `2`) to stay in candidate-dependent, well-determined boxes.

No candidate-dependent relation was found for the ratio object in any scanned prefix box.

- `total degree <= 1`: no hit for ratio-object prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.
- `total degree <= 2`: no hit for ratio-object prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.

## Ratio-Object Self-Quotient Finite-Product Scan

We also checked a simple finite-product self-quotient box for the ratio object:

```text
F(t) / F(t^m) = prod_{r=1}^{m-1} (1 - t^r)^{e_r}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a Mahler-style finite-product functional equation: a hit would give a compact recursive product description, but a miss does not rule out general q-Pochhammer products.

No finite-product self-quotient relation was found in any scanned modulus.

- No hit for moduli `m=2`, `m=3`, `m=4`, `m=5`, `m=6`, `m=12`, `m=20`.

## Ratio-Object Eta-Quotient Scan

We also checked whether the ratio object itself is already a small-level eta-quotient:

```text
F = prod_{d|N} (t^d; t^d)_inf^{e_d}
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- This is a direct closed-form recognition lane rather than another transform-elimination box.

No eta-quotient relation was found in any scanned level.

- No hit for eta levels `N=1`, `N=2`, `N=3`, `N=4`, `N=5`, `N=6`, `N=12`, `N=20`.

## Ratio-Object Multiplicative RR-Tower Scan

We also searched for exact multiplicative corrections built from the benchmark tower:

```text
F = prod_i B_i^e_i
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- `B1 = rogers_ramanujan_q3_normalized`
- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.
- Exponents are solved exactly from the log-series constraints, then verified by exact series re-expansion.

No multiplicative ratio-object relation was found in any scanned prefix box.

- No hit for multiplicative prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.

## Ratio-Object Fractional-Linear RR-Tower Scan

We also searched for low-complexity fractional-linear corrections built from the benchmark tower:

```text
F = (1 + sum a_i*(B_i - 1)) / (1 + sum b_i*(B_i - 1))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- `B1 = rogers_ramanujan_q3_normalized`
- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.
- Each prefix solves an exact linear system for the numerator and denominator correction coefficients.

No fractional-linear ratio-object relation was found in any scanned prefix box.

- No hit for fractional-linear prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.

## Ratio-Object Two-Layer Fractional-Linear RR-Tower Scan

We then expanded to a second-ring nonlinear box built from two single-basis fractional-linear factors:

```text
F = ((1 + a0*(X1 - 1)) / (1 + b0*(Y1 - 1))) * ((1 + a1*(X2 - 1)) / (1 + b1*(Y2 - 1)))
```

- `F = candidate / rogers_ramanujan_q3_normalized`
- `B1 = rogers_ramanujan_q3_normalized`
- `B2 = B1(t^2)`
- `B3 = B1(t^3)`
- `B4 = B1(t^4)`
- `B5 = B1(t^5)`
- `B6 = B1(t^6)`
- `B12 = B1(t^12)`
- `B20 = B1(t^20)`

- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.
- Each prefix scans low-complexity two-factor templates and verifies any candidate hit by exact series re-expansion.

No two-layer fractional-linear ratio-object relation was found in any scanned prefix box.

- No hit for two-layer fractional-linear prefixes ending at `B2`, `B3`, `B4`, `B5`, `B6`, `B12`, `B20`.
