# Identification Note: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized`
- Variable view: `t = q^3`
- Depth: `40`
- Series order: `90`
- Polynomial relation search: total degree `<= 4`

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
underdetermined polynomial relation search: 715 monomials > 90 constraints (increase order, lower max_total_degree, or reduce variables)
```

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
