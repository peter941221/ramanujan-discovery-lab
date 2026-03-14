# Research Note: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized` (7 shared digits)
- Candidate template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=1;ner=6;nek=6;dc=1;ds=1;dr=3;dk=3`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Depth: `40`
- Series order request: `151`
- Raw q-series computed order: `61`

## Ratio Series

- Computed `candidate / rogers_ramanujan_q3_normalized` as a truncated q-series.
- Residue support mod step `3` (excluding constant term): `[0]`

```text
39*q^39 - 14*q^36 - 2*q^33 + 10*q^30 - 9*q^27 + 4*q^24 - q^21 - q^18 + 2*q^15 - q^12 + 1
```

## Step-Reduced View (t-Series)

- Reduced variable: `t = q^3`
- Reduced candidate template: `qcf:top=1;base=1;ns=1;nr=1;nk=1;nes=1;ner=2;nek=2;dc=1;ds=1;dr=1;dk=1`
- Reduced benchmark template: `qcf:top=1;base=1;ns=1;nr=1;nk=1;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`

### Ratio(t)

```text
630109*t^39 + 142254*t^38 - 432394*t^37 + 381487*t^36 - 197523*t^35 + 30915*t^34 + 57631*t^33 - 73508*t^32 + 49694*t^31 - 17910*t^30 - 4022*t^29 + 12247*t^28 - 10840*t^27 + 5633*t^26 - 874*t^25 - 1642*t^24 + 2076*t^23 - 1411*t^22 + 518*t^21 + 117*t^20 - 352*t^19 + 304*t^18 - 159*t^17 + 28*t^16 + 47*t^15 - 61*t^14 + 39*t^13 - 14*t^12 - 2*t^11 + 10*t^10 - 9*t^9 + 4*t^8 - t^7 - t^6 + 2*t^5 - t^4 + 1
```

### Euler Product Exponents

- Using the truncated Euler transform: `Ratio(t) = Π_{n>=1} (1 - t^n)^{c_n}` modulo the visible order.
- Nonzero count in first 30 exponents: `27`
- Max |c_n| among first 30 exponents (rough): `4728`

```text
c_1..c_30:
0, 0, 0, 1, -2, 1, 1, -4, 7, -6, 1, 8, -23, 34, -25, -12, 75, -138, 150, -49, -197, 516, -725, 556, 260, -1720, 3240, -3562, 1168, 4728
```

### Structured Fit Attempts

- Periodic Pochhammer fit (bounded |e_r|<=8, period<=12): `no fit`
- Eta-quotient fit (bounded |e_d|<=8, level<=12): `no fit`

## Direct 1-Step Bauer-Muir Obstruction

- `RR reciprocal`: matching `b0` forces `w0 = 0` and matching `b1` forces `w1 = t`. That leaves the first transformed numerator at `t` instead of the target `t^2 + t`, so no direct 1-step transform exists.
- `cubic reciprocal`: the forced first-stage data `w0 = 0`, `w1 = t` does preserve the first numerator, but matching `b2` also forces `w2 = t^2`. Then the second transformed numerator becomes `t^4 - t^3 + t^2 - t` instead of `t^4 + t^2`, so this direct path is also ruled out.

## Heine hcf2 Specialization Check (c=bz=-1)

- In hcf2, setting `c = b z = -1` forces denominators `b_n = 1 + t^n` and `b0 = 1`.
- Under the standardized coefficient extraction used here, the first numerator term becomes a constant in `t`:

```text
a1_hcf2 = -a/2 - a/(2*b) + 1/2 + 1/(2*b)
```

- But the target reciprocal for this candidate has `a1 = t + t^2`, so this specialization cannot match at the coefficient level.

## Constrained Bauer-Muir Search

- Match target: reciprocal coefficients through depth `4`, checked at exact rational sample points `t = 1/10, 1/7, 1/5`.
- Pattern family per step: `13` low-complexity modifiers (`0`, `±(t^n-1)`, `±(t^(2n)-1)`, `±(t^n-t^(2n))`, with scales `1` or `2`).
- RR source, 1-step search space: `13`; hits: `0`
- RR source, 2-step search space: `169`; hits: `0`
- RR source, 3-step search space: `2197`; hits: `0`
- Cubic source, 1-step search space: `13`; hits: `0`
- Cubic source, 2-step search space: `169`; hits: `0`
- Cubic source, 3-step search space: `2197`; hits: `0`
