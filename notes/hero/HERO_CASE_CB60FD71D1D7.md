# Hero Case Analysis: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized` (7 shared digits)
- Benchmark kind: `classical_family`
- Novelty status: `review`
- Family bucket: `hybrid_perturbed_family::num_scale=1::extra_ratio=2::den_scale=1`
- Equivalence key: `review::closest=rogers_ramanujan_q3_normalized::num_scale=1::main_shift_delta=0::main_step_delta=0::extra_exists=1::extra_scale=1::extra_shift_delta=3::extra_step_delta=3::den_exists=1::den_scale=1::den_shift_delta=3::den_step_delta=3`

## Structural Delta vs Closest Benchmark

- Candidate template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=1;ner=6;nek=6;dc=1;ds=1;dr=3;dk=3`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Main numerator shift delta: `0`
- Main numerator step delta: `0`
- Extra numerator present: `True`
- Denominator perturbation present: `True`

## Symbolic q-Series

- Candidate series (depth `12`, order `31`):
```text
29*q^30 - 18*q^27 + 5*q^24 + 3*q^21 - 5*q^18 + 4*q^15 - 2*q^12 + q^6 - q^3 + 1
```
- Benchmark series `rogers_ramanujan_q3_normalized`:
```text
2*q^30 - q^27 + q^21 - q^18 + q^15 - q^12 + q^6 - q^3 + 1
```
- Candidate minus benchmark:
```text
27*q^30 - 17*q^27 + 5*q^24 + 2*q^21 - 4*q^18 + 3*q^15 - q^12
```
- First divergence order: `12`

## Low-Order Multiplicative Fit

- Fitted modifier multiplying the closest benchmark:
```text
-q^18 + 2*q^15 - q^12 + 1
```
- Truncated residual after this fit:
```text
0
```

## Same-Step Cross-Family Comparison

- Alternative benchmark: `ramanujan_cubic_q3_normalized`
- Description: Ramanujan cubic continued fraction benchmark evaluated at q^3.
```text
-8*q^30 + 8*q^27 - q^24 - 4*q^21 + 4*q^18 - q^15 - 2*q^12 + 2*q^9 - q^3 + 1
```
- Candidate minus alternative benchmark:
```text
37*q^30 - 26*q^27 + 6*q^24 + 7*q^21 - 9*q^18 + 5*q^15 - 2*q^9 + q^6
```
- First divergence order vs `ramanujan_cubic_q3_normalized`: `6`
