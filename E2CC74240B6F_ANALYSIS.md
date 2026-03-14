# Hero Case Analysis: `e2cc74240b6f`

## Snapshot

- Candidate id: `e2cc74240b6f`
- Closest benchmark: `rogers_ramanujan_q3_normalized` (5 shared digits)
- Benchmark kind: `classical_family`
- Novelty status: `review`
- Family bucket: `single_numerator_family::num_scale=1::extra_ratio=0::den_scale=0`
- Equivalence key: `review::closest=rogers_ramanujan_q3_normalized::num_scale=1::main_shift_delta=0::main_step_delta=1::extra_exists=0::extra_scale=0::extra_shift_delta=0::extra_step_delta=0::den_exists=0::den_scale=0::den_shift_delta=0::den_step_delta=0`

## Structural Delta vs Closest Benchmark

- Candidate template: `qcf:top=1;base=1;ns=1;nr=3;nk=4;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Main numerator shift delta: `0`
- Main numerator step delta: `1`
- Extra numerator present: `False`
- Denominator perturbation present: `False`

## Symbolic q-Series

- Candidate series (depth `12`, order `31`):
```text
15*q^30 - 15*q^29 + 9*q^28 - 8*q^27 + 10*q^26 - 6*q^25 + 4*q^24 - 6*q^23 + 5*q^22 - 2*q^21 + 3*q^20 - 4*q^19 + q^18 - q^17 + 3*q^16 - q^15 - 2*q^13 + q^12 + q^10 - q^9 + q^6 - q^3 + 1
```
- Benchmark series `rogers_ramanujan_q3_normalized`:
```text
2*q^30 - q^27 + q^21 - q^18 + q^15 - q^12 + q^6 - q^3 + 1
```
- Candidate minus benchmark:
```text
13*q^30 - 15*q^29 + 9*q^28 - 7*q^27 + 10*q^26 - 6*q^25 + 4*q^24 - 6*q^23 + 5*q^22 - 3*q^21 + 3*q^20 - 4*q^19 + 2*q^18 - q^17 + 3*q^16 - 2*q^15 - 2*q^13 + 2*q^12 + q^10 - q^9
```
- First divergence order: `9`

## Low-Order Multiplicative Fit

- Fitted modifier multiplying the closest benchmark:
```text
q^12 + q^10 - q^9 + 1
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
23*q^30 - 15*q^29 + 9*q^28 - 16*q^27 + 10*q^26 - 6*q^25 + 5*q^24 - 6*q^23 + 5*q^22 + 2*q^21 + 3*q^20 - 4*q^19 - 3*q^18 - q^17 + 3*q^16 - 2*q^13 + 3*q^12 + q^10 - 3*q^9 + q^6
```
- First divergence order vs `ramanujan_cubic_q3_normalized`: `6`
