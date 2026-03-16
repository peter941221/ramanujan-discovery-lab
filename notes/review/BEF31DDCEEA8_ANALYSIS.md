# Hero Case Analysis: `bef31ddceea8`

## Snapshot

- Candidate id: `bef31ddceea8`
- Closest benchmark: `rogers_ramanujan_q4_normalized` (7 shared digits)
- Benchmark kind: `classical_family`
- Novelty status: `review`
- Family bucket: `single_numerator_family::num_scale=1::extra_ratio=0::den_scale=0`
- Equivalence key: `review::closest=rogers_ramanujan_q4_normalized::num_scale=1::main_shift_delta=0::main_step_delta=-1::extra_exists=0::extra_scale=0::extra_shift_delta=0::extra_step_delta=0::den_exists=0::den_scale=0::den_shift_delta=0::den_step_delta=0`

## Structural Delta vs Closest Benchmark

- Candidate template: `qcf:top=1;base=1;ns=1;nr=4;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=4;nk=4;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Main numerator shift delta: `0`
- Main numerator step delta: `-1`
- Extra numerator present: `False`
- Denominator perturbation present: `False`

## Symbolic q-Series

- Candidate series (depth `12`, order `31`):
```text
10*q^30 - 7*q^29 + q^28 + 5*q^27 - 6*q^26 + 3*q^25 + q^24 - 4*q^23 + 3*q^22 - q^21 - q^20 + 3*q^19 - q^18 + q^16 - 2*q^15 - q^12 + q^11 + q^8 - q^4 + 1
```
- Benchmark series `rogers_ramanujan_q4_normalized`:
```text
q^28 - q^24 + q^20 - q^16 + q^8 - q^4 + 1
```
- Candidate minus benchmark:
```text
10*q^30 - 7*q^29 + 5*q^27 - 6*q^26 + 3*q^25 + 2*q^24 - 4*q^23 + 3*q^22 - q^21 - 2*q^20 + 3*q^19 - q^18 + 2*q^16 - 2*q^15 - q^12 + q^11
```
- First divergence order: `11`

## Low-Order Multiplicative Fit

- Fitted modifier multiplying the closest benchmark:
```text
-q^15 - q^12 + q^11 + 1
```
- Truncated residual after this fit:
```text
0
```
