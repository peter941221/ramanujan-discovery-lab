# `e2cc74240b6f` Plain-Step Audit

Audit date: `2026-03-14`

## Snapshot

- Candidate: `e2cc74240b6f`
- Closest benchmark: `rogers_ramanujan_q3_normalized`
- Shared digits: `5`
- First divergence vs `RR(q^3)`: `q^9`
- First divergence vs `cubic(q^3)`: `q^6`

## Structural Read

This candidate is the plain opposite-side step perturbation of `RR(q^3)`.

Benchmark:

```text
RR(q^3) = 1 / (1 + K q^(3 + 3(n-1)) / 1)
```

Candidate:

```text
e2cc74240b6f = 1 / (1 + K q^(3 + 4(n-1)) / 1)
```

So only one structural parameter changes:

```text
numerator_q_step: 3 -> 4
```

What does not change:

- numerator start exponent stays at `q^3`
- no extra numerator branch appears
- no denominator perturbation appears
- no same-step cubic numerator pattern appears

## Low-Order Fit

The current visible correction is

```text
candidate / RR(q^3) = 1 - q^9 + q^10 + q^12.
```

This is a clean small perturbation, but it does not show the deeper delayed
agreement or the pure step-reduced ratio signal seen in `cb60fd71d1d7`.

## Ratio-Series Signal

In the current symbolic window:

- residue support mod `3` is `{0, 1, 2}`
- so the ratio is not a pure function of `t = q^3`

That is the opposite of the current hero case:

```text
cb60fd71d1d7 / RR(q^3)
```

which stays in residue class `0 mod 3` and becomes a pure `t = q^3` series in
the visible window.

## Why This Matters

Right now the simplest reading is also the strongest reading:

```text
e2cc74240b6f
├─ is interesting
├─ is stable
└─ still looks like a plain RR step perturbation
```

That makes it a legitimate second-tier review candidate, but not the best
manual-audit target in the current `q^3` neighborhood.

## Current Conclusion

- keep `e2cc74240b6f` in the public review set
- treat it as the plain step branch of the current `RR(q^3)` neighborhood
- keep `cb60fd71d1d7` ahead of it as the richer hybrid branch
