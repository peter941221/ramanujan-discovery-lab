# `RR(q^4)` Step-Ladder Audit

Audit date: `2026-03-14`

## Scope

This note systematizes the `RR(q^4)` single-ladder perturbation phenomenon that
contains `bef31ddceea8`.

The family visible in the current verified snapshot is:

```text
1 / (1 + K q^(4 + k(n-1)) / 1)
```

with fixed start exponent `4` and step `k` running down the ladder

```text
4 -> 3 -> 2 -> 1.
```

## Ladder Members

| Candidate | Step `k` | Status | Shared digits vs `RR(q^4)` | First divergence |
| --- | ---: | --- | ---: | ---: |
| `0f617e2863fe` | `4` | `known` | `155` | exact |
| `bef31ddceea8` | `3` | `review` | `7` | `q^11` |
| `e42a0d5f2679` | `2` | `review` | `6` | `q^10` |
| `9dbafd59364c` | `1` | `review` | `5` | `q^9` |

So the public `RR(q^4)` review candidates are not isolated accidents. They sit
on a monotone one-parameter step ladder.

## Structural Invariant

All four members share the same shape:

- same top term
- same base denominator
- same numerator start exponent `q^4`
- no extra numerator branch
- no denominator perturbation

Only one parameter moves:

```text
numerator_q_step = 4, 3, 2, 1.
```

That makes this family structurally much simpler than the current hero case
`cb60fd71d1d7`.

## Low-Order Multiplicative Fits

- step `4`:
  - canonical benchmark, so no correction term
- step `3` / `bef31ddceea8`:
  - `candidate / RR(q^4) = 1 + q^11 - q^12 - q^15`
- step `2` / `e42a0d5f2679`:
  - `candidate / RR(q^4) = 1 + q^10 - q^12 - q^14`
- step `1` / `9dbafd59364c`:
  - `candidate / RR(q^4) = 1 + q^9 - q^12 - q^13`

The first divergence order drops monotonically as the step moves farther away
from the benchmark step `4`:

```text
step 4: exact
step 3: q^11
step 2: q^10
step 1: q^9
```

## Ratio-Series Congruence Signals

Residue support for `candidate / RR(q^4)` in the visible symbolic window:

- step `3` / `bef31ddceea8`:
  - residues mod `4`: `{0, 1, 2, 3}`
- step `2` / `e42a0d5f2679`:
  - residues mod `4`: `{0, 2}`
- step `1` / `9dbafd59364c`:
  - residues mod `4`: `{0, 1, 2, 3}`

The step-`2` member is the only one with an extra congruence signal: its ratio
stays even in `q`, so it is a function of `q^2` in the visible window. But it
still remains a plain single-ladder perturbation, not a mixed-family hybrid.

## Interpretation

This ladder currently reads like a catalog-adjacent perturbation family rather
than a standout transform-heavy lead.

Why:

- the whole family is explained by one monotone structural move
- every member keeps the same simple RR reciprocal shape
- there is no extra numerator branch and no denominator deformation
- the visible agreement degrades smoothly as the step drifts away from `4`

So `bef31ddceea8` remains the strongest member of this ladder, but the ladder as
a whole still looks more like a controlled neighborhood phenomenon than a
headline discovery candidate.

## Ranking Impact

This strengthens the current ranking choice:

1. `cb60fd71d1d7`
2. `bef31ddceea8`
3. `e2cc74240b6f`

`bef31ddceea8` stays ahead of the other second-tier candidates because it is the
best-scoring member of the clean `RR(q^4)` step ladder, but it still trails the
hero case because the ladder itself is structurally plain.
