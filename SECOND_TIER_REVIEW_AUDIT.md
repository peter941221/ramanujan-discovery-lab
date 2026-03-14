# Second-Tier Review Audit

Audit date: `2026-03-14`

## Scope

This note starts the next layer of manual audit below the current hero case
`cb60fd71d1d7`.

Current focus:

- `bef31ddceea8`
- `e2cc74240b6f`

These are the two strongest non-hero public `review` candidates that still sit
near classical Rogers-Ramanujan family benchmarks rather than internal fixtures.

## Snapshot

| Candidate | Closest benchmark | Shared digits | First divergence | Structure |
| --- | --- | ---: | ---: | --- |
| `bef31ddceea8` | `rogers_ramanujan_q4_normalized` | 7 | `q^11` | single-ladder step perturbation: start `q^4`, step `4 -> 3` |
| `e2cc74240b6f` | `rogers_ramanujan_q3_normalized` | 5 | `q^9` | single-ladder step perturbation: start `q^3`, step `3 -> 4` |

## Structural Read

- Both candidates are much simpler than `cb60fd71d1d7`.
- Neither candidate introduces:
  - an extra numerator branch
  - a denominator perturbation
- Both are plain single-numerator step perturbations of nearby RR family
  benchmarks.

In template form:

```text
bef31ddceea8: 1 / (1 + K q^(4 + 3(n-1)) / 1)
RR(q^4):      1 / (1 + K q^(4 + 4(n-1)) / 1)

e2cc74240b6f: 1 / (1 + K q^(3 + 4(n-1)) / 1)
RR(q^3):      1 / (1 + K q^(3 + 3(n-1)) / 1)
```

So the leading hypothesis for both is still "catalog-adjacent step ladder"
rather than "structurally hybrid transformed family."

## Low-Order Fits

- `bef31ddceea8`
  - `candidate / RR(q^4) = 1 + q^11 - q^12 - q^15`
- `e2cc74240b6f`
  - `candidate / RR(q^3) = 1 - q^9 + q^10 + q^12`

These fits explain the first visible divergence cleanly, but they do not show
the special step-reduced purity seen in `cb60fd71d1d7`.

## Ratio-Series Signals

From the current `research` notes:

- `bef31ddceea8 / RR(q^4)`
  - residue support mod `4`: `{0, 1, 2, 3}`
  - no clean reduction to a pure function of `t = q^4`
- `e2cc74240b6f / RR(q^3)`
  - residue support mod `3`: `{0, 1, 2}`
  - no clean reduction to a pure function of `t = q^3`

This matters because the current hero case is different in exactly this place:

```text
cb60fd71d1d7 / RR(q^3)
```

stays inside one residue class and becomes a pure `t = q^3` series in the
visible window.

## Cross-Family Signal

Only `e2cc74240b6f` has a nearby same-step cubic comparison in the current
analysis pass:

- divergence vs `RR(q^3)`: `q^9`
- divergence vs `cubic(q^3)`: `q^6`

So even in its own `q^3` neighborhood, it still reads more like a perturbed RR
ladder than like a mixed RR/cubic hybrid.

## Current Ranking Impact

This first second-tier audit keeps the current ordering intact:

1. `cb60fd71d1d7`
2. `bef31ddceea8`
3. `e2cc74240b6f`

Reason:

- `bef31ddceea8` keeps the stronger digit tier, but it looks structurally plain.
- `e2cc74240b6f` is also structurally plain and has shallower agreement.
- neither candidate currently shows the same transform-pressure signal,
  step-reduced purity, or mixed-family structure that keeps `cb60...` in front.

## `RR(q^4)` Ladder Update

The `bef31ddceea8` branch now has a dedicated ladder note:

- `RR_Q4_STEP_LADDER_AUDIT.md`

That audit confirms the current public `RR(q^4)` neighborhood is a monotone
single-ladder family with steps

```text
4 -> 3 -> 2 -> 1
```

and with `bef31ddceea8` as the strongest noncanonical member of that ladder.

## Immediate Next Checks If Promoted

- `bef31ddceea8`
  - extend the current `RR(q^4)` ladder note if more step-family members survive
    future discovery runs
- `e2cc74240b6f`
  - test whether the `3 -> 4` step change is just the symmetric opposite-side
    ladder member of the same RR neighborhood phenomenon
- both
  - only start heavier transform audit if they survive after the hero-case path
    is exhausted or stalls
