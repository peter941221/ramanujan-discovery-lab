# `RR(q^3)` Neighborhood Audit

Audit date: `2026-03-14`

## Scope

This note systematizes the current public neighborhood around the classical
benchmark

```text
RR(q^3) = 1 / (1 + K q^(3 + 3(n-1)) / 1).
```

Unlike the `RR(q^4)` side, the visible `q^3` neighborhood does not collapse to
one plain monotone step ladder.

## Public Members In The Current Snapshot

| Candidate | Structure | Shared digits vs `RR(q^3)` | First divergence |
| --- | --- | ---: | ---: |
| `05fe7657191e` | canonical `RR(q^3)` | `152` | exact |
| `cb60fd71d1d7` | RR main ladder + cubic extra numerator + denominator perturbation | `7` | `q^12` |
| `e2cc74240b6f` | single-ladder step perturbation `3 -> 4` | `5` | `q^9` |

## Branch Split

The current public `RR(q^3)` neighborhood splits into two qualitatively
different branches.

### Branch A: Plain Step Perturbation

`e2cc74240b6f` keeps the pure RR reciprocal shape and changes only one
parameter:

```text
RR(q^3):      1 / (1 + K q^(3 + 3(n-1)) / 1)
e2cc74240b6f: 1 / (1 + K q^(3 + 4(n-1)) / 1)
```

So it is the simple opposite-side step perturbation of the benchmark.

A dedicated focused note for this branch now lives at:

- `E2CC74240B6F_PLAIN_STEP_AUDIT.md`

### Branch B: Hybrid Perturbation

`cb60fd71d1d7` keeps the same main step `3`, but adds two extra structural
features:

- cubic-style extra numerator `q^(6n)`
- denominator perturbation `1 + q^(3n)`

In reduced `t = q^3` form, its reciprocal is

```text
1 + K (t^n + t^(2n)) / (1 + t^n),
```

which is qualitatively different from the plain RR ladder.

## Low-Order Fits

- `cb60fd71d1d7`
  - `candidate / RR(q^3) = 1 - q^12 + 2 q^15 - q^18`
- `e2cc74240b6f`
  - `candidate / RR(q^3) = 1 - q^9 + q^10 + q^12`

So the hybrid branch agrees longer and diverges later.

## Ratio-Series Signals

- `cb60fd71d1d7 / RR(q^3)`
  - residue support mod `3`: `{0}`
  - visible ratio is a pure function of `t = q^3`
- `e2cc74240b6f / RR(q^3)`
  - residue support mod `3`: `{0, 1, 2}`
  - no clean step reduction appears in the visible window

This is the sharpest local separator inside the current `q^3` neighborhood.

## Same-Step Cross-Family Signal

Only the hybrid branch reuses the cubic-style numerator pattern. That is why
`cb60fd71d1d7` naturally triggers transform and contraction questions around
both the RR and cubic neighborhoods.

For `e2cc74240b6f`, the current low-order comparison is simpler:

- divergence vs `RR(q^3)`: `q^9`
- divergence vs `cubic(q^3)`: `q^6`

So even nearby cubic comparisons do not make it look hybrid.

## Interpretation

The `RR(q^3)` side is not currently a one-parameter public ladder analogous to
the `RR(q^4)` side.

Instead it looks like:

```text
RR(q^3) neighborhood
├─ plain step branch
│  └─ e2cc74240b6f
└─ hybrid branch
   └─ cb60fd71d1d7
```

That is exactly why `cb60fd71d1d7` remains the hero case:

- deeper visible agreement
- pure step-reduced ratio signal
- richer mixed-family structure
- stronger transform pressure than the plain step branch

## Ranking Impact

This audit keeps the current local ordering:

1. `cb60fd71d1d7`
2. `e2cc74240b6f`

Inside the `q^3` neighborhood, the hybrid branch is clearly the more valuable
manual-audit target.
