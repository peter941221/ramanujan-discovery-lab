# `q^3 / q^4` Review Comparison

Audit date: `2026-03-13`

## Scope

This note compares the public `review` candidates whose nearest built-in
benchmarks live in the classical `RR(q^3)` or `RR(q^4)` neighborhoods.

Included candidates:

- `cb60fd71d1d7`
- `bef31ddceea8`
- `e42a0d5f2679`
- `9dbafd59364c`
- `e2cc74240b6f`

The fixture-adjacent candidate `1125ffe48b3b` is listed separately at the end,
because its nearest benchmark is the internal `shifted_rr_fixture`, not a
classical public family.

## Snapshot

| Candidate | Closest benchmark | Shared digits | First divergence | Structure |
| --- | --- | ---: | ---: | --- |
| `cb60fd71d1d7` | `rogers_ramanujan_q3_normalized` | 7 | `q^12` | RR main ladder + cubic extra numerator + denominator perturbation |
| `bef31ddceea8` | `rogers_ramanujan_q4_normalized` | 7 | `q^11` | single-numerator step perturbation |
| `e42a0d5f2679` | `rogers_ramanujan_q4_normalized` | 6 | `q^10` | single-numerator step perturbation |
| `9dbafd59364c` | `rogers_ramanujan_q4_normalized` | 5 | `q^9` | single-numerator step perturbation |
| `e2cc74240b6f` | `rogers_ramanujan_q3_normalized` | 5 | `q^9` | single-numerator step perturbation |

## Low-Order Fits

- `cb60fd71d1d7`
  - `candidate / RR(q^3) = 1 - q^12 + 2 q^15 - q^18`
- `bef31ddceea8`
  - `candidate / RR(q^4) = 1 + q^11 - q^12 - q^15`
- `e42a0d5f2679`
  - `candidate / RR(q^4) = 1 + q^10 - q^12 - q^14`
- `9dbafd59364c`
  - `candidate / RR(q^4) = 1 + q^9 - q^12 - q^13`
- `e2cc74240b6f`
  - `candidate / RR(q^3) = 1 - q^9 + q^10 + q^12`

## Ratio-Series Congruence Signals (q-series up to `q^59`)

- `cb60fd71d1d7 / RR(q^3)`:
  - only exponents `q^(3k)` appear (a pure function of `t = q^3`)
  - sample truncation:
    - `1 - q^12 + 2 q^15 - q^18 - q^21 + 4 q^24 - 9 q^27 + ... - 352 q^57`
- `e42a0d5f2679 / RR(q^4)`:
  - only even exponents appear (a function of `q^2`)
  - sample truncation:
    - `1 + q^10 - q^12 - q^14 + 2 q^20 - q^28 + q^30 - 2 q^32 - ... + 2 q^58`
- `bef31ddceea8`, `9dbafd59364c`, `e2cc74240b6f`:
  - exponents appear in all residue classes modulo 3 and modulo 4 in the same truncation window

## Structural Reading

- `cb60fd71d1d7` is the only candidate in this group that is not a plain
  single-ladder step shift.
- `bef31ddceea8`, `e42a0d5f2679`, and `9dbafd59364c` form an obvious `RR(q^4)`
  ladder:
  - same numerator start at `q^4`
  - step changed from `4` to `3`, `2`, `1`
  - no extra numerator
  - no denominator perturbation
- `e2cc74240b6f` is the analogous `RR(q^3)` single-ladder perturbation:
  - same numerator start at `q^3`
  - step changed from `3` to `4`
  - no extra numerator
  - no denominator perturbation

## Why `cb60fd71d1d7` Still Leads

- It ties for the top shared-digit tier in the public `review` set.
- It has the deepest first divergence among the classical `q^3 / q^4` review
  candidates now on the board.
- Its nearest benchmark is a classical public family, not an internal fixture.
- Its structure is qualitatively richer than the other `q^3 / q^4` reviews:
  it mixes the `RR(q^3)` ladder, a cubic-style extra numerator, and a
  denominator perturbation.
- It has now survived the strongest transform audit performed so far:
  direct specialization, arbitrary `n`-dependent equivalence, and the relevant
  `cor2cf` contraction branches were all checked and not matched.

## Why the Others Are Weaker Hero Cases

- The `q^4` trio looks like a clean monotone ladder of step perturbations.
  That makes them interesting, but also makes them look more like catalog
  extensions than like a standout structurally hybrid formula.
- `e2cc74240b6f` is closer in neighborhood to `cb60fd71d1d7`, but it agrees
  with `RR(q^3)` for fewer visible orders and has a much simpler single-ladder
  structure.
- None of the non-`cb60...` candidates currently carries the combination of:
  top-tier digit match, classical benchmark proximity, deeper symbolic agreement,
  and a surviving nontrivial transform audit.

## Fixture-Adjacent Note

- `1125ffe48b3b` remains interesting numerically.
- It is weaker as a public hero case because its nearest benchmark is the
  internal `shifted_rr_fixture`, not a classical public family.
- Its current low-order fit is

```text
candidate / shifted_rr_fixture = 1 - q^11 - q^12 + q^13
```

so it is better treated as a secondary audit target rather than the headline lead.

## Current Ranking

1. `cb60fd71d1d7`
2. `bef31ddceea8`
3. `e2cc74240b6f`
4. `e42a0d5f2679`
5. `9dbafd59364c`
6. `1125ffe48b3b` as a separate fixture-adjacent side lead
