# `cb60fd71d1d7` Phase 1 Sprint Board

Status date: `2026-03-20`

## Phase Definition

Phase 1 here means the document-defined deep-research sprint:

```text
freeze the current local conclusion
+ verify the missing primary sources
+ update the literature artifacts
+ distill one sharper next implementation target
```

## Checklist

### 1. Freeze the current local tail conclusion

- Status: `done`
- Frozen local conclusion:
  first-gap and second-gap tail residuals still do not look like one nearby
  source core times a small eta tail.
- Stable local artifacts:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
  - `notes/hero/CB60FD71D1D7_GG_Q34_OBSTRUCTION_NOTE.md`
  - `notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`

### 2. Verify the missing primary-source spine around nonlinear `GG` relations

- Status: `done`
- Verified spine for this sprint:
  - Chan--Huang, 1997:
    exact `GG` modular-equation and quotient-coordinate source spine around the
    `q^3` and `q^4` lanes
  - Cho--Koo--Park, 2009:
    modular-curve / modular-unit extension of the `GG` story to all odd primes
    via affine models of modular curves
  - Akkarapakam--Morton, 2024:
    Weber-Schlafli / periodic-point reformulation that exposes deeper
    coordinates than the raw `Q_3`, `Q_4` pair
  - Yui--Zagier, 1997:
    classical Weber modular-function singular-value anchor for the
    modular-curve side

### 3. Update bibliography and literature artifacts

- Status: `done`
- Updated artifacts:
  - `notes/hero/CB60FD71D1D7_BIBLIOGRAPHY_MATRIX.md`
  - `notes/hero/CB60FD71D1D7_LITERATURE_LOG.md`
  - `notes/hero/CB60FD71D1D7_AWARD_TRACK_EXECUTION_PLAN.md`
  - `notes/hero/CB60FD71D1D7_AWARD_TRACK_GATE_DASHBOARD.md`

### 4. Distill the next implementation target

- Status: `done`
- Rejected as the immediate next coding move:
  broader anonymous prefix growth
- Rejected as the immediate positive-recognition target:
  the raw `Q_3`, `Q_4` pair by itself
- New exact implementation target:

```text
Weber-Schlafli coordinate lane
├─ start from the `GG` variable v(τ)
├─ pass to the deeper coordinate
│  `P = p(8τ)` with `2P = 1/v - v`
├─ keep the companion Weber coordinate
│  `B = b(4τ)` when needed
└─ test the first exact low-degree source equations on `P`, `P(2τ)`, and `B`
```

## Why This Is The Phase-1 Outcome

Analogy:

- the old plan said “search the right neighborhood more carefully”
- Phase 1 now says “the next door to try is not the front door `Q_3/Q_4`,
  but the side entrance labeled `p(8τ)`”

So Phase 1 does **not** prove the hero identity.

It does something narrower and still valuable:

- it upgrades the nearby literature spine from “general direction” to
  “one concrete next coordinate family”

## Phase-1 Exit Condition

Phase 1 is considered complete for now because it ends with one sharper local
target:

1. keep the `GG/Weber` neighborhood
2. stop broadening anonymous boxes
3. implement the Weber-Schlafli coordinate lane before any larger orbit growth

## Immediate Phase-2 Hand-Off

The first Phase-2 coding pass should try a small exact source-faithful box
built from:

- `v(τ)` or the normalized `GG` coordinate already used locally
- `P = p(8τ)` via `2P = 1/v - v`
- the first Morton polynomial relation in `P` and `P(2τ)`
- optional companion coordinate `B = b(4τ)` only if the one-coordinate `P`
  lane still misses
