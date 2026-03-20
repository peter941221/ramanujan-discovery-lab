# `cb60fd71d1d7` Phase 4 Sprint Board

Status date: `2026-03-20`

## Phase Definition

Phase 4 here means:

```text
take the new Weber companion obstruction
+ package it in Lean
+ splice it into the award-track exact waypoint
+ close the P/B branch as a reusable theorem-shaped handoff
+ open the next named Weber-coordinate phase
```

## Checklist

### 1. Package the Weber companion obstruction in Lean

- Status: `done`
- New Lean shell:
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`
- What it records:
  - the companion coordinate
    `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`
  - the two first exact companion templates
  - the `12` sampled tail-family witnesses
  - the universal constant-term obstruction class
    `(-2, 16)`

### 2. Splice the new shell into the award-track exact waypoint

- Status: `done`
- Updated hub:
  `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current effect:
  the exclusion waypoint now explicitly includes the direct Weber companion
  obstruction layer alongside the page-43, RR/cubic, `cor2cf`, and `GG`
  quotient-coordinate shells.

### 3. Close the P/B branch as a reusable handoff

- Status: `done`
- Closed local branch:

```text
P branch: repeated low-order obstruction classes
B branch: one universal constant-term obstruction class
```

- Stable notes:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
  - `notes/hero/CB60FD71D1D7_WEBER_COMPANION_OBSTRUCTION_NOTE.md`

### 4. Open the next named Weber-coordinate phase

- Status: `done`
- Rejected as the next move:
  - reopening anonymous scan growth
  - re-running the same `P`/`B` box without theorem upgrades
- New phase handoff:

```text
Next named Weber phase
├─ keep the same GG/Weber orbit
├─ treat P/B as exact blocked waypoints
├─ either add the next named Weber coordinate
└─ or upgrade the current companion shell into a stronger exact theorem
```

## Why This Counts As A Phase Advance

Analogy:

- the previous step said “the P/B doors are both locked”
- this phase adds “and now the locks are mounted onto the main proof hallway”

So the project is no longer just collecting negative evidence.

It is converting that evidence into reusable theorem-grade scaffolding.

## Phase-4 Exit Condition

Phase 4 is complete for now because:

1. the direct Weber companion obstruction is now a Lean module
2. `HeroCaseFinalIdentity` imports that module into the current exact waypoint
3. the next phase is no longer “try B”, but “move past P/B to the next named
   Weber coordinate or theorem upgrade”
