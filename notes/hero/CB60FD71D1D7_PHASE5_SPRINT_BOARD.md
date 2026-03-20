# `cb60fd71d1d7` Phase 5 Sprint Board

Status date: `2026-03-20`

## Phase Definition

Phase 5 here means:

```text
upgrade the Weber companion Lean shell
+ expose explicit first-failure theorem families
+ keep the shell wired into HeroCaseFinalIdentity
+ sync the public-facing hero narrative
+ leave the next phase pointed at a new named Weber coordinate
```

## Checklist

### 1. Upgrade the companion shell from data packaging to theorem-family packaging

- Status: `done`
- Updated Lean shell:
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`
- New explicit theorem families:
  - `squareFirstFailure_true`
  - `quarticFirstFailure_true`
  - `pairedFirstFailure_true`

### 2. Preserve the exact waypoint integration

- Status: `done`
- Main hub:
  `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current effect:
  the award-track exact waypoint still imports the Weber companion shell, but
  the imported shell now carries stronger sample-indexed first-failure facts
  instead of only a witness-table summary.

### 3. Sync the public-facing hero narrative

- Status: `done`
- Updated public artifacts:
  - `notes/hero/CB60FD71D1D7_PUBLIC_SUMMARY.md`
  - `notes/hero/CB60FD71D1D7_CASE_STUDY.md`
- New public-facing message:
  the Weber `P/B` branch is still source-faithful and still negative, but it
  is now negative in a sharper, Lean-backed way.

### 4. Set the next post-Phase-5 handoff

- Status: `done`
- Next handoff:

```text
Post-Phase-5 direction
├─ keep the GG/Weber orbit
├─ treat the direct P/B branch as exact-blocked
├─ do not reopen anonymous scan growth
└─ move to the next named Weber coordinate
```

## Why This Counts As A Phase Advance

Analogy:

- Phase 4 mounted the `B` obstruction shell onto the proof hallway
- Phase 5 adds labeled circuit breakers, so each sample-indexed failure is now
  a named reusable theorem rather than just a packaged table

## Phase-5 Exit Condition

Phase 5 is complete for now because:

1. the Weber companion shell now exposes explicit first-failure theorem
   families
2. the public-facing narrative now reflects the same P/B obstruction story
3. the next phase is no longer "strengthen B again", but "move to the next
   named Weber coordinate"
