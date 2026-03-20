# `cb60fd71d1d7` Progress Board

Status date: `2026-03-20`

## Snapshot

```text
Award-track state
├─ Hero object: `cb60fd71d1d7`
├─ Current status: `unexplained candidate`
├─ Proof scaffold status: `exclusion_waypoint`
├─ Main structural asset: exact stationary tail law
└─ Main open need: positive final identity with theorem-grade closure
```

## Gate Board

### 1. Theorem Gate

```text
State
├─ Exact tail-family law extracted
├─ Many nearby source families excluded in bounded or exact lanes
├─ No accepted closed form identified yet
└─ Final positive source object still open
```

- Status: `in progress`
- Strongest current positive structure:
  `T(x) = 1 + x + x*(t + x)/T(t*x)`
- Strongest current negative structure:
  exact and bounded exclusion layers now cover page-43 neighborhoods, exact RR/cubic subsequence lanes, bounded Bauer-Muir neighborhoods, exact `GG` quotient-coordinate witnesses on sampled tail objects, and the first weighted `GG` correction ladder.

### 2. Lean Gate

```text
Lean status
├─ Rational-equivalence layer: present
├─ Exact exclusion waypoint: present
├─ GG quotient-coordinate waypoint: present
├─ GG weighted-correction waypoint: present
├─ Tail-operator waypoint: present
└─ Final positive identity theorem: absent
```

- Status: `in progress`
- Main hub:
  `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current marker:
  `AwardTrackStatus: exclusion_waypoint`

### 3. Literature Gate

```text
Literature status
├─ RR / cubic / GG / S source spine: active
├─ Recent RR modularity and modular-equation papers: logged
├─ Morton periodic-point lane: logged and tested locally
└─ Closure strong enough for prize-caliber novelty claim: not yet
```

- Status: `in progress`
- Current reading:
  the literature still supports staying in eta / modular-unit plus named modular-coordinate lanes, while keeping public language at `unexplained candidate`.

## Current Working Lanes

### A. Tail-family-first recognition

- Exact tail-family note:
  `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Current verdict:
  the sampled `U_t2`, `U_t3`, `U_t4` ladder and their gap-normalized residuals still give `0` hits in the checked one-core source-family eta boxes, direct eta boxes, direct modular-unit / eta boxes, Morton periodic-point templates, and the current `GG/Weber` modular-equation templates.

### B. `GG` quotient-coordinate obstruction

- Obstruction note:
  `notes/hero/CB60FD71D1D7_GG_Q34_OBSTRUCTION_NOTE.md`
- Lean waypoint:
  `proofs/Proofs/HeroCaseGGQuotientCoordinateObstruction.lean`
- Current verdict:
  all `12` sampled tail-family objects miss the exact Chan--Huang `Q_3` and `Q_4` lanes, and the misses compress to a shared leading `3:2` obstruction pattern.

### C. Weighted `GG` correction ladder

- Main coordinate:
  `W_34 = Q_3^3 / Q_4^2`
- Current corrections:
  `F / W_34`, `G_W34`, `G2_W34`
- Lean waypoint:
  `proofs/Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`
- Current verdict:
  the first weighted coordinate and its first two normalized corrections still do not collapse into the checked small eta-quotient, modular-unit / eta, or one-core `RR/GG` source-family correction boxes.

### D. Operator-first proof scaffold

- Operator note:
  `notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Lean waypoint:
  `proofs/Proofs/HeroCaseTailOperatorWaypoint.lean`
- Current verdict:
  the affine q-difference / Mahler box is still a `0`-hit diagnostic lane, but it now gives the project a concrete proof-oriented endgame scaffold.

## Latest Verified Commands

- Full identify refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --benchmark-powers 2,3,4,5,6,12,20 --out notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`
- Latest identify timing:
  `Build elapsed seconds before final render: 2091.88`
- Full tail-family refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Full tail-operator refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-operator-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 36 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Targeted Python validation:
  `$env:PYTHONPATH='src'; pytest tests/test_identification.py -q -k "test_cli_tail_operator_note_writes_tail_operator_note or test_cli_tail_note_writes_tail_family_note or test_detect_reduced_tail_transfer_equation_finds_stationary_hero_tail or test_build_reduced_tail_anchor_builds_stage_three_tail_and_normalization or test_scan_gg_modular_equation_box_reports_weighted_q3q4_coordinate"`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseGGQuotientCoordinateObstruction.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseTailOperatorWaypoint.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseFinalIdentity.lean`

## Priority Queue

1. Keep the next positive-recognition attempt in the modular-function / eta-recognition trunk rather than widening anonymous boxes again.
2. Push the named `GG/Weber` coordinate search deeper than `Q_3`, `Q_4`, `W_34`, `G_W34`, and `G2_W34`.
3. Preserve the operator lane as theorem scaffolding, but do not let it outrank source recognition until a stronger source object appears.
