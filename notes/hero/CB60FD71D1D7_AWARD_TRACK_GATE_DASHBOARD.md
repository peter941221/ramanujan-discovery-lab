# `cb60fd71d1d7` Award-Track Gate Dashboard

Status date: `2026-03-20`

## Purpose

This dashboard is the shortest reusable read for the current award-track
status.

Use it when we need to answer:

```text
Which gate is blocked?
What evidence already exists?
What is the next non-noisy move?
```

## Gate Summary

| Gate | Status | What is already true | Main blocker | Next best move |
| --- | --- | --- | --- | --- |
| Theorem gate | `yellow` | Exact stationary tail law exists; many nearby named or semi-named lanes are excluded | No positive final source object or uniqueness-grade characterization | Stay in source-faithful modular / eta and named `GG/Weber` coordinate lanes |
| Lean gate | `yellow` | `HeroCaseFinalIdentity` now packages the current exclusion waypoint, including the new Weber companion shell | No final positive identity theorem to formalize | Upgrade from waypoint shells to a positive source theorem once recognition lands |
| Literature gate | `yellow` | Core RR / cubic / GG / S spine is logged, plus recent modularity and periodic-point sources | Not enough closure for a prize-grade novelty claim | Close only the winning source orbit after positive recognition, not the whole universe first |

## 1. Theorem Gate

```text
Current evidence
├─ Positive structure
│  └─ `T(x) = 1 + x + x*(t + x)/T(t*x)`
├─ Strong negative structure
│  ├─ page-43 lanes
│  ├─ subsequence lanes
│  ├─ Bauer-Muir neighborhoods
│  ├─ sampled `GG` quotient-coordinate exact lanes
│  └─ first weighted `GG` correction ladder
└─ Missing deliverable
   └─ accepted closed form or equivalent uniqueness theorem
```

- Current status: `blocked on positive identification`
- Best current interpretation:
  we have a good map of nearby cliffs, but not the summit marker yet.
- Immediate practical consequence:
  another generic box widening step is lower value than a source-faithful
  modular-function recognition step.

## 2. Lean Gate

```text
Lean ladder
├─ Rational equivalence: yes
├─ Exact exclusion waypoint: yes
├─ `GG Q_3/Q_4` obstruction waypoint: yes
├─ Weber companion obstruction waypoint: yes
├─ Weighted `GG` correction waypoint: yes
├─ Tail-operator waypoint: yes
└─ Final positive identity theorem: no
```

- Current status: `waiting for the right source object`
- Main file:
  `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current marker:
  `AwardTrackStatus: exclusion_waypoint`
- Practical interpretation:
  the Lean side is no longer missing scaffolding; it is missing the final
  mathematical object worth scaffolding around.

## 3. Literature Gate

```text
Literature closure
├─ Classical spine: usable
├─ Modern RR modularity / modular-equation papers: usable
├─ GG modular-equation and arithmetic papers: usable
├─ Morton periodic-point lane: usable
└─ Prize-grade closure around the final winning object: not yet
```

- Current status: `not ready for novelty claim`
- Public-language consequence:
  keep saying `unexplained candidate`, not `new identity`.
- Practical interpretation:
  literature closure should now trail the likely winning orbit, not lead with a
  giant survey of everything Ramanujan-adjacent.

## Named Working Lanes

### Lane A. Tail-family-first recognition

- Status: `active`
- Current artifact:
  `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Current verdict:
  `12` sampled tail-family objects still give `0` hits in the checked direct
  eta, modular-unit / eta, Morton periodic-point, and current `GG/Weber`
  modular-equation boxes.

### Lane B. Exact `GG` quotient coordinates

- Status: `active as obstruction lane`
- Current artifact:
  `notes/hero/CB60FD71D1D7_GG_Q34_OBSTRUCTION_NOTE.md`
- Current verdict:
  all sampled tail-family objects miss `Q_3 = GG(t^3)/GG(t)` and
  `Q_4 = GG(t^4)/GG(t)`, with a compressed shared leading `3:2` failure shape.

### Lane C. Weighted `GG` correction ladder

- Status: `active as deeper diagnostic`
- Current objects:
  `W_34`, `F / W_34`, `G_W34`, `G2_W34`
- Current verdict:
  the first weighted coordinate and the next two visible corrections still do
  not collapse into the checked small eta / modular-unit / one-core correction
  boxes.

### Lane E. Weber-Schlafli coordinate hand-off

- Status: `implemented P/B first pass; still 0-hit`
- Lean waypoint:
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`
- Literature trigger:
  Akkarapakam--Morton rewrites the relevant Ramanujan continued-fraction object
  through Weber-Schlafli coordinates, with `2*p(8τ) = 1/v(τ) - v(τ)`.
- First coded coordinate:
  `P_ws = (1/F - F) / 2`
- First coded exact template:
  `P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0`
- Companion coded coordinate:
  `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`
- Companion exact templates:
  `B_ws^2 - B_ws,2 - 4 = 0`, `B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0`
- Current practical meaning:
  after the raw `Q_3`, `Q_4` and weighted `W_34` ladders all miss, the first
  exact source-faithful coordinate change and its direct companion closure
  have now both been tested and still give `0` hits across the sampled
  tail-family objects.
- Current verdict:
  the `P_ws` lane closes as a repeated low-order obstruction layer, while the
  `B_ws` companion lane closes even harder as a uniform constant-term
  obstruction recorded in
  `notes/hero/CB60FD71D1D7_WEBER_COMPANION_OBSTRUCTION_NOTE.md`.
  That companion obstruction is now also packaged in Lean, so the next move
  inside the same orbit should therefore be another named Weber coordinate,
  not anonymous box growth.

### Lane D. Operator-first scaffold

- Status: `keep, but do not over-promote`
- Current artifact:
  `notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Current verdict:
  the low-degree affine q-difference / Mahler box remains a `0`-hit diagnostic
  lane, but it gives a plausible proof endgame once the right source object is
  known.

## Latest Verified Runs

- Full identify refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery identify --in results/verified.jsonl --candidate-id cb60fd71d1d7 --benchmark-powers 2,3,4,5,6,12,20 --out notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`
- Full identify note timing:
  `Build elapsed seconds before final render: 2091.88`
- Full tail-family refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Full tail-family note timing:
  about `344s`
- Full tail-operator refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-operator-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 36 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Targeted Weber Python validation:
  `$env:PYTHONPATH='src'; pytest tests/test_identification.py -q -k "weber_schlafli or tail_note or morton_periodic_point"`

## Priority Order

1. Source-faithful modular-function / eta recognition.
2. Deeper named `GG/Weber` coordinates beyond `Q_3`, `Q_4`, `W_34`, `G_W34`,
   and `G2_W34`; the direct `P = p(8τ)` and `B = b(4τ)` passes are now both
   done, so the next source-faithful move is another named Weber coordinate or
   a theorem-grade package of the new companion obstruction.
3. Operator / factorization conversion only after a stronger source object is
   in hand.
