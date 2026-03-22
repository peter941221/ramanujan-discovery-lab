# `cb60fd71d1d7` Award-Track Gate Dashboard

Status date: `2026-03-23`

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
  The same scan now also includes the Weber class-invariant quotient/template
  anchor `G_X_ws` plus its canonical `j`-side lift
  `J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)`.

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

- Status: `active: direct plus nested P/B bridge implemented; first one-coordinate orbit pass also done; still 0-hit`
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
  exact source-faithful coordinate change, its direct companion closure, the
  direct `P/B` normalized bridge, and the first nested quotient layer have now
  all been tested and still give `0` hits across the sampled tail-family
  objects.
- Current verdict:
  the `P_ws` lane closes as a repeated low-order obstruction layer, while the
  `B_ws` companion lane closes even harder as a uniform constant-term
  obstruction recorded in
  `notes/hero/CB60FD71D1D7_WEBER_COMPANION_OBSTRUCTION_NOTE.md`.
  The direct `P/B` bridge then first fails at `t^3` with coefficient `3/4`,
  and the nested bridge first fails at `t^3` with coefficient `35/24`; both
  bridge polynomial boxes light up only once, at degree `3`, while
  `Q_PB_ws`, `K_PB_ws`, `Q_PK_ws`, and `L_PK_ws` all keep the same
  `3 / 18` self-polynomial and `0`-hit side-box profile.
  The first source-faithful one-coordinate orbit pass on that seam is now also
  complete: on true `GG`, the direct `Q_PB_ref_ws` and `Q_PK_ref_ws` ladders do
  register bounded polynomial recognition, but on the hero candidate the same
  focused direct / quotient / mixed boxes remain flat.
  The canonical classical-Weber `j`-side lane is now also explicit through
  `J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)`, and the first focused
  `J_f2_ws ↔ Q_PB_ws` cross-check is also complete:
  on the direct hero series, both the direct and nested bridge boxes stay
  flat, and the follow-up objects `Q_JPB_ws`, `K_JPB_ws`, `Q_JKPB_ws`, and
  `L_JKPB_ws` also keep `0` hits in the focused named-`GG` direct / quotient /
  mixed prefix boxes.
  That means the current seam is structured but still bulky, and its first
  same-orbit compression pass is already exhausted once.
  It also means the cleanest named `j`-side coordinate is not already recovered
  by the first focused `J/PB` bridge box.
  The next move inside the same city should therefore be a different named
  Weber coordinate or only then a wider comparison bridge, not more blind
  growth of the same tiny box.

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
  `Build elapsed seconds before final render: 1809.79`
- Full tail-family refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Full tail-family note timing:
  about `9706s`
- Full tail-operator refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-operator-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 36 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Targeted Weber Python validation:
  `$env:PYTHONPATH='src'; pytest -q tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_matches_true_gg tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_reports_hero_ratio_gap tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_skips_focused_named_gg_lane_in_smoke_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_matches_true_gg_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_reports_hero_gap_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_skips_named_gg_lane_in_smoke_profile tests/test_identification.py::test_cli_tail_note_writes_tail_family_note`

## Priority Order

1. Source-faithful modular-function / eta recognition.
2. Deeper named `GG/Weber` coordinates beyond `Q_3`, `Q_4`, `W_34`, `G_W34`,
   and `G2_W34`; the direct `P = p(8τ)` and `B = b(4τ)` passes are now both
   done, the first focused one-coordinate orbit pass on `Q_PB_ws` / `Q_PK_ws`
   is now also done and flat, and the first focused `J_f2_ws ↔ Q_PB_ws`
   cross-check is also done and flat, so the next source-faithful move is a
   different named Weber coordinate or only then a wider comparison bridge.
3. Operator / factorization conversion only after a stronger source object is
   in hand.
