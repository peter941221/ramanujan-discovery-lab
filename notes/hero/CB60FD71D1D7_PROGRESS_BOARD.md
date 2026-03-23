# `cb60fd71d1d7` Progress Board

Status date: `2026-03-23`

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
  exact and bounded exclusion layers now cover page-43 neighborhoods, exact RR/cubic subsequence lanes, bounded Bauer-Muir neighborhoods, exact `GG` quotient-coordinate witnesses on sampled tail objects, the first weighted `GG` correction ladder, and the direct plus nested `P/B` normalized Weber bridge lane.

### 2. Lean Gate

```text
Lean status
├─ Rational-equivalence layer: present
├─ Exact exclusion waypoint: present
├─ GG quotient-coordinate waypoint: present
├─ Weber companion waypoint: present
├─ GG weighted-correction waypoint: present
├─ Weber residual-bridge waypoint: present
├─ Tail-operator waypoint: present
└─ Final positive identity theorem: absent
```

- Status: `in progress`
- Main hub:
  `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current phase board:
  `notes/hero/CB60FD71D1D7_PHASE5_SPRINT_BOARD.md`
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
  the sampled `U_t2`, `U_t3`, `U_t4` ladder and their gap-normalized residuals still give `0` hits in the checked one-core source-family eta boxes, direct eta boxes, direct modular-unit / eta boxes, Morton periodic-point templates, the explicit squared Morton coordinate `X_mt = F^2`, the first Weber-Schlafli coordinate template, the direct Weber companion templates, and the current `GG/Weber` modular-equation templates.
  The new squared-coordinate lane is already exact-blocked by one shared
  constant-term obstruction:
  `X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 -> (t^0, 4)`.
  Its literature-backed linear-fractional follow-up
  `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)` is now also checked and is
  already exact-blocked by the same kind of uniform constant-term miss:
  `T_mt^2 - (T_mt,2^2 - 4*T_mt,2 + 1)*T_mt + T_mt,2^2 -> (t^0, 8)`.

### A2. Weber companion obstruction

- Obstruction note:
  `notes/hero/CB60FD71D1D7_WEBER_COMPANION_OBSTRUCTION_NOTE.md`
- Lean waypoint:
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`
- Current verdict:
  the direct companion coordinate `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`
  also gives `0` hits across all `12` sampled tail-family objects, and the
  misses compress to one universal constant-term class:
  `B_ws^2 - B_ws,2 - 4 -> (t^0, -2)`,
  `B_ws,2^4 - P_ws^8 - 16*P_ws^4 -> (t^0, 16)`.
  The Lean shell now also exposes this as explicit sample-indexed
  first-failure theorem families rather than only as a packaged witness table.

### A3. Weber class-invariant compression

- Main note:
  `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Lean waypoint:
  `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean`
- Current coordinates:
  `g12_ws = 4*t*(Z_g - 1/Z_g)`,
  `p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`,
  where
  `Z_g = ((1 - t*F^2)^2) / (4*t*F^2)`.
- Current verdict:
  on the true normalized `GG` source these compress exactly to the named
  Ramanujan-Weber objects
  `(t^2; t^4)_inf^12` and `(-t^2; t^4)_inf^12`, so the lane is now a real
  positive-recognition lane rather than only a negative obstruction lane.
  On the hero tail-family samples, both coordinates still give `0` hits after
  the direct template check, the first eta / modular-unit correction check, and
  the first narrower plus-Pochhammer correction check.
  Berndt--Chan--Zhang identifies Ramanujan-Weber `G_n` / `g_n` with classical
  Weber `f` / `f1`, and Yui--Zagier supplies the classical Weber
  `f`, `f1`, `f2` trio, so the same
  `g12_ws / p12_ws / G_f2_ws`
  shell is now source-backed as a named classical Weber trio in project
  normalization rather than as an anonymous product gadget.
  The next focused bridge now keeps `G_g12_ws` as the primary eta-side
  residual, ties `G_p12_ws` to it through the exact coordinate identity
  `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`,
  now also records the classical Weber `f2` tri-product coordinate
  `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2) = G_g12_ws*G_p12_ws`,
  and on the hero sample that coordinate first differs from `1` at `t^1` with
  coefficient `-4`; its normalized follow-up
  `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)` then first differs from `1` at `t^1`
  with coefficient `1/2`, while both layers still give `0` hits in the same
  first self-polynomial, self-fractional-linear, self-quotient finite-product,
  eta / modular-unit, and plus-Pochhammer boxes.
  compresses that bridge further through
  `X_g_ws = 16*t^2 / g12_ws^2` and
  `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0`
  with `Q_gp_ws = p12_ws / g12_ws`.
  The same lane now also keeps the template-normalized coordinate
  `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`,
  which is exact `1` on true normalized `GG` and on the hero sample first
  differs from `1` at `t^1` with coefficient `4`; its first normalized
  follow-up `H_X_ws = (G_X_ws - 1) / (4*t^1)` still gives `0` hits in the same
  first theorem-shaped closure boxes.
  The same template-normalized anchor now also carries its own canonical Weber
  `j`-side lift
  `J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)`, and on the hero sample that
  coordinate first differs from `1` at `t^1` with coefficient `56/17`; its
  normalized follow-up `H_J_X_ws = (J_X_ws - 1) / (56/17*t^1)` then first
  differs from `1` at `t^1` with coefficient `1083/238`, while both layers
  still give `0` hits in the same first theorem-shaped closure boxes.
  The two canonical Weber `j`-side coordinates are now also compared directly
  through a lightweight lift-bridge:
  `D_JX_ws = J_f2_ws - J_X_ws` and `Q_JX_ws = J_f2_ws / J_X_ws`.
  On the hero series, `Q_JX_ws - 1` first fails at `t^1` with coefficient
  `128/85`, and the bounded degree-`<= 3` polynomial bridge plus the one-coordinate
  fractional-linear bridge both still report `0` hits.
  The same lift-bridge quotient is now also checked in a focused source-faithful
  lane (on `Q_JX_ws` and its normalized follow-up `K_JX_ws`):
  the eta and modular-unit boxes stay flat (`0` hits), and the focused named-`GG`
  modular-equation lane also stays flat (`0` hits) in the checked
  direct/quotient/mixed prefix and exact template boxes.
  The same anchor now also carries an alternate signed `j`-lift
  `J_X15_ws = (16*G_X_ws - 1)^3 / (3375*G_X_ws^2)`:
  on the hero sample it first differs from `1` at `t^1` with coefficient `24/5`.
  The corresponding alternate lift-bridge
  `D_JX15_ws = J_f2_ws - J_X15_ws` and `Q_JX15_ws = J_f2_ws / J_X15_ws`
  is materially tighter: `Q_JX15_ws - 1` first fails only at `t^3` with
  coefficient `-576/5`, but the focused eta / modular-unit / named-`GG` boxes on
  `Q_JX15_ws` and its normalized follow-up `K_JX15_ws` still stay flat.
  The first direct bridge between the two hero-side normalized follow-ups is
  now also explicit:
  `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws`, and both first
  fail at `t^2` with coefficient `-24`; even after checking degree-`<= 3`
  polynomial bridges and a one-coordinate fractional-linear bridge, that
  comparison lane still gives `0` hits.
  The quotient side of that comparison bridge now also has its own normalized
  follow-up
  `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, and on the hero sample it first
  differs from `1` at `t^1` with coefficient `2` while still giving `0` hits
  in the first theorem-shaped closure boxes.
  That quotient-follow-up lane is now also compared back to the
  template-normalized branch:
  `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws` both first
  fail at `t^1` with coefficient `-5/2`, and their stripped quotient follow-up
  `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)` still gives `0` hits in the same
  first theorem-shaped closure boxes.
  The Lean shell now also packages that whole local hand-off as the grouped
  theorem-shaped shell
  `quotientFollowupReturnBridgeProp` /
  `quotientFollowupReturnBridge_true`, so this is no longer just a loose list
  of witnesses.
  A later source-faithful literature closure pass then pushed `Q_XK_ws` and
  `L_XK_ws` back through the named Chan--Huang `GG` modular-equation basis and
  still got `0` hits in the checked direct, quotient, and mixed
  quotient-coordinate boxes.
  Current reading:
  these are still best treated as derived return-bridge objects, not as
  recovered named `GG` coordinates.
  After that, the bridge still probes the residual quotient
  `R_gp_ws = G_p12_ws / G_g12_ws`; that quotient still gives `0` hits in the
  first self-polynomial, self-fractional-linear, self-quotient finite-product,
  eta / modular-unit, and plus-Pochhammer boxes, and first differs from `1`
  at `t^3` with coefficient `96`.
  The next normalized follow-up
  `H_gp_ws = (R_gp_ws - 1) / (96*t^3)` is now also tracked explicitly; it still
  gives `0` hits in the same first theorem-shaped closure boxes.

### A4. Weber-Schlafli normalized P/B bridge

- Bridge note:
  `notes/hero/CB60FD71D1D7_PB_BRIDGE_DECISION_NOTE.md`
- Current verdict:
  `N_P_ws` and `N_B_ws` both show `3 / 18` self-polynomial hits and `0` hits in the first surrounding fractional-linear / eta / product boxes.
  The direct bridge layer `D_PB_ws`, `Q_PB_ws`, `K_PB_ws` first fails at `t^3` with coefficient `3/4`; its bridge polynomial box has exactly `1` hit, only at degree `3`.
  The nested bridge layer `D_PK_ws`, `Q_PK_ws`, `L_PK_ws` first fails at `t^3` with coefficient `35/24`; its bridge polynomial box also has exactly `1` hit, again only at degree `3`.
  `Q_PB_ws`, `K_PB_ws`, `Q_PK_ws`, and `L_PK_ws` all keep the same profile: `3 / 18` self-polynomial hits and `0` hits in the surrounding fractional-linear, self-quotient finite-product, eta, modular-unit / eta, and plus-Pochhammer boxes.
  The first source-faithful one-coordinate orbit pass on that same ladder is now also complete:
  on true `GG`, the direct `Q_PB_ref_ws` and `Q_PK_ref_ws` ladders each show `3 / 8` direct polynomial-prefix hits,
  but on the hero candidate the corresponding focused direct / quotient / mixed boxes stay flat for `Q_PB_ws`, `K_PB_ws`, `Q_PK_ws`, and `L_PK_ws`.
  Reading: this is still the cleanest local Weber seam, but the first same-orbit compression pass is now exhausted once and still does not give a clean literature-facing coordinate.

### A5. Canonical Weber `j`-side cross-check

- Decision note:
  `notes/hero/CB60FD71D1D7_WEBER_J_SIDE_DECISION_NOTE.md`
- Current verdict:
  `J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)` is now the cleanest
  literature-facing Weber coordinate in the current city.
  The first focused `J_f2_ws ↔ Q_PB_ws` bridge is now also complete:
  on the direct hero series, the direct bridge
  `D_JPB_ws`, `Q_JPB_ws`, `K_JPB_ws`
  and the nested bridge
  `D_JKPB_ws`, `Q_JKPB_ws`, `L_JKPB_ws`
  both stay flat in the checked bridge polynomial / fractional-linear boxes.
  The follow-up objects
  `Q_JPB_ws`, `K_JPB_ws`, `Q_JKPB_ws`, and `L_JKPB_ws`
  also stay flat in the focused named-`GG` direct / quotient / mixed prefix
  boxes and exact modular-equation template boxes.
  Reading: `J_f2_ws` remains the best literature-facing Weber lane, but the
  easiest "maybe it is just the P/B seam in disguise" explanation has now been
  tested once and stayed flat.
  The new `Q_JX15_ws` pivot rails now also compare that tighter lift-bridge quotient
  back to the strongest local seams:
  `Q_JX15PB_ws` and `Q_JX15JPB_ws` first fail at `t^1` with coefficient `-24/5`,
  while `Q_JX15XKJ_ws` first fails at `t^1` with coefficient `-97/30`.
  All three pivots stay flat in the same focused eta / modular-unit / named-`GG` boxes.

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
  `Build elapsed seconds before final render: 1809.79`
- Full formalize refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md --lean-out proofs/Proofs/Generated/Cb60fd71d1d7.lean`
- Latest formalize runtime:
  about `1092s`
- Quick syntax validation:
  `python -m py_compile src/ramanujan_discovery/identification.py tests/test_identification.py`
- Targeted Weber regression slice:
  `$env:PYTHONPATH='src'; pytest -q tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_matches_true_gg tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_reports_hero_ratio_gap tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_skips_focused_named_gg_lane_in_smoke_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_matches_true_gg_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_reports_hero_gap_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_skips_named_gg_lane_in_smoke_profile tests/test_identification.py::test_cli_tail_note_writes_tail_family_note`
- Full tail-family refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Full tail-family note timing:
  about `9706s`
- Full tail-operator refresh:
  `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-operator-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 36 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_OPERATOR_NOTE.md`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseGGQuotientCoordinateObstruction.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseTailOperatorWaypoint.lean`
- Targeted Lean validation:
  `Set-Location proofs; lake env lean Proofs/HeroCaseFinalIdentity.lean`

## Priority Queue

1. Keep the next positive-recognition attempt in the modular-function / eta-recognition trunk rather than widening anonymous boxes again.
2. Keep the named `GG/Weber` coordinate search deeper than `Q_3`, `Q_4`, `W_34`, `G_W34`, `G2_W34`, `Q_XK_ws`, `L_XK_ws`, and the new classical Weber `f2` tri-product lane `G_f2_ws` / `H_f2_ws`; the direct plus nested `P/B` normalized bridge lane is still the best local seam inside that orbit.
3. Treat the new Weber class-invariant pair `g12_ws` / `p12_ws` as the current best constructive hand-off, but note that the first focused one-coordinate orbit pass on `Q_PB_ws` / `Q_PK_ws` is now complete and flat on the hero candidate.
4. Prefer the new bridge-normalized reading of that hand-off: use `G_g12_ws` as the primary residual, treat `G_p12_ws` as the algebraically constrained companion, keep `X_g_ws = 16*t^2 / g12_ws^2`, `G_X_ws = 1 / G_g12_ws^2` together with its canonical `j`-side lift `J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)`, and the `P_ws -> B_ws -> Q_PB_ws -> Q_PK_ws` seam in the foreground; the canonical `J_f2_ws` lane is now the cleanest literature-facing comparison rail, but the first `J_f2_ws ↔ Q_PB_ws` cross-check is also complete and flat, so the next move should be a different named Weber coordinate or only then a different comparison bridge.
5. Preserve the operator lane as theorem scaffolding, but do not let it outrank source recognition until a stronger source object appears.
