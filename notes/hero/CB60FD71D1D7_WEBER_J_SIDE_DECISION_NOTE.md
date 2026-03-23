# Weber J-Side Decision Note: `cb60fd71d1d7`

## Snapshot

- Date: `2026-03-23`
- Scope: compare the current one-coordinate Weber-facing lanes
  - canonical classical-Weber `j`-side lane
    `J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)`
  - template-normalized class-invariant lane
    `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`
  - anchor-derived canonical `j`-side lift of that template lane
    `J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)`
  - Weber-Schlafli normalized bridge lane
    `Q_PB_ws = N_B_ws / N_P_ws`

## Why This Note Exists

- The `P/B` bridge note already showed that `Q_PB_ws` is the strongest structural bridge seam.
- The new `J_f2_ws` lane adds a more literature-facing coordinate built directly from the existing classical Weber `f2` shell.
- This note is the short comparison board for deciding which lane should lead the next focused Weber pass.
- The first focused `J_f2_ws ↔ Q_PB_ws` cross-check is now also complete, so this note can record not only the candidate lanes but also the first bridge between the cleanest named lane and the strongest local seam.
- Analogy:
  - `G_X_ws` is the clean template-calibrated measuring stick
  - `Q_PB_ws` is the strongest bridge seam inside the local terrain
  - `J_f2_ws` is the most canonical city-map coordinate because it points toward `j`
  - `J_X_ws` is the anchor's own `j`-side lift (a city-map coordinate derived from the ruler)

## Three-Way Comparison

```text
Weber Coordinate Board
├─ J_f2_ws
│  ├─ Type: canonical classical-Weber j-side coordinate
│  ├─ True GG behavior: matches 1 through the checked truncation
│  ├─ Hero first failure: J_f2_ws - 1 at t^1 with coefficient 24/5
│  ├─ Hero normalized follow-up: H_J_f2_ws = (J_f2_ws - 1) / (24/5*t^1)
│  ├─ Follow-up first failure: t^1 with coefficient 409/90
│  └─ First small-box verdict: 0 hits in the checked self-polynomial, fractional-linear,
│     self-quotient, eta, modular-unit, and plus-Pochhammer boxes
├─ G_X_ws
│  ├─ Type: template-normalized positive-recognition coordinate
│  ├─ True GG behavior: exactly 1 in the normalized source lane
│  ├─ Hero first failure: G_X_ws - 1 at t^1 with coefficient 4
│  ├─ Hero normalized follow-up: H_X_ws = (G_X_ws - 1) / (4*t^1)
│  ├─ Follow-up first failure: t^1 with coefficient 9/2
│  ├─ Derived canonical j-lift:
│  │  ├─ J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)
│  │  ├─ Hero first failure: J_X_ws - 1 at t^1 with coefficient 56/17
│  │  ├─ Hero normalized follow-up: H_J_X_ws = (J_X_ws - 1) / (56/17*t^1)
│  │  ├─ Follow-up first failure: t^1 with coefficient 1083/238
│  │  └─ First small-box verdict: 0 hits in the same theorem-shaped closure boxes
│  └─ First small-box verdict: 0 hits in the same theorem-shaped closure boxes
└─ Q_PB_ws
   ├─ Type: Weber-Schlafli normalized bridge quotient
   ├─ True GG behavior: the focused source-faithful one-coordinate orbit is real there
   ├─ Hero first failure: Q_PB_ws - 1 at t^3 with coefficient 3/4
   ├─ Hero normalized follow-up: K_PB_ws = (Q_PB_ws - 1) / (3/4*t^3)
   ├─ Structural signal: Q_PB_ws, K_PB_ws, Q_PK_ws, L_PK_ws each keep 3 / 18 self-polynomial hits
   └─ First source-faithful orbit verdict on hero: direct 0 / 7, quotient 0 / 6,
      mixed 0 / 7, plus 0 multiplicative and 0 fractional-linear hits
```

## Focused `J/PB` Cross-Check

```text
J/PB Bridge Status
├─ Direct hero-series bridge
│  ├─ Objects: D_JPB_ws, Q_JPB_ws, K_JPB_ws
│  ├─ Reading: the quotient already breaks at the first t^1 step
│  ├─ Bridge-box verdict: polynomial 0 / 3, fractional-linear 0 / 1
│  └─ Named-GG verdict: Q_JPB_ws and K_JPB_ws stay flat in the focused
│     direct / quotient / mixed prefix boxes and exact modular-equation boxes
└─ Nested follow-up
   ├─ Objects: D_JKPB_ws, Q_JKPB_ws, L_JKPB_ws
   ├─ Reading: the nested quotient also breaks immediately, again at t^1
   └─ Verdict: the same bridge boxes and named-GG boxes stay flat
```

## Focused `J_f2_ws ↔ J_X_ws` Lift-Bridge Cross-Check

```text
J_f2_ws vs J_X_ws Lift-Bridge
├─ Objects
│  ├─ D_JX_ws = J_f2_ws - J_X_ws
│  ├─ Q_JX_ws = J_f2_ws / J_X_ws
│  └─ K_JX_ws = (Q_JX_ws - 1) / (128/85*t^1)   (hero gap-normalized follow-up)
├─ Hero reading
│  ├─ Q_JX_ws - 1 first fails at t^1 with coefficient 128/85
│  └─ Bridge-box verdict: polynomial 0 / 3, fractional-linear 0 / 1
└─ Named-GG modular-equation reading (focused order 16)
   ├─ Q_JX_ws: 0 exact template hits; 0 direct/quotient/mixed prefix hits
   └─ K_JX_ws: 0 exact template hits; 0 direct/quotient/mixed prefix hits
```

## Alternate Anchor `j`-Lift: `J_X15_ws`

```text
J_f2_ws vs J_X15_ws Alternate Lift-Bridge
├─ Alternate anchor coordinate
│  ├─ J_X15_ws = (16*G_X_ws - 1)^3 / (3375*G_X_ws^2)
│  └─ Hero first failure: J_X15_ws - 1 at t^1 with coefficient 24/5
├─ Bridge objects
│  ├─ D_JX15_ws = J_f2_ws - J_X15_ws
│  ├─ Q_JX15_ws = J_f2_ws / J_X15_ws
│  └─ K_JX15_ws = (Q_JX15_ws - 1) / (-576/5*t^3)   (hero gap-normalized follow-up)
├─ Hero reading
│  ├─ Q_JX15_ws - 1 first fails at t^3 with coefficient -576/5
│  └─ Bridge-box verdict: polynomial 0 / 3, fractional-linear 0 / 1
└─ Source-faithful verdict (focused order 16)
   ├─ eta boxes: 0 hits (Q_JX15_ws and K_JX15_ws)
   ├─ modular-unit / eta boxes: 0 hits (Q_JX15_ws and K_JX15_ws)
   └─ named-GG modular-equation boxes: 0 hits (Q_JX15_ws and K_JX15_ws)
```

## Pivot Rails Around `Q_JX15_ws`

```text
Q_JX15 Pivot Bridges (hero order 24)
├─ J_X15_ws ↔ Q_PB_ws
│  ├─ D_JX15PB_ws = Q_PB_ws - J_X15_ws
│  ├─ Q_JX15PB_ws = Q_PB_ws / J_X15_ws
│  └─ Hero first failure: Q_JX15PB_ws - 1 at t^1 with coefficient -24/5 (0 hits in eta / modular-unit / named-GG)
├─ Q_JX15_ws ↔ Q_JPB_ws
│  ├─ D_JX15JPB_ws = Q_JPB_ws - Q_JX15_ws
│  ├─ Q_JX15JPB_ws = Q_JPB_ws / Q_JX15_ws
│  └─ Hero first failure: Q_JX15JPB_ws - 1 at t^1 with coefficient -24/5 (0 hits in eta / modular-unit / named-GG)
└─ Q_JX15_ws ↔ Q_XKJ_ws
   ├─ D_JX15XKJ_ws = Q_XKJ_ws - Q_JX15_ws
   ├─ Q_JX15XKJ_ws = Q_XKJ_ws / Q_JX15_ws
   └─ Hero first failure: Q_JX15XKJ_ws - 1 at t^1 with coefficient -97/30 (0 hits in eta / modular-unit / named-GG)
```

## Reading

- `G_X_ws` is still the cleanest positive-recognition anchor.
- `J_X_ws` is still the cleanest *template-lifted* `j`-side coordinate of that anchor, but it inherits the same template-driven origin.
- The alternate signed lift `J_X15_ws` is materially closer to `J_f2_ws` on the hero side (the bridge quotient first breaks at `t^3`, not `t^1`).
- `Q_PB_ws` is still the strongest structural bridge seam.
- `J_f2_ws` now becomes the cleanest literature-facing next step.
- The new bridge says `J_f2_ws` is not secretly the same object as the current `P/B` seam inside the first focused bridge and named-`GG` boxes.
- Why:
  - `G_X_ws` is excellent as a ruler, but it is still obviously derived from the eta-side normalization.
  - `Q_PB_ws` has more internal structure on the hero side, but it still looks like a bridge object rather than a canonical named coordinate.
  - `J_f2_ws` sits in the middle in a useful way: less ad hoc than `Q_PB_ws`, less template-bound than `G_X_ws`, and directly attached to the Yui--Zagier Weber cubic.
- In plainer terms:
  - `G_X_ws` tells us whether the source template is lining up.
  - `Q_PB_ws` tells us where the hero bridge has internal stress lines.
  - `J_f2_ws` is the best candidate for the next "clean named coordinate" pass.

## Decision

```text
Decision
├─ Keep
│  ├─ G_X_ws as the source-faithful template anchor
│  ├─ J_X_ws as the anchor-derived canonical j lift (for clean `j`-side comparisons)
│  ├─ Q_PB_ws as the strongest structural bridge seam
│  └─ J_f2_ws as the cleanest literature-facing comparison lane
├─ Mark Complete
│  └─ the first focused J_f2_ws ↔ Q_PB_ws bridge pass and its nested follow-up
   └─ Do Next
   ├─ prefer a different named Weber coordinate or a different comparison bridge
   └─ do not blindly widen the same J/PB box again until a new source-faithful clue appears

- Update (2026-03-23):
  the focused named-`GG` modular-equation scan on the new lift-bridge quotient
  `Q_JX_ws` (and its follow-up `K_JX_ws`) stays flat in the first theorem-shaped
  boxes; the added eta / modular-unit boxes also stay flat, so `J_f2_ws` and
  `J_X_ws` still look genuinely distinct in the source-faithful lane rather
  than being a small named `GG` orbit disguise.

- Update (2026-03-23):
  the alternate anchor lift `J_X15_ws = (16*G_X_ws - 1)^3 / (3375*G_X_ws^2)`
  is now a better comparator against `J_f2_ws` than `J_X_ws` in the
  bridge-first sense: `Q_JX15_ws - 1` first breaks only at `t^3` (not `t^1`),
  but the focused named-`GG` + eta + modular-unit boxes on `Q_JX15_ws` and
  `K_JX15_ws` still stay flat, so this is a strong *alignment* clue without a
  positive recognition hit yet.

- Update (2026-03-23):
  added three pivot rails to compare `Q_JX15_ws` back to the strongest local seams
  (`Q_PB_ws` and `Q_XKJ_ws`) and to the completed `J/PB` quotient `Q_JPB_ws`.
  All three pivots break immediately at `t^1` and stay flat in the same focused
  eta / modular-unit / named-`GG` boxes, so this still reads as a “tighter
  alignment lens” rather than a positive recognition.
```

## Why `J_f2_ws` Still Leads, But Not The Same Box

- The first `Q_PB_ws` orbit pass has already been tried and stayed flat on the hero candidate.
- The first `J_f2_ws ↔ Q_PB_ws` bridge pass has now also been tried and stayed flat.
- `G_X_ws` is already doing its job as the positive-recognition anchor.
- `J_f2_ws` is the lane that most improves the balance between:
  - literature faithfulness
  - one-coordinate simplicity
  - compatibility with the existing constant-1 scan machinery
- So if the next move stays inside the named Weber city, `J_f2_ws` is still the most natural street to walk next.
- But the next box should not just be a wider replay of the same `J/PB` bridge; it should be a different named coordinate or a different comparison rail.

## Backing Artifacts

- Main tail-family artifact:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- P/B bridge note:
  - `notes/hero/CB60FD71D1D7_PB_BRIDGE_DECISION_NOTE.md`
- Formalization waypoint:
  - `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean`
- Research-side frontier hub:
  - `proofs/Proofs/HeroCaseFinalIdentity.lean`

## Latest Validation

- `python -m py_compile src/ramanujan_discovery/identification.py tests/test_identification.py`
- `$env:PYTHONPATH='src'; pytest -q tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_matches_true_gg tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_reports_hero_ratio_gap tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_skips_focused_named_gg_lane_in_smoke_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_matches_true_gg_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_reports_hero_gap_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_skips_named_gg_lane_in_smoke_profile tests/test_identification.py::test_cli_tail_note_writes_tail_family_note`
- Full hero regeneration command:
  - `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Full hero formalize regeneration command:
  - `$env:PYTHONPATH='src'; python -m ramanujan_discovery formalize --in results/verified.jsonl --candidate-id cb60fd71d1d7 --out notes/hero/CB60FD71D1D7_FORMALIZATION_NOTE.md --lean-out proofs/Proofs/Generated/Cb60fd71d1d7.lean`
- Latest successful runtimes:
  - tail-note: about `9706s`
  - formalize: about `1092s`
