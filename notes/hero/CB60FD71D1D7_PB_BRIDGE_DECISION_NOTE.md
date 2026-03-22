# P/B Bridge Decision Note: `cb60fd71d1d7`

## Snapshot

- Date: `2026-03-22`
- Scope: explicit bridge between the Weber-Schlafli normalized lanes
  - `N_P_ws`
  - `N_B_ws`
- Bridge objects:
  - direct layer
    `D_PB_ws = N_B_ws - N_P_ws`
    `Q_PB_ws = N_B_ws / N_P_ws`
    `K_PB_ws = (Q_PB_ws - 1) / (leading term)`
  - nested layer
    `D_PK_ws = K_PB_ws - N_P_ws`
    `Q_PK_ws = K_PB_ws / N_P_ws`
    `L_PK_ws = (Q_PK_ws - 1) / (leading term)`

## Why This Note Exists

- The earlier odd-prime decision note said the next best local seam was likely the normalized `P/B` lane.
- This note records the result of actually building that seam, pushing it one normalized quotient layer deeper, and then testing the first source-faithful one-coordinate compression on the resulting bridge objects.
- It now also records the first focused cross-check from that seam back to the canonical classical-Weber `j`-side coordinate `J_f2_ws`.
- Analogy:
  - before, `P_ws` and `B_ws` were two neighboring ridges on the same map
  - now we have drawn the saddle path between them and checked whether it collapses into one cleaner coordinate

## Hero Result

```text
Hero P/B Bridge
├─ N_P_ws
│  └─ 3 / 18 self-polynomial hits, 0 elsewhere in first micro-box
├─ N_B_ws
│  └─ 3 / 18 self-polynomial hits, 0 elsewhere in first micro-box
├─ Direct bridge layer
│  ├─ D_PB_ws first failure: t^3 with coefficient 3/4
│  ├─ Q_PB_ws - 1 first failure: t^3 with coefficient 3/4
│  ├─ direct bridge polynomial hits: 1 / 3, only at degree 3
│  └─ Q_PB_ws and K_PB_ws: 3 / 18 self-polynomial hits, 0 elsewhere
├─ Nested bridge layer
│  ├─ D_PK_ws first failure: t^3 with coefficient 35/24
│  ├─ Q_PK_ws - 1 first failure: t^3 with coefficient 35/24
│  ├─ nested bridge polynomial hits: 1 / 3, only at degree 3
│  └─ Q_PK_ws and L_PK_ws: 3 / 18 self-polynomial hits, 0 elsewhere
└─ First source-faithful one-coordinate orbit pass
   ├─ true GG reference orbit is real: Q_PB_ref_ws and Q_PK_ref_ws each show 3 / 8 direct polynomial-prefix hits
   ├─ hero Q_PB_ws, K_PB_ws, Q_PK_ws, L_PK_ws all show 0 / 7 direct polynomial-prefix hits
   ├─ the same four hero objects also show 0 / 6 quotient-prefix hits
   └─ and 0 / 7 mixed-prefix hits, with 0 multiplicative and 0 fractional-linear hits throughout
```

## New `J/PB` Cross-Check

```text
P/B vs J
├─ Comparison rails
│  ├─ Canonical named coordinate: J_f2_ws
│  └─ Strongest local seam: Q_PB_ws
├─ Direct bridge
│  ├─ D_JPB_ws, Q_JPB_ws, K_JPB_ws already break at the first t^1 step on the direct hero series
│  ├─ bridge polynomial boxes: 0 / 3
│  └─ bridge fractional-linear box: 0 / 1
└─ Named-GG follow-up
   ├─ Q_JPB_ws, K_JPB_ws, Q_JKPB_ws, L_JKPB_ws all stay flat in the focused
   │  direct / quotient / mixed prefix boxes
   └─ the exact modular-equation template boxes also stay at 0 hits
```

## Reading

- The `P/B` seam is real.
- It is stronger than a vague "both lanes look structured" statement, because the bridge now gives a reproducible quotient object and a reproducible nested quotient follow-up.
- The direct bridge and the nested bridge both only light up once in the bounded polynomial box, and both only at degree `3`.
- That repeated degree-`3` behavior is useful evidence that the seam has internal structure rather than pure noise.
- The first source-faithful one-coordinate Weber orbit pass is also real on the true `GG` source: the direct ladders around `Q_PB_ref_ws` and `Q_PK_ref_ws` do register bounded polynomial recognition there.
- But on the hero candidate, the same orbit stays flat across the direct, quotient, and mixed prefix boxes.
- The new `J/PB` comparison matters because it shows the strongest seam is not already collapsing to the cleanest `j`-side coordinate inside the first focused bridge box.
- Practical interpretation:
  - the lane is not fake; the lock exists on the source side
  - but the hero bridge is still not turning that lock yet
  - so this is still a useful structural waypoint, not yet a final named Weber coordinate

## Decision

```text
Decision
├─ Keep
│  └─ the P/B normalized bridge ladder as a strong structural waypoint
├─ Mark Complete
│  ├─ the first source-faithful one-coordinate Weber orbit pass on Q_PB_ws / Q_PK_ws and their normalized follow-ups
│  └─ the first focused J_f2_ws ↔ Q_PB_ws bridge pass and its nested follow-up
└─ Do Next
   ├─ try a different named Weber coordinate inside the same city
   └─ or, if staying bridge-first, widen a different comparison rail rather than blindly replaying the same J/PB box
```

## Why This Beats Returning To `Q_XR_ws`

- `Q_XR_ws` is still useful as a bridge between the template and residual Weber branches.
- But the current named-`GG` pass on `Q_XR_ws` and `K_XR_ws` is flat.
- The `P/B` lane, by contrast, now has:
  - structure on both source-facing normalized coordinates
  - an explicit quotient bridge
  - a reproducible next quotient-followup object
  - and now one more reproducible nested quotient bridge layer beyond that
  - plus one explicit failed comparison against the canonical `j` lane, which rules out the easiest "same object in light disguise" reading
- In plainer terms:
  - `Q_XR_ws` is still a good diagnostic road
  - `P/B` is now the better excavation site

## Backing Artifacts

- Main tail-family artifact:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Public status summary:
  - `notes/hero/CB60FD71D1D7_PUBLIC_SUMMARY.md`
- Earlier broader decision note:
  - `notes/hero/CB60FD71D1D7_ODD_PRIME_DESCENDANT_DECISION_NOTE.md`
- Short execution checklist:
  - `notes/hero/CB60FD71D1D7_ODD_PRIME_DESCENDANT_CHECKLIST.md`

## Latest Validation

- `python -m py_compile src/ramanujan_discovery/identification.py tests/test_identification.py`
- `$env:PYTHONPATH='src'; pytest -q tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_matches_true_gg tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_reports_hero_ratio_gap tests/test_identification.py::test_scan_weber_class_invariant_bridge_box_skips_focused_named_gg_lane_in_smoke_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_matches_true_gg_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_reports_hero_gap_profile tests/test_identification.py::test_build_weber_j_pb_bridge_scan_skips_named_gg_lane_in_smoke_profile tests/test_identification.py::test_cli_tail_note_writes_tail_family_note`
- Full hero regeneration:
  - `$env:PYTHONPATH='src'; python -m ramanujan_discovery tail-note --in results/verified.jsonl --candidate-id cb60fd71d1d7 --tail-stages 3,4,5 --max-gap-depth 3 --out notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
  - latest successful runtime: about `9706s`

## Follow-up

- The next clean local question is not "is the old quotient box wider?"
- It was:
  - can `Q_PB_ws` or `Q_PK_ws` be compressed into a cleaner Weber one-coordinate function?
- Current answer:
  - the first focused source-faithful one-coordinate orbit pass still says no in the checked direct, quotient, and mixed prefix boxes
  - and the first focused `J_f2_ws ↔ Q_PB_ws` cross-check also says no in the checked bridge and named-`GG` boxes
- So the next question is:
  - which different named Weber coordinate should get the next focused pass, or which different comparison lane should take over now that both the seam-local orbit pass and the first J/PB bridge pass have been exhausted once?
