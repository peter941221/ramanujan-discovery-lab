# Odd-Prime Descendant Checklist: `cb60fd71d1d7`

## Goal

- Keep the current best `GG/Weber` neighborhood in view.
- Avoid spending another full cycle just widening the same tiny odd-prime box.

## Current Read

```text
Objects Checked
├─ H_X_ws  -> 0 / 9 direct, 0 / 9 quotient
├─ H_gp_ws -> 0 / 9 direct, 0 / 9 quotient
├─ K_XR_ws -> 0 / 9 direct, 0 / 9 quotient
└─ L_XK_ws -> 0 / 9 direct, 0 / 9 quotient
```

## Checklist

- [x] Confirm `H_X_ws` tiny odd-prime lane is rendered in tail-note output.
- [x] Confirm `H_gp_ws` tiny odd-prime lane is rendered in tail-note output.
- [x] Confirm `K_XR_ws` tiny odd-prime lane is rendered in tail-note output.
- [x] Confirm `L_XK_ws` tiny odd-prime lane is rendered in tail-note output.
- [x] Add direct `H_X_ws` vs `H_gp_ws` comparison summary.
- [x] Keep public language at `unexplained candidate`.
- [x] Prefer the next pass to test a different named coordinate inside the same `GG/Weber` orbit.
- [x] Execute a first-round `P_ws` micro-scan.
- [x] Execute a first-round `B_ws` companion micro-scan.
- [x] Execute a direct named-`GG` bridge scan on `Q_XR_ws` / `K_XR_ws`.
- [x] Avoid broad anonymous box growth unless a named coordinate starts showing asymmetry.
- [x] Re-check whether the next best move should live on a quotient coordinate, a return-bridge coordinate, or a Morton/Weber compression instead.
- [x] Compare the `P_ws` / `B_ws` normalized lanes against each other before widening the odd-prime box again.
- [x] Build the direct `P/B` normalized bridge objects `D_PB_ws`, `Q_PB_ws`, and `K_PB_ws`.
- [x] Build the first nested `P/B` bridge objects `D_PK_ws`, `Q_PK_ws`, and `L_PK_ws`.
- [x] Decide whether the next bridge object should come from a `P/B` compression or from a new Weber one-coordinate function.
- [x] Run the first focused source-faithful one-coordinate orbit pass on `Q_PB_ws` / `K_PB_ws`.
- [x] Run the first focused source-faithful one-coordinate orbit pass on `Q_PK_ws` / `L_PK_ws`.
- [ ] Decide which different named Weber coordinate should get the next focused pass now that the first `P/B` orbit ladder stayed flat.
- [ ] Only after that, decide whether to return to a wider quotient bridge box beyond the current `Q_XR_ws` lane.

## Practical Decision

- Best current interpretation:
  - the tiny odd-prime lane is useful as a rejection filter, not yet as a recognition lane
  - the first `P_ws` / `B_ws` micro-boxes show some internal self-polynomial structure, but not a small source-faithful closure
- Best next move:
  - keep shifting coordinates before widening the odd-prime descendant box
  - the first preferred lane is now the explicit `P/B` normalized bridge with its first nested quotient layer
  - if that bridge ladder still stays bulky, only then promote another wider quotient box

## Next Coordinates

- next named Weber coordinate
  - preferred, because the first `Q_PB_ws` / `Q_PK_ws` source-faithful orbit pass is now done and still flat
- `Q_XR_ws`
  - keep only as the fallback derived quotient bridge if the next named Weber coordinate also stays bulky
- broader bridge lane
  - only after the current `P/B` orbit ladder and one more named Weber coordinate both miss

## Backing Artifacts

- Full decision note:
  - `notes/hero/CB60FD71D1D7_ODD_PRIME_DESCENDANT_DECISION_NOTE.md`
- Dedicated P/B bridge note:
  - `notes/hero/CB60FD71D1D7_PB_BRIDGE_DECISION_NOTE.md`
- Full regenerated tail note:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
