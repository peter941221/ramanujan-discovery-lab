# Odd-Prime Descendant Decision Note: `cb60fd71d1d7`

## Snapshot

- Date: `2026-03-21`
- Scope: focused tiny odd-prime descendant probes inside the current `GG/Weber` follow-up lane
- Objects compared:
  - `H_X_ws`
  - `H_gp_ws`
  - `K_XR_ws`
  - `L_XK_ws`
- Descendant ladder:
  - direct: `GG5`, `GG7`, `GG11`
  - quotient: `Q_5`, `Q_7`, `Q_11`
- Tiny box used:
  - order `8`
  - degree box `(1,)`
  - multiplicative and fractional-linear prefix checks with max exponent `2`

## Result

```text
Odd-Prime Tiny Lane
├─ H_X_ws
│  ├─ direct hits: 0 / 9
│  └─ quotient hits: 0 / 9
├─ H_gp_ws
│  ├─ direct hits: 0 / 9
│  └─ quotient hits: 0 / 9
├─ K_XR_ws
│  ├─ direct hits: 0 / 9
│  └─ quotient hits: 0 / 9
└─ L_XK_ws
   ├─ direct hits: 0 / 9
   └─ quotient hits: 0 / 9
```

## Reading

- The current odd-prime tiny lane does not distinguish `H_X_ws` from `H_gp_ws`.
- The first stripped return object `K_XR_ws` also stays unresolved in the same tiny lane.
- The deeper return object `L_XK_ws` still shows the same no-hit profile.
- Practical interpretation:
  - this looks less like "the box is just slightly too small"
  - and more like "the coordinate is still not the right literature-facing compression"

## Decision

```text
Decision
├─ Keep
│  └─ the named `GG/Weber` orbit as the current best neighborhood
├─ Do Not Do Next
│  └─ do not spend the next cycle merely widening the same odd-prime tiny box
└─ Prefer Next
   ├─ try a different named coordinate or quotient compression in the same orbit
   ├─ compare against return-bridge objects before broad anonymous scans
   └─ keep public language at `unexplained candidate`
```

## Executed Coordinate Shortlist

```text
Executed Shortlist
├─ 1) P_ws
│  └─ first Weber-Schlafli one-coordinate lane
├─ 2) B_ws
│  └─ explicit Weber-Schlafli companion lane
└─ 3) Q_XR_ws
   └─ current derived quotient bridge between the two hero-side normalized branches
```

- `P_ws` result:
  - on the hero sample, the first leading-term-normalized micro-box is
    `N_P_ws = P_ws / t^3`
  - it gives `3 / 18` self-polynomial hits and `0` hits in the checked
    fractional-linear, self-quotient, eta, modular-unit, and plus-Pochhammer boxes
- `B_ws` result:
  - the companion normalization is `N_B_ws = B_ws / 2`
  - it shows the same first-round profile: `3 / 18` self-polynomial hits and `0`
    hits elsewhere in the current micro-box
- `P/B` bridge result:
  - the direct bridge is now explicit:
    `D_PB_ws = N_B_ws - N_P_ws` and `Q_PB_ws = N_B_ws / N_P_ws` both first fail at
    `t^3` with coefficient `3/4` on the hero candidate
  - the direct bridge polynomial box lights up only once, at degree `3`
  - the bridge quotient `Q_PB_ws` and its normalized follow-up `K_PB_ws` still
    show `3 / 18` self-polynomial hits and `0` hits elsewhere in the first micro-box
  - the same lane now also has a nested quotient bridge:
    `D_PK_ws = K_PB_ws - N_P_ws` and `Q_PK_ws = K_PB_ws / N_P_ws` both first fail at
    `t^3` with coefficient `35/24`; the nested bridge polynomial box also lights up
    only once, again at degree `3`, while `Q_PK_ws` and `L_PK_ws` keep the same
    `3 / 18` self-polynomial and `0`-hit side-box profile
- `Q_XR_ws` result:
  - the direct named `GG` modular-equation pass on `Q_XR_ws` now gives `0` hits
    in the checked direct, quotient, and mixed boxes
  - the normalized follow-up `K_XR_ws` is equally flat in that same named lane
- Reading:
  - `P_ws` and `B_ws` now look more like structured intermediate compressions
    than like immediate closed-form recognitions
  - the new `P/B` bridge is informative, but still too bulky to count as a clean recognition lane
  - `Q_XR_ws` remains a clean comparator, but not yet a recovered literature coordinate

## Evidence Source

- Real CLI confirmation file:
  - `tmp/tail-note-fast-full.md`
- Full regenerated hero note:
  - `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- Relevant rendered sections now include:
  - `H_X_ws` odd-prime descendant preview + micro-scan
  - `H_gp_ws` odd-prime descendant preview + micro-scan
  - `K_XR_ws` odd-prime descendant preview + micro-scan
  - `L_XK_ws` odd-prime descendant preview + micro-scan
  - direct `H_X_ws` vs `H_gp_ws` comparison summary

## Follow-up

- The full hero regeneration path is now working again, but still heavy.
- The latest full note build took about `4453.12s`, so this lane is still reproducible but noticeably heavier after the focused `P/B` orbit pass.
- Until a coordinate starts breaking symmetry inside a named lane, the next local move should be a better coordinate choice, not merely a larger odd-prime box.
- Best current local seam: the `P_ws` / `B_ws` normalized bridge lane now includes a first nested quotient layer and its first focused one-coordinate orbit pass; that pass stayed flat on the hero candidate, so the next local move should try a different named Weber coordinate or only then widen the older `Q_XR_ws` comparison box.
