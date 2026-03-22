# `cb60fd71d1d7` Public Summary

## What It Is

`cb60fd71d1d7` is the current strongest `review` candidate in the Ramanujan Discovery Lab search pipeline.

Its template is

```text
1 / (1 + (q^3 + q^6)/(1 + q^3 + (q^6 + q^12)/(1 + q^6 + ...)))
```

and its nearest built-in benchmark is `rogers_ramanujan_q3_normalized`.

## Why It Matters

- It is not random search noise.
- Symbolically, it matches `RR(q^3)` through `q^9`.
- The first visible divergence from `RR(q^3)` occurs at `q^12`.
- It is materially closer to `RR(q^3)` than to `cubic(q^3)` in the first visible symbolic orders.

## Current Best Description

The candidate looks like a structured Rogers-Ramanujan-adjacent hybrid perturbation:

- it keeps the `RR(q^3)` main numerator progression
- it adds a cubic-style extra numerator term
- it also adds a Rogers-Ramanujan-type denominator perturbation

## What Has Been Ruled Out

- No direct constant-parameter specialization of Ramanujan Entry 6.4.4 / the nearest page-43 ratio family has been found.
- No exact match for the mixed pattern

```text
(q^(3n) + q^(6n)) / (1 + q^(3n))
```

has been found in the primary sources checked so far.
- No source-faithful `GG` modular-equation hit has been found in the first
  explicit box built from `GG(q^3)`, `GG(-q^3)`, `GG(q^6)`, `GG(q^9)`,
  `GG(q^12)`, `GG(q^15)`, `GG(q^21)`, `GG(q^33)` and the corresponding
  quotient coordinates against `GG(q^3)`, even after adding the mixed
  quotient-coordinate pass that keeps `GG(q^3)` explicit while scanning
  low-degree corrections in `GG(-q^3)/GG(q^3)` and `GG(q^{3p})/GG(q^3)`.
- The first exact Chan--Huang modular-equation polynomials in the `q^3` and
  `q^4` lanes also still give `0` hits, both in the direct
  `F` vs. `GG(q^{3m})` form and in the quotient-coordinate
  `F` vs. `GG(q^{3m})/GG(q^3)` form.
- The Morton 2024 tail-family lane now also isolates the explicit squared
  coordinate from Proposition 3.2 / Theorem B:
  `X_mt = F^2`.
  On all `12` sampled tail-family objects, the first exact template
  `X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 = 0`
  still gives `0` hits and already fails in one uniform constant-term class
  with coefficient `4`.
- The same Morton square lane now also carries the paper's linear-fractional
  transformed square coordinate from equation `(3.6)`:
  `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`, `sigma = -1 + sqrt(2)`.
  That transformed square lane also gives `0` hits on all `12` sampled
  tail-family objects, and it likewise collapses immediately to one uniform
  constant-term obstruction, now with coefficient `8`.
- The next source-faithful Weber-Schlafli `P/B` branch also still gives `0`
  hits:
  `P_ws = (1/F - F) / 2` misses in repeated low-order classes, while the
  direct companion
  `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)` already fails in one uniform
  constant-term class across all `12` sampled tail-family objects.
- A deeper Ramanujan-Weber class-invariant compression now also exists on the
  same normalized `GG` variable:
  `g12_ws = 4*t*(Z_g - 1/Z_g)` and
  `p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)`, with
  `Z_g = ((1 - t*F^2)^2) / (4*t*F^2)`.
  On the true `GG` source these recover the explicit named objects
  `(t^2; t^4)_inf^12` and `(-t^2; t^4)_inf^12`.
  Berndt--Chan--Zhang then identifies the Ramanujan-Weber `G_n` / `g_n`
  normalization with classical Weber `f` / `f1`, while Yui--Zagier supplies
  the classical Weber `f`, `f1`, `f2` trio, so the current
  `g12_ws` / `p12_ws` / `G_f2_ws` shell is now being read as that named Weber
  trio in the project's normalization rather than as an anonymous product
  gadget.
  But on the current hero
  tail-family ladder they still give `0` hits after the first direct-template,
  eta / modular-unit, and plus-Pochhammer correction passes.
- That same Weber pair is now tied together by an exact bridge
  `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`.
  This lets the project keep `G_g12_ws` as the current primary eta-side
  residual, treat `G_p12_ws` as the constrained companion, and insert the
  classical Weber `f2` tri-product coordinate
  `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2) = G_g12_ws*G_p12_ws`.
  On the hero sample, that product coordinate first differs from `1` at `t^1`
  with coefficient `-4`, and its normalized follow-up
  `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)` then first differs from `1` at `t^1`
  with coefficient `1/2`; both layers still give `0` hits in the same first
  narrow theorem-shaped closure boxes, so this product lane currently looks
  more like another structured obstruction than like the final RHS.
  The same compression then inserts the
  tighter quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2` together with the
  exact bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0`,
  `Q_gp_ws = p12_ws / g12_ws`.
  The same compression now also exposes the template-normalized coordinate
  `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`, which is
  exact `1` on true normalized `GG`; on the hero sample it first differs from
  `1` at `t^1` with coefficient `4`, and its first follow-up
  `H_X_ws = (G_X_ws - 1) / (4*t^1)` still stays outside the same first
  theorem-shaped closure boxes.
  That same anchor branch also lifts to a canonical Weber `j`-side constant-1
  object
  `J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)`.
  On the hero sample, `J_X_ws - 1` first fails at `t^1` with coefficient
  `56/17`, and the normalized follow-up
  `H_J_X_ws = (J_X_ws - 1) / (56/17*t^1)` then first fails at `t^1` with
  coefficient `1083/238`; both layers still give `0` hits in the same first
  theorem-shaped closure boxes.
  With the classical Weber naming closure in place, this
  `Q_gp_ws -> X_g_ws -> G_X_ws`
  branch is now best read as the source-faithful classical Weber
  quotient/template lane, and it still looks cleaner than the product branch as
  the main constructive trunk.
  The same lane now also compares the two hero-side normalized follow-ups
  directly:
  `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws` both first
  differ at `t^2` with coefficient `-24`, and no degree-`<= 3` polynomial or
  one-coordinate fractional-linear bridge was found in that first comparison
  box.
  The quotient side of that same bridge now also carries a stripped follow-up
  `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, which on the hero sample first
  differs from `1` at `t^1` with coefficient `2` and still gives `0` hits in
  the first self-polynomial, self-fractional-linear, self-quotient
  finite-product, eta / modular-unit, and plus-Pochhammer boxes.
  That stripped quotient lane is now also compared back to the
  template-normalized branch:
  `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws` both first
  differ at `t^1` with coefficient `-5/2`, and the next quotient follow-up
  `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)` still gives `0` hits in the same
  first narrow theorem-shaped closure boxes.
  Even after pushing `Q_XK_ws` and `L_XK_ws` back through the named
  Chan--Huang `GG` modular-equation basis, both objects still give `0` hits in
  the checked direct, quotient, and mixed quotient-coordinate boxes.
  The exact first-failure fingerprints are now also explicit:
  `Q_XK_ws` misses the direct `GG3/GG4` templates with coefficients
  `(-9/2, -6)` and the quotient `Q_3/Q_4` templates with coefficients
  `(-3, -9/2)`;
  `L_XK_ws` misses the same direct templates with
  `(593/10, 1186/15)` and the quotient templates with
  `(593/15, 593/10)`.
  That is a strong hint that the return bridge has not yet settled onto the
  present Chan--Huang basis itself.
  So the current honest reading is:
  these are derived return-bridge coordinates, not recovered named
  literature coordinates.
  In plainer terms: we have a bridge back into the right city, but not yet the
  street address. The next constructive guess should therefore be a deeper
  Weber / modular-function coordinate, not a wider re-scan of the same
  `GG3/GG4/Q_3/Q_4` box.
  After that compression, the next positive-recognition step still stays
  focused on the residual quotient `R_gp_ws = G_p12_ws / G_g12_ws`, which
  still gives `0` hits in the first narrow self-polynomial,
  self-fractional-linear, self-quotient finite-product, eta / modular-unit,
  and plus-Pochhammer boxes.
  The same lane now also keeps a normalized follow-up
  `H_gp_ws = (R_gp_ws - 1) / (96*t^3)` under watch, and that follow-up also
  remains outside the first theorem-shaped closure boxes.
- That direct Weber companion obstruction is now also packaged in Lean under
  `proofs/Proofs/HeroCaseWeberCompanionObstruction.lean`, so this is no longer
  only a notebook-level negative result.
- The new focused Weber residual bridge is likewise packaged under
  `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean` and exposed
  from `Proofs/HeroCaseFinalIdentity.lean` as the current theorem-shaped
  hand-off after the `g/p` compression.

## What Has Not Been Claimed

- This is not a proof of novelty.
- This is not yet a publishable new-formula claim.
- It is still best treated as a high-value audit target.

## Current Direction Snapshot

- A focused odd-prime descendant sanity check still runs on the current
  Weber-side follow-up ladder:
  `H_X_ws`, `H_gp_ws`, `K_XR_ws`, and `L_XK_ws`.
  In the current tiny box, all four objects still give
  `0 / 9` direct hits and `0 / 9` quotient hits.
- The old shortlist has now been executed once rather than merely proposed:
  - `P_ws`: on the hero sample, the leading-term-normalized coordinate
    `N_P_ws = P_ws / t^3` gives `3 / 18` self-polynomial hits, but `0` hits in
    fractional-linear, self-quotient finite-product, eta, modular-unit, and
    plus-Pochhammer boxes
  - `B_ws`: the companion normalization `N_B_ws = B_ws / 2` shows the same
    shape, again with `3 / 18` self-polynomial hits and `0` hits elsewhere in
    the first micro-box
  - the direct `P/B` bridge is now explicit:
    `D_PB_ws = N_B_ws - N_P_ws` and `Q_PB_ws = N_B_ws / N_P_ws` both first
    differ at `t^3` with coefficient `3/4`; the bridge polynomial box lights up
    only once at degree `3`, while `Q_PB_ws` and its normalized follow-up
    `K_PB_ws` again show `3 / 18` self-polynomial hits and `0` hits in the
    surrounding fractional-linear / eta / product boxes
  - the same `P/B` lane now also has a first nested quotient layer:
    `D_PK_ws = K_PB_ws - N_P_ws` and `Q_PK_ws = K_PB_ws / N_P_ws` both first
    differ at `t^3` with coefficient `35/24`; that nested bridge polynomial box
    also lights up only once at degree `3`, while `Q_PK_ws` and its normalized
    follow-up `L_PK_ws` keep the same `3 / 18` self-polynomial and `0`-hit
    side-box profile
  - the first source-faithful one-coordinate Weber orbit pass on that same seam
    is now also executed:
    on the true `GG` source, the direct ladders around `Q_PB_ref_ws` and
    `Q_PK_ref_ws` each show `3 / 8` direct polynomial-prefix hits, so the lane
    itself is not empty.
    But on the hero candidate, `Q_PB_ws`, `K_PB_ws`, `Q_PK_ws`, and `L_PK_ws`
    all still give `0 / 7` direct polynomial-prefix hits, `0 / 6` quotient
    polynomial-prefix hits, and `0 / 7` mixed-prefix hits, with `0`
    multiplicative and `0` fractional-linear hits throughout that focused box.
  - the canonical classical-Weber `j`-side lane is now also explicit:
    `J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)`.
    It is the cleanest literature-facing one-coordinate lane in the current
    Weber city.
  - the first focused `J_f2_ws ↔ Q_PB_ws` comparison is now also executed:
    on the direct hero series, the direct and nested bridge boxes both stay at
    `0` hits in the checked polynomial / fractional-linear bridge boxes, and
    `Q_JPB_ws`, `K_JPB_ws`, `Q_JKPB_ws`, and `L_JKPB_ws` all stay flat in the
    focused named-`GG` direct / quotient / mixed prefix boxes
  - `Q_XR_ws`: the direct named `GG` modular-equation pass now also reports
    `0` hits in the checked direct, quotient, and mixed boxes, and the
    normalized follow-up `K_XR_ws` stays equally flat in that same named lane
- Reading:
  - `P_ws` and `B_ws` are not blank walls anymore; they do carry internal
    low-degree self-polynomial structure
  - the new `P/B` bridge ladder says the two lanes are related, but still
    only by bulky degree-`3` comparisons, not by a clean literature-facing coordinate
  - the first source-faithful one-coordinate orbit pass on that seam is now
    also flat on the hero candidate, even though the same lane does recognize
    the true `GG` source in the bounded direct-prefix box
  - `J_f2_ws` is now the cleanest literature-facing Weber coordinate, but the
    first `J_f2_ws ↔ Q_PB_ws` cross-check also stayed flat, so the seam does
    not simply collapse back to the first canonical `j` lane
  - `Q_XR_ws` remains a useful diagnostic bridge object, not a recovered
    literature coordinate
- In plainer terms:
  the next move is still a better coordinate in the same city, but the first
  `P/B` self-orbit echo test and the first `J/PB` bridge test have now both
  been exhausted once and stayed flat.
  So the repo should next try a different named Weber coordinate or only then
  widen a different comparison bridge.

## Current Working Formula

To low order,

```text
candidate = RR(q^3) * (1 - q^12 + 2 q^15 - q^18) + O(q^31)
```

## Sources

- Bowman, Mc Laughlin, Wyshinski, 2006: https://arxiv.org/abs/1901.00584
- Lee, Mc Laughlin, Sohn, 2020: https://arxiv.org/abs/1906.11991
- Chan, Huang, 1997: https://mrc.sdu.edu.cn/ziliao/8.pdf
- Akkarapakam, Morton, 2024: https://nyjm.albany.edu/j/2024/30-36.html
- Adiga, Kim, et al., 2017: https://www.filomat.org/index.php/filomat/article/download/4026/4026/31-13-2-4026.pdf
- Berndt, Chan, Zhang, 1997: https://mrc.sdu.edu.cn/ziliao/10.pdf
- Yui, Zagier, 1997: https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1090/S0025-5718-97-00854-5/fulltext.pdf
