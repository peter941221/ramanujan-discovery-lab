# Formalization Prep: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized` (7 shared digits)
- Candidate template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=1;ner=6;nek=6;dc=1;ds=1;dr=3;dk=3`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`
- Build profile: `full`

## Current Theorem Status

- No complete source theorem is identified yet.
- This candidate is therefore **not ready** for a full Lean/Coq formalization of a final identity.
- The award-track Lean scaffold now packages an exact intermediate waypoint: finite convergents agree with the reduced-by-factor model over the rational-function reverse-equivalence layer, the page-43 nearest-shift cube is ruled out as a source-family-specific exact lane family, the full zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active is also ruled out exactly in both page-43 families, the nearest RR/cubic arithmetic-subsequence source lanes are excluded exactly, and the direct RR/cubic plus Heine-`cor2cf` low-stage mismatch layers are packaged as named local obstructions.
- The correct near-term target is to formalize exact local lemmas and keep bounded search evidence clearly separated from theorem-grade statements.

## Exact Objects To Formalize

- Reduced variable: `t = q^3`
- Target reciprocal object:

```text
C(t) = 1 + K_(n>=1) a_n / b_n
b0 = 1
a_n = t^n + t^2n
b_n = 1 + t^n
```

- Closest benchmark reciprocal object:

```text
R(t) = 1 + K_(n>=1) a_n / b_n
b0 = 1
a_n = t^n
b_n = 1
```

## Exact Reduction And Equivalence Witness

- Exact convergent gcd factors were checked through stage `8`.
- First common factors:

```text
g1 = t + 1
g2 = t^2 + 1
g3 = t^3 + 1
g4 = t^4 + 1
```

- After cancellation, the induced reduced-by-factor object begins:

```text
b0_red = 1
a1_red = t
b1_red = 1
a2_red = t^2
b2_red = t + 1
a3_red = t^4 + t^3
b3_red = t^2 + 1
a4_red = t^6 + t^4
b4_red = t^3 + 1
```

- Reverse equivalence transform stage scales:

```text
r1 = t + 1
r2 = (t^2 + 1)/(t + 1)
r3 = (t^3 + 1)/(t^2 + 1)
r4 = (t^4 + 1)/(t^3 + 1)
```

- These reverse scales are rational functions in `t`; the finite-stage reverse step now lives in Lean over `RatFunc Rat`, but a final identity still needs an infinite-object bridge and a fuller fraction-field coefficient layer beyond the current exact waypoint.

## Exact Lemma Candidates

### Direct 1-Step Bauer-Muir Obstructions

- RR source: `w0 = 0`, `w1 = t`, so the first transformed numerator is `t` instead of `t^2 + t`.
- Cubic source: `w0 = 0`, `w1 = t`, `w2 = t^2`, so the second transformed numerator is `t^4 - t^3 + t^2 - t` instead of `t^4 + t^2`.

### Simple Cubic Contraction Obstructions

- Odd contraction: initial term is `t^2 + t + 1` instead of the target `1`.
- Even contraction: first numerator does match as `t^2 + t`, but the first denominator is `t^4 + t^2 + 1` instead of `t + 1`.

### Heine `cor2cf` Odd/Even Branch Obstructions

- Lean mirror module: `proofs/Proofs/HeroCaseHeineCor2cf.lean`.
- In the relevant `a = 0` lane, the odd part already has initial term `lambda*t + 1` instead of `1`.
- The even part keeps initial term `1`, but its first numerator is `lambda*t` instead of `t^2 + t`.
- The odd-of-even branch changes the initial term to `(b*t + lambda*t^2 + lambda*t + 1)/(b*t + lambda*t^2 + 1)`, so it fails before the first nontrivial numerator.
- The even-of-even branch keeps initial term `1`, but its first numerator is `b*lambda*t^3 + lambda^2*t^5 + lambda^2*t^4 + lambda*t`. That numerator has no `t^2` term, so it cannot equal the target `t + t^2`.

### Exact `f2` / `gcf3` `n`-Dependent Equivalence Lane

- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.
- Prioritized source-family-specific lane: zero-shift `f2` / `gcf3` under arbitrary `n`-dependent equivalence factors.
- Write `m = t^(n-1)` and enforce the exact necessary identity

```text
alpha_n * (1 + t^(n-1)) = t^n * beta_(n-1) * beta_n
```

- In the zero-shift `f2` / `gcf3` lane this becomes

```text
alpha_n = lambda*m*t - a*b*m^2*t^2
beta_(n-1) = 1 + b*m + a*m*t
beta_n = 1 + b*m*t + a*m*t^2
```

- Residual polynomial:

```text
-a^2*m^3*t^4 - 2*a*b*m^3*t^3 - a*b*m^3*t^2 - a*b*m^2*t^2 - a*m^2*t^3 - a*m^2*t^2 - b^2*m^3*t^2 - b*m^2*t^2 - b*m^2*t + lambda*m^2*t + lambda*m*t - m*t
```

- `m^3` coefficient is `-a^2*t^4 - 2*a*b*t^3 - a*b*t^2 - b^2*t^2`; exact vanishing forces `a = 0`, `b = 0`.
- After that specialization, the `m^1` coefficient is `lambda*t - t`; exact vanishing forces `lambda = 1`.
- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes `t`, still nonzero.
- So no arbitrary `n`-dependent equivalence transformation sends the hero case into this zero-shift `f2` / `gcf3` lane.

### Exact `f4` / `gcf2` `n`-Dependent Equivalence Lane

- The same Lean mirror module also now covers the zero-shift `f4` / `gcf2` lane.
- In that lane the necessary identity becomes

```text
alpha_n = a*t + lambda*m*t
beta_(n-1) = 1 - a*t + b*m
beta_n = 1 - a*t + b*m*t
```

- Residual polynomial:

```text
-a^2*m*t^3 + a*b*m^2*t^3 + a*b*m^2*t^2 + 2*a*m*t^2 + a*m*t + a*t - b^2*m^3*t^2 - b*m^2*t^2 - b*m^2*t + lambda*m^2*t + lambda*m*t - m*t
```

- `m^0` coefficient is `a*t`; exact vanishing forces `a = 0`.
- After that, the `m^3` coefficient is `-b^2*t^2`; exact vanishing forces `b = 0`.
- Then the `m^1` coefficient is `lambda*t - t`; exact vanishing forces `lambda = 1`.
- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes `t`, still nonzero.
- So no arbitrary `n`-dependent equivalence transformation sends the hero case into this zero-shift `f4` / `gcf2` lane either.

### Exact Unit-Shift `lambda` Page-43 Lanes

- The same exact-equivalence layer now also covers the next nearby shift choice `lambda -> lambda*t` with zero `a`/`b` shifts.
- For `f2/gcf3`, the `m^3` coefficient is unchanged from the zero-shift lane, so exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the `m^1` coefficient becomes `lambda*t^2 - t`; no constant `lambda` can make that polynomial vanish identically.
- For `f4/gcf2`, exact vanishing again forces `a = 0`, then `b = 0`, and the surviving `m^1` coefficient becomes `lambda*t^2 - t`.
- So the first nonzero `lambda`-shift nearest lanes already fail before any final `m^2` obstruction is needed.

### Exact Mixed Unit-Shift `a`/`lambda` Page-43 Lanes

- The same exact-equivalence layer now also covers the first mixed nearby shift choice `a -> a*t`, `lambda -> lambda*t` with zero `b` shift.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^6 - 2*a*b*t^4 - a*b*t^3 - b^2*t^2`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the surviving `m^1` coefficient becomes `lambda*t^2 - t`; no constant `lambda` can make that polynomial vanish identically.
- For `f4/gcf2`, the constant coefficient `a*t^2` forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^2` forces `b = 0`, and the surviving `m^1` coefficient becomes `lambda*t^2 - t`.
- So the first mixed unit-`a` / unit-`lambda` lanes already fail before any final `m^2` obstruction is needed.

### Exact Mixed Unit-Shift `b`/`lambda` Page-43 Lanes

- The same exact-equivalence layer now also covers the first mixed nearby shift choice `b -> b*t`, `lambda -> lambda*t` with zero `a` shift.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^4 - 2*a*b*t^4 - a*b*t^3 - b^2*t^4`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the surviving `m^1` coefficient becomes `lambda*t^2 - t`; no constant `lambda` can make that polynomial vanish identically.
- For `f4/gcf2`, the constant coefficient `a*t` still forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^4` forces `b = 0`, and the surviving `m^1` coefficient becomes `lambda*t^2 - t`.
- So the first mixed unit-`b` / unit-`lambda` lanes also fail before any final `m^2` obstruction is needed.

### Exact Mixed Unit-Shift `a`/`b`/`lambda` Page-43 Lanes

- The same exact-equivalence layer now also covers the first full nearby shift choice `a -> a*t`, `b -> b*t`, `lambda -> lambda*t`.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^6 - 2*a*b*t^5 - a*b*t^4 - b^2*t^4`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the surviving `m^1` coefficient becomes `lambda*t^2 - t`; no constant `lambda` can make that polynomial vanish identically.
- For `f4/gcf2`, the constant coefficient `a*t^2` forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^4` forces `b = 0`, and the surviving `m^1` coefficient becomes `lambda*t^2 - t`.
- So the first full three-parameter nearest lane also fails before any final `m^2` obstruction is needed.
- Together with the zero-shift and the other seven nearest-shift cases, this closes the full `{0,1}^3` nearest-shift cube at the current page-43 exact-equivalence level.
- Lean now exposes that cube not only as summary theorems, but also as Bool-parameterized theorems over the shift bits.

### Exact Unit-Shift `a` Page-43 Lanes

- The same exact-equivalence layer now also covers the next nearby shift choice `a -> a*t` with zero `b`/`lambda` shifts.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^6 - 2*a*b*t^4 - a*b*t^3 - b^2*t^2`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the `m^1` coefficient becomes `lambda*t - t`; exact vanishing forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- For `f4/gcf2`, the constant coefficient `a*t^2` already forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^2` forces `b = 0`, and then the `m^1` coefficient `lambda*t - t` forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- So the first nonzero `a`-shift nearest lanes also fail by an exact final obstruction.

### Exact Unit-Shift `b` Page-43 Lanes

- The same exact-equivalence layer now also covers the next nearby shift choice `b -> b*t` with zero `a`/`lambda` shifts.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^4 - 2*a*b*t^4 - a*b*t^3 - b^2*t^4`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the `m^1` coefficient becomes `lambda*t - t`; exact vanishing forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- For `f4/gcf2`, the constant coefficient `a*t` still forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^4` forces `b = 0`, and then the `m^1` coefficient `lambda*t - t` forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- So the first nonzero `b`-shift nearest lanes also fail by an exact final obstruction.

### Exact Mixed Unit-Shift `a`/`b` Page-43 Lanes

- The same exact-equivalence layer now also covers the first mixed nearby shift choice `a -> a*t`, `b -> b*t` with zero `lambda` shift.
- For `f2/gcf3`, the `m^3` coefficient becomes `-a^2*t^6 - 2*a*b*t^5 - a*b*t^4 - b^2*t^4`, and exact vanishing still forces `a = 0`, `b = 0`.
- After that specialization, the `m^1` coefficient becomes `lambda*t - t`; exact vanishing forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- For `f4/gcf2`, the constant coefficient `a*t^2` already forces `a = 0`.
- After that, the `m^3` coefficient `-b^2*t^4` forces `b = 0`, and then the `m^1` coefficient `lambda*t - t` forces `lambda = 1`.
- The surviving `m^2` coefficient is then `t`, still nonzero.
- So the first mixed unit-`a` / unit-`b` lane also fails by an exact final obstruction.

### Exact Zero-Shift Polynomial Single-Prefactor Page-43 Lanes

- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.
- This exact layer now covers the whole polynomial sub-box `phi in {1, 1+t}` with zero shifts and at most one non-plain prefactor active.
- `f2/gcf3`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining stage-1 numerator difference is `lambda*t - t^2 - t`.
- `f2/gcf3`, `phi_a = 1+t`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining stage-1 numerator difference is `lambda*t - t^2 - t`.
- `f2/gcf3`, `phi_a = 1`, `phi_b = 1+t`, `phi_lambda = 1`: denominator matching forces `a = -1, b = 1`, but the remaining stage-1 numerator difference is `lambda*t + t^3 - t`.
- `f2/gcf3`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1+t`: denominator matching forces `a = 0, b = 1`, stage-1 numerator matching then forces `lambda = 1`, but the stage-2 numerator becomes `t^3 + t^2` instead of `t^4 + t^2`.
- `f4/gcf2`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining stage-1 numerator difference is `lambda*t - t^2 - t`.
- `f4/gcf2`, `phi_a = 1+t`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining stage-1 numerator difference is `lambda*t - t^2 - t`.
- `f4/gcf2`, `phi_a = 1`, `phi_b = 1+t`, `phi_lambda = 1`: denominator matching is already incompatible.
- `f4/gcf2`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1+t`: denominator matching forces `a = 0, b = 1`, stage-1 numerator matching then forces `lambda = 1`, but the stage-2 numerator becomes `t^3 + t^2` instead of `t^4 + t^2`.
- So the whole zero-shift polynomial single-prefactor sub-box is excluded exactly in both page-43 families.

### Exact Zero-Shift Reciprocal Single-Prefactor Page-43 Lanes

- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.
- This exact layer now also covers the reciprocal sub-box `phi in {1/(1+t)}` via cross-multiplied coefficient identities.
- `f2/gcf3`, `phi_a = 1/(1+t)`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining cross-multiplied stage-1 numerator difference has numerator `lambda*t - t^2 - t`.
- `f2/gcf3`, `phi_a = 1`, `phi_b = 1/(1+t)`, `phi_lambda = 1`: cross-multiplied denominator matching is already incompatible.
- `f2/gcf3`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1/(1+t)`: denominator matching forces `a = 0, b = 1`, but the remaining cross-multiplied stage-1 numerator difference has numerator `lambda*t - t^3 - 2*t^2 - t`.
- `f4/gcf2`, `phi_a = 1/(1+t)`, `phi_b = 1`, `phi_lambda = 1`: denominator matching forces `a = 0, b = 1`, but the remaining cross-multiplied stage-1 numerator difference has numerator `lambda*t - t^2 - t`.
- `f4/gcf2`, `phi_a = 1`, `phi_b = 1/(1+t)`, `phi_lambda = 1`: cross-multiplied denominator matching is already incompatible.
- `f4/gcf2`, `phi_a = 1`, `phi_b = 1`, `phi_lambda = 1/(1+t)`: denominator matching forces `a = 0, b = 1`, but the remaining cross-multiplied stage-1 numerator difference has numerator `lambda*t - t^3 - 2*t^2 - t`.
- So the whole zero-shift reciprocal single-prefactor sub-box is excluded exactly in both page-43 families.
- Taken together, the current zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` is now fully exactified at theorem grade.

## Bounded Exact Exclusion Results

- Arithmetic subsequence contractions up to stride `4` with stage comparison depth `3`: RR hits `0`, cubic hits `0`.
- These are exact statements for the bounded class being checked, but they do not identify a final source theorem.
- Page-43 monomial substitutions in the current `[-3,3]` shift box with `3` matched stages: `f2/gcf3` hits `0`, `f4/gcf2` hits `0`.
- Page-43 low-complexity rational-prefactor box: `phi in {1, 1+t, 1/(1+t)}` with at most `3` non-plain prefactors active, shift box `[-1,1]`, and `3` matched stages: `f2/gcf3` hits `0`, `f4/gcf2` hits `0`.
- These are bounded symbolic searches, useful for narrowing the theorem statement but not substitutes for a full origin proof.
- Within that bounded prefactor box, the full zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active is now theorem-grade exact.

## Formalization Order

1. Formalize generalized continued fractions and convergent recurrence for finite truncations.
2. Reuse the exact convergent-factor reduction theorem for the candidate-side local model.
3. Bridge the current rational-function reverse-equivalence layer from finite convergents to the infinite object.
4. Extend `Proofs/HeroCasePage43Equivalence.lean` beyond the currently formalized nearest-shift cube plus the exact zero-shift single-prefactor box, especially if shifted or multi-prefactor rational lanes become theorem-relevant.
5. Formalize the direct 1-step Bauer-Muir obstruction lemmas against the reduced target.
6. Formalize odd/even contraction reconstruction together with the cubic and Heine-`cor2cf` low-stage mismatch lemmas.
7. Defer the bounded search exclusions until a final theorem statement makes them clearly necessary.
8. Do not start a full Lean/Coq origin theorem until a unique source family or exact identity is identified.

## Why This Is Still Not A Full Theorem

- The current exact lemmas only rule out nearby transforms and simple contraction sources.
- They do not prove what the candidate *is*.
- A full formal proof needs a final theorem statement of the form `C(t) = known_object(t)` or a uniquely characterizing theorem that has not been found yet.

## Award-Track Endgame Hook

- This candidate matches the current hero-case structural signature in reduced variable `t`.
- Award-track target module (Lean scaffold): `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current state: the module compiles, carries status marker `exclusion_waypoint`, and now exposes the current exact waypoint as a certificate object `currentExactWaypointCertificate`, whose two fields are proved by `finiteConvergentReductionWaypoint_true` and `knownSourceOrbitExclusionWaypoint_true` and then repackaged by `exactWaypointStatement_true`. The lower-level exact ingredients now also explicitly include `mortonSquaredCoordinateExcluded`, `weberSchlafliCoordinateExcluded`, and the older `directLocalObstructions` / `simpleCor2cfBranchesExcluded` layer, while the research-side hand-off now threads through `mortonNamedCoordinateResearchWaypoint_true`, `currentNamedWeberOrbitResearchCertificate`, and `namedWeberOrbitResearchWaypoint_true` before the combined exact-plus-research frontier is exposed as `currentRecognitionFrontierCertificate` / `currentRecognitionFrontierWaypoint_true`.
- Current Weber-Schlafli exact shell: `proofs/Proofs/HeroCaseWeberSchlafliCoordinateObstruction.lean` packages `P_ws = (1/F - F) / 2` and the first Weber-Schlafli template `P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0` as a seven-class repeated witness table `(2,-1)`, `(2,3)`, `(2,8)`, `(4,3)`, `(6,3)`, `(8,3)`, `(10,3)` rather than a universal one-class obstruction.
- Current Morton named-coordinate waypoint: `proofs/Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean` packages the unified source-faithful stack `X_mt = F^2`, `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`, `P_ws = (1/F - F) / 2`, and `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`, with universal obstruction classes already visible in the `X_mt`, `T_mt`, and `B_ws` lanes and an exact witness-table shell on `P_ws`.
- Current Morton square-coordinate exact shell: `proofs/Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean` still packages the explicit source-faithful lane `X_mt = F^2`, the Proposition 3.2 template `X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 = 0`, and the current universal constant-term obstruction `(t^0, 4)` across all sampled tail-family objects.
- Current modular-coordinate waypoint: `proofs/Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean` records the first failure of `F / W_34`, the normalized probe `G_W34`, the deeper follow-up `G2_W34`, and the current empty small correction boxes at both normalized layers.
- Current Weber residual waypoint: `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean` records the primary residual choice `G_g12_ws`, the exact coordinate bridge tying `G_p12_ws` back to it, the classical Weber `f2` tri-product coordinate `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2) = G_g12_ws*G_p12_ws` together with its normalized follow-up `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)`, and now also records the source-backed reading that Berndt--Chan--Zhang identify Ramanujan-Weber `G_n` / `g_n` with classical Weber `f` / `f1` while Yui--Zagier supplies the classical Weber `f`, `f1`, `f2` trio, so this `g12_ws` / `p12_ws` / `G_f2_ws` shell should be read as that named Weber trio in the project's normalization. The same waypoint then keeps the derived quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2`, the exact quotient-coordinate bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0` with `Q_gp_ws = p12_ws / g12_ws`, the template-normalized coordinate `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`, its normalized follow-up `H_X_ws = (G_X_ws - 1) / (4*t^1)`, the direct follow-up bridge `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws`, the quotient-bridge follow-up `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, the next stripped comparison back to the template-normalized lane `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws`, its quotient follow-up `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)`, the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws`, the normalized follow-up `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`, the theorem-shaped return-bridge closure shell saying `Q_XK_ws` and `L_XK_ws` still have `0` hits in the checked direct / quotient / mixed named `GG` modular-equation boxes, and now also packages the exact Chan--Huang obstruction quartets on the later return bridge: for `Q_XK_ws`, direct failures `(-9/2, -6)` against `GG3/GG4` and quotient failures `(-3, -9/2)` against `Q_3/Q_4`; for `L_XK_ws`, direct failures `(593/10, 1186/15)` and quotient failures `(593/15, 593/10)`, alongside the current empty theorem-shaped uniqueness / closure boxes on the hero-side follow-up objects.
- `finalIdentityStatement` is still a placeholder because no final closed form is identified yet.
- Only replace the placeholder after a concrete closed form is identified.
