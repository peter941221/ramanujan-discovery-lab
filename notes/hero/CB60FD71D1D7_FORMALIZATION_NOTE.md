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

- These reverse scales are rational functions in `t`, so a full formalization of this step likely needs a fraction-field coefficient layer in addition to the current polynomial one.

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

## Bounded Exact Exclusion Results

- Arithmetic subsequence contractions up to stride `4` with stage comparison depth `3`: RR hits `0`, cubic hits `0`.
- These are exact statements for the bounded class being checked, but they do not identify a final source theorem.
- Page-43 monomial substitutions in the current `[-3,3]` shift box with `3` matched stages: `f2/gcf3` hits `0`, `f4/gcf2` hits `0`.
- Page-43 low-complexity rational-prefactor box: `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active, shift box `[-0,0]`, and `2` matched stages: `f2/gcf3` hits `0`, `f4/gcf2` hits `0`.
- These are bounded symbolic searches, useful for narrowing the theorem statement but not substitutes for a full origin proof.

## Formalization Order

1. Formalize generalized continued fractions and convergent recurrence for finite truncations.
2. Reuse the exact convergent-factor reduction theorem for the candidate-side local model.
3. Add a rational-function or fraction-field coefficient layer for the reverse equivalence transform.
4. Formalize the direct 1-step Bauer-Muir obstruction lemmas against the reduced target.
5. Formalize odd/even contraction reconstruction together with the cubic and Heine-`cor2cf` low-stage mismatch lemmas.
6. Defer the bounded search exclusions until a final theorem statement makes them clearly necessary.
7. Do not start a full Lean/Coq origin theorem until a unique source family or exact identity is identified.

## Why This Is Still Not A Full Theorem

- The current exact lemmas only rule out nearby transforms and simple contraction sources.
- They do not prove what the candidate *is*.
- A full formal proof needs a final theorem statement of the form `C(t) = known_object(t)` or a uniquely characterizing theorem that has not been found yet.

## Award-Track Endgame Hook

- This candidate matches the current hero-case structural signature in reduced variable `t`.
- Award-track target module (Lean scaffold): `proofs/Proofs/HeroCaseFinalIdentity.lean`
- Current state: the module compiles, but `finalIdentityStatement` is still a placeholder.
- Only replace the placeholder after a concrete closed form is identified.
