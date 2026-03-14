# Formalization Prep: `cb60fd71d1d7`

## Snapshot

- Candidate id: `cb60fd71d1d7`
- Closest benchmark: `rogers_ramanujan_q3_normalized` (7 shared digits)
- Candidate template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=1;ner=6;nek=6;dc=1;ds=1;dr=3;dk=3`
- Benchmark template: `qcf:top=1;base=1;ns=1;nr=3;nk=3;nes=0;ner=0;nek=0;dc=1;ds=0;dr=0;dk=0`

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

## Exact Lemma Candidates

### Direct 1-Step Bauer-Muir Obstructions

- RR source: `w0 = 0`, `w1 = t`, so the first transformed numerator is `t` instead of `t^2 + t`.
- Cubic source: `w0 = 0`, `w1 = t`, `w2 = t^2`, so the second transformed numerator is `t^4 - t^3 + t^2 - t` instead of `t^4 + t^2`.

### Simple Cubic Contraction Obstructions

- Odd contraction: initial term is `t^2 + t + 1` instead of the target `1`.
- Even contraction: first numerator does match as `t^2 + t`, but the first denominator is `t^4 + t^2 + 1` instead of `t + 1`.

## Bounded Exact Exclusion Results

- Arithmetic subsequence contractions up to stride `4`: RR hits `0`, cubic hits `0`.
- These are exact statements for the bounded class being checked, but they do not identify a final source theorem.
- Page-43 monomial substitutions in the current `[-3,3]` shift box: `f2/gcf3` hits `0`, `f4/gcf2` hits `0`.
- These are bounded symbolic searches, useful for narrowing the theorem statement but not substitutes for a full origin proof.

## Formalization Order

1. Formalize generalized continued fractions and convergent recurrence for finite truncations.
2. Formalize the direct 1-step Bauer-Muir obstruction lemmas against the reduced target.
3. Formalize odd/even contraction reconstruction and the cubic denominator mismatch lemma.
4. Defer the bounded search exclusions until a final theorem statement makes them clearly necessary.
5. Do not start a full Lean/Coq origin theorem until a unique source family or exact identity is identified.

## Why This Is Still Not A Full Theorem

- The current exact lemmas only rule out nearby transforms and simple contraction sources.
- They do not prove what the candidate *is*.
- A full formal proof needs a final theorem statement of the form `C(t) = known_object(t)` or a uniquely characterizing theorem that has not been found yet.
