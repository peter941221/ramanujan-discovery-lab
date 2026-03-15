# Lean Proofs

This subproject is the first machine-checked proof layer for the Ramanujan
Discovery Lab.

Current scope:

- `Proofs/GeneralizedCF.lean` formalizes finite-truncation continuants and the
  basic convergent recurrence
- `Proofs/HeroCaseLocal.lean` proves local hero-case obstruction lemmas and
  stronger polynomial no-match theorems
- `Proofs/HeroCaseLocal.lean` also proves an exact convergent-factor reduction:
  the `n`th hero-case convergent numerator and denominator both contain the
  common factor `1 + t^n`, and the reduced convergents come from a second,
  simpler continued fraction
- `Proofs/RationalEquivalence.lean` is a first rational-function testbed for the
  reverse equivalence witness: it lifts the stage scales into `RatFunc Rat` and
  checks the first few retransformed stages exactly
- the reverse step back to the original reciprocal object is currently tracked
  as an explicit equivalence transform witness in the Python-generated research
  and formalization notes; because its stage scales are rational functions, a
  fuller Lean treatment will likely need a fraction-field coefficient layer
- `Proofs/Generated/Cb60fd71d1d7.lean` is an auto-generated, fully checked
  candidate-specific proof module driven by the Python `formalize` command
- does not yet formalize a final source theorem

Commands:

```powershell
Set-Location proofs
lake build
lake env lean Proofs/Generated/Cb60fd71d1d7.lean
lake env lean Proofs/RationalEquivalence.lean
```

Primary file:

- `Proofs/GeneralizedCF.lean`
- `Proofs/HeroCaseLocal.lean`
- `Proofs/RationalEquivalence.lean`
- `Proofs/Generated/Cb60fd71d1d7.lean`
