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
- `Proofs/Generated/Cb60fd71d1d7.lean` is an auto-generated, fully checked
  candidate-specific proof module driven by the Python `formalize` command
- does not yet formalize a final source theorem

Commands:

```powershell
Set-Location proofs
lake build
lake env lean Proofs/Generated/Cb60fd71d1d7.lean
```

Primary file:

- `Proofs/GeneralizedCF.lean`
- `Proofs/HeroCaseLocal.lean`
- `Proofs/Generated/Cb60fd71d1d7.lean`
