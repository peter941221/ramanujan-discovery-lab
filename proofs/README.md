# Lean Proofs

This subproject is the first machine-checked proof layer for the Ramanujan
Discovery Lab.

Current scope:

- `Proofs/GeneralizedCF.lean` formalizes finite-truncation continuants and the
  basic convergent recurrence
- `Proofs/HeroCaseObjects.lean` defines the canonical hero-case objects used by
  the local obstruction layer and the rational-function equivalence layer
- `Proofs/HeroCaseLocal.lean` proves local hero-case obstruction lemmas and
  stronger polynomial no-match theorems
- `Proofs/HeroCaseLocal.lean` also proves an exact convergent-factor reduction:
  the `n`th hero-case convergent numerator and denominator both contain the
  common factor `1 + t^n`, and the reduced convergents come from a second,
  simpler continued fraction
- `Proofs/HeroCaseHeineCor2cf.lean` mirrors the decisive low-stage Heine `cor2cf`
  obstruction formulas in the forced `a = 0` lane and proves exact odd/even
  branch mismatches against the hero target
- `Proofs/HeroCasePage43.lean` formalizes a stage-2 exclusion theorem for the
  page-43 monomial-substitution families (`f2`, `f4`) in Laurent polynomials
- `Proofs/HeroCasePage43Equivalence.lean` formalizes the zero-shift
  `n`-dependent equivalence obstructions for the page-43 `f2/gcf3` and
  `f4/gcf2` lanes
- `Proofs/HeroCaseSubsequence.lean` mirrors the bounded arithmetic subsequence
  contraction scan (stride ≤ 4) as a computation-checked theorem over exact
  rational sample points
- `Proofs/HeroCaseSubsequenceExact.lean` strengthens the subsequence exclusion
  for the RR and cubic source lanes without a stride bound, by proving forced
  low-degree coefficient gaps in the source convergents
- `Proofs/HeroCaseFinalIdentity.lean` is the award-track target scaffold where the
  final closed-form identity (once identified) will be stated and proved
- `Proofs/HeroCaseBauerMuir.lean` mirrors the bounded 1/2/3-step Bauer–Muir scan
  (fixed low-complexity `wₙ` templates) as a computation-checked theorem over
  exact rational sample points
- `Proofs/RationalEquivalence.lean` proves the reverse equivalence witness over
  `RatFunc Rat`: the retransformed coefficient data agrees with the canonical
  hero-case data at all stages, and convergents agree for every finite truncation
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
- `Proofs/HeroCaseObjects.lean`
- `Proofs/HeroCaseLocal.lean`
- `Proofs/HeroCaseHeineCor2cf.lean`
- `Proofs/HeroCasePage43.lean`
- `Proofs/HeroCasePage43Equivalence.lean`
- `Proofs/HeroCaseSubsequence.lean`
- `Proofs/HeroCaseSubsequenceExact.lean`
- `Proofs/HeroCaseFinalIdentity.lean`
- `Proofs/HeroCaseBauerMuir.lean`
- `Proofs/RationalEquivalence.lean`
- `Proofs/Generated/Cb60fd71d1d7.lean`
