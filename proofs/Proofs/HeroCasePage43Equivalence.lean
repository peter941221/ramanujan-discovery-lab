import Mathlib
import Proofs.HeroCaseObjects

open Polynomial

namespace Proofs
namespace HeroCase
namespace Page43Equivalence

/-!
Exact zero-shift `n`-dependent equivalence obstructions for the page-43 lanes.

These theorems formalize the current source-family-specific exact lanes produced by
the Python research pipeline:

* `f2 / gcf3` at zero shift
* `f4 / gcf2` at zero shift

For each family, the necessary residual identity is expanded as a polynomial in
`m = t^(n-1)` with coefficients in `Rat[X]`.  Vanishing of all coefficients would
be necessary for an arbitrary `n`-dependent equivalence transformation to match
the hero-case reciprocal.  We prove those coefficient conditions are inconsistent.
-/

abbrev QPoly := Polynomial Rat

noncomputable section

def f2ZeroShiftM1 (lam : Rat) : QPoly :=
  C (lam - 1) * X

def f2ZeroShiftM2 (a b lam : Rat) : QPoly :=
  C (-a) * X ^ 3 + C (-(a * b + a + b)) * X ^ 2 + C (lam - b) * X

def f2ZeroShiftM3 (a b : Rat) : QPoly :=
  C (-(a * a)) * X ^ 4 + C (-(2 * a * b)) * X ^ 3 + C (-(a * b + b * b)) * X ^ 2

def f4ZeroShiftM0 (a : Rat) : QPoly :=
  C a * X

def f4ZeroShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 3 + C (2 * a) * X ^ 2 + C (a + lam - 1) * X

def f4ZeroShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 3 + C (a * b - b) * X ^ 2 + C (lam - b) * X

def f4ZeroShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 2

theorem f2ZeroShiftM3_forces_a_zero (a b : Rat) (h : f2ZeroShiftM3 a b = 0) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have haSq : a * a = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2ZeroShiftM3_forces_b_zero (a b : Rat) (h : f2ZeroShiftM3 a b = 0) :
    b = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have hbSq : b * b = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2ZeroShiftM1_forces_lambda_one (lam : Rat) (h : f2ZeroShiftM1 lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2ZeroShiftM1] at hCoeff
  linarith

theorem f2ZeroShiftM2_specialized :
    f2ZeroShiftM2 0 0 1 = X := by
  simp [f2ZeroShiftM2]

theorem f2ZeroShiftM2_specialized_nonzero :
    f2ZeroShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2ZeroShiftM2] at hCoeff

theorem noZeroShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftM1 lam = 0 ∧ f2ZeroShiftM2 a b lam = 0 ∧ f2ZeroShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  have ha : a = 0 := f2ZeroShiftM3_forces_a_zero a b hM3
  have hb : b = 0 := f2ZeroShiftM3_forces_b_zero a b hM3
  have hlam : lam = 1 := f2ZeroShiftM1_forces_lambda_one lam hM1
  have hSpec : f2ZeroShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f2ZeroShiftM2_specialized_nonzero hSpec

theorem f4ZeroShiftM0_forces_a_zero (a : Rat) (h : f4ZeroShiftM0 a = 0) :
    a = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftM0] at hEval
  exact hEval

theorem f4ZeroShiftM3_forces_b_zero (b : Rat) (h : f4ZeroShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftM3] at hEval
  exact hEval

theorem f4ZeroShiftM1_forces_lambda_one (a lam : Rat)
    (ha : a = 0) (h : f4ZeroShiftM1 a lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4ZeroShiftM1, ha] at hCoeff
  linarith

theorem f4ZeroShiftM2_specialized :
    f4ZeroShiftM2 0 0 1 = X := by
  simp [f4ZeroShiftM2]

theorem f4ZeroShiftM2_specialized_nonzero :
    f4ZeroShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4ZeroShiftM2] at hCoeff

theorem noZeroShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4ZeroShiftM1 a lam = 0 ∧
      f4ZeroShiftM2 a b lam = 0 ∧
      f4ZeroShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4ZeroShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4ZeroShiftM3_forces_b_zero b hM3
  have hlam : lam = 1 := f4ZeroShiftM1_forces_lambda_one a lam ha hM1
  have hSpec : f4ZeroShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f4ZeroShiftM2_specialized_nonzero hSpec

end

end Page43Equivalence
end HeroCase
end Proofs
