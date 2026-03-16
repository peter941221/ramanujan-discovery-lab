import Proofs.HeroCaseObjects

open Polynomial

namespace Proofs
namespace HeroCase
namespace HeineCor2cf

/-!
Exact low-stage obstruction lemmas for the nearby Heine `cor2cf` lane after the
`a = 0` specialization forced by the hero-case initial term.

The Python research layer reconstructs the odd/even contraction formulas exactly.
This Lean module records the decisive polynomial identities behind that output:

* odd branch initial term: `1 + lambda * t`
* even branch first numerator: `lambda * t`
* odd-of-even branch initial-term gap: `lambda * t`
* even-of-even branch first numerator:
  `lambda * t * (1 + b * t^2 + lambda * t^3 + lambda * t^4)`

Those formulas are enough to prove exact low-stage mismatches against the hero
target `t + t^2` / `1`.
-/

abbrev ZPoly := Proofs.HeroCase.ZPoly

noncomputable section

def aZeroOddInitialDenominatorPoly (lam : Int) : ZPoly :=
  1 + C lam * X

def aZeroEvenFirstNumeratorPoly (lam : Int) : ZPoly :=
  C lam * X

def aZeroEvenOddInitialNumeratorPoly (b lam : Int) : ZPoly :=
  1 + C b * X + C lam * X + C lam * X ^ 2

def aZeroEvenOddInitialDenominatorPoly (b lam : Int) : ZPoly :=
  1 + C b * X + C lam * X ^ 2

def aZeroEvenEvenFirstNumeratorPoly (b lam : Int) : ZPoly :=
  C lam * X * (1 + C b * X ^ 2 + C lam * X ^ 3 + C lam * X ^ 4)

private lemma aZeroEvenFirstNumerator_coeff_two (lam : Int) :
    (aZeroEvenFirstNumeratorPoly lam).coeff 2 = 0 := by
  rw [aZeroEvenFirstNumeratorPoly, Polynomial.coeff_mul_X]
  simpa using (Polynomial.coeff_C (a := lam) (n := 1))

private lemma c_mul_x_pow_coeff_eq_zero_of_ne (a : Int) {k n : Nat} (h : n ≠ k) :
    ((C a) * X ^ k : ZPoly).coeff n = 0 := by
  rw [Polynomial.coeff_C_mul]
  simp [Polynomial.coeff_X_pow, h]

private lemma heroFirstNumerator_coeff_two :
    (heroTargetNumeratorPoly 1).coeff 2 = 1 := by
  simp [heroTargetNumeratorPoly, Polynomial.coeff_X, Polynomial.coeff_X_pow]

private lemma aZeroEvenEvenFirstNumerator_expanded (b lam : Int) :
    aZeroEvenEvenFirstNumeratorPoly b lam =
      C lam * X
        + C (lam * b) * X ^ 3
        + C (lam * lam) * X ^ 4
        + C (lam * lam) * X ^ 5 := by
  simp [aZeroEvenEvenFirstNumeratorPoly]
  ring

private lemma aZeroEvenEvenFirstNumerator_coeff_two (b lam : Int) :
    (aZeroEvenEvenFirstNumeratorPoly b lam).coeff 2 = 0 := by
  have h1 : ((C lam) * X : ZPoly).coeff 2 = 0 := by
    simpa using (c_mul_x_pow_coeff_eq_zero_of_ne lam (k := 1) (n := 2) (by decide))
  have h2 : ((C (lam * b)) * X ^ 3 : ZPoly).coeff 2 = 0 := by
    exact c_mul_x_pow_coeff_eq_zero_of_ne (lam * b) (k := 3) (n := 2) (by decide)
  have h3 : ((C (lam * lam)) * X ^ 4 : ZPoly).coeff 2 = 0 := by
    exact c_mul_x_pow_coeff_eq_zero_of_ne (lam * lam) (k := 4) (n := 2) (by decide)
  have h4 : ((C (lam * lam)) * X ^ 5 : ZPoly).coeff 2 = 0 := by
    exact c_mul_x_pow_coeff_eq_zero_of_ne (lam * lam) (k := 5) (n := 2) (by decide)
  rw [aZeroEvenEvenFirstNumerator_expanded]
  rw [Polynomial.coeff_add, Polynomial.coeff_add, Polynomial.coeff_add, h1, h2, h3, h4]
  norm_num

theorem aZeroOddInitialDenominator_mismatch (lam : Int) (h : lam ≠ 0) :
    aZeroOddInitialDenominatorPoly lam ≠ heroInitialDenominatorPoly := by
  intro hEq
  have hCoeff := congrArg (fun p : ZPoly => p.coeff 1) hEq
  have hLam : lam = 0 := by
    simpa [aZeroOddInitialDenominatorPoly, heroInitialDenominatorPoly] using hCoeff
  exact h hLam

theorem aZeroEvenFirstNumerator_mismatch (lam : Int) :
    aZeroEvenFirstNumeratorPoly lam ≠ heroTargetNumeratorPoly 1 := by
  intro hEq
  have hCoeff := congrArg (fun p : ZPoly => p.coeff 2) hEq
  change (aZeroEvenFirstNumeratorPoly lam).coeff 2 = (heroTargetNumeratorPoly 1).coeff 2 at hCoeff
  rw [aZeroEvenFirstNumerator_coeff_two, heroFirstNumerator_coeff_two] at hCoeff
  have hZero : (0 : Int) = 1 := hCoeff
  norm_num at hZero

theorem aZeroEvenOddInitialDenominator_mismatch (b lam : Int) (h : lam ≠ 0) :
    aZeroEvenOddInitialNumeratorPoly b lam ≠ aZeroEvenOddInitialDenominatorPoly b lam := by
  intro hEq
  have hCoeff := congrArg (fun p : ZPoly => p.coeff 1) hEq
  have hLam : lam = 0 := by
    simpa [aZeroEvenOddInitialNumeratorPoly, aZeroEvenOddInitialDenominatorPoly] using hCoeff
  exact h hLam

theorem aZeroEvenEvenFirstNumerator_mismatch (b lam : Int) :
    aZeroEvenEvenFirstNumeratorPoly b lam ≠ heroTargetNumeratorPoly 1 := by
  intro hEq
  have hCoeff := congrArg (fun p : ZPoly => p.coeff 2) hEq
  change (aZeroEvenEvenFirstNumeratorPoly b lam).coeff 2 = (heroTargetNumeratorPoly 1).coeff 2 at hCoeff
  rw [aZeroEvenEvenFirstNumerator_coeff_two, heroFirstNumerator_coeff_two] at hCoeff
  have hZero : (0 : Int) = 1 := hCoeff
  norm_num at hZero

theorem noSimpleCor2cfOddBranch (lam : Int) (h : lam ≠ 0) :
    aZeroOddInitialDenominatorPoly lam ≠ heroData.b0 := by
  simpa [heroData] using aZeroOddInitialDenominator_mismatch lam h

theorem noSimpleCor2cfEvenBranch (lam : Int) :
    aZeroEvenFirstNumeratorPoly lam ≠ heroData.a 0 := by
  simpa [heroData] using aZeroEvenFirstNumerator_mismatch lam

theorem noSimpleCor2cfOddOfEvenBranch (b lam : Int) (h : lam ≠ 0) :
    aZeroEvenOddInitialNumeratorPoly b lam ≠ aZeroEvenOddInitialDenominatorPoly b lam := by
  exact aZeroEvenOddInitialDenominator_mismatch b lam h

theorem noSimpleCor2cfEvenOfEvenBranch (b lam : Int) :
    aZeroEvenEvenFirstNumeratorPoly b lam ≠ heroData.a 0 := by
  simpa [heroData] using aZeroEvenEvenFirstNumerator_mismatch b lam

end

end HeineCor2cf
end HeroCase
end Proofs
