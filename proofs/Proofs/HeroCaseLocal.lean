import Mathlib
import Proofs.GeneralizedCF

open Polynomial

namespace Proofs
namespace HeroCase

abbrev ZPoly := Polynomial Int

def heroTargetNumerator (n : Nat) (t : Int) : Int :=
  t ^ n + t ^ (2 * n)

def heroTargetDenominator (n : Nat) (t : Int) : Int :=
  1 + t ^ n

def heroInitialDenominator (_t : Int) : Int :=
  1

def rrDirectFirstNumerator (t : Int) : Int :=
  t

def cubicDirectSecondNumerator (t : Int) : Int :=
  t ^ 4 - t ^ 3 + t ^ 2 - t

def cubicOddInitialDenominator (t : Int) : Int :=
  t ^ 2 + t + 1

def cubicEvenFirstDenominator (t : Int) : Int :=
  t ^ 4 + t ^ 2 + 1

noncomputable section

def heroTargetNumeratorPoly (n : Nat) : ZPoly :=
  X ^ n + X ^ (2 * n)

def heroTargetDenominatorPoly (n : Nat) : ZPoly :=
  1 + X ^ n

def heroInitialDenominatorPoly : ZPoly :=
  1

def rrDirectFirstNumeratorPoly : ZPoly :=
  X

def cubicDirectSecondNumeratorPoly : ZPoly :=
  X ^ 4 - X ^ 3 + X ^ 2 - X

def cubicOddInitialDenominatorPoly : ZPoly :=
  X ^ 2 + X + 1

def cubicEvenFirstDenominatorPoly : ZPoly :=
  X ^ 4 + X ^ 2 + 1

def heroData : GCFData ZPoly where
  b0 := 1
  a := fun k => heroTargetNumeratorPoly (k + 1)
  b := fun k => heroTargetDenominatorPoly (k + 1)

def rrReciprocalData : GCFData ZPoly where
  b0 := 1
  a := fun k => X ^ (k + 1)
  b := fun _ => 1

def cubicReciprocalData : GCFData ZPoly where
  b0 := 1
  a := fun k => X ^ (k + 1) + X ^ (2 * (k + 1))
  b := fun _ => 1

theorem heroTargetNumerator_one (t : Int) :
    heroTargetNumerator 1 t = t + t ^ 2 := by
  simp [heroTargetNumerator]

theorem heroTargetNumerator_two (t : Int) :
    heroTargetNumerator 2 t = t ^ 2 + t ^ 4 := by
  simp [heroTargetNumerator]

theorem heroTargetDenominator_one (t : Int) :
    heroTargetDenominator 1 t = 1 + t := by
  simp [heroTargetDenominator]

theorem heroData_stage0_a :
    heroData.a 0 = heroTargetNumeratorPoly 1 := rfl

theorem heroData_stage1_a :
    heroData.a 1 = heroTargetNumeratorPoly 2 := rfl

theorem heroData_stage0_b :
    heroData.b 0 = heroTargetDenominatorPoly 1 := rfl

theorem heroFirstConvergentNumerator :
    continuantNum heroData 1 = 1 + 2 * X + X ^ 2 := by
  simp [heroData, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
  ring_nf

theorem heroFirstConvergentDenominator :
    continuantDen heroData 1 = 1 + X := by
  simp [heroData, heroTargetDenominatorPoly]

theorem heroSecondConvergentNumerator :
    continuantNum heroData 2 = 1 + 2 * X + 3 * X ^ 2 + 2 * X ^ 3 + 2 * X ^ 4 := by
  simp [continuantNum, heroData, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
  ring_nf

theorem heroSecondConvergentDenominator :
    continuantDen heroData 2 = 1 + X + 2 * X ^ 2 + X ^ 3 + X ^ 4 := by
  simp [continuantDen, heroData, heroTargetDenominatorPoly]
  simp [heroTargetNumeratorPoly]
  ring_nf

theorem rrDirectFirstNumerator_mismatch_at_one :
    rrDirectFirstNumerator 1 ≠ heroTargetNumerator 1 1 := by
  norm_num [rrDirectFirstNumerator, heroTargetNumerator]

theorem rrDirectFirstNumerator_mismatch :
    rrDirectFirstNumerator ≠ fun t => heroTargetNumerator 1 t := by
  intro h
  exact rrDirectFirstNumerator_mismatch_at_one (by simpa using congrFun h 1)

theorem cubicDirectSecondNumerator_mismatch_at_one :
    cubicDirectSecondNumerator 1 ≠ heroTargetNumerator 2 1 := by
  norm_num [cubicDirectSecondNumerator, heroTargetNumerator]

theorem cubicDirectSecondNumerator_mismatch :
    cubicDirectSecondNumerator ≠ fun t => heroTargetNumerator 2 t := by
  intro h
  exact cubicDirectSecondNumerator_mismatch_at_one (by simpa using congrFun h 1)

theorem cubicOddInitialDenominator_mismatch_at_one :
    cubicOddInitialDenominator 1 ≠ heroInitialDenominator 1 := by
  norm_num [cubicOddInitialDenominator, heroInitialDenominator]

theorem cubicOddInitialDenominator_mismatch :
    cubicOddInitialDenominator ≠ heroInitialDenominator := by
  intro h
  exact cubicOddInitialDenominator_mismatch_at_one (by simpa using congrFun h 1)

theorem cubicEvenFirstDenominator_mismatch_at_one :
    cubicEvenFirstDenominator 1 ≠ heroTargetDenominator 1 1 := by
  norm_num [cubicEvenFirstDenominator, heroTargetDenominator]

theorem cubicEvenFirstDenominator_mismatch :
    cubicEvenFirstDenominator ≠ fun t => heroTargetDenominator 1 t := by
  intro h
  exact cubicEvenFirstDenominator_mismatch_at_one (by simpa using congrFun h 1)

theorem rrDirectFirstNumeratorPoly_mismatch :
    rrDirectFirstNumeratorPoly ≠ heroTargetNumeratorPoly 1 := by
  intro h
  have hEval := congrArg (fun p : ZPoly => Polynomial.eval 1 p) h
  norm_num [rrDirectFirstNumeratorPoly, heroTargetNumeratorPoly] at hEval

theorem cubicDirectSecondNumeratorPoly_mismatch :
    cubicDirectSecondNumeratorPoly ≠ heroTargetNumeratorPoly 2 := by
  intro h
  have hEval := congrArg (fun p : ZPoly => Polynomial.eval 1 p) h
  norm_num [cubicDirectSecondNumeratorPoly, heroTargetNumeratorPoly] at hEval

theorem cubicOddInitialDenominatorPoly_mismatch :
    cubicOddInitialDenominatorPoly ≠ heroInitialDenominatorPoly := by
  intro h
  have hEval := congrArg (fun p : ZPoly => Polynomial.eval 1 p) h
  norm_num [cubicOddInitialDenominatorPoly, heroInitialDenominatorPoly] at hEval

theorem cubicEvenFirstDenominatorPoly_mismatch :
    cubicEvenFirstDenominatorPoly ≠ heroTargetDenominatorPoly 1 := by
  intro h
  have hEval := congrArg (fun p : ZPoly => Polynomial.eval 1 p) h
  norm_num [cubicEvenFirstDenominatorPoly, heroTargetDenominatorPoly] at hEval

theorem noDirectRRBauerMuirMatchAtStage1 :
    rrDirectFirstNumeratorPoly ≠ heroData.a 0 := by
  simpa [heroData] using rrDirectFirstNumeratorPoly_mismatch

theorem noDirectCubicBauerMuirMatchAtStage2 :
    cubicDirectSecondNumeratorPoly ≠ heroData.a 1 := by
  simpa [heroData] using cubicDirectSecondNumeratorPoly_mismatch

theorem noSimpleCubicOddContraction :
    cubicOddInitialDenominatorPoly ≠ heroData.b0 := by
  simpa [heroData] using cubicOddInitialDenominatorPoly_mismatch

theorem noSimpleCubicEvenContraction :
    cubicEvenFirstDenominatorPoly ≠ heroData.b 0 := by
  simpa [heroData] using cubicEvenFirstDenominatorPoly_mismatch

end

end HeroCase
end Proofs
