import Proofs.HeroCaseObjects

open Polynomial

namespace Proofs
namespace HeroCase

theorem heroTargetNumerator_one (t : Int) :
    heroTargetNumerator 1 t = t + t ^ 2 := by
  simp [heroTargetNumerator]

theorem heroTargetNumerator_two (t : Int) :
    heroTargetNumerator 2 t = t ^ 2 + t ^ 4 := by
  simp [heroTargetNumerator]

theorem heroTargetDenominator_one (t : Int) :
    heroTargetDenominator 1 t = 1 + t := by
  simp [heroTargetDenominator]

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

theorem reducedHeroData_stage2_a :
    reducedHeroData.a 2 = X ^ 3 * (1 + X) := by
  simp [reducedHeroData, reducedHeroA]

theorem reducedHeroData_stage0_b :
    reducedHeroData.b 0 = 1 := rfl

theorem reducedHeroData_stage1_b :
    reducedHeroData.b 1 = 1 + X := by
  simp [reducedHeroData, reducedHeroB]

theorem reducedHeroData_stage2_b :
    reducedHeroData.b 2 = 1 + X ^ 2 := by
  simp [reducedHeroData, reducedHeroB]

theorem reducedHeroFirstConvergentNumerator :
    continuantNum reducedHeroData 1 = 1 + X := by
  simp [reducedHeroData, reducedHeroA, reducedHeroB]

theorem reducedHeroFirstConvergentDenominator :
    continuantDen reducedHeroData 1 = 1 := by
  simp [reducedHeroData, reducedHeroB]

theorem reducedHeroSecondConvergentNumerator :
    continuantNum reducedHeroData 2 = 1 + 2 * X + 2 * X ^ 2 := by
  simp [continuantNum, reducedHeroData, reducedHeroA, reducedHeroB]
  ring_nf

theorem reducedHeroSecondConvergentDenominator :
    continuantDen reducedHeroData 2 = 1 + X + X ^ 2 := by
  simp [continuantDen, reducedHeroData, reducedHeroA, reducedHeroB]

theorem heroConvergentFactor_pair (n : Nat) :
    (continuantNum heroData (n + 1) = (1 + X ^ (n + 1)) * continuantNum reducedHeroData (n + 1)) ∧
      (continuantDen heroData (n + 1) = (1 + X ^ (n + 1)) * continuantDen reducedHeroData (n + 1)) := by
  induction' n using Nat.twoStepInduction with n ih0 ih1
  · constructor
    · simp [heroData, reducedHeroData, reducedHeroA, reducedHeroB, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
      ring_nf
    · simp [heroData, reducedHeroData, reducedHeroB, heroTargetDenominatorPoly]
  · constructor
    · simp [continuantNum, heroData, reducedHeroData, reducedHeroA, reducedHeroB, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
      ring_nf
    · simp [continuantDen, heroData, reducedHeroData, reducedHeroA, reducedHeroB, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
      ring_nf
  · rcases ih0 with ⟨ihNum0, ihDen0⟩
    rcases ih1 with ⟨ihNum1, ihDen1⟩
    constructor
    · calc
        continuantNum heroData (n + 3)
            = heroData.b (n + 2) * continuantNum heroData (n + 2) + heroData.a (n + 2) * continuantNum heroData (n + 1) := by
                simp [continuantNum]
        _ = (1 + X ^ (n + 3)) *
              ((1 + X ^ (n + 2)) * continuantNum reducedHeroData (n + 2) +
                (X ^ (n + 3) * (1 + X ^ (n + 1))) * continuantNum reducedHeroData (n + 1)) := by
              rw [ihNum1, ihNum0]
              simp [heroData, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
              ring_nf
        _ = (1 + X ^ (n + 3)) * continuantNum reducedHeroData (n + 3) := by
              simp [continuantNum, reducedHeroData, reducedHeroA, reducedHeroB]
    · calc
        continuantDen heroData (n + 3)
            = heroData.b (n + 2) * continuantDen heroData (n + 2) + heroData.a (n + 2) * continuantDen heroData (n + 1) := by
                simp [continuantDen]
        _ = (1 + X ^ (n + 3)) *
              ((1 + X ^ (n + 2)) * continuantDen reducedHeroData (n + 2) +
                (X ^ (n + 3) * (1 + X ^ (n + 1))) * continuantDen reducedHeroData (n + 1)) := by
              rw [ihDen1, ihDen0]
              simp [heroData, heroTargetNumeratorPoly, heroTargetDenominatorPoly]
              ring_nf
        _ = (1 + X ^ (n + 3)) * continuantDen reducedHeroData (n + 3) := by
              simp [continuantDen, reducedHeroData, reducedHeroA, reducedHeroB]

theorem heroConvergentNumerator_factor (n : Nat) :
    continuantNum heroData (n + 1) = (1 + X ^ (n + 1)) * continuantNum reducedHeroData (n + 1) :=
  (heroConvergentFactor_pair n).1

theorem heroConvergentDenominator_factor (n : Nat) :
    continuantDen heroData (n + 1) = (1 + X ^ (n + 1)) * continuantDen reducedHeroData (n + 1) :=
  (heroConvergentFactor_pair n).2

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

end HeroCase
end Proofs
