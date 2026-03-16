import Proofs.HeroCaseObjects

namespace Proofs
namespace HeroCase

/-!
Exact arithmetic-subsequence obstruction for the two nearest source lanes:

- reduced RR reciprocal
- reduced cubic reciprocal

This file removes the old `stride <= 4 / sample points / bounded stages` limitations
for these two specific sources by proving forced low-degree coefficient gaps in the
source convergents.
-/

namespace SubsequenceExact

open Proofs
open Polynomial

abbrev ZPoly := Proofs.HeroCase.ZPoly

noncomputable section

def rrDelta (n : Nat) : ZPoly :=
  continuantNum rrReciprocalData n - continuantDen rrReciprocalData n

def rrStepGap (n : Nat) : ZPoly :=
  continuantNum rrReciprocalData n - (1 + X) * continuantDen rrReciprocalData n

def cubicDelta (n : Nat) : ZPoly :=
  continuantNum cubicReciprocalData n - continuantDen cubicReciprocalData n

def cubicStepGap (n : Nat) : ZPoly :=
  continuantNum cubicReciprocalData n - (1 + X) * continuantDen cubicReciprocalData n

private lemma coeff_X_pow_mul_eq_zero_of_lt (p : ZPoly) (d k : Nat) (h : d < k) :
    (X ^ k * p).coeff d = 0 := by
  rw [coeff_X_pow_mul']
  have hk : ¬ k ≤ d := not_le_of_gt h
  simp [hk]

private lemma rrDelta_rec (n : Nat) :
    rrDelta (n + 2) = rrDelta (n + 1) + X ^ (n + 2) * rrDelta n := by
  -- Unfold the delta and expand the `(n+2)` continuants once.
  -- We avoid `rw` here because `dsimp` can already unfold the definitional equations.
  dsimp [rrDelta]
  simp [rrReciprocalData]
  ring_nf

private lemma rrStepGap_rec (n : Nat) :
    rrStepGap (n + 4) = rrStepGap (n + 3) + X ^ (n + 4) * rrStepGap (n + 2) := by
  dsimp [rrStepGap]
  simp [rrReciprocalData]
  ring_nf

private lemma cubicDelta_rec (n : Nat) :
    cubicDelta (n + 2) =
      cubicDelta (n + 1) + (X ^ (n + 2) + X ^ (2 * (n + 2))) * cubicDelta n := by
  dsimp [cubicDelta]
  simp [cubicReciprocalData]
  ring_nf

private lemma cubicStepGap_rec (n : Nat) :
    cubicStepGap (n + 4) =
      cubicStepGap (n + 3) + (X ^ (n + 4) + X ^ (2 * (n + 4))) * cubicStepGap (n + 2) := by
  dsimp [cubicStepGap]
  simp [cubicReciprocalData]
  ring_nf

theorem rrDelta_coeff_one (n : Nat) :
    (rrDelta (n + 1)).coeff 1 = 1 := by
  induction n with
  | zero =>
      -- rrDelta 1 = X
      norm_num [rrDelta, rrReciprocalData, continuantNum, continuantDen]
  | succ n ih =>
      have hrec := rrDelta_rec n
      have hz : (X ^ (n + 2) * rrDelta n).coeff 1 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      have : (rrDelta (n + 2)).coeff 1 = (rrDelta (n + 1)).coeff 1 := by
        rw [hrec]
        simp [coeff_add, hz]
      -- `ih` is the statement for `n`
      simpa [Nat.add_assoc] using (this.trans ih)

theorem cubicDelta_coeff_one (n : Nat) :
    (cubicDelta (n + 1)).coeff 1 = 1 := by
  induction n with
  | zero =>
      norm_num [cubicDelta, cubicReciprocalData, continuantNum, continuantDen]
  | succ n ih =>
      have hrec := cubicDelta_rec n
      have hz1 : ((X ^ (n + 2) : ZPoly) * cubicDelta n).coeff 1 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      -- `simp` will normalize `2 * (n + 2)` to `2 * n + 4`, so we state the vanishing in that form.
      have hz2 : ((X ^ (2 * n + 4) : ZPoly) * cubicDelta n).coeff 1 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      have : (cubicDelta (n + 2)).coeff 1 = (cubicDelta (n + 1)).coeff 1 := by
        rw [hrec]
        -- distribute the product before taking coefficients; both pieces vanish at degree 1.
        simp [coeff_add, mul_add, add_mul, hz1, hz2]
      simpa [Nat.add_assoc] using (this.trans ih)

private lemma rrStepGap_coeff_three_base2 :
    (rrStepGap 2).coeff 3 = -1 := by
  have h : rrStepGap 2 = -X ^ 3 := by
    simp [rrStepGap, rrReciprocalData, continuantNum, continuantDen]
    ring_nf
  simp [h, coeff_X_pow]

private lemma rrStepGap_coeff_three_base3 :
    (rrStepGap 3).coeff 3 = -1 := by
  have h : rrStepGap 3 = -X ^ 3 := by
    simp [rrStepGap, rrReciprocalData, continuantNum, continuantDen]
    ring_nf
  simp [h, coeff_X_pow]

theorem rrStepGap_coeff_three (n : Nat) :
    (rrStepGap (n + 2)).coeff 3 = -1 := by
  induction' n using Nat.twoStepInduction with n ih0 ih1
  · simpa using rrStepGap_coeff_three_base2
  · simpa using rrStepGap_coeff_three_base3
  · have hrec := rrStepGap_rec n
    have : (rrStepGap (n + 4)).coeff 3 = (rrStepGap (n + 3)).coeff 3 := by
      rw [hrec]
      -- `3 < n+4`, so the product term cannot affect `coeff 3`.
      have hz : (X ^ (n + 4) * rrStepGap (n + 2)).coeff 3 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      simp [coeff_add, hz]
    -- `ih1` gives the value at `n+3`
    simpa [Nat.add_assoc] using (this.trans ih1)

private lemma cubicStepGap_coeff_three_base2 :
    (cubicStepGap 2).coeff 3 = -1 := by
  have h : cubicStepGap 2 = X ^ 2 - X ^ 3 - X ^ 5 := by
    simp [cubicStepGap, cubicReciprocalData, continuantNum, continuantDen]
    ring_nf
  -- Only the `-X^3` term contributes to `coeff 3`.
  simp [h, sub_eq_add_neg, coeff_add, coeff_X_pow]

private lemma cubicStepGap_coeff_three_base3 :
    (cubicStepGap 3).coeff 3 = -1 := by
  have h : cubicStepGap 3 = X ^ 2 - X ^ 3 + X ^ 8 := by
    simp [cubicStepGap, cubicReciprocalData, continuantNum, continuantDen]
    ring_nf
  simp [h, sub_eq_add_neg, coeff_add, coeff_X_pow]

theorem cubicStepGap_coeff_three (n : Nat) :
    (cubicStepGap (n + 2)).coeff 3 = -1 := by
  induction' n using Nat.twoStepInduction with n ih0 ih1
  · simpa using cubicStepGap_coeff_three_base2
  · simpa using cubicStepGap_coeff_three_base3
  · have hrec := cubicStepGap_rec n
    have : (cubicStepGap (n + 4)).coeff 3 = (cubicStepGap (n + 3)).coeff 3 := by
      rw [hrec]
      have hz1 : ((X ^ (n + 4) : ZPoly) * cubicStepGap (n + 2)).coeff 3 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      -- As above, `simp` normalizes `2 * (n + 4)` to `2 * n + 8`.
      have hz2 : ((X ^ (2 * n + 8) : ZPoly) * cubicStepGap (n + 2)).coeff 3 = 0 := by
        exact coeff_X_pow_mul_eq_zero_of_lt _ _ _ (by omega)
      -- distribute first: (X^k + X^m) * p = X^k*p + X^m*p
      simp [coeff_add, mul_add, add_mul, hz1, hz2]
    simpa [Nat.add_assoc] using (this.trans ih1)

theorem noRRArithmeticSubsequenceContraction
    (stride offset : Nat) (hstride : 2 ≤ stride) (hoff : offset < stride) :
    ¬ (continuantNum rrReciprocalData offset = continuantDen rrReciprocalData offset ∧
        continuantNum rrReciprocalData (offset + stride) =
          (1 + X) * continuantDen rrReciprocalData (offset + stride)) := by
  intro h
  rcases h with ⟨h0, h1⟩
  cases offset with
  | zero =>
      -- stride >= 2, so stride = (stride-2) + 2
      have hs : stride = (stride - 2) + 2 := by omega
      have hCoeff : (rrStepGap stride).coeff 3 = 0 := by
        -- rrStepGap stride = 0 under the assumed equality
        have h1' : continuantNum rrReciprocalData stride = (1 + X) * continuantDen rrReciprocalData stride := by
          simpa [Nat.zero_add] using h1
        have : rrStepGap stride = 0 := by
          simp [rrStepGap, h1']
        simp [this]
      -- But the forced value is -1 for any stride >= 2.
      have hForced : (rrStepGap stride).coeff 3 = -1 := by
        rw [hs]
        simpa using rrStepGap_coeff_three (stride - 2)
      linarith
  | succ n =>
      have hCoeff : (rrDelta (n + 1)).coeff 1 = 0 := by
        have : rrDelta (n + 1) = 0 := by
          simp [rrDelta, h0]
        simp [this]
      have hForced : (rrDelta (n + 1)).coeff 1 = 1 := by
        simpa using rrDelta_coeff_one n
      linarith

theorem noCubicArithmeticSubsequenceContraction
    (stride offset : Nat) (hstride : 2 ≤ stride) (hoff : offset < stride) :
    ¬ (continuantNum cubicReciprocalData offset = continuantDen cubicReciprocalData offset ∧
        continuantNum cubicReciprocalData (offset + stride) =
          (1 + X) * continuantDen cubicReciprocalData (offset + stride)) := by
  intro h
  rcases h with ⟨h0, h1⟩
  cases offset with
  | zero =>
      have hs : stride = (stride - 2) + 2 := by omega
      have hCoeff : (cubicStepGap stride).coeff 3 = 0 := by
        have h1' : continuantNum cubicReciprocalData stride = (1 + X) * continuantDen cubicReciprocalData stride := by
          simpa [Nat.zero_add] using h1
        have : cubicStepGap stride = 0 := by
          simp [cubicStepGap, h1']
        simp [this]
      have hForced : (cubicStepGap stride).coeff 3 = -1 := by
        rw [hs]
        simpa using cubicStepGap_coeff_three (stride - 2)
      linarith
  | succ n =>
      have hCoeff : (cubicDelta (n + 1)).coeff 1 = 0 := by
        have : cubicDelta (n + 1) = 0 := by
          simp [cubicDelta, h0]
        simp [this]
      have hForced : (cubicDelta (n + 1)).coeff 1 = 1 := by
        simpa using cubicDelta_coeff_one n
      linarith

end

end SubsequenceExact

end HeroCase
end Proofs
