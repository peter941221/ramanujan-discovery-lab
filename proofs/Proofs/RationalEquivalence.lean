import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Proofs.GeneralizedCF

open Polynomial

namespace Proofs
namespace RationalEquivalence

abbrev QPoly := Polynomial Rat
abbrev QRatFunc := RatFunc Rat

noncomputable section

def heroA : Nat → QRatFunc
  | k => RatFunc.X ^ (k + 1) + RatFunc.X ^ (2 * (k + 1))

def heroB : Nat → QRatFunc
  | k => 1 + RatFunc.X ^ (k + 1)

def reducedHeroA : Nat → QRatFunc
  | 0 => RatFunc.X
  | 1 => RatFunc.X ^ 2
  | k + 2 => RatFunc.X ^ (k + 3) * (1 + RatFunc.X ^ (k + 1))

def reducedHeroB : Nat → QRatFunc
  | 0 => 1
  | k + 1 => 1 + RatFunc.X ^ (k + 1)

def reverseScale : Nat → QRatFunc
  | 0 => 1 + RatFunc.X
  | n + 1 => (1 + RatFunc.X ^ (n + 2)) / (1 + RatFunc.X ^ (n + 1))

def heroDataRat : Proofs.GCFData QRatFunc where
  b0 := 1
  a := heroA
  b := heroB

def reducedHeroDataRat : Proofs.GCFData QRatFunc where
  b0 := 1
  a := reducedHeroA
  b := reducedHeroB

def retransformedHeroDataRat : Proofs.GCFData QRatFunc :=
  equivalenceTransform reducedHeroDataRat reverseScale

def reverseTransformedA : Nat → QRatFunc :=
  retransformedHeroDataRat.a

def reverseTransformedB : Nat → QRatFunc :=
  retransformedHeroDataRat.b

theorem gcfData_ext {R : Type*} {g h : Proofs.GCFData R}
    (hb0 : g.b0 = h.b0) (ha : g.a = h.a) (hb : g.b = h.b) : g = h := by
  cases g
  cases h
  cases hb0
  cases ha
  cases hb
  rfl

lemma one_add_X_pow_ne_zero (n : Nat) : (1 + RatFunc.X ^ (n + 1) : QRatFunc) ≠ 0 := by
  have hpoly : ((1 : QPoly) + X ^ (n + 1)) ≠ 0 := by
    intro h
    have hEval := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
    simp at hEval
  simpa only [RatFunc.X, map_add, map_pow, map_one] using (RatFunc.algebraMap_ne_zero hpoly)

theorem reverseScale_mul_denominator (n : Nat) :
    reverseScale (n + 1) * (1 + RatFunc.X ^ (n + 1)) = 1 + RatFunc.X ^ (n + 2) := by
  simp [reverseScale]
  field_simp [one_add_X_pow_ne_zero n]

theorem reverseScale_mul_denominator_stage1 :
    reverseScale 1 * (1 + RatFunc.X) = 1 + RatFunc.X ^ 2 := by
  simpa using reverseScale_mul_denominator 0

theorem reverseScale_mul_denominator_stage2 :
    reverseScale 2 * (1 + RatFunc.X ^ 2) = 1 + RatFunc.X ^ 3 := by
  simpa using reverseScale_mul_denominator 1

theorem reverseScale_stage0 :
    reverseScale 0 = 1 + RatFunc.X := rfl

theorem reverseScale_stage1 :
    reverseScale 1 = (1 + RatFunc.X ^ 2) / (1 + RatFunc.X) := by
  simp [reverseScale]

theorem reverseScale_stage2 :
    reverseScale 2 = (1 + RatFunc.X ^ 3) / (1 + RatFunc.X ^ 2) := by
  simp [reverseScale]

theorem reverseTransform_stage0_a :
    reverseTransformedA 0 = heroA 0 := by
  simp [reverseTransformedA, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroA, heroA]
  ring_nf

theorem reverseTransform_stage1_a :
    reverseTransformedA 1 = heroA 1 := by
  calc
    reverseTransformedA 1
        = reverseScale 1 * (reverseScale 0 * RatFunc.X ^ 2) := by
            simp [reverseTransformedA, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroA, mul_left_comm, mul_comm]
    _ = reverseScale 1 * ((1 + RatFunc.X) * RatFunc.X ^ 2) := by
          simp [reverseScale, mul_left_comm, mul_comm]
    _ = (reverseScale 1 * (1 + RatFunc.X)) * RatFunc.X ^ 2 := by
          ring
    _ = (1 + RatFunc.X ^ 2) * RatFunc.X ^ 2 := by
          rw [reverseScale_mul_denominator_stage1]
    _ = heroA 1 := by
          simp [heroA]
          ring

theorem reverseTransform_stage2_a :
    reverseTransformedA 2 = heroA 2 := by
  calc
    reverseTransformedA 2
        = reverseScale 2 * (reverseScale 1 * (RatFunc.X ^ 3 * (1 + RatFunc.X))) := by
            simp [reverseTransformedA, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroA, mul_left_comm, mul_comm]
    _ = reverseScale 2 * ((reverseScale 1 * (1 + RatFunc.X)) * RatFunc.X ^ 3) := by
          ring
    _ = reverseScale 2 * ((1 + RatFunc.X ^ 2) * RatFunc.X ^ 3) := by
          rw [reverseScale_mul_denominator_stage1]
    _ = (reverseScale 2 * (1 + RatFunc.X ^ 2)) * RatFunc.X ^ 3 := by
          ring
    _ = (1 + RatFunc.X ^ 3) * RatFunc.X ^ 3 := by
          rw [reverseScale_mul_denominator_stage2]
    _ = heroA 2 := by
          simp [heroA]
          ring

theorem reverseTransform_stage0_b :
    reverseTransformedB 0 = heroB 0 := by
  simp [reverseTransformedB, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroB, heroB]

theorem reverseTransform_stage1_b :
    reverseTransformedB 1 = heroB 1 := by
  calc
    reverseTransformedB 1 = reverseScale 1 * (1 + RatFunc.X) := by
      simp [reverseTransformedB, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reducedHeroB]
    _ = 1 + RatFunc.X ^ 2 := by
      rw [reverseScale_mul_denominator_stage1]
    _ = heroB 1 := by
      simp [heroB]

theorem reverseTransform_stage2_b :
    reverseTransformedB 2 = heroB 2 := by
  calc
    reverseTransformedB 2 = reverseScale 2 * (1 + RatFunc.X ^ 2) := by
      simp [reverseTransformedB, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reducedHeroB]
    _ = 1 + RatFunc.X ^ 3 := by
      rw [reverseScale_mul_denominator_stage2]
    _ = heroB 2 := by
      simp [heroB]

theorem reverseTransformA_all (k : Nat) :
    reverseTransformedA k = heroA k := by
  cases k with
  | zero =>
      simpa using reverseTransform_stage0_a
  | succ k =>
      cases k with
      | zero =>
          simpa using reverseTransform_stage1_a
      | succ k =>
          -- General stage: the (1 + X^(k+1)) and (1 + X^(k+2)) factors cancel.
          simp [reverseTransformedA, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroA, heroA]
          field_simp [one_add_X_pow_ne_zero k, one_add_X_pow_ne_zero (k + 1)]
          ring_nf

theorem reverseTransformB_all (k : Nat) :
    reverseTransformedB k = heroB k := by
  cases k with
  | zero =>
      simp [reverseTransformedB, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reverseScale, reducedHeroB, heroB]
  | succ k =>
      simpa [reverseTransformedB, retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, reducedHeroB, heroB] using reverseScale_mul_denominator k

theorem retransformedHeroDataRat_eq :
    retransformedHeroDataRat = heroDataRat := by
  apply gcfData_ext (g := retransformedHeroDataRat) (h := heroDataRat)
  · simp [retransformedHeroDataRat, equivalenceTransform, reducedHeroDataRat, heroDataRat]
  · funext k
    exact reverseTransformA_all k
  · funext k
    exact reverseTransformB_all k

lemma reverseScale_ne_zero : ∀ k : Nat, reverseScale k ≠ 0 := by
  intro k
  cases k with
  | zero =>
      simpa [reverseScale] using one_add_X_pow_ne_zero 0
  | succ k =>
      have hnum : (1 + RatFunc.X ^ (k + 2) : QRatFunc) ≠ 0 := by
        simpa [Nat.add_assoc] using one_add_X_pow_ne_zero (k + 1)
      have hden : (1 + RatFunc.X ^ (k + 1) : QRatFunc) ≠ 0 := by
        simpa using one_add_X_pow_ne_zero k
      simpa [reverseScale] using (div_ne_zero hnum hden)

theorem hero_convergent_eq_reduced_convergent (n : Nat) :
    Proofs.convergent heroDataRat n = Proofs.convergent reducedHeroDataRat n := by
  have h :=
    convergent_equivalenceTransform
      (g := reducedHeroDataRat)
      (r := reverseScale)
      reverseScale_ne_zero
      n
  have h' : Proofs.convergent retransformedHeroDataRat n = Proofs.convergent reducedHeroDataRat n := by
    simpa [retransformedHeroDataRat] using h
  calc
    Proofs.convergent heroDataRat n = Proofs.convergent retransformedHeroDataRat n := by
      simp [retransformedHeroDataRat_eq]
    _ = Proofs.convergent reducedHeroDataRat n := h'

end

end RationalEquivalence
end Proofs
