import Mathlib
import Mathlib.FieldTheory.RatFunc.AsPolynomial
import Proofs.GeneralizedCF

open Polynomial

namespace Proofs
namespace HeroCase

/-!
Core objects for the hero case `cb60fd71d1d7` in the step-reduced variable `t`.

This file is intentionally definition-heavy and proof-light: it provides the
canonical Lean definitions that other modules (local obstructions, rational
equivalence, etc.) build on.
-/

abbrev ZPoly := Polynomial Int
abbrev QPoly := Polynomial Rat
abbrev QRatFunc := RatFunc Rat

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

def reducedHeroA : Nat → ZPoly
  | 0 => X
  | 1 => X ^ 2
  | k + 2 => X ^ (k + 3) * (1 + X ^ (k + 1))

def reducedHeroB : Nat → ZPoly
  | 0 => 1
  | k + 1 => 1 + X ^ (k + 1)

def reducedHeroData : GCFData ZPoly where
  b0 := 1
  a := reducedHeroA
  b := reducedHeroB

def rrReciprocalData : GCFData ZPoly where
  b0 := 1
  a := fun k => X ^ (k + 1)
  b := fun _ => 1

def cubicReciprocalData : GCFData ZPoly where
  b0 := 1
  a := fun k => X ^ (k + 1) + X ^ (2 * (k + 1))
  b := fun _ => 1

section Generic

variable {R : Type*} [Semiring R]

def heroA (t : R) : Nat → R
  | k => t ^ (k + 1) + t ^ (2 * (k + 1))

def heroB (t : R) : Nat → R
  | k => 1 + t ^ (k + 1)

def heroDataGeneric (t : R) : GCFData R where
  b0 := 1
  a := heroA t
  b := heroB t

def reducedHeroAGeneric (t : R) : Nat → R
  | 0 => t
  | 1 => t ^ 2
  | k + 2 => t ^ (k + 3) * (1 + t ^ (k + 1))

def reducedHeroBGeneric (t : R) : Nat → R
  | 0 => 1
  | k + 1 => 1 + t ^ (k + 1)

def reducedHeroDataGeneric (t : R) : GCFData R where
  b0 := 1
  a := reducedHeroAGeneric t
  b := reducedHeroBGeneric t

def rrReciprocalDataGeneric (t : R) : GCFData R where
  b0 := 1
  a := fun k => t ^ (k + 1)
  b := fun _ => 1

def cubicReciprocalDataGeneric (t : R) : GCFData R where
  b0 := 1
  a := fun k => t ^ (k + 1) + t ^ (2 * (k + 1))
  b := fun _ => 1

end Generic

section ReverseScales

variable {F : Type*} [Field F]

def reverseScaleGeneric (t : F) : Nat → F
  | 0 => 1 + t
  | n + 1 => (1 + t ^ (n + 2)) / (1 + t ^ (n + 1))

def retransformedHeroDataGeneric (t : F) : GCFData F :=
  equivalenceTransform (reducedHeroDataGeneric (t := t)) (reverseScaleGeneric (t := t))

end ReverseScales

def heroDataQ : GCFData QPoly :=
  heroDataGeneric (t := (X : QPoly))

def reducedHeroDataQ : GCFData QPoly :=
  reducedHeroDataGeneric (t := (X : QPoly))

def rrReciprocalDataQ : GCFData QPoly :=
  rrReciprocalDataGeneric (t := (X : QPoly))

def heroDataRatFunc : GCFData QRatFunc :=
  heroDataGeneric (t := RatFunc.X)

def reducedHeroDataRatFunc : GCFData QRatFunc :=
  reducedHeroDataGeneric (t := RatFunc.X)

def rrReciprocalDataRatFunc : GCFData QRatFunc :=
  rrReciprocalDataGeneric (t := RatFunc.X)

def reverseScaleRatFunc : Nat → QRatFunc :=
  reverseScaleGeneric (t := RatFunc.X)

theorem heroData_stage0_a :
    heroData.a 0 = heroTargetNumeratorPoly 1 := rfl

theorem heroData_stage1_a :
    heroData.a 1 = heroTargetNumeratorPoly 2 := rfl

theorem heroData_stage0_b :
    heroData.b 0 = heroTargetDenominatorPoly 1 := rfl

theorem reducedHeroData_stage0_a :
    reducedHeroData.a 0 = X := rfl

theorem reducedHeroData_stage1_a :
    reducedHeroData.a 1 = X ^ 2 := rfl

end

end HeroCase
end Proofs
