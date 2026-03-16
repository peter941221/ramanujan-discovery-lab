import Proofs.HeroCaseObjects

namespace Proofs
namespace HeroCase

/-!
Bounded arithmetic subsequence contraction exclusions.

The Python research pipeline checks whether the step-reduced target reciprocal continued fraction
can be recovered as an arithmetic subsequence contraction of nearby sources (RR reciprocal and
the cubic reciprocal) for `stride ≤ 4`.

Here we formalize a matching *necessary condition* using exact rational sample points:
if a contraction matched the target coefficient data, then the corresponding subsequence of
source convergents would agree with the target convergents at all sample points.

We compute these finite checks and prove that the hit list is empty.
-/

namespace Subsequence

open Proofs

def samples : List Rat :=
  [(1 : Rat) / 10, (1 : Rat) / 7, (1 : Rat) / 5]

def convergentEq (g₁ g₂ : GCFData Rat) (n₁ n₂ : Nat) : Bool :=
  let p₁ := continuantNum g₁ n₁
  let q₁ := continuantDen g₁ n₁
  let p₂ := continuantNum g₂ n₂
  let q₂ := continuantDen g₂ n₂
  decide (q₁ ≠ 0 ∧ q₂ ≠ 0 ∧ p₁ * q₂ = p₂ * q₁)

def subsequenceMatches (gSource gTarget : GCFData Rat) (stride offset stages : Nat) : Bool :=
  if stride < 2 then
    false
  else if offset ≥ stride then
    false
  else
    (List.range (stages + 1)).all (fun m =>
      convergentEq gSource gTarget (offset + stride * m) m)

def strideOffsetPairs (maxStride : Nat) : List (Nat × Nat) :=
  List.flatMap
    (fun stride => (List.range stride).map (fun offset => (stride, offset)))
    ((List.range (maxStride + 1)).drop 2)

def hitsAtSamples (source target : Rat → GCFData Rat) (maxStride stages : Nat) : List (Nat × Nat) :=
  (strideOffsetPairs maxStride).filter (fun p =>
    samples.all (fun t => subsequenceMatches (source t) (target t) p.1 p.2 stages))

def rrHits : List (Nat × Nat) :=
  hitsAtSamples (source := rrReciprocalDataGeneric) (target := heroDataGeneric) (maxStride := 4) (stages := 3)

def cubicHits : List (Nat × Nat) :=
  hitsAtSamples (source := cubicReciprocalDataGeneric) (target := heroDataGeneric) (maxStride := 4) (stages := 3)

theorem rrHits_eq_nil : rrHits = [] := by
  native_decide

theorem cubicHits_eq_nil : cubicHits = [] := by
  native_decide

end Subsequence

end HeroCase
end Proofs
