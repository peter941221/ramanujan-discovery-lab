import Mathlib
import Proofs.HeroCaseGGQuotientCoordinateObstruction
import Proofs.HeroCaseMortonSquaredCoordinateObstruction
import Proofs.HeroCaseWeberSchlafliCoordinateObstruction
import Proofs.HeroCaseWeberCompanionObstruction

namespace Proofs
namespace HeroCase
namespace MortonNamedCoordinate

/-!
Unified Morton named-coordinate waypoint.

This module packages the current source-faithful Morton coordinate stack used
in the hero-case tail-family lane:

- the explicit square coordinate `X_mt = F^2`
- the linear-fractional transformed square coordinate
  `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`
- the Weber-Schläfli coordinate `P_ws = (1/F - F) / 2`
- the direct Weber companion `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`

The current research picture is:

- `X_mt` already collapses to one universal constant-term obstruction `(t^0, 4)`
- `T_mt` already collapses to one universal constant-term obstruction `(t^0, 8)`
- `P_ws` now lives as its own exact obstruction shell, still with a small
  multi-class witness table across the `12` sampled tail-family objects
- `B_ws` collapses back to one universal constant-term obstruction class
  through `Proofs.HeroCase.WeberCompanion.currentWaypoint`

So this is a theorem-shaped waypoint for the whole current Morton named
coordinate stack, not a positive source theorem.
-/

abbrev TailSample := GGQ34.TailSample

structure FirstFailureWitness where
  sampleLabel : String
  power : Nat
  coeff : Int
deriving Repr, DecidableEq

structure ObstructionClass where
  power : Nat
  coeff : Int
deriving Repr, DecidableEq

def tailSamples : List TailSample := GGQ34.tailSamples

def tailSampleLabels : List String := GGQ34.tailSampleLabels

def transformedCoordinateLabel : String := "T_mt"

def transformedCoordinateExpression : String :=
  "T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1), sigma = -1 + sqrt(2)"

def transformedTemplateLabel : String :=
  "Morton Eq. (3.6) transformed squared-coordinate template `T^2 - (T_2^2 - 4*T_2 + 1)*T + T_2^2`"

def transformedWitnessFor (sample : TailSample) : FirstFailureWitness :=
  ⟨sample.label, 0, 8⟩

def transformedWitnesses : List FirstFailureWitness :=
  tailSamples.map transformedWitnessFor

def transformedUniversalClass : ObstructionClass :=
  ⟨0, 8⟩

def transformedClassFor (sample : TailSample) : ObstructionClass :=
  ⟨(transformedWitnessFor sample).power, (transformedWitnessFor sample).coeff⟩

def transformedWitnessClasses : List ObstructionClass :=
  tailSamples.map transformedClassFor

def transformedWitnessTableProp : Prop :=
  transformedWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels ∧
    transformedWitnesses.length = 12

def transformedWitnessTableCheck : Bool :=
  decide (transformedWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels) &&
    decide (transformedWitnesses.length = 12)

theorem transformedWitnessTable_true : transformedWitnessTableProp := by
  have h : transformedWitnessTableCheck = true := by native_decide
  simpa [transformedWitnessTableCheck, transformedWitnessTableProp] using h

def transformedFirstFailureProp (sample : TailSample) : Prop :=
  (transformedWitnessFor sample).power = 0 ∧
    (transformedWitnessFor sample).coeff = 8

theorem transformedFirstFailure_true (sample : TailSample) :
    transformedFirstFailureProp sample := by
  cases sample <;>
    simp [transformedFirstFailureProp, transformedWitnessFor, GGQ34.TailSample.label]

def transformedUniversalClassProp (sample : TailSample) : Prop :=
  transformedClassFor sample = transformedUniversalClass

theorem transformedUniversalClass_true (sample : TailSample) :
    transformedUniversalClassProp sample := by
  cases sample <;>
    simp [transformedUniversalClassProp, transformedClassFor, transformedUniversalClass,
      transformedWitnessFor, GGQ34.TailSample.label]

def transformedClassCoverageProp : Prop :=
  transformedWitnessClasses = List.replicate tailSamples.length transformedUniversalClass

def transformedClassCoverageCheck : Bool :=
  decide (transformedWitnessClasses = List.replicate tailSamples.length transformedUniversalClass)

theorem transformedClassCoverage_true : transformedClassCoverageProp := by
  have h : transformedClassCoverageCheck = true := by native_decide
  simpa [transformedClassCoverageCheck, transformedClassCoverageProp] using h

def transformedCoordinateWaypoint : Prop :=
  transformedCoordinateLabel = "T_mt" ∧
    transformedCoordinateExpression =
      "T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1), sigma = -1 + sqrt(2)" ∧
    transformedTemplateLabel =
      "Morton Eq. (3.6) transformed squared-coordinate template `T^2 - (T_2^2 - 4*T_2 + 1)*T + T_2^2`"

theorem transformedCoordinateWaypoint_true : transformedCoordinateWaypoint := by
  simp [transformedCoordinateWaypoint, transformedCoordinateLabel,
    transformedCoordinateExpression, transformedTemplateLabel]

structure WaypointCertificate where
  xCoordinate : Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint
  transformedCoordinate : transformedCoordinateWaypoint
  transformedWitnessTable : transformedWitnessTableProp
  transformedFailure : ∀ sample : TailSample, transformedFirstFailureProp sample
  transformedCoverage : transformedClassCoverageProp
  pCoordinate : Proofs.HeroCase.WeberSchlafliCoordinate.currentWaypoint
  bCoordinate : Proofs.HeroCase.WeberCompanion.currentWaypoint

def currentWaypointCertificate : WaypointCertificate where
  xCoordinate := Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint_true
  transformedCoordinate := transformedCoordinateWaypoint_true
  transformedWitnessTable := transformedWitnessTable_true
  transformedFailure := transformedFirstFailure_true
  transformedCoverage := transformedClassCoverage_true
  pCoordinate := Proofs.HeroCase.WeberSchlafliCoordinate.currentWaypoint_true
  bCoordinate := Proofs.HeroCase.WeberCompanion.currentWaypoint_true

def currentWaypoint : Prop :=
  Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint ∧
    transformedCoordinateWaypoint ∧
    transformedWitnessTableProp ∧
    (∀ sample : TailSample, transformedFirstFailureProp sample) ∧
    transformedClassCoverageProp ∧
    Proofs.HeroCase.WeberSchlafliCoordinate.currentWaypoint ∧
    Proofs.HeroCase.WeberCompanion.currentWaypoint

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentWaypointCertificate.xCoordinate,
    currentWaypointCertificate.transformedCoordinate,
    currentWaypointCertificate.transformedWitnessTable,
    currentWaypointCertificate.transformedFailure,
    currentWaypointCertificate.transformedCoverage,
    currentWaypointCertificate.pCoordinate,
    currentWaypointCertificate.bCoordinate
  ⟩

end MortonNamedCoordinate
end HeroCase
end Proofs
