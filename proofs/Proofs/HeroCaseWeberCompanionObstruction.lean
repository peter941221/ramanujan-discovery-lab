import Mathlib
import Proofs.HeroCaseGGQuotientCoordinateObstruction

namespace Proofs
namespace HeroCase
namespace WeberCompanion

/-!
Weber companion-coordinate obstruction scaffold.

This file packages the current Python-side exact obstruction data for the
direct Weber companion lane

  `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`

where

  `P_ws = (1 / Y - Y) / 2`

for the same `12` sampled tail-family objects used in the `GG` quotient-
coordinate obstruction shell.

The current computational picture is especially rigid:

- every sampled tail object gives `0` exact hits in the first two companion
  templates
- both templates fail immediately at the constant term
- the whole witness table collapses to one universal class

So this module is again a theorem-shaped research shell, not yet a full
modular-function exclusion theorem.
-/

abbrev TailSample := GGQ34.TailSample

structure FirstFailureWitness where
  sampleLabel : String
  power : Nat
  coeff : Int
deriving Repr, DecidableEq

structure CompanionObstructionClass where
  squarePower : Nat
  squareCoeff : Int
  quarticPower : Nat
  quarticCoeff : Int
deriving Repr, DecidableEq

def companionCoordinateLabel : String := "B_ws"

def companionCoordinateExpression : String :=
  "B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)"

def squareTemplateLabel : String :=
  "Morton Weber companion template `B^2 - B_2 - 4`"

def quarticTemplateLabel : String :=
  "Morton Weber companion template `B_2^4 - P^8 - 16*P^4`"

def tailSamples : List TailSample := GGQ34.tailSamples

def tailSampleLabels : List String := GGQ34.tailSampleLabels

def squareWitnessFor (sample : TailSample) : FirstFailureWitness :=
  ⟨sample.label, 0, -2⟩

def quarticWitnessFor (sample : TailSample) : FirstFailureWitness :=
  ⟨sample.label, 0, 16⟩

def squareWitnesses : List FirstFailureWitness :=
  tailSamples.map squareWitnessFor

def quarticWitnesses : List FirstFailureWitness :=
  tailSamples.map quarticWitnessFor

def universalClass : CompanionObstructionClass :=
  ⟨0, -2, 0, 16⟩

def classFor (sample : TailSample) : CompanionObstructionClass :=
  ⟨(squareWitnessFor sample).power, (squareWitnessFor sample).coeff,
    (quarticWitnessFor sample).power, (quarticWitnessFor sample).coeff⟩

def witnessClasses : List CompanionObstructionClass :=
  tailSamples.map classFor

def sampleWitnessProp (sample : TailSample) : Prop :=
  (squareWitnessFor sample).sampleLabel = sample.label ∧
    (quarticWitnessFor sample).sampleLabel = sample.label ∧
    (squareWitnessFor sample).power = 0 ∧
    (squareWitnessFor sample).coeff = -2 ∧
    (quarticWitnessFor sample).power = 0 ∧
    (quarticWitnessFor sample).coeff = 16

theorem sampleWitness_true (sample : TailSample) : sampleWitnessProp sample := by
  cases sample <;>
    simp [sampleWitnessProp, squareWitnessFor, quarticWitnessFor, GGQ34.TailSample.label]

def witnessTableProp : Prop :=
  squareWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels ∧
    quarticWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels ∧
    squareWitnesses.length = 12 ∧
    quarticWitnesses.length = 12

def witnessTableCheck : Bool :=
  decide (squareWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels) &&
    decide (quarticWitnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels) &&
    decide (squareWitnesses.length = 12) &&
    decide (quarticWitnesses.length = 12)

theorem witnessTable_true : witnessTableProp := by
  have h : witnessTableCheck = true := by native_decide
  simpa [witnessTableCheck, witnessTableProp, and_assoc] using h

def squareFirstFailureProp (sample : TailSample) : Prop :=
  (squareWitnessFor sample).power = 0 ∧
    (squareWitnessFor sample).coeff = -2

theorem squareFirstFailure_true (sample : TailSample) : squareFirstFailureProp sample := by
  cases sample <;>
    simp [squareFirstFailureProp, squareWitnessFor, GGQ34.TailSample.label]

def quarticFirstFailureProp (sample : TailSample) : Prop :=
  (quarticWitnessFor sample).power = 0 ∧
    (quarticWitnessFor sample).coeff = 16

theorem quarticFirstFailure_true (sample : TailSample) : quarticFirstFailureProp sample := by
  cases sample <;>
    simp [quarticFirstFailureProp, quarticWitnessFor, GGQ34.TailSample.label]

def pairedFirstFailureProp (sample : TailSample) : Prop :=
  squareFirstFailureProp sample ∧
    quarticFirstFailureProp sample

theorem pairedFirstFailure_true (sample : TailSample) : pairedFirstFailureProp sample := by
  exact ⟨squareFirstFailure_true sample, quarticFirstFailure_true sample⟩

def universalClassProp (sample : TailSample) : Prop :=
  classFor sample = universalClass

theorem universalClass_true (sample : TailSample) : universalClassProp sample := by
  cases sample <;>
    simp [universalClassProp, classFor, universalClass, squareWitnessFor, quarticWitnessFor,
      GGQ34.TailSample.label]

def classCoverageProp : Prop :=
  witnessClasses = List.replicate tailSamples.length universalClass

def classCoverageCheck : Bool :=
  decide (witnessClasses = List.replicate tailSamples.length universalClass)

theorem classCoverage_true : classCoverageProp := by
  have h : classCoverageCheck = true := by native_decide
  simpa [classCoverageCheck, classCoverageProp] using h

def coordinateWaypoint : Prop :=
  companionCoordinateLabel = "B_ws" ∧
    companionCoordinateExpression = "B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)" ∧
    squareTemplateLabel = "Morton Weber companion template `B^2 - B_2 - 4`" ∧
    quarticTemplateLabel = "Morton Weber companion template `B_2^4 - P^8 - 16*P^4`"

theorem coordinateWaypoint_true : coordinateWaypoint := by
  simp [
    coordinateWaypoint,
    companionCoordinateLabel,
    companionCoordinateExpression,
    squareTemplateLabel,
    quarticTemplateLabel
  ]

structure WaypointCertificate where
  coordinate : coordinateWaypoint
  witnessTable : witnessTableProp
  pairFailure : ∀ sample : TailSample, pairedFirstFailureProp sample
  classCoverage : classCoverageProp
  sampleUniversal : ∀ sample : TailSample, universalClassProp sample

def currentWaypointCertificate : WaypointCertificate where
  coordinate := coordinateWaypoint_true
  witnessTable := witnessTable_true
  pairFailure := pairedFirstFailure_true
  classCoverage := classCoverage_true
  sampleUniversal := universalClass_true

def currentWaypoint : Prop :=
  coordinateWaypoint ∧
    witnessTableProp ∧
    (∀ sample : TailSample, pairedFirstFailureProp sample) ∧
    classCoverageProp ∧
    (∀ sample : TailSample, universalClassProp sample)

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentWaypointCertificate.coordinate,
    currentWaypointCertificate.witnessTable,
    currentWaypointCertificate.pairFailure,
    currentWaypointCertificate.classCoverage,
    currentWaypointCertificate.sampleUniversal
  ⟩

end WeberCompanion
end HeroCase
end Proofs
