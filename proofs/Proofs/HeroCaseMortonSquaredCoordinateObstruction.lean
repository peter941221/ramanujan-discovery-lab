import Mathlib
import Proofs.HeroCaseGGQuotientCoordinateObstruction

namespace Proofs
namespace HeroCase
namespace MortonSquaredCoordinate

/-!
Morton square-coordinate obstruction scaffold.

This module packages the current Python-side exact obstruction data for the
explicit square-level coordinate from Akkarapakam--Morton Proposition 3.2 /
Theorem B:

  `X_mt = F^2`

with template

  `X_2^2 - (X^2 - 4*X + 1)*X_2 + X^2 = 0`

on the same `12` sampled tail-family objects used elsewhere in the current
hero-case research trunk.

At the current exact waypoint, every sampled tail object gives `0` hits in this
named square lane and the whole witness table collapses to one universal
constant-term obstruction class `(t^0, 4)`.

So this file is a theorem-shaped obstruction shell, not a positive source
identification theorem.
-/

abbrev TailSample := GGQ34.TailSample

structure FirstFailureWitness where
  sampleLabel : String
  power : Nat
  coeff : Int
deriving Repr, DecidableEq

structure SquaredCoordinateObstructionClass where
  power : Nat
  coeff : Int
deriving Repr, DecidableEq

def coordinateLabel : String := "X_mt"

def coordinateExpression : String :=
  "X_mt = F^2"

def templateLabel : String :=
  "Morton Prop. 3.2 squared-coordinate template `X_2^2 - (X^2 - 4*X + 1)*X_2 + X^2`"

def tailSamples : List TailSample := GGQ34.tailSamples

def tailSampleLabels : List String := GGQ34.tailSampleLabels

def witnessFor (sample : TailSample) : FirstFailureWitness :=
  ⟨sample.label, 0, 4⟩

def witnesses : List FirstFailureWitness :=
  tailSamples.map witnessFor

def universalClass : SquaredCoordinateObstructionClass :=
  ⟨0, 4⟩

def classFor (sample : TailSample) : SquaredCoordinateObstructionClass :=
  ⟨(witnessFor sample).power, (witnessFor sample).coeff⟩

def witnessClasses : List SquaredCoordinateObstructionClass :=
  tailSamples.map classFor

def sampleWitnessProp (sample : TailSample) : Prop :=
  (witnessFor sample).sampleLabel = sample.label ∧
    (witnessFor sample).power = 0 ∧
    (witnessFor sample).coeff = 4

theorem sampleWitness_true (sample : TailSample) : sampleWitnessProp sample := by
  cases sample <;>
    simp [sampleWitnessProp, witnessFor, GGQ34.TailSample.label]

def witnessTableProp : Prop :=
  witnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels ∧
    witnesses.length = 12

def witnessTableCheck : Bool :=
  decide (witnesses.map FirstFailureWitness.sampleLabel = tailSampleLabels) &&
    decide (witnesses.length = 12)

theorem witnessTable_true : witnessTableProp := by
  have h : witnessTableCheck = true := by native_decide
  simpa [witnessTableCheck, witnessTableProp] using h

def firstFailureProp (sample : TailSample) : Prop :=
  (witnessFor sample).power = 0 ∧
    (witnessFor sample).coeff = 4

theorem firstFailure_true (sample : TailSample) : firstFailureProp sample := by
  cases sample <;>
    simp [firstFailureProp, witnessFor, GGQ34.TailSample.label]

def universalClassProp (sample : TailSample) : Prop :=
  classFor sample = universalClass

theorem universalClass_true (sample : TailSample) : universalClassProp sample := by
  cases sample <;>
    simp [universalClassProp, classFor, universalClass, witnessFor, GGQ34.TailSample.label]

def classCoverageProp : Prop :=
  witnessClasses = List.replicate tailSamples.length universalClass

def classCoverageCheck : Bool :=
  decide (witnessClasses = List.replicate tailSamples.length universalClass)

theorem classCoverage_true : classCoverageProp := by
  have h : classCoverageCheck = true := by native_decide
  simpa [classCoverageCheck, classCoverageProp] using h

def coordinateWaypoint : Prop :=
  coordinateLabel = "X_mt" ∧
    coordinateExpression = "X_mt = F^2" ∧
    templateLabel =
      "Morton Prop. 3.2 squared-coordinate template `X_2^2 - (X^2 - 4*X + 1)*X_2 + X^2`"

theorem coordinateWaypoint_true : coordinateWaypoint := by
  simp [coordinateWaypoint, coordinateLabel, coordinateExpression, templateLabel]

structure WaypointCertificate where
  coordinate : coordinateWaypoint
  witnessTable : witnessTableProp
  sampleFailure : ∀ sample : TailSample, firstFailureProp sample
  classCoverage : classCoverageProp
  sampleUniversal : ∀ sample : TailSample, universalClassProp sample

def currentWaypointCertificate : WaypointCertificate where
  coordinate := coordinateWaypoint_true
  witnessTable := witnessTable_true
  sampleFailure := firstFailure_true
  classCoverage := classCoverage_true
  sampleUniversal := universalClass_true

def currentWaypoint : Prop :=
  coordinateWaypoint ∧
    witnessTableProp ∧
    (∀ sample : TailSample, firstFailureProp sample) ∧
    classCoverageProp ∧
    (∀ sample : TailSample, universalClassProp sample)

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentWaypointCertificate.coordinate,
    currentWaypointCertificate.witnessTable,
    currentWaypointCertificate.sampleFailure,
    currentWaypointCertificate.classCoverage,
    currentWaypointCertificate.sampleUniversal
  ⟩

end MortonSquaredCoordinate
end HeroCase
end Proofs
