import Mathlib
import Proofs.HeroCaseGGQuotientCoordinateObstruction

namespace Proofs
namespace HeroCase
namespace WeberSchlafliCoordinate

/-!
Weber-Schlafli coordinate obstruction scaffold.

This file packages the current Python-side exact obstruction data for the first
deeper Weber-Schlafli coordinate

  `P_ws = (1 / Y - Y) / 2`

on the same `12` sampled tail-family objects used elsewhere in the current
hero-case tail-family lane.

Unlike the square-level `X_mt`, transported-square `T_mt`, and companion
`B_ws` lanes, this `P_ws` lane does not collapse to one universal obstruction
class. Instead, the first exact failure of the Morton Weber-Schlafli template

  `P^2*P_2^2 + P^2 - 2*P_2`

falls into a small repeated witness table with seven currently realized classes.

That seven-class table also compresses one step further:

- the coefficient support is exactly `{-1, 3, 8}`
- the coefficient `3` lives on the even power ladder `2, 4, 6, 8, 10`
- the exceptional coefficients `-1` and `8` only occur at power `2`

So this module is a theorem-shaped exact obstruction shell, not yet a positive
source identification theorem.
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

def coordinateLabel : String := "P_ws"

def coordinateExpression : String :=
  "P_ws = (1/F - F) / 2"

def templateLabel : String :=
  "Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`"

def tailSamples : List TailSample := GGQ34.tailSamples

def tailSampleLabels : List String := GGQ34.tailSampleLabels

def witnessFor : TailSample → FirstFailureWitness
  | .u_t2 => ⟨"U_t2", 6, 3⟩
  | .u_t2_g1 => ⟨"U_t2_g1", 2, 3⟩
  | .u_t2_g2 => ⟨"U_t2_g2", 2, -1⟩
  | .u_t2_g3 => ⟨"U_t2_g3", 2, 8⟩
  | .u_t3 => ⟨"U_t3", 8, 3⟩
  | .u_t3_g1 => ⟨"U_t3_g1", 4, 3⟩
  | .u_t3_g2 => ⟨"U_t3_g2", 2, -1⟩
  | .u_t3_g3 => ⟨"U_t3_g3", 2, 3⟩
  | .u_t4 => ⟨"U_t4", 10, 3⟩
  | .u_t4_g1 => ⟨"U_t4_g1", 6, 3⟩
  | .u_t4_g2 => ⟨"U_t4_g2", 2, -1⟩
  | .u_t4_g3 => ⟨"U_t4_g3", 2, 3⟩

def witnesses : List FirstFailureWitness :=
  tailSamples.map witnessFor

def classFor (sample : TailSample) : ObstructionClass :=
  ⟨(witnessFor sample).power, (witnessFor sample).coeff⟩

def witnessClasses : List ObstructionClass :=
  tailSamples.map classFor

def expectedClasses : List ObstructionClass :=
  [
    ⟨2, -1⟩,
    ⟨2, 3⟩,
    ⟨2, 8⟩,
    ⟨4, 3⟩,
    ⟨6, 3⟩,
    ⟨8, 3⟩,
    ⟨10, 3⟩
  ]

def expectedCoefficients : List Int := [-1, 3, 8]

def coefficientSupport : List Int :=
  witnesses.map FirstFailureWitness.coeff

def coeffThreePowerLadder : List Nat := [2, 4, 6, 8, 10]

def sampleWitnessProp (sample : TailSample) : Prop :=
  (witnessFor sample).sampleLabel = sample.label ∧
    (classFor sample).power = (witnessFor sample).power ∧
    (classFor sample).coeff = (witnessFor sample).coeff

theorem sampleWitness_true (sample : TailSample) : sampleWitnessProp sample := by
  cases sample <;>
    simp [sampleWitnessProp, witnessFor, classFor, GGQ34.TailSample.label]

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
  (witnessFor sample).power = (classFor sample).power ∧
    (witnessFor sample).coeff = (classFor sample).coeff

theorem firstFailure_true (sample : TailSample) : firstFailureProp sample := by
  cases sample <;>
    simp [firstFailureProp, witnessFor, classFor]

def classCoverageProp : Prop :=
  witnessClasses.all (fun obstruction => decide (obstruction ∈ expectedClasses)) = true ∧
    expectedClasses.all (fun obstruction => decide (obstruction ∈ witnessClasses)) = true

def classCoverageCheck : Bool :=
  witnessClasses.all (fun obstruction => decide (obstruction ∈ expectedClasses)) &&
    expectedClasses.all (fun obstruction => decide (obstruction ∈ witnessClasses))

theorem classCoverage_true : classCoverageProp := by
  have h : classCoverageCheck = true := by native_decide
  simpa [classCoverageCheck, classCoverageProp] using h

def coefficientCoverageProp : Prop :=
  coefficientSupport.all (fun coeff => decide (coeff ∈ expectedCoefficients)) = true ∧
    expectedCoefficients.all (fun coeff => decide (coeff ∈ coefficientSupport)) = true

def coefficientCoverageCheck : Bool :=
  coefficientSupport.all (fun coeff => decide (coeff ∈ expectedCoefficients)) &&
    expectedCoefficients.all (fun coeff => decide (coeff ∈ coefficientSupport))

theorem coefficientCoverage_true : coefficientCoverageProp := by
  have h : coefficientCoverageCheck = true := by native_decide
  simpa [coefficientCoverageCheck, coefficientCoverageProp] using h

def coeffThreeLadderProp : Prop :=
  coeffThreePowerLadder.all (fun power => decide (⟨power, 3⟩ ∈ witnessClasses)) = true ∧
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ 3 ∨ obstruction.power ∈ coeffThreePowerLadder)) = true

def coeffThreeLadderCheck : Bool :=
  coeffThreePowerLadder.all (fun power => decide (⟨power, 3⟩ ∈ witnessClasses)) &&
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ 3 ∨ obstruction.power ∈ coeffThreePowerLadder))

theorem coeffThreeLadder_true : coeffThreeLadderProp := by
  have h : coeffThreeLadderCheck = true := by native_decide
  simpa [coeffThreeLadderCheck, coeffThreeLadderProp] using h

def exceptionalCoeffPowerProp : Prop :=
  ⟨2, -1⟩ ∈ witnessClasses ∧
    ⟨2, 8⟩ ∈ witnessClasses ∧
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ -1 ∨ obstruction.power = 2)) = true ∧
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ 8 ∨ obstruction.power = 2)) = true

def exceptionalCoeffPowerCheck : Bool :=
  decide (⟨2, -1⟩ ∈ witnessClasses) &&
    decide (⟨2, 8⟩ ∈ witnessClasses) &&
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ -1 ∨ obstruction.power = 2)) &&
    witnessClasses.all (fun obstruction =>
      decide (obstruction.coeff ≠ 8 ∨ obstruction.power = 2))

theorem exceptionalCoeffPower_true : exceptionalCoeffPowerProp := by
  have h : exceptionalCoeffPowerCheck = true := by native_decide
  simpa [exceptionalCoeffPowerCheck, exceptionalCoeffPowerProp, and_assoc] using h

def coordinateWaypoint : Prop :=
  coordinateLabel = "P_ws" ∧
    coordinateExpression = "P_ws = (1/F - F) / 2" ∧
    templateLabel = "Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`"

theorem coordinateWaypoint_true : coordinateWaypoint := by
  simp [coordinateWaypoint, coordinateLabel, coordinateExpression, templateLabel]

structure WaypointCertificate where
  coordinate : coordinateWaypoint
  witnessTable : witnessTableProp
  sampleWitness : ∀ sample : TailSample, sampleWitnessProp sample
  firstFailure : ∀ sample : TailSample, firstFailureProp sample
  classCoverage : classCoverageProp
  coefficientCoverage : coefficientCoverageProp
  coeffThreeLadder : coeffThreeLadderProp
  exceptionalCoeffPower : exceptionalCoeffPowerProp

def currentWaypointCertificate : WaypointCertificate where
  coordinate := coordinateWaypoint_true
  witnessTable := witnessTable_true
  sampleWitness := sampleWitness_true
  firstFailure := firstFailure_true
  classCoverage := classCoverage_true
  coefficientCoverage := coefficientCoverage_true
  coeffThreeLadder := coeffThreeLadder_true
  exceptionalCoeffPower := exceptionalCoeffPower_true

def currentWaypoint : Prop :=
  coordinateWaypoint ∧
    witnessTableProp ∧
    (∀ sample : TailSample, sampleWitnessProp sample) ∧
    (∀ sample : TailSample, firstFailureProp sample) ∧
    classCoverageProp ∧
    coefficientCoverageProp ∧
    coeffThreeLadderProp ∧
    exceptionalCoeffPowerProp

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentWaypointCertificate.coordinate,
    currentWaypointCertificate.witnessTable,
    currentWaypointCertificate.sampleWitness,
    currentWaypointCertificate.firstFailure,
    currentWaypointCertificate.classCoverage,
    currentWaypointCertificate.coefficientCoverage,
    currentWaypointCertificate.coeffThreeLadder,
    currentWaypointCertificate.exceptionalCoeffPower
  ⟩

end WeberSchlafliCoordinate
end HeroCase
end Proofs
