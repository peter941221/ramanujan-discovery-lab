import Mathlib

namespace Proofs
namespace HeroCase
namespace GGQ34

/-!
Tail-family-first `GG` quotient-coordinate obstruction scaffold.

This file does not yet prove the full modular-equation exclusion theorem.
Instead, it packages the current Python-side exact obstruction witnesses into
Lean data and small proposition shells so later exact proofs have a stable
landing spot.
-/

inductive TailSample where
  | u_t2
  | u_t2_g1
  | u_t2_g2
  | u_t2_g3
  | u_t3
  | u_t3_g1
  | u_t3_g2
  | u_t3_g3
  | u_t4
  | u_t4_g1
  | u_t4_g2
  | u_t4_g3
deriving Repr, DecidableEq

structure ObstructionWitness where
  sampleLabel : String
  q3Power : Nat
  q3Coeff : Int
  q4Power : Nat
  q4Coeff : Int
deriving Repr, DecidableEq

structure ObstructionClass where
  q3Power : Nat
  q3Coeff : Int
  q4Power : Nat
  q4Coeff : Int
deriving Repr, DecidableEq

structure ResidualPrefix where
  q3Coeff0 : Int
  q3Coeff1 : Int
  q3Coeff2 : Int
  q4Coeff0 : Int
  q4Coeff1 : Int
  q4Coeff2 : Int
deriving Repr, DecidableEq

structure TailDerivedCoordinate where
  failurePower : Nat
  obstructionScalar : Int
deriving Repr, DecidableEq

def TailSample.label : TailSample → String
  | .u_t2 => "U_t2"
  | .u_t2_g1 => "U_t2_g1"
  | .u_t2_g2 => "U_t2_g2"
  | .u_t2_g3 => "U_t2_g3"
  | .u_t3 => "U_t3"
  | .u_t3_g1 => "U_t3_g1"
  | .u_t3_g2 => "U_t3_g2"
  | .u_t3_g3 => "U_t3_g3"
  | .u_t4 => "U_t4"
  | .u_t4_g1 => "U_t4_g1"
  | .u_t4_g2 => "U_t4_g2"
  | .u_t4_g3 => "U_t4_g3"

def tailSamples : List TailSample :=
  [
    .u_t2,
    .u_t2_g1,
    .u_t2_g2,
    .u_t2_g3,
    .u_t3,
    .u_t3_g1,
    .u_t3_g2,
    .u_t3_g3,
    .u_t4,
    .u_t4_g1,
    .u_t4_g2,
    .u_t4_g3
  ]

def tailSampleLabels : List String :=
  [
    "U_t2",
    "U_t2_g1",
    "U_t2_g2",
    "U_t2_g3",
    "U_t3",
    "U_t3_g1",
    "U_t3_g2",
    "U_t3_g3",
    "U_t4",
    "U_t4_g1",
    "U_t4_g2",
    "U_t4_g3"
  ]

def witnessFor : TailSample → ObstructionWitness
  | .u_t2 => ⟨"U_t2", 1, 2, 1, 3⟩
  | .u_t2_g1 => ⟨"U_t2_g1", 1, 4, 1, 6⟩
  | .u_t2_g2 => ⟨"U_t2_g2", 2, -4, 2, -6⟩
  | .u_t2_g3 => ⟨"U_t2_g3", 1, 6, 1, 9⟩
  | .u_t3 => ⟨"U_t3", 1, 2, 1, 3⟩
  | .u_t3_g1 => ⟨"U_t3_g1", 1, 2, 1, 3⟩
  | .u_t3_g2 => ⟨"U_t3_g2", 2, -2, 2, -3⟩
  | .u_t3_g3 => ⟨"U_t3_g3", 1, 4, 1, 6⟩
  | .u_t4 => ⟨"U_t4", 1, 2, 1, 3⟩
  | .u_t4_g1 => ⟨"U_t4_g1", 1, 2, 1, 3⟩
  | .u_t4_g2 => ⟨"U_t4_g2", 2, -2, 2, -3⟩
  | .u_t4_g3 => ⟨"U_t4_g3", 1, 4, 1, 6⟩

def tailWitnesses : List ObstructionWitness :=
  [
    ⟨"U_t2", 1, 2, 1, 3⟩,
    ⟨"U_t2_g1", 1, 4, 1, 6⟩,
    ⟨"U_t2_g2", 2, -4, 2, -6⟩,
    ⟨"U_t2_g3", 1, 6, 1, 9⟩,
    ⟨"U_t3", 1, 2, 1, 3⟩,
    ⟨"U_t3_g1", 1, 2, 1, 3⟩,
    ⟨"U_t3_g2", 2, -2, 2, -3⟩,
    ⟨"U_t3_g3", 1, 4, 1, 6⟩,
    ⟨"U_t4", 1, 2, 1, 3⟩,
    ⟨"U_t4_g1", 1, 2, 1, 3⟩,
    ⟨"U_t4_g2", 2, -2, 2, -3⟩,
    ⟨"U_t4_g3", 1, 4, 1, 6⟩
  ]

def witnessClasses : List ObstructionClass :=
  tailWitnesses.map fun witness =>
    ⟨witness.q3Power, witness.q3Coeff, witness.q4Power, witness.q4Coeff⟩

def expectedClasses : List ObstructionClass :=
  [
    ⟨1, 2, 1, 3⟩,
    ⟨1, 4, 1, 6⟩,
    ⟨2, -2, 2, -3⟩,
    ⟨2, -4, 2, -6⟩,
    ⟨1, 6, 1, 9⟩
  ]

def classFor (sample : TailSample) : ObstructionClass :=
  let witness := witnessFor sample
  ⟨witness.q3Power, witness.q3Coeff, witness.q4Power, witness.q4Coeff⟩

def residualPrefixFor : TailSample → ResidualPrefix
  | .u_t2 => ⟨0, 2, 0, 0, 3, 0⟩
  | .u_t2_g1 => ⟨0, 4, 0, 0, 6, 0⟩
  | .u_t2_g2 => ⟨0, 0, -4, 0, 0, -6⟩
  | .u_t2_g3 => ⟨0, 6, 0, 0, 9, 0⟩
  | .u_t3 => ⟨0, 2, 0, 0, 3, 0⟩
  | .u_t3_g1 => ⟨0, 2, 0, 0, 3, 0⟩
  | .u_t3_g2 => ⟨0, 0, -2, 0, 0, -3⟩
  | .u_t3_g3 => ⟨0, 4, 0, 0, 6, 0⟩
  | .u_t4 => ⟨0, 2, 0, 0, 3, 0⟩
  | .u_t4_g1 => ⟨0, 2, 0, 0, 3, 0⟩
  | .u_t4_g2 => ⟨0, 0, -2, 0, 0, -3⟩
  | .u_t4_g3 => ⟨0, 4, 0, 0, 6, 0⟩

def tailDerivedCoordinateFor : TailSample → TailDerivedCoordinate
  | .u_t2 => ⟨1, 1⟩
  | .u_t2_g1 => ⟨1, 2⟩
  | .u_t2_g2 => ⟨2, -2⟩
  | .u_t2_g3 => ⟨1, 3⟩
  | .u_t3 => ⟨1, 1⟩
  | .u_t3_g1 => ⟨1, 1⟩
  | .u_t3_g2 => ⟨2, -1⟩
  | .u_t3_g3 => ⟨1, 2⟩
  | .u_t4 => ⟨1, 1⟩
  | .u_t4_g1 => ⟨1, 1⟩
  | .u_t4_g2 => ⟨2, -1⟩
  | .u_t4_g3 => ⟨1, 2⟩

def expectedTailDerivedCoordinates : List TailDerivedCoordinate :=
  [
    ⟨1, 1⟩,
    ⟨1, 2⟩,
    ⟨2, -1⟩,
    ⟨2, -2⟩,
    ⟨1, 3⟩
  ]

def firstFailureAtPrefix (coeff0 coeff1 coeff2 : Int) (power : Nat) (coeff : Int) : Prop :=
  coeff0 = 0 ∧
    ((power = 1 ∧ coeff1 = coeff ∧ coeff ≠ 0) ∨
      (power = 2 ∧ coeff1 = 0 ∧ coeff2 = coeff ∧ coeff ≠ 0))

def tailWitnessTableCheck : Bool :=
  decide (tailWitnesses.map ObstructionWitness.sampleLabel = tailSampleLabels) &&
  decide (tailWitnesses.length = 12) &&
  tailWitnesses.all fun witness =>
    decide (witness.q3Coeff ≠ 0 ∧ witness.q4Coeff ≠ 0)

def tailWitnessTableProp : Prop :=
  tailWitnessTableCheck = true

theorem tailWitnessTable_true : tailWitnessTableProp := by
  simpa [tailWitnessTableProp] using (show tailWitnessTableCheck = true by native_decide)

def tailWitnessClassCoverageCheck : Bool :=
  witnessClasses.all fun obstruction =>
    decide (obstruction ∈ expectedClasses) &&
  expectedClasses.all fun obstruction =>
    decide (obstruction ∈ witnessClasses)

def tailWitnessClassCoverageProp : Prop :=
  tailWitnessClassCoverageCheck = true

theorem tailWitnessClassCoverage_true : tailWitnessClassCoverageProp := by
  simpa [tailWitnessClassCoverageProp] using
    (show tailWitnessClassCoverageCheck = true by native_decide)

def sampleWitnessCheck (sample : TailSample) : Bool :=
  let witness := witnessFor sample
  let obstruction := classFor sample
  decide (witness.sampleLabel = sample.label) &&
    decide (witness ∈ tailWitnesses) &&
    decide (obstruction ∈ expectedClasses) &&
    decide (witness.q3Coeff ≠ 0) &&
    decide (witness.q4Coeff ≠ 0)

def sampleWitnessProp (sample : TailSample) : Prop :=
  sampleWitnessCheck sample = true

theorem sampleWitness_true (sample : TailSample) : sampleWitnessProp sample := by
  cases sample <;>
    simp [sampleWitnessProp, sampleWitnessCheck, witnessFor, classFor, TailSample.label,
      tailWitnesses, expectedClasses]

def q3FirstFailureProp (sample : TailSample) : Prop :=
  firstFailureAtPrefix
    (residualPrefixFor sample).q3Coeff0
    (residualPrefixFor sample).q3Coeff1
    (residualPrefixFor sample).q3Coeff2
    (witnessFor sample).q3Power
    (witnessFor sample).q3Coeff

theorem q3FirstFailure_true (sample : TailSample) : q3FirstFailureProp sample := by
  cases sample <;>
    simp [q3FirstFailureProp, residualPrefixFor, witnessFor, firstFailureAtPrefix]

def q4FirstFailureProp (sample : TailSample) : Prop :=
  firstFailureAtPrefix
    (residualPrefixFor sample).q4Coeff0
    (residualPrefixFor sample).q4Coeff1
    (residualPrefixFor sample).q4Coeff2
    (witnessFor sample).q4Power
    (witnessFor sample).q4Coeff

theorem q4FirstFailure_true (sample : TailSample) : q4FirstFailureProp sample := by
  cases sample <;>
    simp [q4FirstFailureProp, residualPrefixFor, witnessFor, firstFailureAtPrefix]

def pairedFirstFailureProp (sample : TailSample) : Prop :=
  q3FirstFailureProp sample ∧
    q4FirstFailureProp sample ∧
    (witnessFor sample).q3Power = (witnessFor sample).q4Power

theorem pairedFirstFailure_true (sample : TailSample) : pairedFirstFailureProp sample := by
  refine ⟨q3FirstFailure_true sample, q4FirstFailure_true sample, ?_⟩
  cases sample <;> simp [witnessFor]

def q3q4ExactObstructionCheck (sample : TailSample) : Bool :=
  let witness := witnessFor sample
  decide (witness.q3Coeff ≠ 0) && decide (witness.q4Coeff ≠ 0)

def q3q4ExactObstructionProp (sample : TailSample) : Prop :=
  q3q4ExactObstructionCheck sample = true

theorem q3q4ExactObstruction_true (sample : TailSample) :
    q3q4ExactObstructionProp sample := by
  cases sample <;>
    simp [q3q4ExactObstructionProp, q3q4ExactObstructionCheck, witnessFor]

def obstructionClassCheck (sample : TailSample) : Bool :=
  decide (classFor sample ∈ expectedClasses)

def obstructionClassProp (sample : TailSample) : Prop :=
  obstructionClassCheck sample = true

theorem obstructionClass_true (sample : TailSample) :
    obstructionClassProp sample := by
  cases sample <;>
    simp [obstructionClassProp, obstructionClassCheck, classFor, witnessFor, expectedClasses]

def tailDerivedCoordinateCheck (sample : TailSample) : Bool :=
  decide (tailDerivedCoordinateFor sample ∈ expectedTailDerivedCoordinates)

def tailDerivedCoordinateProp (sample : TailSample) : Prop :=
  tailDerivedCoordinateCheck sample = true

theorem tailDerivedCoordinate_true (sample : TailSample) :
    tailDerivedCoordinateProp sample := by
  cases sample <;>
    simp [tailDerivedCoordinateProp, tailDerivedCoordinateCheck, tailDerivedCoordinateFor,
      expectedTailDerivedCoordinates]

def sharedLeadingObstructionProp (sample : TailSample) : Prop :=
  (witnessFor sample).q3Power = (tailDerivedCoordinateFor sample).failurePower ∧
    (witnessFor sample).q4Power = (tailDerivedCoordinateFor sample).failurePower ∧
    (witnessFor sample).q3Coeff = 2 * (tailDerivedCoordinateFor sample).obstructionScalar ∧
    (witnessFor sample).q4Coeff = 3 * (tailDerivedCoordinateFor sample).obstructionScalar ∧
    (tailDerivedCoordinateFor sample).obstructionScalar ≠ 0

theorem sharedLeadingObstruction_true (sample : TailSample) :
    sharedLeadingObstructionProp sample := by
  cases sample <;>
    simp [sharedLeadingObstructionProp, witnessFor, tailDerivedCoordinateFor]

def everyTailSampleHasExactObstructionProp : Prop :=
  ∀ sample : TailSample, q3q4ExactObstructionProp sample

theorem everyTailSampleHasExactObstruction_true :
    everyTailSampleHasExactObstructionProp := by
  intro sample
  exact q3q4ExactObstruction_true sample

def everyTailSampleHasFirstFailureProp : Prop :=
  ∀ sample : TailSample, pairedFirstFailureProp sample

theorem everyTailSampleHasFirstFailure_true :
    everyTailSampleHasFirstFailureProp := by
  intro sample
  exact pairedFirstFailure_true sample

def everyTailSampleHasSharedLeadingObstructionProp : Prop :=
  ∀ sample : TailSample, sharedLeadingObstructionProp sample

theorem everyTailSampleHasSharedLeadingObstruction_true :
    everyTailSampleHasSharedLeadingObstructionProp := by
  intro sample
  exact sharedLeadingObstruction_true sample

structure ObstructionCertificate where
  table : tailWitnessTableProp
  classes : tailWitnessClassCoverageProp
  firstFailure : everyTailSampleHasFirstFailureProp
  samplewise : everyTailSampleHasExactObstructionProp
  sharedLeading : everyTailSampleHasSharedLeadingObstructionProp

def currentObstructionCertificate : ObstructionCertificate where
  table := tailWitnessTable_true
  classes := tailWitnessClassCoverage_true
  firstFailure := everyTailSampleHasFirstFailure_true
  samplewise := everyTailSampleHasExactObstruction_true
  sharedLeading := everyTailSampleHasSharedLeadingObstruction_true

def currentWaypoint : Prop :=
  tailWitnessTableProp ∧
  tailWitnessClassCoverageProp ∧
  everyTailSampleHasFirstFailureProp ∧
  everyTailSampleHasExactObstructionProp ∧
  everyTailSampleHasSharedLeadingObstructionProp

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentObstructionCertificate.table,
    currentObstructionCertificate.classes,
    currentObstructionCertificate.firstFailure,
    currentObstructionCertificate.samplewise,
    currentObstructionCertificate.sharedLeading
  ⟩

end GGQ34
end HeroCase
end Proofs
