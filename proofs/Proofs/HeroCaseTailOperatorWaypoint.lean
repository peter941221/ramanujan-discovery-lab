import Mathlib

namespace Proofs
namespace HeroCase
namespace TailOperator

/-!
Tail-family operator-lane waypoint shell.

This file does not prove the final q-difference / operator statement.
It records the current sampled operator lane so the proof workspace has a
stable landing spot once an exact recurrence or factorization theorem appears.
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

def checkedModuli : List Nat := [2, 3]

def checkedRecurrenceDepths : List Nat := [2, 3]

def checkedTDegrees : List Nat := [1, 2, 3]

def currentWaypoint : Prop :=
  tailSamples.length = 12 ∧
    checkedModuli = [2, 3] ∧
    checkedRecurrenceDepths = [2, 3] ∧
    checkedTDegrees = [1, 2, 3]

theorem currentWaypoint_true : currentWaypoint := by
  simp [currentWaypoint, tailSamples, checkedModuli, checkedRecurrenceDepths, checkedTDegrees]

end TailOperator
end HeroCase
end Proofs
