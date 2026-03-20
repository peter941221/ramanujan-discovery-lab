import Mathlib

namespace Proofs
namespace HeroCase
namespace WeberClassInvariantBridge

/-!
Focused Weber class-invariant residual bridge shell.

This module records the current Python-side hand-off after the Ramanujan-Weber
class-invariant compression lane:

- `g12_ws` and `p12_ws` are the two named Weber coordinates
- the eta-side residual `G_g12_ws` is treated as the current primary residual
- the plus-side residual `G_p12_ws` is kept as an algebraically constrained
  companion through the exact coordinate bridge
- the derived quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2` compresses the
  bridge into a smaller algebraic lane next to the quotient
- the exact quotient-coordinate bridge
  `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0` with
  `Q_gp_ws = p12_ws / g12_ws` is the current tighter elimination relation
- the template-normalized coordinate
  `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`
  is the current source-faithful positive-recognition object
- the normalized follow-up `H_X_ws = (G_X_ws - 1) / (4*t^1)` is the next
  stripped version of that template-normalized coordinate on the hero side
- the first direct bridge between the two normalized hero-side follow-ups is
  now also recorded:
  `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws`
- the quotient-bridge follow-up
  `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)` is the next stripped constant-1 object
  on that comparison lane
- the same quotient-follow-up lane is now also compared back to the
  template-normalized branch:
  `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws`
- the quotient-follow-up bridge now also has its own stripped quotient object
  `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)`
- the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws` is the next narrow
  diagnostic object
- the normalized follow-up `H_gp_ws = (R_gp_ws - 1) / (96*t^3)` is the current
  next object after removing the first hero-specific quotient failure term
- the first theorem-shaped uniqueness / closure boxes on that quotient are
  still empty: low-degree self-polynomial, self-fractional-linear, and
  self-quotient finite-product scans, in addition to the first eta /
  modular-unit / plus-Pochhammer boxes
- the same first theorem-shaped boxes are also still empty on the normalized
  follow-up object

This is again a theorem-shaped research shell, not yet a final positive source
identity theorem.
-/

structure FirstFailureWitness where
  power : Nat
  coeff : Rat
deriving Repr, DecidableEq

structure SmallBoxStatus where
  etaHits : Nat
  modularUnitEtaHits : Nat
  selfPolynomialHits : Nat
  selfFractionalLinearHits : Nat
  selfQuotientProductHits : Nat
  plusPochhammerHits : Nat
  plusPochhammerEtaHits : Nat
deriving Repr, DecidableEq

def primaryResidualLabel : String := "G_g12_ws"

def primaryResidualExpression : String :=
  "G_g12_ws = g12_ws / (t^2; t^4)_inf^12"

def primaryResidualReason : String :=
  "The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its algebraically constrained companion."

def companionResidualLabel : String := "G_p12_ws"

def companionResidualExpression : String :=
  "G_p12_ws = p12_ws / (-t^2; t^4)_inf^12"

def exactCoordinateBridgeExpression : String :=
  "g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0"

def exactResidualBridgeExpression : String :=
  "(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2 - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4 + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2 + 4096*t^6 = 0"

def focusedQuotientCoordinateLabel : String := "X_g_ws"

def focusedQuotientCoordinateExpression : String :=
  "X_g_ws = 16*t^2 / g12_ws^2"

def exactQuotientCoordinateBridgeExpression : String :=
  "Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws"

def templateNormalizedCoordinateLabel : String := "G_X_ws"

def templateNormalizedCoordinateExpression : String :=
  "G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2"

def templateNormalizedCoordinateBridgeExpression : String :=
  "G_X_ws*G_g12_ws^2 - 1 = 0"

def templateNormalizedCoordinateFirstFailure : FirstFailureWitness := ⟨1, 4⟩

def templateNormalizedFollowupLabel : String := "H_X_ws"

def templateNormalizedFollowupExpression : String :=
  "H_X_ws = (G_X_ws - 1) / (4*t^1)"

def templateNormalizedFollowupFirstFailure : FirstFailureWitness := ⟨1, (9 : Rat) / 2⟩

def followupBridgeDifferenceLabel : String := "D_XR_ws"

def followupBridgeDifferenceExpression : String :=
  "D_XR_ws = H_gp_ws - H_X_ws"

def followupBridgeDifferenceFirstFailure : FirstFailureWitness := ⟨2, -24⟩

def followupBridgeQuotientLabel : String := "Q_XR_ws"

def followupBridgeQuotientExpression : String :=
  "Q_XR_ws = H_gp_ws / H_X_ws"

def followupBridgeQuotientFirstFailure : FirstFailureWitness := ⟨2, -24⟩

def followupBridgeQuotientFollowupLabel : String := "K_XR_ws"

def followupBridgeQuotientFollowupExpression : String :=
  "K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)"

def followupBridgeQuotientFollowupFirstFailure : FirstFailureWitness := ⟨1, 2⟩

def quotientFollowupBridgeDifferenceLabel : String := "D_XK_ws"

def quotientFollowupBridgeDifferenceExpression : String :=
  "D_XK_ws = K_XR_ws - H_X_ws"

def quotientFollowupBridgeDifferenceFirstFailure : FirstFailureWitness := ⟨1, (-5 : Rat) / 2⟩

def quotientFollowupBridgeQuotientLabel : String := "Q_XK_ws"

def quotientFollowupBridgeQuotientExpression : String :=
  "Q_XK_ws = K_XR_ws / H_X_ws"

def quotientFollowupBridgeQuotientFirstFailure : FirstFailureWitness := ⟨1, (-5 : Rat) / 2⟩

def quotientFollowupBridgeQuotientFollowupLabel : String := "L_XK_ws"

def quotientFollowupBridgeQuotientFollowupExpression : String :=
  "L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)"

def quotientFollowupBridgeQuotientFollowupFirstFailure : FirstFailureWitness := ⟨1, (563 : Rat) / 30⟩

def followupBridgePolynomialDegrees : List Nat := [1, 2, 3]

def followupBridgePolynomialHitCount : Nat := 0

def followupBridgeFractionalLinearHit : Bool := false

def followupBridgeQuotientSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def followupBridgeQuotientFollowupSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def quotientFollowupBridgePolynomialDegrees : List Nat := [1, 2, 3]

def quotientFollowupBridgePolynomialHitCount : Nat := 0

def quotientFollowupBridgeFractionalLinearHit : Bool := false

def quotientFollowupBridgeQuotientSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def quotientFollowupBridgeQuotientFollowupSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def focusedQuotientLabel : String := "R_gp_ws"

def focusedQuotientExpression : String := "R_gp_ws = G_p12_ws / G_g12_ws"

def focusedQuotientFirstFailure : FirstFailureWitness := ⟨3, 96⟩

def normalizedFollowupLabel : String := "H_gp_ws"

def normalizedFollowupExpression : String :=
  "H_gp_ws = (R_gp_ws - 1) / (96*t^3)"

def checkedEtaLevels : List Nat := [1, 2, 4]

def checkedModuli : List Nat := [2, 3, 4]

def templateNormalizedCoordinateSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def templateNormalizedFollowupSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def focusedQuotientSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def normalizedFollowupSmallBoxStatus : SmallBoxStatus := ⟨0, 0, 0, 0, 0, 0, 0⟩

def exactCoordinateBridgeProp : Prop :=
  exactCoordinateBridgeExpression =
    "g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0"

theorem exactCoordinateBridge_true : exactCoordinateBridgeProp := by
  simp [exactCoordinateBridgeProp, exactCoordinateBridgeExpression]

def exactQuotientCoordinateBridgeProp : Prop :=
  focusedQuotientCoordinateLabel = "X_g_ws" ∧
    focusedQuotientCoordinateExpression = "X_g_ws = 16*t^2 / g12_ws^2" ∧
    exactQuotientCoordinateBridgeExpression =
      "Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, Q_gp_ws = p12_ws / g12_ws"

theorem exactQuotientCoordinateBridge_true : exactQuotientCoordinateBridgeProp := by
  simp [
    exactQuotientCoordinateBridgeProp,
    focusedQuotientCoordinateLabel,
    focusedQuotientCoordinateExpression,
    exactQuotientCoordinateBridgeExpression
  ]

def templateNormalizedCoordinateProp : Prop :=
  templateNormalizedCoordinateLabel = "G_X_ws" ∧
    templateNormalizedCoordinateExpression =
      "G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2" ∧
    templateNormalizedCoordinateBridgeExpression = "G_X_ws*G_g12_ws^2 - 1 = 0"

theorem templateNormalizedCoordinate_true : templateNormalizedCoordinateProp := by
  simp [
    templateNormalizedCoordinateProp,
    templateNormalizedCoordinateLabel,
    templateNormalizedCoordinateExpression,
    templateNormalizedCoordinateBridgeExpression
  ]

def primaryResidualSelectionProp : Prop :=
  primaryResidualLabel = "G_g12_ws" ∧
    primaryResidualExpression = "G_g12_ws = g12_ws / (t^2; t^4)_inf^12" ∧
    companionResidualLabel = "G_p12_ws" ∧
    companionResidualExpression = "G_p12_ws = p12_ws / (-t^2; t^4)_inf^12" ∧
    focusedQuotientCoordinateLabel = "X_g_ws" ∧
    focusedQuotientCoordinateExpression = "X_g_ws = 16*t^2 / g12_ws^2" ∧
    templateNormalizedCoordinateLabel = "G_X_ws" ∧
    templateNormalizedCoordinateExpression =
      "G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2" ∧
    templateNormalizedFollowupLabel = "H_X_ws" ∧
    templateNormalizedFollowupExpression = "H_X_ws = (G_X_ws - 1) / (4*t^1)" ∧
    followupBridgeDifferenceLabel = "D_XR_ws" ∧
    followupBridgeDifferenceExpression = "D_XR_ws = H_gp_ws - H_X_ws" ∧
    followupBridgeQuotientLabel = "Q_XR_ws" ∧
    followupBridgeQuotientExpression = "Q_XR_ws = H_gp_ws / H_X_ws" ∧
    followupBridgeQuotientFollowupLabel = "K_XR_ws" ∧
    followupBridgeQuotientFollowupExpression = "K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)" ∧
    quotientFollowupBridgeDifferenceLabel = "D_XK_ws" ∧
    quotientFollowupBridgeDifferenceExpression = "D_XK_ws = K_XR_ws - H_X_ws" ∧
    quotientFollowupBridgeQuotientLabel = "Q_XK_ws" ∧
    quotientFollowupBridgeQuotientExpression = "Q_XK_ws = K_XR_ws / H_X_ws" ∧
    quotientFollowupBridgeQuotientFollowupLabel = "L_XK_ws" ∧
    quotientFollowupBridgeQuotientFollowupExpression = "L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)" ∧
    focusedQuotientLabel = "R_gp_ws" ∧
    focusedQuotientExpression = "R_gp_ws = G_p12_ws / G_g12_ws" ∧
    normalizedFollowupLabel = "H_gp_ws" ∧
    normalizedFollowupExpression = "H_gp_ws = (R_gp_ws - 1) / (96*t^3)"

theorem primaryResidualSelection_true : primaryResidualSelectionProp := by
  simp [
    primaryResidualSelectionProp,
    primaryResidualLabel,
    primaryResidualExpression,
    companionResidualLabel,
    companionResidualExpression,
    focusedQuotientCoordinateLabel,
    focusedQuotientCoordinateExpression,
    templateNormalizedCoordinateLabel,
    templateNormalizedCoordinateExpression,
    templateNormalizedFollowupLabel,
    templateNormalizedFollowupExpression,
    followupBridgeDifferenceLabel,
    followupBridgeDifferenceExpression,
    followupBridgeQuotientLabel,
    followupBridgeQuotientExpression,
    followupBridgeQuotientFollowupLabel,
    followupBridgeQuotientFollowupExpression,
    quotientFollowupBridgeDifferenceLabel,
    quotientFollowupBridgeDifferenceExpression,
    quotientFollowupBridgeQuotientLabel,
    quotientFollowupBridgeQuotientExpression,
    quotientFollowupBridgeQuotientFollowupLabel,
    quotientFollowupBridgeQuotientFollowupExpression,
    focusedQuotientLabel,
    focusedQuotientExpression,
    normalizedFollowupLabel,
    normalizedFollowupExpression
  ]

def templateNormalizedCoordinateFirstFailureProp : Prop :=
  templateNormalizedCoordinateFirstFailure.power = 1 ∧
    templateNormalizedCoordinateFirstFailure.coeff = 4

theorem templateNormalizedCoordinateFirstFailure_true :
    templateNormalizedCoordinateFirstFailureProp := by
  simp [templateNormalizedCoordinateFirstFailureProp, templateNormalizedCoordinateFirstFailure]

def templateNormalizedFollowupFirstFailureProp : Prop :=
  templateNormalizedFollowupFirstFailure.power = 1 ∧
    templateNormalizedFollowupFirstFailure.coeff = (9 : Rat) / 2

theorem templateNormalizedFollowupFirstFailure_true :
    templateNormalizedFollowupFirstFailureProp := by
  simp [templateNormalizedFollowupFirstFailureProp, templateNormalizedFollowupFirstFailure]

def followupBridgeDifferenceFirstFailureProp : Prop :=
  followupBridgeDifferenceFirstFailure.power = 2 ∧
    followupBridgeDifferenceFirstFailure.coeff = -24

theorem followupBridgeDifferenceFirstFailure_true :
    followupBridgeDifferenceFirstFailureProp := by
  simp [followupBridgeDifferenceFirstFailureProp, followupBridgeDifferenceFirstFailure]

def followupBridgeQuotientFirstFailureProp : Prop :=
  followupBridgeQuotientFirstFailure.power = 2 ∧
    followupBridgeQuotientFirstFailure.coeff = -24

theorem followupBridgeQuotientFirstFailure_true :
    followupBridgeQuotientFirstFailureProp := by
  simp [followupBridgeQuotientFirstFailureProp, followupBridgeQuotientFirstFailure]

def followupBridgeQuotientFollowupFirstFailureProp : Prop :=
  followupBridgeQuotientFollowupFirstFailure.power = 1 ∧
    followupBridgeQuotientFollowupFirstFailure.coeff = 2

theorem followupBridgeQuotientFollowupFirstFailure_true :
    followupBridgeQuotientFollowupFirstFailureProp := by
  simp [
    followupBridgeQuotientFollowupFirstFailureProp,
    followupBridgeQuotientFollowupFirstFailure
  ]

def quotientFollowupBridgeDifferenceFirstFailureProp : Prop :=
  quotientFollowupBridgeDifferenceFirstFailure.power = 1 ∧
    quotientFollowupBridgeDifferenceFirstFailure.coeff = (-5 : Rat) / 2

theorem quotientFollowupBridgeDifferenceFirstFailure_true :
    quotientFollowupBridgeDifferenceFirstFailureProp := by
  simp [
    quotientFollowupBridgeDifferenceFirstFailureProp,
    quotientFollowupBridgeDifferenceFirstFailure
  ]

def quotientFollowupBridgeQuotientFirstFailureProp : Prop :=
  quotientFollowupBridgeQuotientFirstFailure.power = 1 ∧
    quotientFollowupBridgeQuotientFirstFailure.coeff = (-5 : Rat) / 2

theorem quotientFollowupBridgeQuotientFirstFailure_true :
    quotientFollowupBridgeQuotientFirstFailureProp := by
  simp [
    quotientFollowupBridgeQuotientFirstFailureProp,
    quotientFollowupBridgeQuotientFirstFailure
  ]

def quotientFollowupBridgeQuotientFollowupFirstFailureProp : Prop :=
  quotientFollowupBridgeQuotientFollowupFirstFailure.power = 1 ∧
    quotientFollowupBridgeQuotientFollowupFirstFailure.coeff = (563 : Rat) / 30

theorem quotientFollowupBridgeQuotientFollowupFirstFailure_true :
    quotientFollowupBridgeQuotientFollowupFirstFailureProp := by
  simp [
    quotientFollowupBridgeQuotientFollowupFirstFailureProp,
    quotientFollowupBridgeQuotientFollowupFirstFailure
  ]

def focusedQuotientFirstFailureProp : Prop :=
  focusedQuotientFirstFailure.power = 3 ∧
    focusedQuotientFirstFailure.coeff = 96

theorem focusedQuotientFirstFailure_true : focusedQuotientFirstFailureProp := by
  simp [focusedQuotientFirstFailureProp, focusedQuotientFirstFailure]

def checkedSearchParametersProp : Prop :=
  checkedEtaLevels = [1, 2, 4] ∧
    checkedModuli = [2, 3, 4]

theorem checkedSearchParameters_true : checkedSearchParametersProp := by
  simp [checkedSearchParametersProp, checkedEtaLevels, checkedModuli]

def templateNormalizedCoordinateSmallBoxStatusProp : Prop :=
  templateNormalizedCoordinateSmallBoxStatus.etaHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.modularUnitEtaHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.selfPolynomialHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.selfQuotientProductHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.plusPochhammerHits = 0 ∧
    templateNormalizedCoordinateSmallBoxStatus.plusPochhammerEtaHits = 0

theorem templateNormalizedCoordinateSmallBoxStatus_true :
    templateNormalizedCoordinateSmallBoxStatusProp := by
  simp [templateNormalizedCoordinateSmallBoxStatusProp, templateNormalizedCoordinateSmallBoxStatus]

def templateNormalizedFollowupSmallBoxStatusProp : Prop :=
  templateNormalizedFollowupSmallBoxStatus.etaHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.modularUnitEtaHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.selfPolynomialHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.selfQuotientProductHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.plusPochhammerHits = 0 ∧
    templateNormalizedFollowupSmallBoxStatus.plusPochhammerEtaHits = 0

theorem templateNormalizedFollowupSmallBoxStatus_true :
    templateNormalizedFollowupSmallBoxStatusProp := by
  simp [templateNormalizedFollowupSmallBoxStatusProp, templateNormalizedFollowupSmallBoxStatus]

def followupBridgeStatusProp : Prop :=
  followupBridgePolynomialDegrees = [1, 2, 3] ∧
    followupBridgePolynomialHitCount = 0 ∧
    followupBridgeFractionalLinearHit = false

theorem followupBridgeStatus_true : followupBridgeStatusProp := by
  simp [
    followupBridgeStatusProp,
    followupBridgePolynomialDegrees,
    followupBridgePolynomialHitCount,
    followupBridgeFractionalLinearHit
  ]

def followupBridgeQuotientSmallBoxStatusProp : Prop :=
  followupBridgeQuotientSmallBoxStatus.etaHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.modularUnitEtaHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.selfPolynomialHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.selfQuotientProductHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.plusPochhammerHits = 0 ∧
    followupBridgeQuotientSmallBoxStatus.plusPochhammerEtaHits = 0

theorem followupBridgeQuotientSmallBoxStatus_true :
    followupBridgeQuotientSmallBoxStatusProp := by
  simp [followupBridgeQuotientSmallBoxStatusProp, followupBridgeQuotientSmallBoxStatus]

def followupBridgeQuotientFollowupSmallBoxStatusProp : Prop :=
  followupBridgeQuotientFollowupSmallBoxStatus.etaHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.modularUnitEtaHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.selfPolynomialHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.selfQuotientProductHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.plusPochhammerHits = 0 ∧
    followupBridgeQuotientFollowupSmallBoxStatus.plusPochhammerEtaHits = 0

theorem followupBridgeQuotientFollowupSmallBoxStatus_true :
    followupBridgeQuotientFollowupSmallBoxStatusProp := by
  simp [
    followupBridgeQuotientFollowupSmallBoxStatusProp,
    followupBridgeQuotientFollowupSmallBoxStatus
  ]

def quotientFollowupBridgeStatusProp : Prop :=
  quotientFollowupBridgePolynomialDegrees = [1, 2, 3] ∧
    quotientFollowupBridgePolynomialHitCount = 0 ∧
    quotientFollowupBridgeFractionalLinearHit = false

theorem quotientFollowupBridgeStatus_true : quotientFollowupBridgeStatusProp := by
  simp [
    quotientFollowupBridgeStatusProp,
    quotientFollowupBridgePolynomialDegrees,
    quotientFollowupBridgePolynomialHitCount,
    quotientFollowupBridgeFractionalLinearHit
  ]

def quotientFollowupBridgeQuotientSmallBoxStatusProp : Prop :=
  quotientFollowupBridgeQuotientSmallBoxStatus.etaHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.modularUnitEtaHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.selfPolynomialHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.selfQuotientProductHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.plusPochhammerHits = 0 ∧
    quotientFollowupBridgeQuotientSmallBoxStatus.plusPochhammerEtaHits = 0

theorem quotientFollowupBridgeQuotientSmallBoxStatus_true :
    quotientFollowupBridgeQuotientSmallBoxStatusProp := by
  simp [
    quotientFollowupBridgeQuotientSmallBoxStatusProp,
    quotientFollowupBridgeQuotientSmallBoxStatus
  ]

def quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp : Prop :=
  quotientFollowupBridgeQuotientFollowupSmallBoxStatus.etaHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.modularUnitEtaHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.selfPolynomialHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.selfQuotientProductHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.plusPochhammerHits = 0 ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus.plusPochhammerEtaHits = 0

theorem quotientFollowupBridgeQuotientFollowupSmallBoxStatus_true :
    quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp := by
  simp [
    quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp,
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus
  ]

def quotientFollowupReturnBridgeProp : Prop :=
  quotientFollowupBridgeDifferenceLabel = "D_XK_ws" ∧
    quotientFollowupBridgeDifferenceExpression = "D_XK_ws = K_XR_ws - H_X_ws" ∧
    quotientFollowupBridgeDifferenceFirstFailureProp ∧
    quotientFollowupBridgeQuotientLabel = "Q_XK_ws" ∧
    quotientFollowupBridgeQuotientExpression = "Q_XK_ws = K_XR_ws / H_X_ws" ∧
    quotientFollowupBridgeQuotientFirstFailureProp ∧
    quotientFollowupBridgeQuotientFollowupLabel = "L_XK_ws" ∧
    quotientFollowupBridgeQuotientFollowupExpression =
      "L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)" ∧
    quotientFollowupBridgeQuotientFollowupFirstFailureProp ∧
    quotientFollowupBridgeStatusProp ∧
    quotientFollowupBridgeQuotientSmallBoxStatusProp ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp

theorem quotientFollowupReturnBridge_true : quotientFollowupReturnBridgeProp := by
  simp [
    quotientFollowupReturnBridgeProp,
    quotientFollowupBridgeDifferenceLabel,
    quotientFollowupBridgeDifferenceExpression,
    quotientFollowupBridgeQuotientLabel,
    quotientFollowupBridgeQuotientExpression,
    quotientFollowupBridgeQuotientFollowupLabel,
    quotientFollowupBridgeQuotientFollowupExpression,
    quotientFollowupBridgeDifferenceFirstFailure_true,
    quotientFollowupBridgeQuotientFirstFailure_true,
    quotientFollowupBridgeQuotientFollowupFirstFailure_true,
    quotientFollowupBridgeStatus_true,
    quotientFollowupBridgeQuotientSmallBoxStatus_true,
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus_true
  ]

def focusedQuotientSmallBoxStatusProp : Prop :=
  focusedQuotientSmallBoxStatus.etaHits = 0 ∧
    focusedQuotientSmallBoxStatus.modularUnitEtaHits = 0 ∧
    focusedQuotientSmallBoxStatus.selfPolynomialHits = 0 ∧
    focusedQuotientSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    focusedQuotientSmallBoxStatus.selfQuotientProductHits = 0 ∧
    focusedQuotientSmallBoxStatus.plusPochhammerHits = 0 ∧
    focusedQuotientSmallBoxStatus.plusPochhammerEtaHits = 0

theorem focusedQuotientSmallBoxStatus_true : focusedQuotientSmallBoxStatusProp := by
  simp [focusedQuotientSmallBoxStatusProp, focusedQuotientSmallBoxStatus]

def normalizedFollowupSmallBoxStatusProp : Prop :=
  normalizedFollowupSmallBoxStatus.etaHits = 0 ∧
    normalizedFollowupSmallBoxStatus.modularUnitEtaHits = 0 ∧
    normalizedFollowupSmallBoxStatus.selfPolynomialHits = 0 ∧
    normalizedFollowupSmallBoxStatus.selfFractionalLinearHits = 0 ∧
    normalizedFollowupSmallBoxStatus.selfQuotientProductHits = 0 ∧
    normalizedFollowupSmallBoxStatus.plusPochhammerHits = 0 ∧
    normalizedFollowupSmallBoxStatus.plusPochhammerEtaHits = 0

theorem normalizedFollowupSmallBoxStatus_true : normalizedFollowupSmallBoxStatusProp := by
  simp [normalizedFollowupSmallBoxStatusProp, normalizedFollowupSmallBoxStatus]

structure WaypointCertificate where
  bridge : exactCoordinateBridgeProp
  quotientCoordinateBridge : exactQuotientCoordinateBridgeProp
  templateCoordinate : templateNormalizedCoordinateProp
  selection : primaryResidualSelectionProp
  templateFirstFailure : templateNormalizedCoordinateFirstFailureProp
  templateFollowupFirstFailure : templateNormalizedFollowupFirstFailureProp
  followupBridgeDifference : followupBridgeDifferenceFirstFailureProp
  followupBridgeQuotient : followupBridgeQuotientFirstFailureProp
  followupBridgeQuotientFollowup : followupBridgeQuotientFollowupFirstFailureProp
  quotientFollowupReturnBridge : quotientFollowupReturnBridgeProp
  quotientFollowupBridgeDifference : quotientFollowupBridgeDifferenceFirstFailureProp
  quotientFollowupBridgeQuotient : quotientFollowupBridgeQuotientFirstFailureProp
  quotientFollowupBridgeQuotientFollowup : quotientFollowupBridgeQuotientFollowupFirstFailureProp
  firstFailure : focusedQuotientFirstFailureProp
  searchParameters : checkedSearchParametersProp
  templateSmallBoxes : templateNormalizedCoordinateSmallBoxStatusProp
  templateFollowupSmallBoxes : templateNormalizedFollowupSmallBoxStatusProp
  followupBridgeStatus : followupBridgeStatusProp
  followupBridgeQuotientSmallBoxes : followupBridgeQuotientSmallBoxStatusProp
  followupBridgeQuotientFollowupSmallBoxes : followupBridgeQuotientFollowupSmallBoxStatusProp
  quotientFollowupBridgeStatus : quotientFollowupBridgeStatusProp
  quotientFollowupBridgeQuotientSmallBoxes : quotientFollowupBridgeQuotientSmallBoxStatusProp
  quotientFollowupBridgeQuotientFollowupSmallBoxes :
    quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp
  smallBoxes : focusedQuotientSmallBoxStatusProp
  normalizedSmallBoxes : normalizedFollowupSmallBoxStatusProp

def currentWaypointCertificate : WaypointCertificate where
  bridge := exactCoordinateBridge_true
  quotientCoordinateBridge := exactQuotientCoordinateBridge_true
  templateCoordinate := templateNormalizedCoordinate_true
  selection := primaryResidualSelection_true
  templateFirstFailure := templateNormalizedCoordinateFirstFailure_true
  templateFollowupFirstFailure := templateNormalizedFollowupFirstFailure_true
  followupBridgeDifference := followupBridgeDifferenceFirstFailure_true
  followupBridgeQuotient := followupBridgeQuotientFirstFailure_true
  followupBridgeQuotientFollowup := followupBridgeQuotientFollowupFirstFailure_true
  quotientFollowupReturnBridge := quotientFollowupReturnBridge_true
  quotientFollowupBridgeDifference := quotientFollowupBridgeDifferenceFirstFailure_true
  quotientFollowupBridgeQuotient := quotientFollowupBridgeQuotientFirstFailure_true
  quotientFollowupBridgeQuotientFollowup := quotientFollowupBridgeQuotientFollowupFirstFailure_true
  firstFailure := focusedQuotientFirstFailure_true
  searchParameters := checkedSearchParameters_true
  templateSmallBoxes := templateNormalizedCoordinateSmallBoxStatus_true
  templateFollowupSmallBoxes := templateNormalizedFollowupSmallBoxStatus_true
  followupBridgeStatus := followupBridgeStatus_true
  followupBridgeQuotientSmallBoxes := followupBridgeQuotientSmallBoxStatus_true
  followupBridgeQuotientFollowupSmallBoxes := followupBridgeQuotientFollowupSmallBoxStatus_true
  quotientFollowupBridgeStatus := quotientFollowupBridgeStatus_true
  quotientFollowupBridgeQuotientSmallBoxes := quotientFollowupBridgeQuotientSmallBoxStatus_true
  quotientFollowupBridgeQuotientFollowupSmallBoxes :=
    quotientFollowupBridgeQuotientFollowupSmallBoxStatus_true
  smallBoxes := focusedQuotientSmallBoxStatus_true
  normalizedSmallBoxes := normalizedFollowupSmallBoxStatus_true

def currentWaypoint : Prop :=
  exactCoordinateBridgeProp ∧
    exactQuotientCoordinateBridgeProp ∧
    templateNormalizedCoordinateProp ∧
    primaryResidualSelectionProp ∧
    templateNormalizedCoordinateFirstFailureProp ∧
    templateNormalizedFollowupFirstFailureProp ∧
    followupBridgeDifferenceFirstFailureProp ∧
    followupBridgeQuotientFirstFailureProp ∧
    followupBridgeQuotientFollowupFirstFailureProp ∧
    quotientFollowupReturnBridgeProp ∧
    quotientFollowupBridgeDifferenceFirstFailureProp ∧
    quotientFollowupBridgeQuotientFirstFailureProp ∧
    quotientFollowupBridgeQuotientFollowupFirstFailureProp ∧
    focusedQuotientFirstFailureProp ∧
    checkedSearchParametersProp ∧
    templateNormalizedCoordinateSmallBoxStatusProp ∧
    templateNormalizedFollowupSmallBoxStatusProp ∧
    followupBridgeStatusProp ∧
    followupBridgeQuotientSmallBoxStatusProp ∧
    followupBridgeQuotientFollowupSmallBoxStatusProp ∧
    quotientFollowupBridgeStatusProp ∧
    quotientFollowupBridgeQuotientSmallBoxStatusProp ∧
    quotientFollowupBridgeQuotientFollowupSmallBoxStatusProp ∧
    focusedQuotientSmallBoxStatusProp ∧
    normalizedFollowupSmallBoxStatusProp

theorem currentWaypoint_true : currentWaypoint := by
  exact ⟨
    currentWaypointCertificate.bridge,
    currentWaypointCertificate.quotientCoordinateBridge,
    currentWaypointCertificate.templateCoordinate,
    currentWaypointCertificate.selection,
    currentWaypointCertificate.templateFirstFailure,
    currentWaypointCertificate.templateFollowupFirstFailure,
    currentWaypointCertificate.followupBridgeDifference,
    currentWaypointCertificate.followupBridgeQuotient,
    currentWaypointCertificate.followupBridgeQuotientFollowup,
    currentWaypointCertificate.quotientFollowupReturnBridge,
    currentWaypointCertificate.quotientFollowupBridgeDifference,
    currentWaypointCertificate.quotientFollowupBridgeQuotient,
    currentWaypointCertificate.quotientFollowupBridgeQuotientFollowup,
    currentWaypointCertificate.firstFailure,
    currentWaypointCertificate.searchParameters,
    currentWaypointCertificate.templateSmallBoxes,
    currentWaypointCertificate.templateFollowupSmallBoxes,
    currentWaypointCertificate.followupBridgeStatus,
    currentWaypointCertificate.followupBridgeQuotientSmallBoxes,
    currentWaypointCertificate.followupBridgeQuotientFollowupSmallBoxes,
    currentWaypointCertificate.quotientFollowupBridgeStatus,
    currentWaypointCertificate.quotientFollowupBridgeQuotientSmallBoxes,
    currentWaypointCertificate.quotientFollowupBridgeQuotientFollowupSmallBoxes,
    currentWaypointCertificate.smallBoxes,
    currentWaypointCertificate.normalizedSmallBoxes
  ⟩

end WeberClassInvariantBridge
end HeroCase
end Proofs
