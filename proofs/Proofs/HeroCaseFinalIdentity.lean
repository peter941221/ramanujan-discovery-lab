import Proofs.HeroCaseObjects
import Proofs.HeroCaseGGQuotientCoordinateObstruction
import Proofs.HeroCaseMortonSquaredCoordinateObstruction
import Proofs.HeroCaseMortonNamedCoordinateWaypoint
import Proofs.HeroCaseWeberCompanionObstruction
import Proofs.HeroCaseGGWeightedCorrectionWaypoint
import Proofs.HeroCaseWeberClassInvariantBridgeWaypoint
import Proofs.HeroCaseHeineCor2cf
import Proofs.HeroCaseLocal
import Proofs.HeroCasePage43Equivalence
import Proofs.HeroCaseSubsequenceExact
import Proofs.HeroCaseTailOperatorWaypoint
import Proofs.RationalEquivalence

open Polynomial

namespace Proofs
namespace HeroCase
namespace AwardTrack

/-!
Award-track target module.

This file is intentionally a *scaffold*:

- It is the canonical place where the final closed-form identity for the hero case
  (`cb60fd71d1d7`, reduced variable `t = q^3`) will be stated and proved.
- It should compile without `sorry` so `lake build` stays meaningful.

Once a concrete closed form is identified (q-product / eta-quotient / theta ratio, etc.),
replace `finalIdentityStatement` with the actual theorem statement and build out the proof.

Current source-family-specific exact lane:

- `Proofs/HeroCasePage43Equivalence.lean` now formalizes the zero-shift
  `f2/gcf3` and `f4/gcf2` `n`-dependent equivalence obstructions
- the same module also now covers the first nearby unit-`a`-shift,
  unit-`b`-shift, mixed unit-`a` / unit-`b`-shift, mixed unit-`a` /
  unit-`lambda`-shift, mixed unit-`b` / unit-`lambda`-shift, mixed
  unit-`a` / unit-`b` / unit-`lambda`-shift, and unit-`lambda`-shift lanes
- those eight nearby cases are now also packaged there as the nearest-shift
  cube summary theorems `noNearestShiftCubeF2ExactEquivalence`,
  `noNearestShiftCubeF4ExactEquivalence`, and
  `noNearestShiftCubeExactEquivalence`
- the same file now also exposes that exact layer as Bool-parameterized
  theorems over the shift bits:
  `noNearestShiftCubeF2ExactEquivalenceFor`,
  `noNearestShiftCubeF4ExactEquivalenceFor`, and
  `noNearestShiftCubeExactEquivalenceFor`
- the same module also now excludes the zero-shift polynomial
  single-prefactor sub-box `phi in {1, 1 + t}` with at most one non-plain
  prefactor active in both page-43 families via
  `noZeroShiftPolynomialSinglePrefactorF2DirectMatches`,
  `noZeroShiftPolynomialSinglePrefactorF4DirectMatches`, and
  `noZeroShiftPolynomialSinglePrefactorDirectMatches`
- the same module also now excludes the zero-shift reciprocal
  single-prefactor sub-box `phi in {1 / (1 + t)}` in both page-43 families via
  the cross-multiplied theorems
  `noZeroShiftReciprocalSinglePrefactorF2CrossMatches`,
  `noZeroShiftReciprocalSinglePrefactorF4CrossMatches`, and
  `noZeroShiftReciprocalSinglePrefactorCrossMatches`
- the unit-`a` shifts still end with the same surviving nonzero `m^2`
  coefficient after the forced parameter specializations
- the unit-`b` shifts also end with that same surviving nonzero `m^2`
  coefficient after the forced parameter specializations
- the first mixed unit-`a` / unit-`b` shifts also end with that same surviving
  nonzero `m^2` coefficient after the forced parameter specializations
- the first mixed unit-`a` / unit-`lambda` shifts fail earlier, where the
  surviving `m^1` coefficient already has no constant-parameter solution
- the first mixed unit-`b` / unit-`lambda` shifts also fail earlier, with the
  same surviving `m^1` coefficient obstruction
- the first mixed unit-`a` / unit-`b` / unit-`lambda` shifts also fail earlier,
  with that same surviving `m^1` coefficient obstruction
- the unit-`lambda` shifts fail earlier, where the forced `m^1` coefficient
  already has no constant-parameter solution
- the zero-shift lanes still end with a surviving nonzero `m^2` coefficient
  after the forced parameter specializations
- so this scaffold should currently treat those lanes as exclusion theorems,
  not as positive final-identity targets
- `Proofs/RationalEquivalence.lean` already proves the finite convergents of the
  reduced-by-factor model agree exactly with the unreduced hero-case convergents
  over `RatFunc Rat`
- `Proofs/HeroCaseSubsequenceExact.lean` already proves the nearest RR/cubic
  arithmetic-subsequence source lanes are excluded exactly, without the old
  stride/sample-point bounds
- `Proofs/HeroCaseLocal.lean` already proves the direct RR/cubic Bauer-Muir and
  simple cubic contraction mismatches at the decisive early stages
- `Proofs/HeroCaseHeineCor2cf.lean` already proves exact low-stage Heine
  `cor2cf` branch mismatches in the forced `a = 0` lane
- `Proofs/HeroCaseGGQuotientCoordinateObstruction.lean` now packages the
  current tail-family-first `GG` quotient-coordinate witness table for
  `Q_3 = GG(t^3)/GG(t)` and `Q_4 = GG(t^4)/GG(t)` as a Lean waypoint shell,
  including a sample-indexed theorem family asserting each currently sampled
  tail object already fails at the recorded first exact `Q_3` and `Q_4`
  residual term, together with the shared leading-obstruction compression
  `q3Coeff = 2*rho`, `q4Coeff = 3*rho`
- `Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean` now packages the next
  hero-ratio diagnostic layer after `W_34 = Q_3^3 / Q_4^2`:
  the recorded first failure of `F / W_34`, the current normalized follow-up
  `G_W34`, the deeper follow-up `G2_W34`, and the checked small eta /
  modular-unit / RR-GG source-family correction boxes at both normalized
  layers
- `Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean` now packages the current
  unified Morton named-coordinate stack:
  the square lane `X_mt = F^2`, the transported square follow-up
  `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`,
  the Weber-Schlafli coordinate `P_ws = (1/F - F) / 2`,
  and the direct Weber companion
  `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`;
  at the current waypoint `X_mt`, `T_mt`, and `B_ws` collapse to universal
  obstruction classes, while `P_ws` remains a small exact witness table across
  the `12` sampled tail-family objects
- `Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean` still packages the
  theorem-shaped exact shell for the explicit square lane `X_mt = F^2`,
  recording that the Proposition 3.2 template
  `X_2^2 - (X^2 - 4*X + 1)*X_2 + X^2`
  still misses across all `12` sampled tail-family objects and collapses to
  one universal constant-term obstruction class `(t^0, 4)`
- `Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean` now packages the
  next focused hand-off after the Ramanujan-Weber class-invariant compression:
  the eta-side residual `G_g12_ws` is the current primary residual, the
  plus-side residual `G_p12_ws` is kept as its algebraically constrained
  companion through the exact coordinate bridge, the classical Weber `f2`
  tri-product coordinate
  `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2)` plus its
  normalized follow-up `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)` also stay
  outside the first checked theorem-shaped closure boxes, the same shell now
  records the source-backed classical Weber trio reading coming from
  Berndt--Chan--Zhang plus Yui--Zagier, and the focused quotient
  `R_gp_ws = G_p12_ws / G_g12_ws` plus the normalized follow-up
  `H_gp_ws = (R_gp_ws - 1) / (96*t^3)` still have `0` hits in the first
  checked self-polynomial / self-fractional-linear / self-quotient-product /
  eta / modular-unit / plus-Pochhammer boxes;
  the same shell now also upgrades the later return-bridge `Q_XK_ws` /
  `L_XK_ws` lane from a pure hit-count shell to an exact Chan--Huang
  obstruction package, recording the first direct and quotient
  modular-equation failure witnesses against `GG3`, `GG4`, `Q_3`, and `Q_4`
- `Proofs/HeroCaseTailOperatorWaypoint.lean` now records the current sampled
  tail-operator lane (`12` sampled tail objects, Mahler-style moduli `2,3`,
  recurrence depths `2,3`, and `t`-degrees `1,2,3`) as a proof-workspace
  landing spot for a later exact operator / factorization theorem
- this scaffold now packages that rational-function equivalence layer together
  with the page-43 nearest-shift-cube exclusions, the exact zero-shift
  single-prefactor exclusions, and the exact RR/cubic
  arithmetic-subsequence obstructions, plus these exact local mismatch layers,
  as the current theorem-grade award-track waypoint
-/

-- AwardTrackStatus: exclusion_waypoint

abbrev QPoly := Polynomial Rat
abbrev QRatFunc := RatFunc Rat

noncomputable section

def heroConvergentRatFunc (n : Nat) : QRatFunc :=
  (continuantNum heroDataRatFunc n) / (continuantDen heroDataRatFunc n)

def reducedHeroConvergentRatFunc (n : Nat) : QRatFunc :=
  Proofs.convergent reducedHeroDataRatFunc n

theorem heroDataRat_eq_heroDataRatFunc :
    Proofs.RationalEquivalence.heroDataRat = heroDataRatFunc := rfl

theorem reducedHeroDataRat_eq_reducedHeroDataRatFunc :
    Proofs.RationalEquivalence.reducedHeroDataRat = reducedHeroDataRatFunc := rfl

theorem reverseScale_eq_reverseScaleRatFunc :
    Proofs.RationalEquivalence.reverseScale = reverseScaleRatFunc := by
  funext k
  cases k <;> rfl

theorem heroConvergentRatFunc_eq_reducedHeroConvergentRatFunc (n : Nat) :
    heroConvergentRatFunc n = reducedHeroConvergentRatFunc n := by
  simpa [heroConvergentRatFunc, reducedHeroConvergentRatFunc,
    heroDataRat_eq_heroDataRatFunc, reducedHeroDataRat_eq_reducedHeroDataRatFunc]
    using Proofs.RationalEquivalence.hero_convergent_eq_reduced_convergent n

theorem reverseEquivalenceRecoversHeroData :
    Proofs.equivalenceTransform reducedHeroDataRatFunc reverseScaleRatFunc = heroDataRatFunc := by
  simpa [heroDataRat_eq_heroDataRatFunc, reducedHeroDataRat_eq_reducedHeroDataRatFunc,
    reverseScale_eq_reverseScaleRatFunc]
    using Proofs.RationalEquivalence.equivalenceTransform_reducedHeroDataRat_eq

def page43NearestShiftCubeExcludedProp
    (shiftA shiftB shiftLambda : Bool) : Prop :=
    (¬ ∃ a b lam : Rat,
      Proofs.HeroCase.Page43Equivalence.f2NearestShiftCubeM1 shiftA shiftB shiftLambda lam = 0 ∧
      Proofs.HeroCase.Page43Equivalence.f2NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      Proofs.HeroCase.Page43Equivalence.f2NearestShiftCubeM3 shiftA shiftB shiftLambda a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      Proofs.HeroCase.Page43Equivalence.f4NearestShiftCubeM0 shiftA shiftB shiftLambda a = 0 ∧
      Proofs.HeroCase.Page43Equivalence.f4NearestShiftCubeM1 shiftA shiftB shiftLambda a lam = 0 ∧
      Proofs.HeroCase.Page43Equivalence.f4NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      Proofs.HeroCase.Page43Equivalence.f4NearestShiftCubeM3 shiftA shiftB shiftLambda b = 0)

theorem page43NearestShiftCubeExcluded
    (shiftA shiftB shiftLambda : Bool) :
    page43NearestShiftCubeExcludedProp shiftA shiftB shiftLambda := by
  exact Proofs.HeroCase.Page43Equivalence.noNearestShiftCubeExactEquivalenceFor
    shiftA shiftB shiftLambda

def page43PolynomialPrefactorExcludedProp : Prop :=
    ((¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlainPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlainPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlainPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlainPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAPlusPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAPlusPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAPlusPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAPlusPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftBPlusPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftBPlusPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftBPlusPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftBPlusPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlusLambdaPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlusLambdaPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlusLambdaPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftPlusLambdaPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ)) ∧
    ((¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlainPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlainPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlainPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlainPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAPlusPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAPlusPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAPlusPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAPlusPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftBPlusPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftBPlusPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftBPlusPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftBPlusPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlusLambdaPrefactorStage1A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage1NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlusLambdaPrefactorStage1B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage1DenominatorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlusLambdaPrefactorStage2A a b lam =
          Proofs.HeroCase.Page43Equivalence.heroStage2NumeratorQ ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftPlusLambdaPrefactorStage2B a b =
          Proofs.HeroCase.Page43Equivalence.heroStage2DenominatorQ))

theorem page43PolynomialPrefactorExcluded :
    page43PolynomialPrefactorExcludedProp := by
  exact Proofs.HeroCase.Page43Equivalence.noZeroShiftPolynomialSinglePrefactorDirectMatches

def page43ReciprocalPrefactorExcludedProp : Prop :=
    ((¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftAReciprocalPrefactorNumResidual a b lam = 0) ∧
      (¬ ∃ a b : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
        Proofs.HeroCase.Page43Equivalence.f2ZeroShiftLambdaReciprocalPrefactorNumResidual a b lam = 0)) ∧
    ((¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftAReciprocalPrefactorNumResidual a lam = 0) ∧
      (¬ ∃ a b : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
      (¬ ∃ a b lam : Rat,
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
        Proofs.HeroCase.Page43Equivalence.f4ZeroShiftLambdaReciprocalPrefactorNumResidual a lam = 0))

theorem page43ReciprocalPrefactorExcluded :
    page43ReciprocalPrefactorExcludedProp := by
  exact Proofs.HeroCase.Page43Equivalence.noZeroShiftReciprocalSinglePrefactorCrossMatches

def nearestArithmeticSubsequenceSourcesExcludedProp
    (stride offset : Nat) : Prop :=
  (¬ (continuantNum rrReciprocalData offset = continuantDen rrReciprocalData offset ∧
      continuantNum rrReciprocalData (offset + stride) =
        (1 + X) * continuantDen rrReciprocalData (offset + stride))) ∧
  (¬ (continuantNum cubicReciprocalData offset = continuantDen cubicReciprocalData offset ∧
      continuantNum cubicReciprocalData (offset + stride) =
        (1 + X) * continuantDen cubicReciprocalData (offset + stride)))

theorem nearestArithmeticSubsequenceSourcesExcluded
    (stride offset : Nat) (hstride : 2 ≤ stride) (hoff : offset < stride) :
    nearestArithmeticSubsequenceSourcesExcludedProp stride offset := by
  exact ⟨
    Proofs.HeroCase.SubsequenceExact.noRRArithmeticSubsequenceContraction
      stride offset hstride hoff,
    Proofs.HeroCase.SubsequenceExact.noCubicArithmeticSubsequenceContraction
      stride offset hstride hoff
  ⟩

def directLocalObstructionsProp : Prop :=
  rrDirectFirstNumeratorPoly ≠ heroData.a 0 ∧
  cubicDirectSecondNumeratorPoly ≠ heroData.a 1 ∧
  cubicOddInitialDenominatorPoly ≠ heroData.b0 ∧
  cubicEvenFirstDenominatorPoly ≠ heroData.b 0

theorem directLocalObstructions : directLocalObstructionsProp := by
  exact ⟨
    Proofs.HeroCase.noDirectRRBauerMuirMatchAtStage1,
    Proofs.HeroCase.noDirectCubicBauerMuirMatchAtStage2,
    Proofs.HeroCase.noSimpleCubicOddContraction,
    Proofs.HeroCase.noSimpleCubicEvenContraction
  ⟩

def simpleCor2cfBranchesExcludedProp (b lam : Int) : Prop :=
  Proofs.HeroCase.HeineCor2cf.aZeroOddInitialDenominatorPoly lam ≠ heroData.b0 ∧
  Proofs.HeroCase.HeineCor2cf.aZeroEvenFirstNumeratorPoly lam ≠ heroData.a 0 ∧
  Proofs.HeroCase.HeineCor2cf.aZeroEvenOddInitialNumeratorPoly b lam ≠
    Proofs.HeroCase.HeineCor2cf.aZeroEvenOddInitialDenominatorPoly b lam ∧
  Proofs.HeroCase.HeineCor2cf.aZeroEvenEvenFirstNumeratorPoly b lam ≠ heroData.a 0

theorem simpleCor2cfBranchesExcluded (b lam : Int) (h : lam ≠ 0) :
    simpleCor2cfBranchesExcludedProp b lam := by
  exact ⟨
    Proofs.HeroCase.HeineCor2cf.noSimpleCor2cfOddBranch lam h,
    Proofs.HeroCase.HeineCor2cf.noSimpleCor2cfEvenBranch lam,
    Proofs.HeroCase.HeineCor2cf.noSimpleCor2cfOddOfEvenBranch b lam h,
    Proofs.HeroCase.HeineCor2cf.noSimpleCor2cfEvenOfEvenBranch b lam
  ⟩

def finiteConvergentReductionWaypoint : Prop :=
  (∀ n : Nat, heroConvergentRatFunc n = reducedHeroConvergentRatFunc n) ∧
  (Proofs.equivalenceTransform reducedHeroDataRatFunc reverseScaleRatFunc = heroDataRatFunc)

theorem finiteConvergentReductionWaypoint_true : finiteConvergentReductionWaypoint := by
  refine ⟨?_, ?_⟩
  · intro n
    exact heroConvergentRatFunc_eq_reducedHeroConvergentRatFunc n
  · exact reverseEquivalenceRecoversHeroData

def knownSourceOrbitExclusionWaypoint : Prop :=
  (∀ shiftA shiftB shiftLambda : Bool, page43NearestShiftCubeExcludedProp shiftA shiftB shiftLambda) ∧
  page43PolynomialPrefactorExcludedProp ∧
  page43ReciprocalPrefactorExcludedProp ∧
  (∀ stride offset : Nat, 2 ≤ stride → offset < stride →
    nearestArithmeticSubsequenceSourcesExcludedProp stride offset) ∧
  directLocalObstructionsProp ∧
  (∀ b lam : Int, lam ≠ 0 → simpleCor2cfBranchesExcludedProp b lam) ∧
  Proofs.HeroCase.GGQ34.currentWaypoint ∧
  Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint ∧
  Proofs.HeroCase.MortonNamedCoordinate.currentWaypoint ∧
  Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint

theorem knownSourceOrbitExclusionWaypoint_true : knownSourceOrbitExclusionWaypoint := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro shiftA shiftB shiftLambda
    exact page43NearestShiftCubeExcluded shiftA shiftB shiftLambda
  · exact page43PolynomialPrefactorExcluded
  · exact page43ReciprocalPrefactorExcluded
  · intro stride offset hstride hoff
    exact nearestArithmeticSubsequenceSourcesExcluded stride offset hstride hoff
  · exact directLocalObstructions
  · intro b lam h
    exact simpleCor2cfBranchesExcluded b lam h
  · exact Proofs.HeroCase.GGQ34.currentWaypoint_true
  · exact Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint_true
  · exact Proofs.HeroCase.MortonNamedCoordinate.currentWaypoint_true
  · exact Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint_true

def exactWaypointStatement : Prop :=
  finiteConvergentReductionWaypoint ∧ knownSourceOrbitExclusionWaypoint

structure ExactWaypointCertificate where
  reduction : finiteConvergentReductionWaypoint
  exclusions : knownSourceOrbitExclusionWaypoint

def currentExactWaypointCertificate : ExactWaypointCertificate where
  reduction := finiteConvergentReductionWaypoint_true
  exclusions := knownSourceOrbitExclusionWaypoint_true

theorem currentExactWaypointCertificate_sound : exactWaypointStatement := by
  exact ⟨
    currentExactWaypointCertificate.reduction,
    currentExactWaypointCertificate.exclusions
  ⟩

theorem exactWaypointStatement_true : exactWaypointStatement := by
  exact currentExactWaypointCertificate_sound

def tailOperatorResearchWaypoint : Prop :=
  Proofs.HeroCase.TailOperator.currentWaypoint

theorem tailOperatorResearchWaypoint_true : tailOperatorResearchWaypoint := by
  exact Proofs.HeroCase.TailOperator.currentWaypoint_true

def ggWeightedCorrectionResearchWaypoint : Prop :=
  Proofs.HeroCase.GGWeightedCorrection.currentWaypoint

theorem ggWeightedCorrectionResearchWaypoint_true : ggWeightedCorrectionResearchWaypoint := by
  exact Proofs.HeroCase.GGWeightedCorrection.currentWaypoint_true

def mortonSquaredCoordinateExcludedProp : Prop :=
  Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint

theorem mortonSquaredCoordinateExcluded : mortonSquaredCoordinateExcludedProp := by
  exact Proofs.HeroCase.MortonSquaredCoordinate.currentWaypoint_true

def mortonNamedCoordinateResearchWaypoint : Prop :=
  Proofs.HeroCase.MortonNamedCoordinate.currentWaypoint

theorem mortonNamedCoordinateResearchWaypoint_true :
    mortonNamedCoordinateResearchWaypoint := by
  exact Proofs.HeroCase.MortonNamedCoordinate.currentWaypoint_true

def weberClassInvariantBridgeResearchWaypoint : Prop :=
  Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint

theorem weberClassInvariantBridgeResearchWaypoint_true :
    weberClassInvariantBridgeResearchWaypoint := by
  exact Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint_true

def namedWeberOrbitResearchWaypoint : Prop :=
  Proofs.HeroCase.GGQ34.currentWaypoint ∧
    mortonNamedCoordinateResearchWaypoint ∧
    Proofs.HeroCase.GGWeightedCorrection.currentWaypoint ∧
    Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint

structure NamedWeberOrbitResearchCertificate where
  ggQ34 : Proofs.HeroCase.GGQ34.currentWaypoint
  mortonNamedCoordinate : mortonNamedCoordinateResearchWaypoint
  ggWeightedCorrection : Proofs.HeroCase.GGWeightedCorrection.currentWaypoint
  weberClassInvariantBridge : Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint

def currentNamedWeberOrbitResearchCertificate : NamedWeberOrbitResearchCertificate where
  ggQ34 := Proofs.HeroCase.GGQ34.currentWaypoint_true
  mortonNamedCoordinate := mortonNamedCoordinateResearchWaypoint_true
  ggWeightedCorrection := Proofs.HeroCase.GGWeightedCorrection.currentWaypoint_true
  weberClassInvariantBridge := Proofs.HeroCase.WeberClassInvariantBridge.currentWaypoint_true

theorem currentNamedWeberOrbitResearchCertificate_sound :
    namedWeberOrbitResearchWaypoint := by
  exact ⟨
    currentNamedWeberOrbitResearchCertificate.ggQ34,
    currentNamedWeberOrbitResearchCertificate.mortonNamedCoordinate,
    currentNamedWeberOrbitResearchCertificate.ggWeightedCorrection,
    currentNamedWeberOrbitResearchCertificate.weberClassInvariantBridge
  ⟩

theorem namedWeberOrbitResearchWaypoint_true :
    namedWeberOrbitResearchWaypoint := by
  exact currentNamedWeberOrbitResearchCertificate_sound

def currentRecognitionFrontierWaypoint : Prop :=
  exactWaypointStatement ∧ namedWeberOrbitResearchWaypoint

structure RecognitionFrontierCertificate where
  exactWaypoint : exactWaypointStatement
  namedWeberOrbit : namedWeberOrbitResearchWaypoint

def currentRecognitionFrontierCertificate : RecognitionFrontierCertificate where
  exactWaypoint := exactWaypointStatement_true
  namedWeberOrbit := namedWeberOrbitResearchWaypoint_true

theorem currentRecognitionFrontierCertificate_sound :
    currentRecognitionFrontierWaypoint := by
  exact ⟨
    currentRecognitionFrontierCertificate.exactWaypoint,
    currentRecognitionFrontierCertificate.namedWeberOrbit
  ⟩

theorem currentRecognitionFrontierWaypoint_true :
    currentRecognitionFrontierWaypoint := by
  exact currentRecognitionFrontierCertificate_sound

/-!
TODO (award track, endgame):

1. Define the closed form target in a formal object (likely a `PowerSeries Rat` / analytic function).
2. Prove the closed form satisfies the same characterization as the continued fraction value
   (functional equation / modular equation / uniqueness lemma).
3. Bridge the current exact waypoint from finite convergents and exact exclusion
   theorems to the infinite object.
4. Generalize the new page-43 exact obstruction layer beyond
   the currently formalized nearest-shift cube plus the exact zero-shift
   single-prefactor box for `f2/gcf3` / `f4/gcf2`,
   especially if shifted or multi-prefactor rational lanes become relevant.
5. Upgrade the current `GG` quotient-coordinate waypoint shell into actual
   exact modular-equation obstruction theorems, and either absorb or replace
   the current weighted-correction research waypoint with an exact theorem.
6. Upgrade the current Weber companion shell and the remaining non-exact parts
   of the Weber class-invariant bridge shell into exact modular-function
   obstruction theorems or a positive source bridge, depending on which side
   of the same named orbit yields the first true uniqueness theorem.
7. Prove the final identity in Lean and remove this placeholder.
-/

def finalIdentityStatement : Prop :=
  False

end

end AwardTrack
end HeroCase
end Proofs
