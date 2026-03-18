import Mathlib
import Proofs.HeroCaseObjects

open Polynomial

namespace Proofs
namespace HeroCase
namespace Page43Equivalence

/-!
Exact page-43 obstruction layers around the current hero-case source lanes.

These theorems formalize the current source-family-specific exact lanes produced by
the Python research pipeline:

* `f2 / gcf3` at zero shift
* `f4 / gcf2` at zero shift
* the zero-shift polynomial single-prefactor sub-box `phi in {1, 1 + t}`
* the zero-shift reciprocal single-prefactor sub-box `phi in {1 / (1 + t)}`
* the nearest unit-`a` shift lanes
* the nearest unit-`b` shift lanes
* the first mixed unit-`a` / unit-`b` shift lanes
* the first mixed unit-`a` / unit-`lambda` shift lanes
* the first mixed unit-`b` / unit-`lambda` shift lanes
* the first mixed unit-`a` / unit-`b` / unit-`lambda` shift lanes
* the nearest unit-`lambda` shift lanes

For the `n`-dependent equivalence lanes, the necessary residual identity is expanded
as a polynomial in `m = t^(n-1)` with coefficients in `Rat[X]`.  Vanishing of all
coefficients would be necessary for an arbitrary `n`-dependent equivalence
transformation to match the hero-case reciprocal.  We prove those coefficient
conditions are inconsistent.  The zero-shift single-prefactor box is
handled one level earlier: denominator matching either already fails or forces
the parameters, then stage-1 or stage-2 numerators disagree exactly; the reciprocal
lanes are recorded via cross-multiplied coefficient identities.  The final `noNearestShiftCube...`
theorems package the full `{0,1}^3` nearest-shift cube for the page-43 audit, and
the `...For` variants expose the same layer as a Bool-parameterized theorem family
over the shift bits.
-/

abbrev QPoly := Polynomial Rat

noncomputable section

def heroStage1NumeratorQ : QPoly :=
  C (1 : Rat) * X + C (1 : Rat) * X ^ 2

def heroStage2NumeratorQ : QPoly :=
  C (1 : Rat) * X ^ 2 + C (1 : Rat) * X ^ 4

def heroStage1DenominatorQ : QPoly :=
  1 + C (1 : Rat) * X

def heroStage2DenominatorQ : QPoly :=
  1 + C (1 : Rat) * X ^ 2

def f2ZeroShiftPlainPrefactorStage1A (a b lam : Rat) : QPoly :=
  C lam * X + C (-(a * b)) * X ^ 2

def f2ZeroShiftPlainPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C b * X + C a * X ^ 2

def f2ZeroShiftPlainPrefactorStage2A (a b lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-(a * b)) * X ^ 4

def f2ZeroShiftPlainPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C b * X ^ 2 + C a * X ^ 3

def f2ZeroShiftAPlusPrefactorStage1A (a b lam : Rat) : QPoly :=
  C lam * X + C (-(a * b)) * X ^ 2 + C (-(a * b)) * X ^ 3

def f2ZeroShiftAPlusPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C b * X + C a * X ^ 2 + C a * X ^ 3

def f2ZeroShiftAPlusPrefactorStage2A (a b lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-(a * b)) * X ^ 4 + C (-(a * b)) * X ^ 5

def f2ZeroShiftAPlusPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C b * X ^ 2 + C a * X ^ 3 + C a * X ^ 4

def f2ZeroShiftBPlusPrefactorStage1A (a b lam : Rat) : QPoly :=
  C lam * X + C (-(a * b)) * X ^ 2 + C (-(a * b)) * X ^ 3

def f2ZeroShiftBPlusPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C b * X + C (a + b) * X ^ 2

def f2ZeroShiftBPlusPrefactorStage2A (a b lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-(a * b)) * X ^ 4 + C (-(a * b)) * X ^ 5

def f2ZeroShiftBPlusPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C (a + b) * X ^ 3 + C b * X ^ 2

def f2ZeroShiftPlusLambdaPrefactorStage1A (a b lam : Rat) : QPoly :=
  C lam * X + C (lam - a * b) * X ^ 2

def f2ZeroShiftPlusLambdaPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C b * X + C a * X ^ 2

def f2ZeroShiftPlusLambdaPrefactorStage2A (a b lam : Rat) : QPoly :=
  C lam * X ^ 2 + C lam * X ^ 3 + C (-(a * b)) * X ^ 4

def f2ZeroShiftPlusLambdaPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C b * X ^ 2 + C a * X ^ 3

def f4ZeroShiftPlainPrefactorStage1A (a _b lam : Rat) : QPoly :=
  C (a + lam) * X

def f4ZeroShiftPlainPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C (b - a) * X

def f4ZeroShiftPlainPrefactorStage2A (a _b lam : Rat) : QPoly :=
  C a * X + C lam * X ^ 2

def f4ZeroShiftPlainPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C (-a) * X + C b * X ^ 2

def f4ZeroShiftAPlusPrefactorStage1A (a _b lam : Rat) : QPoly :=
  C (a + lam) * X + C a * X ^ 2

def f4ZeroShiftAPlusPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C (b - a) * X + C (-a) * X ^ 2

def f4ZeroShiftAPlusPrefactorStage2A (a _b lam : Rat) : QPoly :=
  C a * X + C (a + lam) * X ^ 2

def f4ZeroShiftAPlusPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C (-a) * X + C (b - a) * X ^ 2

def f4ZeroShiftBPlusPrefactorStage1A (a _b lam : Rat) : QPoly :=
  C (a + lam) * X

def f4ZeroShiftBPlusPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C (b - a) * X + C b * X ^ 2

def f4ZeroShiftBPlusPrefactorStage2A (a _b lam : Rat) : QPoly :=
  C a * X + C lam * X ^ 2

def f4ZeroShiftBPlusPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C (-a) * X + C b * X ^ 2 + C b * X ^ 3

def f4ZeroShiftPlusLambdaPrefactorStage1A (a _b lam : Rat) : QPoly :=
  C (a + lam) * X + C lam * X ^ 2

def f4ZeroShiftPlusLambdaPrefactorStage1B (a b : Rat) : QPoly :=
  1 + C (b - a) * X

def f4ZeroShiftPlusLambdaPrefactorStage2A (a _b lam : Rat) : QPoly :=
  C a * X + C lam * X ^ 2 + C lam * X ^ 3

def f4ZeroShiftPlusLambdaPrefactorStage2B (a b : Rat) : QPoly :=
  1 + C (-a) * X + C b * X ^ 2

def f2ZeroShiftAReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C (b - 1) + C (a + b - 1) * X

def f2ZeroShiftAReciprocalPrefactorNumResidual (a b lam : Rat) : QPoly :=
  C (lam - 1) + C (lam - a * b - 2) * X + C (-1 : Rat) * X ^ 2

def f2ZeroShiftBReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C (b - 1) + C (a - 1) * X + C a * X ^ 2

def f2ZeroShiftLambdaReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C (b - 1) + C a * X

def f2ZeroShiftLambdaReciprocalPrefactorNumResidual (a b lam : Rat) : QPoly :=
  C (lam - 1) + C (-(a * b) - 2) * X + C (-(a * b) - 1) * X ^ 2

def f4ZeroShiftAReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C (-a) + C (b - 1) * X + C (b - 1) * X ^ 2

def f4ZeroShiftAReciprocalPrefactorNumResidual (a lam : Rat) : QPoly :=
  C (a + lam - 1) + C (lam - 2) * X + C (-1 : Rat) * X ^ 2

def f4ZeroShiftBReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C a + C (a - b + 1) * X + X ^ 2

def f4ZeroShiftLambdaReciprocalPrefactorDenResidual (a b : Rat) : QPoly :=
  C (-a) + C (b - 1) * X

def f4ZeroShiftLambdaReciprocalPrefactorNumResidual (a lam : Rat) : QPoly :=
  C (a + lam - 1) + C (a - 2) * X + C (-1 : Rat) * X ^ 2

def f2ZeroShiftM1 (lam : Rat) : QPoly :=
  C (lam - 1) * X

def f2ZeroShiftM2 (a b lam : Rat) : QPoly :=
  C (-a) * X ^ 3 + C (-(a * b + a + b)) * X ^ 2 + C (lam - b) * X

def f2ZeroShiftM3 (a b : Rat) : QPoly :=
  C (-(a * a)) * X ^ 4 + C (-(2 * a * b)) * X ^ 3 + C (-(a * b + b * b)) * X ^ 2

def f4ZeroShiftM0 (a : Rat) : QPoly :=
  C a * X

def f4ZeroShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 3 + C (2 * a) * X ^ 2 + C (a + lam - 1) * X

def f4ZeroShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 3 + C (a * b - b) * X ^ 2 + C (lam - b) * X

def f4ZeroShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 2

def f2UnitLambdaShiftM1 (lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-1) * X

def f2UnitLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (-a) * X ^ 3 + C (-(a * b + a + b) + lam) * X ^ 2 + C (-b) * X

def f2UnitLambdaShiftM3 (a b : Rat) : QPoly :=
  f2ZeroShiftM3 a b

def f2UnitAShiftM1 (lam : Rat) : QPoly :=
  C (lam - 1) * X

def f2UnitAShiftM2 (a b lam : Rat) : QPoly :=
  C (-a) * X ^ 4 + C (-(a * b + a)) * X ^ 3 + C (-b) * X ^ 2 + C (lam - b) * X

def f2UnitAShiftM3 (a b : Rat) : QPoly :=
  C (-(a * a)) * X ^ 6 + C (-(2 * a * b)) * X ^ 4 + C (-(a * b)) * X ^ 3 + C (-(b * b)) * X ^ 2

def f2UnitBShiftM1 (lam : Rat) : QPoly :=
  f2ZeroShiftM1 lam

def f2UnitBShiftM2 (a b lam : Rat) : QPoly :=
  C (-(a * b + a + b)) * X ^ 3 + C (-(a + b)) * X ^ 2 + C lam * X

def f2UnitBShiftM3 (a b : Rat) : QPoly :=
  C (-(a * b)) * X ^ 3 + C (-(a * a + 2 * a * b + b * b)) * X ^ 4

def f4UnitLambdaShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 3 + C (2 * a + lam) * X ^ 2 + C (a - 1) * X

def f4UnitLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 3 + C (a * b - b + lam) * X ^ 2 + C (-b) * X

def f4UnitLambdaShiftM3 (b : Rat) : QPoly :=
  f4ZeroShiftM3 b

def f4UnitAShiftM0 (a : Rat) : QPoly :=
  C a * X ^ 2

def f4UnitAShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 5 + C (2 * a) * X ^ 3 + C a * X ^ 2 + C (lam - 1) * X

def f4UnitAShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 4 + C (a * b) * X ^ 3 + C (-b) * X ^ 2 + C (lam - b) * X

def f4UnitAShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 2

def f4UnitBShiftM0 (a : Rat) : QPoly :=
  f4ZeroShiftM0 a

def f4UnitBShiftM1 (a lam : Rat) : QPoly :=
  f4ZeroShiftM1 a lam

def f4UnitBShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 4 + C (a * b - b) * X ^ 3 + C (-b) * X ^ 2 + C lam * X

def f4UnitBShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 4

def f2UnitABShiftM1 (lam : Rat) : QPoly :=
  f2ZeroShiftM1 lam

def f2UnitABShiftM2 (a b lam : Rat) : QPoly :=
  C (-(a * b + a)) * X ^ 4 + C (-(a + b)) * X ^ 3 + C (-b) * X ^ 2 + C lam * X

def f2UnitABShiftM3 (a b : Rat) : QPoly :=
  C (-(a * a)) * X ^ 6 + C (-(2 * a * b)) * X ^ 5 + C (-(a * b + b * b)) * X ^ 4

def f4UnitABShiftM0 (a : Rat) : QPoly :=
  C a * X ^ 2

def f4UnitABShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 5 + C (2 * a) * X ^ 3 + C a * X ^ 2 + C (lam - 1) * X

def f4UnitABShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 5 + C (a * b) * X ^ 4 + C (-b) * X ^ 3 + C (-b) * X ^ 2 + C lam * X

def f4UnitABShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 4

def f2UnitALambdaShiftM1 (lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-1) * X

def f2UnitALambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (-a) * X ^ 4 + C (-(a * b + a)) * X ^ 3 + C (-b + lam) * X ^ 2 + C (-b) * X

def f2UnitALambdaShiftM3 (a b : Rat) : QPoly :=
  C (-(a * a)) * X ^ 6 + C (-(2 * a * b)) * X ^ 4 + C (-(a * b)) * X ^ 3 + C (-(b * b)) * X ^ 2

def f4UnitALambdaShiftM0 (a : Rat) : QPoly :=
  C a * X ^ 2

def f4UnitALambdaShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 5 + C (2 * a + lam) * X ^ 3 + C a * X ^ 2 + C (-1) * X

def f4UnitALambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 4 + C (a * b + lam) * X ^ 3 + C (-b) * X ^ 2 + C (-b) * X

def f4UnitALambdaShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 2

def f2UnitBLambdaShiftM1 (lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-1) * X

def f2UnitBLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (-(a * b + a + b)) * X ^ 3 + C (-(a + b) + lam) * X ^ 2

def f2UnitBLambdaShiftM3 (a b : Rat) : QPoly :=
  C (-(a * b)) * X ^ 3 + C (-(a * a + 2 * a * b + b * b)) * X ^ 4

def f4UnitBLambdaShiftM0 (a : Rat) : QPoly :=
  C a * X

def f4UnitBLambdaShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 3 + C (2 * a + lam) * X ^ 2 + C (a - 1) * X

def f4UnitBLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 4 + C (a * b - b + lam) * X ^ 3 + C (-b) * X ^ 2

def f4UnitBLambdaShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 4

def f2UnitABLambdaShiftM1 (lam : Rat) : QPoly :=
  C lam * X ^ 2 + C (-1) * X

def f2UnitABLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (-(a * b + a)) * X ^ 4 + C (-(a + b)) * X ^ 3 + C (-b + lam) * X ^ 2

def f2UnitABLambdaShiftM3 (a b : Rat) : QPoly :=
  f2UnitABShiftM3 a b

def f4UnitABLambdaShiftM0 (a : Rat) : QPoly :=
  C a * X ^ 2

def f4UnitABLambdaShiftM1 (a lam : Rat) : QPoly :=
  C (-(a * a)) * X ^ 5 + C (2 * a) * X ^ 3 + C (a + lam) * X ^ 2 + C (-1) * X

def f4UnitABLambdaShiftM2 (a b lam : Rat) : QPoly :=
  C (a * b) * X ^ 5 + C (a * b) * X ^ 4 + C (-b) * X ^ 3 + C (-b + lam) * X ^ 2

def f4UnitABLambdaShiftM3 (b : Rat) : QPoly :=
  C (-(b * b)) * X ^ 4

def f2NearestShiftCubeM1 (shiftA shiftB shiftLambda : Bool) (lam : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f2ZeroShiftM1 lam
  | false, false, true => f2UnitLambdaShiftM1 lam
  | false, true, false => f2UnitBShiftM1 lam
  | false, true, true => f2UnitBLambdaShiftM1 lam
  | true, false, false => f2UnitAShiftM1 lam
  | true, false, true => f2UnitALambdaShiftM1 lam
  | true, true, false => f2UnitABShiftM1 lam
  | true, true, true => f2UnitABLambdaShiftM1 lam

def f2NearestShiftCubeM2 (shiftA shiftB shiftLambda : Bool) (a b lam : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f2ZeroShiftM2 a b lam
  | false, false, true => f2UnitLambdaShiftM2 a b lam
  | false, true, false => f2UnitBShiftM2 a b lam
  | false, true, true => f2UnitBLambdaShiftM2 a b lam
  | true, false, false => f2UnitAShiftM2 a b lam
  | true, false, true => f2UnitALambdaShiftM2 a b lam
  | true, true, false => f2UnitABShiftM2 a b lam
  | true, true, true => f2UnitABLambdaShiftM2 a b lam

def f2NearestShiftCubeM3 (shiftA shiftB shiftLambda : Bool) (a b : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f2ZeroShiftM3 a b
  | false, false, true => f2UnitLambdaShiftM3 a b
  | false, true, false => f2UnitBShiftM3 a b
  | false, true, true => f2UnitBLambdaShiftM3 a b
  | true, false, false => f2UnitAShiftM3 a b
  | true, false, true => f2UnitALambdaShiftM3 a b
  | true, true, false => f2UnitABShiftM3 a b
  | true, true, true => f2UnitABLambdaShiftM3 a b

def f4NearestShiftCubeM0 (shiftA shiftB shiftLambda : Bool) (a : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f4ZeroShiftM0 a
  | false, false, true => f4ZeroShiftM0 a
  | false, true, false => f4UnitBShiftM0 a
  | false, true, true => f4UnitBLambdaShiftM0 a
  | true, false, false => f4UnitAShiftM0 a
  | true, false, true => f4UnitALambdaShiftM0 a
  | true, true, false => f4UnitABShiftM0 a
  | true, true, true => f4UnitABLambdaShiftM0 a

def f4NearestShiftCubeM1 (shiftA shiftB shiftLambda : Bool) (a lam : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f4ZeroShiftM1 a lam
  | false, false, true => f4UnitLambdaShiftM1 a lam
  | false, true, false => f4UnitBShiftM1 a lam
  | false, true, true => f4UnitBLambdaShiftM1 a lam
  | true, false, false => f4UnitAShiftM1 a lam
  | true, false, true => f4UnitALambdaShiftM1 a lam
  | true, true, false => f4UnitABShiftM1 a lam
  | true, true, true => f4UnitABLambdaShiftM1 a lam

def f4NearestShiftCubeM2 (shiftA shiftB shiftLambda : Bool) (a b lam : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f4ZeroShiftM2 a b lam
  | false, false, true => f4UnitLambdaShiftM2 a b lam
  | false, true, false => f4UnitBShiftM2 a b lam
  | false, true, true => f4UnitBLambdaShiftM2 a b lam
  | true, false, false => f4UnitAShiftM2 a b lam
  | true, false, true => f4UnitALambdaShiftM2 a b lam
  | true, true, false => f4UnitABShiftM2 a b lam
  | true, true, true => f4UnitABLambdaShiftM2 a b lam

def f4NearestShiftCubeM3 (shiftA shiftB shiftLambda : Bool) (b : Rat) : QPoly :=
  match shiftA, shiftB, shiftLambda with
  | false, false, false => f4ZeroShiftM3 b
  | false, false, true => f4UnitLambdaShiftM3 b
  | false, true, false => f4UnitBShiftM3 b
  | false, true, true => f4UnitBLambdaShiftM3 b
  | true, false, false => f4UnitAShiftM3 b
  | true, false, true => f4UnitALambdaShiftM3 b
  | true, true, false => f4UnitABShiftM3 b
  | true, true, true => f4UnitABLambdaShiftM3 b

theorem f2ZeroShiftM3_forces_a_zero (a b : Rat) (h : f2ZeroShiftM3 a b = 0) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have haSq : a * a = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2ZeroShiftM3_forces_b_zero (a b : Rat) (h : f2ZeroShiftM3 a b = 0) :
    b = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have hbSq : b * b = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2ZeroShiftM1_forces_lambda_one (lam : Rat) (h : f2ZeroShiftM1 lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2ZeroShiftM1] at hCoeff
  linarith

theorem f2ZeroShiftM2_specialized :
    f2ZeroShiftM2 0 0 1 = X := by
  simp [f2ZeroShiftM2]

theorem f2ZeroShiftM2_specialized_nonzero :
    f2ZeroShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2ZeroShiftM2] at hCoeff

theorem noZeroShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftM1 lam = 0 ∧ f2ZeroShiftM2 a b lam = 0 ∧ f2ZeroShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  have ha : a = 0 := f2ZeroShiftM3_forces_a_zero a b hM3
  have hb : b = 0 := f2ZeroShiftM3_forces_b_zero a b hM3
  have hlam : lam = 1 := f2ZeroShiftM1_forces_lambda_one lam hM1
  have hSpec : f2ZeroShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f2ZeroShiftM2_specialized_nonzero hSpec

theorem f4ZeroShiftM0_forces_a_zero (a : Rat) (h : f4ZeroShiftM0 a = 0) :
    a = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftM0] at hEval
  exact hEval

theorem f4ZeroShiftM3_forces_b_zero (b : Rat) (h : f4ZeroShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftM3] at hEval
  exact hEval

theorem f4ZeroShiftM1_forces_lambda_one (a lam : Rat)
    (ha : a = 0) (h : f4ZeroShiftM1 a lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4ZeroShiftM1, ha] at hCoeff
  linarith

theorem f4ZeroShiftM2_specialized :
    f4ZeroShiftM2 0 0 1 = X := by
  simp [f4ZeroShiftM2]

theorem f4ZeroShiftM2_specialized_nonzero :
    f4ZeroShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4ZeroShiftM2] at hCoeff

theorem noZeroShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4ZeroShiftM1 a lam = 0 ∧
      f4ZeroShiftM2 a b lam = 0 ∧
      f4ZeroShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4ZeroShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4ZeroShiftM3_forces_b_zero b hM3
  have hlam : lam = 1 := f4ZeroShiftM1_forces_lambda_one a lam ha hM1
  have hSpec : f4ZeroShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f4ZeroShiftM2_specialized_nonzero hSpec

theorem f2ZeroShiftPlainPrefactorStage1B_forces_a_zero
    (a b : Rat) (h : f2ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftPlainPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftPlainPrefactorStage1B_forces_b_one
    (a b : Rat) (h : f2ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftPlainPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftPlainPrefactorStage1A_specialized_nonmatch (lam : Rat) :
    f2ZeroShiftPlainPrefactorStage1A 0 1 lam ≠ heroStage1NumeratorQ := by
  intro h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f2ZeroShiftPlainPrefactorStage1A, heroStage1NumeratorQ] at hEval1 hEval2
  nlinarith [hEval1, hEval2]

theorem noZeroShiftPlainPrefactorF2DirectMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, hB1, _, _⟩
  have ha : a = 0 := f2ZeroShiftPlainPrefactorStage1B_forces_a_zero a b hB1
  have hb : b = 1 := f2ZeroShiftPlainPrefactorStage1B_forces_b_one a b hB1
  have hSpec : f2ZeroShiftPlainPrefactorStage1A 0 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using hA1
  exact f2ZeroShiftPlainPrefactorStage1A_specialized_nonmatch lam hSpec

theorem f2ZeroShiftAPlusPrefactorStage1B_forces_a_zero
    (a b : Rat) (h : f2ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftAPlusPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftAPlusPrefactorStage1B_forces_b_one
    (a b : Rat) (h : f2ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ) :
    b = 1 := by
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftAPlusPrefactorStage1B, heroStage1DenominatorQ] at hEvalNeg1
  linarith

theorem f2ZeroShiftAPlusPrefactorStage1A_specialized_nonmatch (lam : Rat) :
    f2ZeroShiftAPlusPrefactorStage1A 0 1 lam ≠ heroStage1NumeratorQ := by
  intro h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f2ZeroShiftAPlusPrefactorStage1A, heroStage1NumeratorQ] at hEval1 hEval2
  nlinarith [hEval1, hEval2]

theorem noZeroShiftAPlusPrefactorF2DirectMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, hB1, _, _⟩
  have ha : a = 0 := f2ZeroShiftAPlusPrefactorStage1B_forces_a_zero a b hB1
  have hb : b = 1 := f2ZeroShiftAPlusPrefactorStage1B_forces_b_one a b hB1
  have hSpec : f2ZeroShiftAPlusPrefactorStage1A 0 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using hA1
  exact f2ZeroShiftAPlusPrefactorStage1A_specialized_nonmatch lam hSpec

theorem f2ZeroShiftBPlusPrefactorStage1B_forces_b_one
    (a b : Rat) (h : f2ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftBPlusPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftBPlusPrefactorStage1B_forces_a_neg_one
    (a b : Rat) (h : f2ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ) :
    a = -1 := by
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftBPlusPrefactorStage1B, heroStage1DenominatorQ] at hEvalNeg1
  linarith

theorem f2ZeroShiftBPlusPrefactorStage1A_specialized_nonmatch (lam : Rat) :
    f2ZeroShiftBPlusPrefactorStage1A (-1) 1 lam ≠ heroStage1NumeratorQ := by
  intro h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f2ZeroShiftBPlusPrefactorStage1A, heroStage1NumeratorQ] at hEval1 hEval2
  nlinarith [hEval1, hEval2]

theorem noZeroShiftBPlusPrefactorF2DirectMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, hB1, _, _⟩
  have ha : a = -1 := f2ZeroShiftBPlusPrefactorStage1B_forces_a_neg_one a b hB1
  have hb : b = 1 := f2ZeroShiftBPlusPrefactorStage1B_forces_b_one a b hB1
  have hSpec : f2ZeroShiftBPlusPrefactorStage1A (-1) 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using hA1
  exact f2ZeroShiftBPlusPrefactorStage1A_specialized_nonmatch lam hSpec

theorem f2ZeroShiftPlusLambdaPrefactorStage1B_forces_a_zero
    (a b : Rat) (h : f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftPlusLambdaPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftPlusLambdaPrefactorStage1B_forces_b_one
    (a b : Rat) (h : f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2ZeroShiftPlusLambdaPrefactorStage1B, heroStage1DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f2ZeroShiftPlusLambdaPrefactorStage1A_forces_lambda_one
    (a b lam : Rat)
    (ha : a = 0) (hb : b = 1)
    (h : f2ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ) :
    lam = 1 := by
  have hSpec : f2ZeroShiftPlusLambdaPrefactorStage1A 0 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using h
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) hSpec
  norm_num [f2ZeroShiftPlusLambdaPrefactorStage1A, heroStage1NumeratorQ] at hEval
  linarith

theorem f2ZeroShiftPlusLambdaPrefactorStage2A_specialized :
    f2ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 = X ^ 2 + X ^ 3 := by
  simp [f2ZeroShiftPlusLambdaPrefactorStage2A]

theorem f2ZeroShiftPlusLambdaPrefactorStage2A_specialized_nonmatch :
    f2ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 ≠ heroStage2NumeratorQ := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 3) h
  norm_num [f2ZeroShiftPlusLambdaPrefactorStage2A, heroStage2NumeratorQ] at hCoeff

theorem noZeroShiftPlusLambdaPrefactorF2DirectMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, hB1, hA2, hB2⟩
  have ha : a = 0 := f2ZeroShiftPlusLambdaPrefactorStage1B_forces_a_zero a b hB1
  have hb : b = 1 := f2ZeroShiftPlusLambdaPrefactorStage1B_forces_b_one a b hB1
  have hlam : lam = 1 := f2ZeroShiftPlusLambdaPrefactorStage1A_forces_lambda_one a b lam ha hb hA1
  have hSpec : f2ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 = heroStage2NumeratorQ := by
    simpa [ha, hb, hlam] using hA2
  exact f2ZeroShiftPlusLambdaPrefactorStage2A_specialized_nonmatch hSpec

theorem f4ZeroShiftPlainPrefactorStage2B_forces_a_zero
    (a b : Rat) (h : f4ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftPlainPrefactorStage2B, heroStage2DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f4ZeroShiftPlainPrefactorStage2B_forces_b_one
    (a b : Rat) (h : f4ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftPlainPrefactorStage2B, heroStage2DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f4ZeroShiftPlainPrefactorStage1A_specialized_nonmatch (lam : Rat) :
    f4ZeroShiftPlainPrefactorStage1A 0 1 lam ≠ heroStage1NumeratorQ := by
  intro h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f4ZeroShiftPlainPrefactorStage1A, heroStage1NumeratorQ] at hEval1 hEval2
  nlinarith [hEval1, hEval2]

theorem noZeroShiftPlainPrefactorF4DirectMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, _, _, hB2⟩
  have ha : a = 0 := f4ZeroShiftPlainPrefactorStage2B_forces_a_zero a b hB2
  have hb : b = 1 := f4ZeroShiftPlainPrefactorStage2B_forces_b_one a b hB2
  have hSpec : f4ZeroShiftPlainPrefactorStage1A 0 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using hA1
  exact f4ZeroShiftPlainPrefactorStage1A_specialized_nonmatch lam hSpec

theorem f4ZeroShiftAPlusPrefactorStage2B_forces_a_zero
    (a b : Rat) (h : f4ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftAPlusPrefactorStage2B, heroStage2DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f4ZeroShiftAPlusPrefactorStage2B_forces_b_one
    (a b : Rat) (h : f4ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) :
    b = 1 := by
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftAPlusPrefactorStage2B, heroStage2DenominatorQ] at hEvalNeg1
  linarith

theorem f4ZeroShiftAPlusPrefactorStage1A_specialized_nonmatch (lam : Rat) :
    f4ZeroShiftAPlusPrefactorStage1A 0 1 lam ≠ heroStage1NumeratorQ := by
  intro h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f4ZeroShiftAPlusPrefactorStage1A, heroStage1NumeratorQ] at hEval1 hEval2
  nlinarith [hEval1, hEval2]

theorem noZeroShiftAPlusPrefactorF4DirectMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, _, _, hB2⟩
  have ha : a = 0 := f4ZeroShiftAPlusPrefactorStage2B_forces_a_zero a b hB2
  have hb : b = 1 := f4ZeroShiftAPlusPrefactorStage2B_forces_b_one a b hB2
  have hSpec : f4ZeroShiftAPlusPrefactorStage1A 0 1 lam = heroStage1NumeratorQ := by
    simpa [ha, hb] using hA1
  exact f4ZeroShiftAPlusPrefactorStage1A_specialized_nonmatch lam hSpec

theorem f4ZeroShiftBPlusPrefactorStage2B_nonmatch (a b : Rat) :
    f4ZeroShiftBPlusPrefactorStage2B a b ≠ heroStage2DenominatorQ := by
  intro h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f4ZeroShiftBPlusPrefactorStage2B, heroStage2DenominatorQ] at hEvalNeg1 hEval1 hEval2
  nlinarith [hEvalNeg1, hEval1, hEval2]

theorem noZeroShiftBPlusPrefactorF4DirectMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, _, _, _, hB2⟩
  exact f4ZeroShiftBPlusPrefactorStage2B_nonmatch a b hB2

theorem f4ZeroShiftPlusLambdaPrefactorStage2B_forces_a_zero
    (a b : Rat) (h : f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftPlusLambdaPrefactorStage2B, heroStage2DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f4ZeroShiftPlusLambdaPrefactorStage2B_forces_b_one
    (a b : Rat) (h : f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f4ZeroShiftPlusLambdaPrefactorStage2B, heroStage2DenominatorQ] at hEval1 hEvalNeg1
  nlinarith [hEval1, hEvalNeg1]

theorem f4ZeroShiftPlusLambdaPrefactorStage1A_forces_lambda_one
    (a b lam : Rat) (ha : a = 0)
    (h : f4ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ) :
    lam = 1 := by
  have hSpec : f4ZeroShiftPlusLambdaPrefactorStage1A 0 b lam = heroStage1NumeratorQ := by
    simpa [ha] using h
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) hSpec
  norm_num [f4ZeroShiftPlusLambdaPrefactorStage1A, heroStage1NumeratorQ] at hEval
  linarith

theorem f4ZeroShiftPlusLambdaPrefactorStage2A_specialized :
    f4ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 = X ^ 2 + X ^ 3 := by
  simp [f4ZeroShiftPlusLambdaPrefactorStage2A]

theorem f4ZeroShiftPlusLambdaPrefactorStage2A_specialized_nonmatch :
    f4ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 ≠ heroStage2NumeratorQ := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 3) h
  norm_num [f4ZeroShiftPlusLambdaPrefactorStage2A, heroStage2NumeratorQ] at hCoeff

theorem noZeroShiftPlusLambdaPrefactorF4DirectMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ := by
  rintro ⟨a, b, lam, hA1, hB1, hA2, hB2⟩
  have ha : a = 0 := f4ZeroShiftPlusLambdaPrefactorStage2B_forces_a_zero a b hB2
  have hb : b = 1 := f4ZeroShiftPlusLambdaPrefactorStage2B_forces_b_one a b hB2
  have hlam : lam = 1 := f4ZeroShiftPlusLambdaPrefactorStage1A_forces_lambda_one a b lam ha hA1
  have hSpec : f4ZeroShiftPlusLambdaPrefactorStage2A 0 1 1 = heroStage2NumeratorQ := by
    simpa [ha, hb, hlam] using hA2
  exact f4ZeroShiftPlusLambdaPrefactorStage2A_specialized_nonmatch hSpec

theorem noZeroShiftPlusLambdaPrefactorDirectMatches :
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) := by
  exact ⟨
    noZeroShiftPlusLambdaPrefactorF2DirectMatch,
    noZeroShiftPlusLambdaPrefactorF4DirectMatch
  ⟩

theorem f2ZeroShiftAReciprocalPrefactorDenResidual_forces_b_one
    (a b : Rat) (h : f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0) :
    b = 1 := by
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  norm_num [f2ZeroShiftAReciprocalPrefactorDenResidual] at hEval0
  linarith

theorem f2ZeroShiftAReciprocalPrefactorDenResidual_forces_a_zero
    (a b : Rat) (h : f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0) :
    a = 0 := by
  have hb : b = 1 := f2ZeroShiftAReciprocalPrefactorDenResidual_forces_b_one a b h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f2ZeroShiftAReciprocalPrefactorDenResidual, hb] at hEval1
  linarith

theorem f2ZeroShiftAReciprocalPrefactorNumResidual_specialized_nonzero (lam : Rat) :
    f2ZeroShiftAReciprocalPrefactorNumResidual 0 1 lam ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f2ZeroShiftAReciprocalPrefactorNumResidual] at hEval0 hEval1
  nlinarith [hEval0, hEval1]

theorem noZeroShiftAReciprocalPrefactorF2CrossMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
      f2ZeroShiftAReciprocalPrefactorNumResidual a b lam = 0 := by
  rintro ⟨a, b, lam, hDen, hNum⟩
  have ha : a = 0 := f2ZeroShiftAReciprocalPrefactorDenResidual_forces_a_zero a b hDen
  have hb : b = 1 := f2ZeroShiftAReciprocalPrefactorDenResidual_forces_b_one a b hDen
  have hSpec : f2ZeroShiftAReciprocalPrefactorNumResidual 0 1 lam = 0 := by
    simpa [ha, hb] using hNum
  exact f2ZeroShiftAReciprocalPrefactorNumResidual_specialized_nonzero lam hSpec

theorem f2ZeroShiftBReciprocalPrefactorDenResidual_nonzero (a b : Rat) :
    f2ZeroShiftBReciprocalPrefactorDenResidual a b ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f2ZeroShiftBReciprocalPrefactorDenResidual] at hEval0 hEval1 hEval2
  nlinarith [hEval0, hEval1, hEval2]

theorem noZeroShiftBReciprocalPrefactorF2CrossMatch :
    ¬ ∃ a b : Rat,
      f2ZeroShiftBReciprocalPrefactorDenResidual a b = 0 := by
  rintro ⟨a, b, hDen⟩
  exact f2ZeroShiftBReciprocalPrefactorDenResidual_nonzero a b hDen

theorem f2ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_b_one
    (a b : Rat) (h : f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0) :
    b = 1 := by
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  norm_num [f2ZeroShiftLambdaReciprocalPrefactorDenResidual] at hEval0
  linarith

theorem f2ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_a_zero
    (a b : Rat) (h : f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0) :
    a = 0 := by
  have hb : b = 1 := f2ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_b_one a b h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f2ZeroShiftLambdaReciprocalPrefactorDenResidual, hb] at hEval1
  linarith

theorem f2ZeroShiftLambdaReciprocalPrefactorNumResidual_specialized_nonzero (lam : Rat) :
    f2ZeroShiftLambdaReciprocalPrefactorNumResidual 0 1 lam ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f2ZeroShiftLambdaReciprocalPrefactorNumResidual] at hEval0 hEval1
  nlinarith [hEval0, hEval1]

theorem noZeroShiftLambdaReciprocalPrefactorF2CrossMatch :
    ¬ ∃ a b lam : Rat,
      f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
      f2ZeroShiftLambdaReciprocalPrefactorNumResidual a b lam = 0 := by
  rintro ⟨a, b, lam, hDen, hNum⟩
  have ha : a = 0 := f2ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_a_zero a b hDen
  have hb : b = 1 := f2ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_b_one a b hDen
  have hSpec : f2ZeroShiftLambdaReciprocalPrefactorNumResidual 0 1 lam = 0 := by
    simpa [ha, hb] using hNum
  exact f2ZeroShiftLambdaReciprocalPrefactorNumResidual_specialized_nonzero lam hSpec

theorem f4ZeroShiftAReciprocalPrefactorDenResidual_forces_a_zero
    (a b : Rat) (h : f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0) :
    a = 0 := by
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  norm_num [f4ZeroShiftAReciprocalPrefactorDenResidual] at hEval0
  linarith

theorem f4ZeroShiftAReciprocalPrefactorDenResidual_forces_b_one
    (a b : Rat) (h : f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0) :
    b = 1 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have ha : a = 0 := f4ZeroShiftAReciprocalPrefactorDenResidual_forces_a_zero a b h
  norm_num [f4ZeroShiftAReciprocalPrefactorDenResidual, ha] at hEval1
  linarith

theorem f4ZeroShiftAReciprocalPrefactorNumResidual_specialized_nonzero (lam : Rat) :
    f4ZeroShiftAReciprocalPrefactorNumResidual 0 lam ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftAReciprocalPrefactorNumResidual] at hEval0 hEval1
  nlinarith [hEval0, hEval1]

theorem noZeroShiftAReciprocalPrefactorF4CrossMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
      f4ZeroShiftAReciprocalPrefactorNumResidual a lam = 0 := by
  rintro ⟨a, b, lam, hDen, hNum⟩
  have ha : a = 0 := f4ZeroShiftAReciprocalPrefactorDenResidual_forces_a_zero a b hDen
  have hSpec : f4ZeroShiftAReciprocalPrefactorNumResidual 0 lam = 0 := by
    simpa [ha] using hNum
  exact f4ZeroShiftAReciprocalPrefactorNumResidual_specialized_nonzero lam hSpec

theorem f4ZeroShiftBReciprocalPrefactorDenResidual_nonzero (a b : Rat) :
    f4ZeroShiftBReciprocalPrefactorDenResidual a b ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEval2 := congrArg (fun p : QPoly => Polynomial.eval (2 : Rat) p) h
  norm_num [f4ZeroShiftBReciprocalPrefactorDenResidual] at hEval0 hEval1 hEval2
  nlinarith [hEval0, hEval1, hEval2]

theorem noZeroShiftBReciprocalPrefactorF4CrossMatch :
    ¬ ∃ a b : Rat,
      f4ZeroShiftBReciprocalPrefactorDenResidual a b = 0 := by
  rintro ⟨a, b, hDen⟩
  exact f4ZeroShiftBReciprocalPrefactorDenResidual_nonzero a b hDen

theorem f4ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_a_zero
    (a b : Rat) (h : f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0) :
    a = 0 := by
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  norm_num [f4ZeroShiftLambdaReciprocalPrefactorDenResidual] at hEval0
  linarith

theorem f4ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_b_one
    (a b : Rat) (h : f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0) :
    b = 1 := by
  have ha : a = 0 := f4ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_a_zero a b h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftLambdaReciprocalPrefactorDenResidual, ha] at hEval1
  linarith

theorem f4ZeroShiftLambdaReciprocalPrefactorNumResidual_specialized_nonzero (lam : Rat) :
    f4ZeroShiftLambdaReciprocalPrefactorNumResidual 0 lam ≠ 0 := by
  intro h
  have hEval0 := congrArg (fun p : QPoly => Polynomial.eval (0 : Rat) p) h
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4ZeroShiftLambdaReciprocalPrefactorNumResidual] at hEval0 hEval1
  nlinarith [hEval0, hEval1]

theorem noZeroShiftLambdaReciprocalPrefactorF4CrossMatch :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
      f4ZeroShiftLambdaReciprocalPrefactorNumResidual a lam = 0 := by
  rintro ⟨a, b, lam, hDen, hNum⟩
  have ha : a = 0 := f4ZeroShiftLambdaReciprocalPrefactorDenResidual_forces_a_zero a b hDen
  have hSpec : f4ZeroShiftLambdaReciprocalPrefactorNumResidual 0 lam = 0 := by
    simpa [ha] using hNum
  exact f4ZeroShiftLambdaReciprocalPrefactorNumResidual_specialized_nonzero lam hSpec

theorem noZeroShiftReciprocalSinglePrefactorF2CrossMatches :
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
      f2ZeroShiftAReciprocalPrefactorNumResidual a b lam = 0) ∧
    (¬ ∃ a b : Rat,
      f2ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
      f2ZeroShiftLambdaReciprocalPrefactorNumResidual a b lam = 0) := by
  exact ⟨
    noZeroShiftAReciprocalPrefactorF2CrossMatch,
    noZeroShiftBReciprocalPrefactorF2CrossMatch,
    noZeroShiftLambdaReciprocalPrefactorF2CrossMatch
  ⟩

theorem noZeroShiftReciprocalSinglePrefactorF4CrossMatches :
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
      f4ZeroShiftAReciprocalPrefactorNumResidual a lam = 0) ∧
    (¬ ∃ a b : Rat,
      f4ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
      f4ZeroShiftLambdaReciprocalPrefactorNumResidual a lam = 0) := by
  exact ⟨
    noZeroShiftAReciprocalPrefactorF4CrossMatch,
    noZeroShiftBReciprocalPrefactorF4CrossMatch,
    noZeroShiftLambdaReciprocalPrefactorF4CrossMatch
  ⟩

theorem noZeroShiftReciprocalSinglePrefactorCrossMatches :
    ((¬ ∃ a b lam : Rat,
        f2ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
        f2ZeroShiftAReciprocalPrefactorNumResidual a b lam = 0) ∧
      (¬ ∃ a b : Rat,
        f2ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
      (¬ ∃ a b lam : Rat,
        f2ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
        f2ZeroShiftLambdaReciprocalPrefactorNumResidual a b lam = 0)) ∧
    ((¬ ∃ a b lam : Rat,
        f4ZeroShiftAReciprocalPrefactorDenResidual a b = 0 ∧
        f4ZeroShiftAReciprocalPrefactorNumResidual a lam = 0) ∧
      (¬ ∃ a b : Rat,
        f4ZeroShiftBReciprocalPrefactorDenResidual a b = 0) ∧
      (¬ ∃ a b lam : Rat,
        f4ZeroShiftLambdaReciprocalPrefactorDenResidual a b = 0 ∧
        f4ZeroShiftLambdaReciprocalPrefactorNumResidual a lam = 0)) := by
  exact ⟨
    noZeroShiftReciprocalSinglePrefactorF2CrossMatches,
    noZeroShiftReciprocalSinglePrefactorF4CrossMatches
  ⟩

theorem noZeroShiftPolynomialSinglePrefactorF2DirectMatches :
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f2ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) := by
  refine ⟨noZeroShiftPlainPrefactorF2DirectMatch, ?_⟩
  refine ⟨noZeroShiftAPlusPrefactorF2DirectMatch, ?_⟩
  refine ⟨noZeroShiftBPlusPrefactorF2DirectMatch, ?_⟩
  exact noZeroShiftPlusLambdaPrefactorF2DirectMatch

theorem noZeroShiftPolynomialSinglePrefactorF4DirectMatches :
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
      f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ) := by
  refine ⟨noZeroShiftPlainPrefactorF4DirectMatch, ?_⟩
  refine ⟨noZeroShiftAPlusPrefactorF4DirectMatch, ?_⟩
  refine ⟨noZeroShiftBPlusPrefactorF4DirectMatch, ?_⟩
  exact noZeroShiftPlusLambdaPrefactorF4DirectMatch

theorem noZeroShiftPolynomialSinglePrefactorDirectMatches :
    ((¬ ∃ a b lam : Rat,
        f2ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f2ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f2ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f2ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f2ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f2ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f2ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f2ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f2ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f2ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f2ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f2ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f2ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f2ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f2ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f2ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ)) ∧
    ((¬ ∃ a b lam : Rat,
        f4ZeroShiftPlainPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f4ZeroShiftPlainPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f4ZeroShiftPlainPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f4ZeroShiftPlainPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f4ZeroShiftAPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f4ZeroShiftAPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f4ZeroShiftAPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f4ZeroShiftAPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f4ZeroShiftBPlusPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f4ZeroShiftBPlusPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f4ZeroShiftBPlusPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f4ZeroShiftBPlusPrefactorStage2B a b = heroStage2DenominatorQ) ∧
      (¬ ∃ a b lam : Rat,
        f4ZeroShiftPlusLambdaPrefactorStage1A a b lam = heroStage1NumeratorQ ∧
        f4ZeroShiftPlusLambdaPrefactorStage1B a b = heroStage1DenominatorQ ∧
        f4ZeroShiftPlusLambdaPrefactorStage2A a b lam = heroStage2NumeratorQ ∧
        f4ZeroShiftPlusLambdaPrefactorStage2B a b = heroStage2DenominatorQ)) := by
  exact ⟨
    noZeroShiftPolynomialSinglePrefactorF2DirectMatches,
    noZeroShiftPolynomialSinglePrefactorF4DirectMatches
  ⟩

theorem f2UnitAShiftM3_forces_a_zero (a b : Rat) (h : f2UnitAShiftM3 a b = 0) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitAShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a + a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have haSq : a * a = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitAShiftM3_forces_b_zero (a b : Rat) (h : f2UnitAShiftM3 a b = 0) :
    b = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitAShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a + a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have hbSq : b * b = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitAShiftM1_forces_lambda_one (lam : Rat) (h : f2UnitAShiftM1 lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitAShiftM1] at hCoeff
  linarith

theorem f2UnitAShiftM2_specialized :
    f2UnitAShiftM2 0 0 1 = X := by
  simp [f2UnitAShiftM2]

theorem f2UnitAShiftM2_specialized_nonzero :
    f2UnitAShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitAShiftM2] at hCoeff

theorem noUnitAShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitAShiftM1 lam = 0 ∧
      f2UnitAShiftM2 a b lam = 0 ∧
      f2UnitAShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  have ha : a = 0 := f2UnitAShiftM3_forces_a_zero a b hM3
  have hb : b = 0 := f2UnitAShiftM3_forces_b_zero a b hM3
  have hlam : lam = 1 := f2UnitAShiftM1_forces_lambda_one lam hM1
  have hSpec : f2UnitAShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f2UnitAShiftM2_specialized_nonzero hSpec

theorem f4UnitAShiftM0_forces_a_zero (a : Rat) (h : f4UnitAShiftM0 a = 0) :
    a = 0 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 2) h
  norm_num [f4UnitAShiftM0] at hCoeff
  exact hCoeff

theorem f4UnitAShiftM3_forces_b_zero (b : Rat) (h : f4UnitAShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitAShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitAShiftM1_forces_lambda_one (a lam : Rat)
    (ha : a = 0) (h : f4UnitAShiftM1 a lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitAShiftM1, ha] at hCoeff
  linarith

theorem f4UnitAShiftM2_specialized :
    f4UnitAShiftM2 0 0 1 = X := by
  simp [f4UnitAShiftM2]

theorem f4UnitAShiftM2_specialized_nonzero :
    f4UnitAShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitAShiftM2] at hCoeff

theorem noUnitAShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitAShiftM0 a = 0 ∧
      f4UnitAShiftM1 a lam = 0 ∧
      f4UnitAShiftM2 a b lam = 0 ∧
      f4UnitAShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4UnitAShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitAShiftM3_forces_b_zero b hM3
  have hlam : lam = 1 := f4UnitAShiftM1_forces_lambda_one a lam ha hM1
  have hSpec : f4UnitAShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f4UnitAShiftM2_specialized_nonzero hSpec

theorem f2UnitBShiftM3_forces_a_zero (a b : Rat) (h : f2UnitBShiftM3 a b = 0) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitBShiftM3] at hEval1 hEvalNeg1
  have hab : a * b = 0 := by
    nlinarith [hEval1, hEvalNeg1]
  have hsum : a * a + 2 * a * b + b * b = 0 := by
    nlinarith [hEval1, hEvalNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [hab, hsum]
  have haSq : a * a = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitBShiftM3_forces_b_zero (a b : Rat) (h : f2UnitBShiftM3 a b = 0) :
    b = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitBShiftM3] at hEval1 hEvalNeg1
  have hab : a * b = 0 := by
    nlinarith [hEval1, hEvalNeg1]
  have hsum : a * a + 2 * a * b + b * b = 0 := by
    nlinarith [hEval1, hEvalNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [hab, hsum]
  have hbSq : b * b = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitBShiftM2_specialized :
    f2UnitBShiftM2 0 0 1 = X := by
  simp [f2UnitBShiftM2]

theorem f2UnitBShiftM2_specialized_nonzero :
    f2UnitBShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitBShiftM2] at hCoeff

theorem noUnitBShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitBShiftM1 lam = 0 ∧
      f2UnitBShiftM2 a b lam = 0 ∧
      f2UnitBShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  have ha : a = 0 := f2UnitBShiftM3_forces_a_zero a b hM3
  have hb : b = 0 := f2UnitBShiftM3_forces_b_zero a b hM3
  have hlam : lam = 1 := f2ZeroShiftM1_forces_lambda_one lam hM1
  have hSpec : f2UnitBShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f2UnitBShiftM2_specialized_nonzero hSpec

theorem f4UnitBShiftM3_forces_b_zero (b : Rat) (h : f4UnitBShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitBShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitBShiftM2_specialized :
    f4UnitBShiftM2 0 0 1 = X := by
  simp [f4UnitBShiftM2]

theorem f4UnitBShiftM2_specialized_nonzero :
    f4UnitBShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitBShiftM2] at hCoeff

theorem noUnitBShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitBShiftM0 a = 0 ∧
      f4UnitBShiftM1 a lam = 0 ∧
      f4UnitBShiftM2 a b lam = 0 ∧
      f4UnitBShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4ZeroShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitBShiftM3_forces_b_zero b hM3
  have hlam : lam = 1 := f4ZeroShiftM1_forces_lambda_one a lam ha hM1
  have hSpec : f4UnitBShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f4UnitBShiftM2_specialized_nonzero hSpec

theorem f2UnitABShiftM3_forces_a_zero (a b : Rat) (h : f2UnitABShiftM3 a b = 0) :
    a = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitABShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have haSq : a * a = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitABShiftM3_forces_b_zero (a b : Rat) (h : f2UnitABShiftM3 a b = 0) :
    b = 0 := by
  have hEval1 := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  have hEvalNeg1 := congrArg (fun p : QPoly => Polynomial.eval (-1 : Rat) p) h
  norm_num [f2UnitABShiftM3] at hEval1 hEvalNeg1
  have h1 : a * a + 3 * a * b + b * b = 0 := by
    nlinarith [hEval1]
  have hNeg1 : a * a - a * b + b * b = 0 := by
    nlinarith [hEvalNeg1]
  have hab : a * b = 0 := by
    nlinarith [h1, hNeg1]
  have hsq : a * a + b * b = 0 := by
    nlinarith [h1, hNeg1]
  have hbSq : b * b = 0 := by
    nlinarith [hab, hsq]
  nlinarith

theorem f2UnitABShiftM2_specialized :
    f2UnitABShiftM2 0 0 1 = X := by
  simp [f2UnitABShiftM2]

theorem f2UnitABShiftM2_specialized_nonzero :
    f2UnitABShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitABShiftM2] at hCoeff

theorem noUnitABShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitABShiftM1 lam = 0 ∧
      f2UnitABShiftM2 a b lam = 0 ∧
      f2UnitABShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  have ha : a = 0 := f2UnitABShiftM3_forces_a_zero a b hM3
  have hb : b = 0 := f2UnitABShiftM3_forces_b_zero a b hM3
  have hlam : lam = 1 := f2ZeroShiftM1_forces_lambda_one lam hM1
  have hSpec : f2UnitABShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f2UnitABShiftM2_specialized_nonzero hSpec

theorem f4UnitABShiftM0_forces_a_zero (a : Rat) (h : f4UnitABShiftM0 a = 0) :
    a = 0 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 2) h
  norm_num [f4UnitABShiftM0] at hCoeff
  exact hCoeff

theorem f4UnitABShiftM3_forces_b_zero (b : Rat) (h : f4UnitABShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitABShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitABShiftM1_forces_lambda_one (a lam : Rat)
    (ha : a = 0) (h : f4UnitABShiftM1 a lam = 0) :
    lam = 1 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitABShiftM1, ha] at hCoeff
  linarith

theorem f4UnitABShiftM2_specialized :
    f4UnitABShiftM2 0 0 1 = X := by
  simp [f4UnitABShiftM2]

theorem f4UnitABShiftM2_specialized_nonzero :
    f4UnitABShiftM2 0 0 1 ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitABShiftM2] at hCoeff

theorem noUnitABShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitABShiftM0 a = 0 ∧
      f4UnitABShiftM1 a lam = 0 ∧
      f4UnitABShiftM2 a b lam = 0 ∧
      f4UnitABShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4UnitABShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitABShiftM3_forces_b_zero b hM3
  have hlam : lam = 1 := f4UnitABShiftM1_forces_lambda_one a lam ha hM1
  have hSpec : f4UnitABShiftM2 0 0 1 = 0 := by
    simpa [ha, hb, hlam] using hM2
  exact f4UnitABShiftM2_specialized_nonzero hSpec

theorem f2UnitALambdaShiftM1_nonzero (lam : Rat) :
    f2UnitALambdaShiftM1 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitALambdaShiftM1] at hCoeff

theorem noUnitALambdaShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitALambdaShiftM1 lam = 0 ∧
      f2UnitALambdaShiftM2 a b lam = 0 ∧
      f2UnitALambdaShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  exact f2UnitALambdaShiftM1_nonzero lam hM1

theorem f4UnitALambdaShiftM0_forces_a_zero (a : Rat) (h : f4UnitALambdaShiftM0 a = 0) :
    a = 0 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 2) h
  norm_num [f4UnitALambdaShiftM0] at hCoeff
  exact hCoeff

theorem f4UnitALambdaShiftM3_forces_b_zero (b : Rat) (h : f4UnitALambdaShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitALambdaShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitALambdaShiftM1_specialized_nonzero (lam : Rat) :
    f4UnitALambdaShiftM1 0 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitALambdaShiftM1] at hCoeff

theorem noUnitALambdaShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitALambdaShiftM0 a = 0 ∧
      f4UnitALambdaShiftM1 a lam = 0 ∧
      f4UnitALambdaShiftM2 a b lam = 0 ∧
      f4UnitALambdaShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4UnitALambdaShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitALambdaShiftM3_forces_b_zero b hM3
  have hSpec : f4UnitALambdaShiftM1 0 lam = 0 := by
    simpa [ha] using hM1
  exact f4UnitALambdaShiftM1_specialized_nonzero lam hSpec

theorem f2UnitBLambdaShiftM1_nonzero (lam : Rat) :
    f2UnitBLambdaShiftM1 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitBLambdaShiftM1] at hCoeff

theorem noUnitBLambdaShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitBLambdaShiftM1 lam = 0 ∧
      f2UnitBLambdaShiftM2 a b lam = 0 ∧
      f2UnitBLambdaShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  exact f2UnitBLambdaShiftM1_nonzero lam hM1

theorem f4UnitBLambdaShiftM0_forces_a_zero (a : Rat) (h : f4UnitBLambdaShiftM0 a = 0) :
    a = 0 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitBLambdaShiftM0] at hCoeff
  exact hCoeff

theorem f4UnitBLambdaShiftM3_forces_b_zero (b : Rat) (h : f4UnitBLambdaShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitBLambdaShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitBLambdaShiftM1_specialized_nonzero (lam : Rat) :
    f4UnitBLambdaShiftM1 0 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitBLambdaShiftM1] at hCoeff

theorem noUnitBLambdaShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitBLambdaShiftM0 a = 0 ∧
      f4UnitBLambdaShiftM1 a lam = 0 ∧
      f4UnitBLambdaShiftM2 a b lam = 0 ∧
      f4UnitBLambdaShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4UnitBLambdaShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitBLambdaShiftM3_forces_b_zero b hM3
  have hSpec : f4UnitBLambdaShiftM1 0 lam = 0 := by
    simpa [ha] using hM1
  exact f4UnitBLambdaShiftM1_specialized_nonzero lam hSpec

theorem f2UnitABLambdaShiftM1_nonzero (lam : Rat) :
    f2UnitABLambdaShiftM1 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitABLambdaShiftM1] at hCoeff

theorem noUnitABLambdaShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitABLambdaShiftM1 lam = 0 ∧
      f2UnitABLambdaShiftM2 a b lam = 0 ∧
      f2UnitABLambdaShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  exact f2UnitABLambdaShiftM1_nonzero lam hM1

theorem f4UnitABLambdaShiftM0_forces_a_zero (a : Rat) (h : f4UnitABLambdaShiftM0 a = 0) :
    a = 0 := by
  have hCoeff := congrArg (fun p : QPoly => p.coeff 2) h
  norm_num [f4UnitABLambdaShiftM0] at hCoeff
  exact hCoeff

theorem f4UnitABLambdaShiftM3_forces_b_zero (b : Rat) (h : f4UnitABLambdaShiftM3 b = 0) :
    b = 0 := by
  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h
  norm_num [f4UnitABLambdaShiftM3] at hEval
  nlinarith [hEval]

theorem f4UnitABLambdaShiftM1_specialized_nonzero (lam : Rat) :
    f4UnitABLambdaShiftM1 0 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitABLambdaShiftM1] at hCoeff

theorem noUnitABLambdaShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4UnitABLambdaShiftM0 a = 0 ∧
      f4UnitABLambdaShiftM1 a lam = 0 ∧
      f4UnitABLambdaShiftM2 a b lam = 0 ∧
      f4UnitABLambdaShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4UnitABLambdaShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4UnitABLambdaShiftM3_forces_b_zero b hM3
  have hSpec : f4UnitABLambdaShiftM1 0 lam = 0 := by
    simpa [ha] using hM1
  exact f4UnitABLambdaShiftM1_specialized_nonzero lam hSpec

theorem f2UnitLambdaShiftM1_nonzero (lam : Rat) :
    f2UnitLambdaShiftM1 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f2UnitLambdaShiftM1] at hCoeff

theorem noUnitLambdaShiftF2ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f2UnitLambdaShiftM1 lam = 0 ∧
      f2UnitLambdaShiftM2 a b lam = 0 ∧
      f2UnitLambdaShiftM3 a b = 0 := by
  rintro ⟨a, b, lam, hM1, hM2, hM3⟩
  exact f2UnitLambdaShiftM1_nonzero lam hM1

theorem f4UnitLambdaShiftM1_specialized_nonzero (lam : Rat) :
    f4UnitLambdaShiftM1 0 lam ≠ 0 := by
  intro h
  have hCoeff := congrArg (fun p : QPoly => p.coeff 1) h
  norm_num [f4UnitLambdaShiftM1] at hCoeff

theorem noUnitLambdaShiftF4ExactEquivalence :
    ¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4UnitLambdaShiftM1 a lam = 0 ∧
      f4UnitLambdaShiftM2 a b lam = 0 ∧
      f4UnitLambdaShiftM3 b = 0 := by
  rintro ⟨a, b, lam, hM0, hM1, hM2, hM3⟩
  have ha : a = 0 := f4ZeroShiftM0_forces_a_zero a hM0
  have hb : b = 0 := f4ZeroShiftM3_forces_b_zero b hM3
  have hSpec : f4UnitLambdaShiftM1 0 lam = 0 := by
    simpa [ha] using hM1
  exact f4UnitLambdaShiftM1_specialized_nonzero lam hSpec

theorem noNearestShiftCubeF2ExactEquivalence :
    (¬ ∃ a b lam : Rat,
      f2ZeroShiftM1 lam = 0 ∧
      f2ZeroShiftM2 a b lam = 0 ∧
      f2ZeroShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitAShiftM1 lam = 0 ∧
      f2UnitAShiftM2 a b lam = 0 ∧
      f2UnitAShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitBShiftM1 lam = 0 ∧
      f2UnitBShiftM2 a b lam = 0 ∧
      f2UnitBShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitABShiftM1 lam = 0 ∧
      f2UnitABShiftM2 a b lam = 0 ∧
      f2UnitABShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitALambdaShiftM1 lam = 0 ∧
      f2UnitALambdaShiftM2 a b lam = 0 ∧
      f2UnitALambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitBLambdaShiftM1 lam = 0 ∧
      f2UnitBLambdaShiftM2 a b lam = 0 ∧
      f2UnitBLambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitABLambdaShiftM1 lam = 0 ∧
      f2UnitABLambdaShiftM2 a b lam = 0 ∧
      f2UnitABLambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitLambdaShiftM1 lam = 0 ∧
      f2UnitLambdaShiftM2 a b lam = 0 ∧
      f2UnitLambdaShiftM3 a b = 0) := by
  exact ⟨
    noZeroShiftF2ExactEquivalence,
    ⟨
      noUnitAShiftF2ExactEquivalence,
      ⟨
        noUnitBShiftF2ExactEquivalence,
        ⟨
          noUnitABShiftF2ExactEquivalence,
          ⟨
            noUnitALambdaShiftF2ExactEquivalence,
            ⟨
              noUnitBLambdaShiftF2ExactEquivalence,
              ⟨
                noUnitABLambdaShiftF2ExactEquivalence,
                noUnitLambdaShiftF2ExactEquivalence
              ⟩
            ⟩
          ⟩
        ⟩
      ⟩
    ⟩
  ⟩

theorem noNearestShiftCubeF4ExactEquivalence :
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4ZeroShiftM1 a lam = 0 ∧
      f4ZeroShiftM2 a b lam = 0 ∧
      f4ZeroShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitAShiftM0 a = 0 ∧
      f4UnitAShiftM1 a lam = 0 ∧
      f4UnitAShiftM2 a b lam = 0 ∧
      f4UnitAShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitBShiftM0 a = 0 ∧
      f4UnitBShiftM1 a lam = 0 ∧
      f4UnitBShiftM2 a b lam = 0 ∧
      f4UnitBShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitABShiftM0 a = 0 ∧
      f4UnitABShiftM1 a lam = 0 ∧
      f4UnitABShiftM2 a b lam = 0 ∧
      f4UnitABShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitALambdaShiftM0 a = 0 ∧
      f4UnitALambdaShiftM1 a lam = 0 ∧
      f4UnitALambdaShiftM2 a b lam = 0 ∧
      f4UnitALambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitBLambdaShiftM0 a = 0 ∧
      f4UnitBLambdaShiftM1 a lam = 0 ∧
      f4UnitBLambdaShiftM2 a b lam = 0 ∧
      f4UnitBLambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitABLambdaShiftM0 a = 0 ∧
      f4UnitABLambdaShiftM1 a lam = 0 ∧
      f4UnitABLambdaShiftM2 a b lam = 0 ∧
      f4UnitABLambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4UnitLambdaShiftM1 a lam = 0 ∧
      f4UnitLambdaShiftM2 a b lam = 0 ∧
      f4UnitLambdaShiftM3 b = 0) := by
  exact ⟨
    noZeroShiftF4ExactEquivalence,
    ⟨
      noUnitAShiftF4ExactEquivalence,
      ⟨
        noUnitBShiftF4ExactEquivalence,
        ⟨
          noUnitABShiftF4ExactEquivalence,
          ⟨
            noUnitALambdaShiftF4ExactEquivalence,
            ⟨
              noUnitBLambdaShiftF4ExactEquivalence,
              ⟨
                noUnitABLambdaShiftF4ExactEquivalence,
                noUnitLambdaShiftF4ExactEquivalence
              ⟩
            ⟩
          ⟩
        ⟩
      ⟩
    ⟩
  ⟩

theorem noNearestShiftCubeExactEquivalence :
    ((¬ ∃ a b lam : Rat,
      f2ZeroShiftM1 lam = 0 ∧
      f2ZeroShiftM2 a b lam = 0 ∧
      f2ZeroShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitAShiftM1 lam = 0 ∧
      f2UnitAShiftM2 a b lam = 0 ∧
      f2UnitAShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitBShiftM1 lam = 0 ∧
      f2UnitBShiftM2 a b lam = 0 ∧
      f2UnitBShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitABShiftM1 lam = 0 ∧
      f2UnitABShiftM2 a b lam = 0 ∧
      f2UnitABShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitALambdaShiftM1 lam = 0 ∧
      f2UnitALambdaShiftM2 a b lam = 0 ∧
      f2UnitALambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitBLambdaShiftM1 lam = 0 ∧
      f2UnitBLambdaShiftM2 a b lam = 0 ∧
      f2UnitBLambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitABLambdaShiftM1 lam = 0 ∧
      f2UnitABLambdaShiftM2 a b lam = 0 ∧
      f2UnitABLambdaShiftM3 a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f2UnitLambdaShiftM1 lam = 0 ∧
      f2UnitLambdaShiftM2 a b lam = 0 ∧
      f2UnitLambdaShiftM3 a b = 0)) ∧
    ((¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4ZeroShiftM1 a lam = 0 ∧
      f4ZeroShiftM2 a b lam = 0 ∧
      f4ZeroShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitAShiftM0 a = 0 ∧
      f4UnitAShiftM1 a lam = 0 ∧
      f4UnitAShiftM2 a b lam = 0 ∧
      f4UnitAShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitBShiftM0 a = 0 ∧
      f4UnitBShiftM1 a lam = 0 ∧
      f4UnitBShiftM2 a b lam = 0 ∧
      f4UnitBShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitABShiftM0 a = 0 ∧
      f4UnitABShiftM1 a lam = 0 ∧
      f4UnitABShiftM2 a b lam = 0 ∧
      f4UnitABShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitALambdaShiftM0 a = 0 ∧
      f4UnitALambdaShiftM1 a lam = 0 ∧
      f4UnitALambdaShiftM2 a b lam = 0 ∧
      f4UnitALambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitBLambdaShiftM0 a = 0 ∧
      f4UnitBLambdaShiftM1 a lam = 0 ∧
      f4UnitBLambdaShiftM2 a b lam = 0 ∧
      f4UnitBLambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4UnitABLambdaShiftM0 a = 0 ∧
      f4UnitABLambdaShiftM1 a lam = 0 ∧
      f4UnitABLambdaShiftM2 a b lam = 0 ∧
      f4UnitABLambdaShiftM3 b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4ZeroShiftM0 a = 0 ∧
      f4UnitLambdaShiftM1 a lam = 0 ∧
      f4UnitLambdaShiftM2 a b lam = 0 ∧
      f4UnitLambdaShiftM3 b = 0)) := by
  exact ⟨noNearestShiftCubeF2ExactEquivalence, noNearestShiftCubeF4ExactEquivalence⟩

theorem noNearestShiftCubeF2ExactEquivalenceFor
    (shiftA shiftB shiftLambda : Bool) :
    ¬ ∃ a b lam : Rat,
      f2NearestShiftCubeM1 shiftA shiftB shiftLambda lam = 0 ∧
      f2NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      f2NearestShiftCubeM3 shiftA shiftB shiftLambda a b = 0 := by
  cases shiftA <;> cases shiftB <;> cases shiftLambda
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noZeroShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitLambdaShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitBShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitBLambdaShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitAShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitALambdaShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitABShiftF2ExactEquivalence
  · simpa [f2NearestShiftCubeM1, f2NearestShiftCubeM2, f2NearestShiftCubeM3] using noUnitABLambdaShiftF2ExactEquivalence

theorem noNearestShiftCubeF4ExactEquivalenceFor
    (shiftA shiftB shiftLambda : Bool) :
    ¬ ∃ a b lam : Rat,
      f4NearestShiftCubeM0 shiftA shiftB shiftLambda a = 0 ∧
      f4NearestShiftCubeM1 shiftA shiftB shiftLambda a lam = 0 ∧
      f4NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      f4NearestShiftCubeM3 shiftA shiftB shiftLambda b = 0 := by
  cases shiftA <;> cases shiftB <;> cases shiftLambda
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noZeroShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitLambdaShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitBShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitBLambdaShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitAShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitALambdaShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitABShiftF4ExactEquivalence
  · simpa [f4NearestShiftCubeM0, f4NearestShiftCubeM1, f4NearestShiftCubeM2, f4NearestShiftCubeM3] using noUnitABLambdaShiftF4ExactEquivalence

theorem noNearestShiftCubeExactEquivalenceFor
    (shiftA shiftB shiftLambda : Bool) :
    (¬ ∃ a b lam : Rat,
      f2NearestShiftCubeM1 shiftA shiftB shiftLambda lam = 0 ∧
      f2NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      f2NearestShiftCubeM3 shiftA shiftB shiftLambda a b = 0) ∧
    (¬ ∃ a b lam : Rat,
      f4NearestShiftCubeM0 shiftA shiftB shiftLambda a = 0 ∧
      f4NearestShiftCubeM1 shiftA shiftB shiftLambda a lam = 0 ∧
      f4NearestShiftCubeM2 shiftA shiftB shiftLambda a b lam = 0 ∧
      f4NearestShiftCubeM3 shiftA shiftB shiftLambda b = 0) := by
  exact ⟨
    noNearestShiftCubeF2ExactEquivalenceFor shiftA shiftB shiftLambda,
    noNearestShiftCubeF4ExactEquivalenceFor shiftA shiftB shiftLambda
  ⟩

end

end Page43Equivalence
end HeroCase
end Proofs
