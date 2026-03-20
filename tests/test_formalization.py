from pathlib import Path

from ramanujan_discovery.cli import main
from ramanujan_discovery.formalization import _formalization_build_profile
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.storage import write_candidates


def test_formalization_build_profile_uses_stronger_full_rational_prefactor_box():
    smoke = _formalization_build_profile(smoke=True)
    full = _formalization_build_profile(smoke=False)

    assert smoke.page43_rational_max_nontrivial_profiles == 1
    assert smoke.page43_stages == 2
    assert full.page43_rational_max_nontrivial_profiles == 3
    assert full.page43_stages == 3


def test_cli_formalize_writes_note(tmp_path: Path):
    verified = tmp_path / "verified.jsonl"
    output_path = tmp_path / "formalization.md"
    lean_output_path = tmp_path / "formalization.lean"

    record = CandidateRecord(
        id="hero",
        template=QCFTemplate(
            numerator_scale=1,
            numerator_q_shift=3,
            numerator_q_step=3,
            numerator_extra_scale=1,
            numerator_extra_q_shift=6,
            numerator_extra_q_step=6,
            denominator_constant=1,
            denominator_scale=1,
            denominator_q_shift=3,
            denominator_q_step=3,
        ),
        q_values=[0.04],
        value_estimates=["0.0"],
        matched_target="unmatched",
        closest_benchmark="rogers_ramanujan_q3_normalized",
        closest_benchmark_digits=7,
        family_bucket="hybrid_perturbed_family",
        equivalence_key="review::demo",
        benchmark_kind="exploratory",
        digits_agree=7,
        stability_score=50,
        novelty_status="review",
        notes="demo",
    )
    write_candidates(verified, [record])

    assert (
        main(
            [
                "formalize",
                "--in",
                str(verified),
                "--candidate-id",
                "hero",
                "--smoke",
                "--out",
                str(output_path),
                "--lean-out",
                str(lean_output_path),
            ]
        )
        == 0
    )
    text = output_path.read_text(encoding="utf-8")
    assert "Formalization Prep: `hero`" in text
    assert "Build profile: `smoke`" in text
    assert "Current Theorem Status" in text
    assert "exact intermediate waypoint" in text
    assert "Exact Lemma Candidates" in text
    assert "Heine `cor2cf` Odd/Even Branch Obstructions" in text
    assert "Exact `f2` / `gcf3` `n`-Dependent Equivalence Lane" in text
    assert "Exact `f4` / `gcf2` `n`-Dependent Equivalence Lane" in text
    assert "Exact Unit-Shift `a` Page-43 Lanes" in text
    assert "Exact Unit-Shift `b` Page-43 Lanes" in text
    assert "Exact Mixed Unit-Shift `a`/`b` Page-43 Lanes" in text
    assert "Exact Mixed Unit-Shift `a`/`lambda` Page-43 Lanes" in text
    assert "Exact Mixed Unit-Shift `b`/`lambda` Page-43 Lanes" in text
    assert "Exact Mixed Unit-Shift `a`/`b`/`lambda` Page-43 Lanes" in text
    assert "Exact Unit-Shift `lambda` Page-43 Lanes" in text
    assert "Exact Zero-Shift Polynomial Single-Prefactor Page-43 Lanes" in text
    assert "Exact Zero-Shift Reciprocal Single-Prefactor Page-43 Lanes" in text
    assert "cross-multiplied coefficient identities" in text
    assert "fully exactified at theorem grade" in text
    assert "nearest-shift cube" in text
    assert "Bool-parameterized" in text
    assert "Bounded Exact Exclusion Results" in text
    assert "Page-43 low-complexity rational-prefactor box" in text
    assert "shift box `[-1,1]`" in text
    assert "`f2/gcf3` hits `0`" in text
    assert "`f4/gcf2` hits `0`" in text
    assert "Exact Reduction And Equivalence Witness" in text
    assert "fraction-field coefficient layer" in text
    assert "Formalization Order" in text
    assert "exclusion_waypoint" in text
    assert "currentExactWaypointCertificate" in text
    assert "finiteConvergentReductionWaypoint_true" in text
    assert "knownSourceOrbitExclusionWaypoint_true" in text
    assert "exactWaypointStatement_true" in text
    assert "mortonSquaredCoordinateExcluded" in text
    assert "weberSchlafliCoordinateExcluded" in text
    assert "mortonNamedCoordinateResearchWaypoint_true" in text
    assert "weberSchlafliBridgeResearchWaypoint_true" in text
    assert "stores that tighter `P_ws -> Weber` bridge directly" in text
    assert "Proofs/HeroCaseWeberSchlafliCoordinateObstruction.lean" in text
    assert "Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean" in text
    assert "Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean" in text
    assert "X_mt = F^2" in text
    assert "T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)" in text
    assert "P_ws = (1/F - F) / 2" in text
    assert "B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)" in text
    assert "(2,-1)" in text
    assert "(10,3)" in text
    assert "{-1, 3, 8}" in text
    assert "2,4,6,8,10" in text
    assert "currentNamedWeberOrbitResearchCertificate" in text
    assert "namedWeberOrbitResearchWaypoint_true" in text
    assert "currentRecognitionFrontierCertificate" in text
    assert "currentRecognitionFrontierWaypoint_true" in text
    assert "exact Chan--Huang obstruction quartets" in text
    assert "(-9/2, -6)" in text
    assert "(593/10, 1186/15)" in text
    assert "not ready" in text

    lean_text = lean_output_path.read_text(encoding="utf-8")
    assert "namespace Proofs" in lean_text
    assert "def candidateData" in lean_text
    assert "theorem candidate_stage0_a" in lean_text
    assert "theorem rr_direct_obstruction" in lean_text
    assert "theorem candidate_second_convergent_num" in lean_text
    assert "reverse equivalence transform" in lean_text
    assert "fraction-field coefficient layer" in lean_text
    assert "currentExactWaypointCertificate" in lean_text
    assert "finiteConvergentReductionWaypoint_true" in lean_text
    assert "knownSourceOrbitExclusionWaypoint_true" in lean_text
    assert "exactWaypointStatement_true" in lean_text
    assert "mortonSquaredCoordinateExcluded" in lean_text
    assert "weberSchlafliCoordinateExcluded" in lean_text
    assert "mortonNamedCoordinateResearchWaypoint_true" in lean_text
    assert "weberSchlafliBridgeResearchWaypoint_true" in lean_text
    assert "stores that tighter `P_ws -> Weber` bridge directly" in lean_text
    assert "Proofs/HeroCaseWeberSchlafliCoordinateObstruction.lean" in lean_text
    assert "Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean" in lean_text
    assert "Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean" in lean_text
    assert "X_mt = F^2" in lean_text
    assert "T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)" in lean_text
    assert "P_ws = (1/F - F) / 2" in lean_text
    assert "B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)" in lean_text
    assert "(2,-1)" in lean_text
    assert "(10,3)" in lean_text
    assert "{-1, 3, 8}" in lean_text
    assert "2,4,6,8,10" in lean_text
    assert "currentNamedWeberOrbitResearchCertificate" in lean_text
    assert "namedWeberOrbitResearchWaypoint_true" in lean_text
    assert "currentRecognitionFrontierCertificate" in lean_text
    assert "currentRecognitionFrontierWaypoint_true" in lean_text
    assert "reverseEquivalenceRecoversHeroData" in lean_text
    assert "nearestArithmeticSubsequenceSourcesExcluded" in lean_text
    assert "directLocalObstructions" in lean_text
    assert "simpleCor2cfBranchesExcluded" in lean_text
    assert "Current source-family-specific exact lanes" in lean_text
    assert "mixed unit-a/unit-b-shift" in lean_text
    assert "mixed unit-a/unit-lambda-shift" in lean_text
    assert "mixed unit-b/unit-lambda-shift" in lean_text
    assert "mixed unit-a/unit-b/unit-lambda-shift" in lean_text
    assert "page43PolynomialPrefactorExcluded" in lean_text
    assert "page43ReciprocalPrefactorExcluded" in lean_text
    assert "noZeroShiftPolynomialSinglePrefactorF2DirectMatches" in lean_text
    assert "noZeroShiftPolynomialSinglePrefactorF4DirectMatches" in lean_text
    assert "noZeroShiftPolynomialSinglePrefactorDirectMatches" in lean_text
    assert "noZeroShiftReciprocalSinglePrefactorF2CrossMatches" in lean_text
    assert "noZeroShiftReciprocalSinglePrefactorF4CrossMatches" in lean_text
    assert "noZeroShiftReciprocalSinglePrefactorCrossMatches" in lean_text
    assert "noNearestShiftCubeF2ExactEquivalence" in lean_text
    assert "noNearestShiftCubeF2ExactEquivalenceFor" in lean_text
    assert "noNearestShiftCubeF4ExactEquivalenceFor" in lean_text
    assert "noNearestShiftCubeExactEquivalenceFor" in lean_text
    assert "Proofs/HeroCasePage43Equivalence.lean" in lean_text
    assert "Suggested next theorem extensions" in lean_text
    assert "rr_direct_obstruction" in lean_text
    assert "sorry" not in lean_text
