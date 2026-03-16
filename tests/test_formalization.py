from pathlib import Path

from ramanujan_discovery.cli import main
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.storage import write_candidates


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
    assert "Exact Lemma Candidates" in text
    assert "Heine `cor2cf` Odd/Even Branch Obstructions" in text
    assert "Exact `f2` / `gcf3` `n`-Dependent Equivalence Lane" in text
    assert "Exact `f4` / `gcf2` `n`-Dependent Equivalence Lane" in text
    assert "Exact Unit-Shift `a` Page-43 Lanes" in text
    assert "Exact Unit-Shift `lambda` Page-43 Lanes" in text
    assert "Bounded Exact Exclusion Results" in text
    assert "Page-43 low-complexity rational-prefactor box" in text
    assert "Exact Reduction And Equivalence Witness" in text
    assert "fraction-field coefficient layer" in text
    assert "Formalization Order" in text
    assert "not ready" in text

    lean_text = lean_output_path.read_text(encoding="utf-8")
    assert "namespace Proofs" in lean_text
    assert "def candidateData" in lean_text
    assert "theorem candidate_stage0_a" in lean_text
    assert "theorem rr_direct_obstruction" in lean_text
    assert "theorem candidate_second_convergent_num" in lean_text
    assert "reverse equivalence transform" in lean_text
    assert "fraction-field coefficient layer" in lean_text
    assert "Current source-family-specific exact lanes" in lean_text
    assert "Proofs/HeroCasePage43Equivalence.lean" in lean_text
    assert "Suggested next theorem extensions" in lean_text
    assert "rr_direct_obstruction" in lean_text
    assert "sorry" not in lean_text
