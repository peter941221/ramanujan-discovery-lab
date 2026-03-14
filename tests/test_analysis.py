from pathlib import Path

from ramanujan_discovery.analysis import build_candidate_analysis_note, build_candidate_terminal_summary
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.storage import write_candidates


def _write_demo_record(input_path: Path) -> None:
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
        q_values=[0.04, 0.11],
        value_estimates=["0.0", "0.0"],
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
    write_candidates(input_path, [record])


def test_build_candidate_analysis_note(tmp_path: Path):
    input_path = tmp_path / "verified.jsonl"
    output_path = tmp_path / "hero.md"
    _write_demo_record(input_path)

    build_candidate_analysis_note(
        input_path=str(input_path),
        candidate_id="hero",
        output_path=str(output_path),
        depth=8,
        series_order=25,
    )

    text = output_path.read_text(encoding="utf-8")
    assert "Hero Case Analysis: `hero`" in text
    assert "Closest benchmark: `rogers_ramanujan_q3_normalized`" in text
    assert "First divergence order: `12`" in text
    assert "Low-Order Multiplicative Fit" in text
    assert "Alternative benchmark: `ramanujan_cubic_q3_normalized`" in text

    summary = build_candidate_terminal_summary(
        input_path=str(input_path),
        candidate_id="hero",
        depth=8,
        series_order=25,
        math_format="unicode",
    )
    assert "closest benchmark: rogers_ramanujan_q3_normalized" in summary
    assert "low-order multiplicative fit:" in summary
