from ramanujan_discovery.equivalence import equivalence_key, family_bucket, select_unique_reviews
from ramanujan_discovery.models import CandidateRecord, QCFTemplate


def _review_record(candidate_id: str, template: QCFTemplate, closest_benchmark: str, digits: int) -> CandidateRecord:
    return CandidateRecord(
        id=candidate_id,
        template=template,
        q_values=[0.05, 0.09],
        value_estimates=["0.0", "0.0"],
        matched_target="unmatched",
        closest_benchmark=closest_benchmark,
        closest_benchmark_digits=digits,
        family_bucket=family_bucket(
            template=template,
            matched_target="unmatched",
            closest_benchmark=closest_benchmark,
            novelty_status="review",
        ),
        equivalence_key=equivalence_key(
            template=template,
            matched_target="unmatched",
            closest_benchmark=closest_benchmark,
            novelty_status="review",
        ),
        benchmark_kind="exploratory",
        digits_agree=digits,
        stability_score=90,
        novelty_status="review",
        notes="synthetic review candidate",
    )


def test_equivalence_key_ignores_inactive_fields():
    template_a = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=4,
        numerator_q_step=3,
        numerator_extra_scale=0,
        numerator_extra_q_shift=8,
        numerator_extra_q_step=9,
        denominator_constant=1,
        denominator_scale=0,
        denominator_q_shift=7,
        denominator_q_step=5,
    )
    template_b = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=4,
        numerator_q_step=3,
        numerator_extra_scale=0,
        numerator_extra_q_shift=1,
        numerator_extra_q_step=1,
        denominator_constant=1,
        denominator_scale=0,
        denominator_q_shift=1,
        denominator_q_step=1,
    )

    assert (
        equivalence_key(template_a, "unmatched", "rogers_ramanujan_q4_normalized", "review")
        == equivalence_key(template_b, "unmatched", "rogers_ramanujan_q4_normalized", "review")
    )


def test_equivalence_key_depends_on_closest_benchmark():
    template = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=2,
        numerator_q_step=2,
        numerator_extra_scale=1,
        numerator_extra_q_shift=6,
        numerator_extra_q_step=6,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=3,
        denominator_q_step=2,
    )

    assert equivalence_key(template, "unmatched", "shifted_rr_fixture", "review") != equivalence_key(
        template, "unmatched", "rogers_ramanujan_q2_normalized", "review"
    )


def test_select_unique_reviews_keeps_distinct_equivalence_keys_in_same_family():
    first = _review_record(
        "review-a",
        QCFTemplate(
            numerator_scale=1,
            numerator_q_shift=4,
            numerator_q_step=3,
            denominator_constant=1,
            denominator_scale=0,
            denominator_q_shift=1,
            denominator_q_step=1,
        ),
        "rogers_ramanujan_q4_normalized",
        digits=8,
    )
    second = _review_record(
        "review-b",
        QCFTemplate(
            numerator_scale=1,
            numerator_q_shift=3,
            numerator_q_step=4,
            denominator_constant=1,
            denominator_scale=0,
            denominator_q_shift=1,
            denominator_q_step=1,
        ),
        "rogers_ramanujan_q4_normalized",
        digits=7,
    )

    assert first.family_bucket == second.family_bucket
    assert first.equivalence_key != second.equivalence_key
    assert [record.id for record in select_unique_reviews([first, second], [], max_review_candidates=6)] == [
        "review-a",
        "review-b",
    ]
