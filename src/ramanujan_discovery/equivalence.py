from __future__ import annotations

from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate


def family_bucket(
    template: QCFTemplate,
    matched_target: str,
    closest_benchmark: str,
    novelty_status: str,
) -> str:
    if novelty_status == "review":
        return template.exploratory_family_key(closest_benchmark)
    if matched_target != "unmatched":
        return matched_target
    return template.signature()


def equivalence_key(
    template: QCFTemplate,
    matched_target: str,
    closest_benchmark: str,
    novelty_status: str,
) -> str:
    if novelty_status != "review":
        if matched_target != "unmatched":
            return matched_target
        return template.signature()

    if not closest_benchmark:
        return template.signature()

    normalized = template.normalized()
    benchmark_template = get_benchmark(closest_benchmark).canonical_template.normalized()

    if normalized.numerator_extra_scale == 0:
        extra_exists = 0
        extra_scale = 0
        extra_shift_delta = 0
        extra_step_delta = 0
    else:
        extra_exists = 1
        extra_scale = normalized.numerator_extra_scale
        extra_shift_delta = normalized.numerator_extra_q_shift - normalized.numerator_q_shift
        extra_step_delta = normalized.numerator_extra_q_step - normalized.numerator_q_step

    if normalized.denominator_scale == 0:
        denominator_exists = 0
        denominator_scale = 0
        denominator_shift_delta = 0
        denominator_step_delta = 0
    else:
        denominator_exists = 1
        denominator_scale = normalized.denominator_scale
        denominator_shift_delta = normalized.denominator_q_shift - benchmark_template.denominator_q_shift
        denominator_step_delta = normalized.denominator_q_step - benchmark_template.denominator_q_step

    return (
        "review::"
        f"closest={closest_benchmark}::"
        f"num_scale={normalized.numerator_scale}::"
        f"main_shift_delta={normalized.numerator_q_shift - benchmark_template.numerator_q_shift}::"
        f"main_step_delta={normalized.numerator_q_step - benchmark_template.numerator_q_step}::"
        f"extra_exists={extra_exists}::"
        f"extra_scale={extra_scale}::"
        f"extra_shift_delta={extra_shift_delta}::"
        f"extra_step_delta={extra_step_delta}::"
        f"den_exists={denominator_exists}::"
        f"den_scale={denominator_scale}::"
        f"den_shift_delta={denominator_shift_delta}::"
        f"den_step_delta={denominator_step_delta}"
    )


def select_unique_reviews(
    ordered_reviews: list[CandidateRecord],
    admitted_records: list[CandidateRecord],
    max_review_candidates: int,
) -> list[CandidateRecord]:
    seen_signatures = {record.template.signature() for record in admitted_records}
    seen_equivalence_keys = {
        record.equivalence_key for record in admitted_records if record.novelty_status == "review"
    }
    unique_reviews: list[CandidateRecord] = []

    for record in ordered_reviews:
        signature = record.template.signature()
        if signature in seen_signatures or record.equivalence_key in seen_equivalence_keys:
            continue
        seen_signatures.add(signature)
        seen_equivalence_keys.add(record.equivalence_key)
        unique_reviews.append(record)
        if len(unique_reviews) >= max_review_candidates:
            break

    return unique_reviews
