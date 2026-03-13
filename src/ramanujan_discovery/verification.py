from __future__ import annotations

from ramanujan_discovery.benchmarks import benchmark_names, get_benchmark, target_value
from ramanujan_discovery.config import VerificationConfig
from ramanujan_discovery.continued_fraction import agreement_digits, evaluate_qcf, format_mpf
from ramanujan_discovery.models import CandidateRecord
from ramanujan_discovery.storage import read_candidates


def _classify(record: CandidateRecord, matched_target: str) -> tuple[str, str]:
    benchmark = get_benchmark(matched_target)
    if record.template.signature() == benchmark.canonical_template.signature():
        if benchmark.kind.startswith("classical"):
            return "known", "Exact canonical match to a classical benchmark."
        return "fixture", "Exact canonical match to an internal regression fixture."

    if benchmark.kind.startswith("classical"):
        return "known_variant", "Strong numerical match to a classical benchmark; review for equivalence."

    return "fixture", "Strong numerical match to an internal regression fixture."


def _verification_stability(record: CandidateRecord, q_values: tuple[float, ...], depth: int, precision: int) -> int:
    scores = []
    for q_value in q_values:
        coarse = evaluate_qcf(record.template, q=q_value, depth=depth, precision=precision)
        refined = evaluate_qcf(record.template, q=q_value, depth=depth + 10, precision=precision)
        scores.append(agreement_digits(coarse, refined, cap=precision - 2))
    return min(scores)


def _rank_key(record: CandidateRecord) -> tuple[int, int, int, str]:
    return (
        record.digits_agree,
        record.stability_score,
        -record.template.complexity_score(),
        record.template.signature(),
    )


def verify_candidates(input_path: str, config: VerificationConfig) -> list[CandidateRecord]:
    verified: list[CandidateRecord] = []
    review_candidates: list[CandidateRecord] = []

    for record in read_candidates(input_path):
        try:
            values = [
                evaluate_qcf(record.template, q=q_value, depth=config.depth, precision=config.precision)
                for q_value in config.q_values
            ]
            stability = _verification_stability(record, config.q_values, config.depth, config.precision)
        except ZeroDivisionError:
            continue
        except ValueError:
            continue

        best_name = ""
        best_digits = -1
        best_kind = ""
        for target_name in benchmark_names():
            digits = []
            for q_value, candidate_value in zip(config.q_values, values):
                target = target_value(target_name, q=q_value, precision=config.precision, depth=config.depth + 10)
                digits.append(agreement_digits(candidate_value, target, cap=config.precision - 2))
            score = min(digits)
            if score > best_digits:
                best_name = target_name
                best_digits = score
                best_kind = get_benchmark(target_name).kind

        if best_digits >= config.min_verified_digits:
            novelty_status, notes = _classify(record, best_name)
            verified.append(
                CandidateRecord(
                    id=record.id,
                    template=record.template,
                    q_values=list(config.q_values),
                    value_estimates=[format_mpf(value) for value in values],
                    matched_target=best_name,
                    closest_benchmark=best_name,
                    closest_benchmark_digits=best_digits,
                    benchmark_kind=best_kind,
                    digits_agree=best_digits,
                    stability_score=stability,
                    novelty_status=novelty_status,
                    notes=notes,
                )
            )
            continue

        if stability >= config.min_review_stability:
            review_candidates.append(
                CandidateRecord(
                    id=record.id,
                    template=record.template,
                    q_values=list(config.q_values),
                    value_estimates=[format_mpf(value) for value in values],
                    matched_target="unmatched",
                    closest_benchmark=best_name,
                    closest_benchmark_digits=best_digits,
                    benchmark_kind="exploratory",
                    digits_agree=best_digits,
                    stability_score=stability,
                    novelty_status="review",
                    notes=(
                        "Stable at higher precision but still unmatched against the benchmark catalog. "
                        f"Closest benchmark: {best_name} ({best_digits} shared digits)."
                    ),
                )
            )

    ordered_reviews = sorted(review_candidates, key=_rank_key, reverse=True)
    seen_signatures = {record.template.signature() for record in verified}
    for record in ordered_reviews:
        signature = record.template.signature()
        if signature in seen_signatures:
            continue
        seen_signatures.add(signature)
        verified.append(record)
        if sum(item.novelty_status == "review" for item in verified) >= config.max_review_candidates:
            break

    return sorted(verified, key=_rank_key, reverse=True)
