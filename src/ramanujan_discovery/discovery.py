from __future__ import annotations

import hashlib
import time
from collections import defaultdict
from ramanujan_discovery.benchmarks import BENCHMARKS, benchmark_names, target_value
from ramanujan_discovery.config import SearchConfig
from ramanujan_discovery.continued_fraction import agreement_digits, evaluate_qcf, format_mpf
from ramanujan_discovery.models import CandidateRecord, QCFTemplate


def iter_search_templates():
    yielded: set[str] = set()

    def emit(template: QCFTemplate):
        signature = template.signature()
        if signature in yielded:
            return
        yielded.add(signature)
        return template

    for shift in (1, 2, 3, 4):
        for step in (1, 2, 3, 4):
            template = emit(
                QCFTemplate(
                    numerator_scale=1,
                    numerator_q_shift=shift,
                    numerator_q_step=step,
                    denominator_constant=1,
                    denominator_scale=0,
                    denominator_q_shift=1,
                    denominator_q_step=1,
                )
            )
            if template is not None:
                yield template

    for shift in (1, 2, 3, 4):
        for step in (1, 2, 3, 4):
            for denominator_shift in (1, 2, 3, 4):
                for denominator_step in (0, 1, 2, 3):
                    template = emit(
                        QCFTemplate(
                            numerator_scale=1,
                            numerator_q_shift=shift,
                            numerator_q_step=step,
                            denominator_constant=1,
                            denominator_scale=1,
                            denominator_q_shift=denominator_shift,
                            denominator_q_step=denominator_step,
                        )
                    )
                    if template is not None:
                        yield template

    for step in (1, 2, 3):
        template = emit(
            QCFTemplate(
                numerator_scale=1,
                numerator_q_shift=step,
                numerator_q_step=step,
                numerator_extra_scale=1,
                numerator_extra_q_shift=2 * step,
                numerator_extra_q_step=2 * step,
                denominator_constant=1,
                denominator_scale=0,
                denominator_q_shift=1,
                denominator_q_step=1,
            )
        )
        if template is not None:
            yield template

    for step in (1, 2, 3):
        for extra_multiplier in (2, 3):
            for denominator_shift in (1, 2, 3):
                template = emit(
                    QCFTemplate(
                        numerator_scale=1,
                        numerator_q_shift=step,
                        numerator_q_step=step,
                        numerator_extra_scale=1,
                        numerator_extra_q_shift=extra_multiplier * step,
                        numerator_extra_q_step=extra_multiplier * step,
                        denominator_constant=1,
                        denominator_scale=1,
                        denominator_q_shift=denominator_shift,
                        denominator_q_step=step,
                    )
                )
                if template is not None:
                    yield template


def _candidate_id(template: QCFTemplate, target_name: str) -> str:
    digest = hashlib.sha1(f"{template.signature()}::{target_name}".encode("utf-8")).hexdigest()
    return digest[:12]


def _stability_score(template: QCFTemplate, q_values: tuple[float, ...], depth: int, precision: int) -> int:
    scores = []
    for q_value in q_values:
        coarse = evaluate_qcf(template, q=q_value, depth=depth, precision=precision)
        refined = evaluate_qcf(template, q=q_value, depth=depth + 8, precision=precision)
        scores.append(agreement_digits(coarse, refined, cap=precision - 2))
    return min(scores)


def _best_target(template: QCFTemplate, config: SearchConfig):
    values = [evaluate_qcf(template, q=value, depth=config.depth, precision=config.precision) for value in config.q_values]
    best_name = ""
    best_kind = ""
    best_digits = -1

    for target_name in benchmark_names():
        digits = []
        for q_value, candidate_value in zip(config.q_values, values):
            target = target_value(target_name, q=q_value, precision=config.precision, depth=config.depth + 8)
            digits.append(agreement_digits(candidate_value, target, cap=config.precision - 2))

        min_digits = min(digits)
        if min_digits > best_digits:
            best_name = target_name
            best_kind = BENCHMARKS[target_name].kind
            best_digits = min_digits

    return values, best_name, best_kind, best_digits


def _rank_key(record: CandidateRecord) -> tuple[int, int, int, str]:
    return (
        record.digits_agree,
        record.stability_score,
        -record.template.complexity_score(),
        record.template.signature(),
    )


def discover_candidates(config: SearchConfig) -> list[CandidateRecord]:
    started_at = time.monotonic()
    grouped: dict[str, list[CandidateRecord]] = defaultdict(list)
    review_records: list[CandidateRecord] = []

    for template in iter_search_templates():
        if time.monotonic() - started_at > config.budget_hours * 3600:
            break

        try:
            values, target_name, benchmark_kind, digits = _best_target(template, config)
            stability = _stability_score(template, config.q_values, config.depth, config.precision)
        except ZeroDivisionError:
            continue
        except ValueError:
            continue

        if digits >= config.min_discovery_digits and stability >= config.min_stability_digits:
            grouped[target_name].append(
                CandidateRecord(
                    id=_candidate_id(template, target_name),
                    template=template,
                    q_values=list(config.q_values),
                    value_estimates=[format_mpf(value) for value in values],
                    matched_target=target_name,
                    closest_benchmark=target_name,
                    closest_benchmark_digits=digits,
                    benchmark_kind=benchmark_kind,
                    digits_agree=digits,
                    stability_score=stability,
                    novelty_status="unreviewed",
                    notes="Matched a built-in benchmark from the template grid.",
                )
            )
            continue

        if stability >= config.min_review_stability:
            review_records.append(
                CandidateRecord(
                    id=_candidate_id(template, "review"),
                    template=template,
                    q_values=list(config.q_values),
                    value_estimates=[format_mpf(value) for value in values],
                    matched_target="unmatched",
                    closest_benchmark=target_name,
                    closest_benchmark_digits=digits,
                    benchmark_kind="exploratory",
                    digits_agree=digits,
                    stability_score=stability,
                    novelty_status="review",
                    notes=(
                        "Stable under deeper truncation but not close to the built-in benchmark catalog. "
                        f"Closest benchmark: {target_name} ({digits} shared digits)."
                    ),
                )
            )

    results: list[CandidateRecord] = []
    for target_name, records in grouped.items():
        ordered = sorted(records, key=_rank_key, reverse=True)
        unique_records: list[CandidateRecord] = []
        seen_signatures: set[str] = set()
        for record in ordered:
            signature = record.template.signature()
            if signature in seen_signatures:
                continue
            seen_signatures.add(signature)
            unique_records.append(record)
        results.extend(unique_records[: config.max_per_target])

    ordered_reviews = sorted(review_records, key=_rank_key, reverse=True)
    seen_signatures = {record.template.signature() for record in results}
    unique_reviews: list[CandidateRecord] = []
    for record in ordered_reviews:
        signature = record.template.signature()
        if signature in seen_signatures:
            continue
        seen_signatures.add(signature)
        unique_reviews.append(record)
        if len(unique_reviews) >= config.max_review_candidates:
            break

    results.extend(unique_reviews)
    return sorted(results, key=_rank_key, reverse=True)
