from ramanujan_discovery.benchmarks import CUBIC_Q3_TEMPLATE, CUBIC_TEMPLATE, RR_Q4_TEMPLATE, RR_TEMPLATE
from ramanujan_discovery.config import SearchConfig
from ramanujan_discovery.discovery import discover_candidates


def test_discovery_rediscovers_classical_benchmark():
    records = discover_candidates(
        SearchConfig(
            depth=36,
            precision=80,
            budget_hours=0.1,
            min_discovery_digits=18,
            min_stability_digits=18,
            min_review_stability=28,
            max_per_target=4,
            max_review_candidates=6,
        )
    )

    assert any(
        record.matched_target == "rogers_ramanujan_normalized"
        and record.template.signature() == RR_TEMPLATE.signature()
        for record in records
    )
    assert any(
        record.matched_target == "ramanujan_cubic_normalized"
        and record.template.signature() == CUBIC_TEMPLATE.signature()
        for record in records
    )
    assert any(
        record.matched_target == "rogers_ramanujan_q4_normalized"
        and record.template.signature() == RR_Q4_TEMPLATE.signature()
        for record in records
    )
    assert any(
        record.matched_target == "ramanujan_cubic_q3_normalized"
        and record.template.signature() == CUBIC_Q3_TEMPLATE.signature()
        for record in records
    )
    review_records = [record for record in records if record.novelty_status == "review"]
    assert review_records
    assert len(review_records) == 4
    assert len(review_records) == len({record.family_bucket for record in review_records})
