from ramanujan_discovery.benchmarks import (
    CUBIC_Q2_TEMPLATE,
    CUBIC_Q3_TEMPLATE,
    CUBIC_TEMPLATE,
    RR_Q2_TEMPLATE,
    RR_Q3_TEMPLATE,
    RR_Q4_TEMPLATE,
    RR_TEMPLATE,
    target_value,
)
from ramanujan_discovery.continued_fraction import agreement_digits, evaluate_qcf


def test_rogers_ramanujan_matches_product_formula():
    for q_value in (0.05, 0.1, 0.15):
        continued_fraction_value = evaluate_qcf(RR_TEMPLATE, q=q_value, depth=48, precision=120)
        product_value = target_value("rogers_ramanujan_normalized", q=q_value, precision=120, depth=60)
        assert agreement_digits(continued_fraction_value, product_value, cap=100) >= 40


def test_extended_classical_benchmarks_match_product_formula():
    checks = [
        ("rogers_ramanujan_q2_normalized", RR_Q2_TEMPLATE),
        ("rogers_ramanujan_q3_normalized", RR_Q3_TEMPLATE),
        ("rogers_ramanujan_q4_normalized", RR_Q4_TEMPLATE),
        ("ramanujan_cubic_normalized", CUBIC_TEMPLATE),
        ("ramanujan_cubic_q2_normalized", CUBIC_Q2_TEMPLATE),
        ("ramanujan_cubic_q3_normalized", CUBIC_Q3_TEMPLATE),
    ]
    for benchmark_name, template in checks:
        for q_value in (0.04, 0.09, 0.15):
            continued_fraction_value = evaluate_qcf(template, q=q_value, depth=72, precision=140)
            product_value = target_value(benchmark_name, q=q_value, precision=140, depth=96)
            assert agreement_digits(continued_fraction_value, product_value, cap=120) >= 50
