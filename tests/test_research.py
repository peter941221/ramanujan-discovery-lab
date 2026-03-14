from pathlib import Path

import sympy as sp

from ramanujan_discovery.cli import main
from ramanujan_discovery.benchmarks import CUBIC_TEMPLATE, RR_TEMPLATE
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.research import (
    arithmetic_subsequence_contraction_coeffs,
    arithmetic_subsequence_contraction_search,
    bauer_muir_pattern_search,
    continued_fraction_convergents,
    convergent_common_factor_reduction,
    direct_bauer_muir_obstruction,
    euler_product_exponents,
    heine_hcf2_standardized_coeffs,
    parity_contraction_coeffs,
    page43_monomial_parameter_search,
    reduce_template_by_step,
    _template_reciprocal_coeffs,
)
from ramanujan_discovery.storage import write_candidates


def test_euler_product_exponents_recovers_simple_factor():
    q = sp.Symbol("q")
    expr = sp.series(1 / (1 - q), q, 0, 10).removeO()
    exponents = euler_product_exponents(expr, q=q, order=10)
    assert exponents[0] == -1
    assert all(value == 0 for value in exponents[1:])


def test_heine_hcf2_c_equals_bz_minus_one_has_constant_first_numerator():
    q = sp.Symbol("q")
    a, b = sp.symbols("a b")
    coeffs = heine_hcf2_standardized_coeffs(a=a, b=b, c=-1, z=-1 / b, q=q, depth=2)
    assert sp.simplify(coeffs.b0 - 1) == 0
    assert sp.simplify(coeffs.b_terms[1] - (1 + q)) == 0
    assert not coeffs.a_terms[1].has(q)


def test_reduce_template_by_step_divides_exponents():
    template = QCFTemplate(
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
    )
    reduced = reduce_template_by_step(template, step=3)
    assert reduced is not None
    assert reduced.numerator_q_shift == 1
    assert reduced.numerator_q_step == 1
    assert reduced.numerator_extra_q_shift == 2
    assert reduced.numerator_extra_q_step == 2
    assert reduced.denominator_q_shift == 1
    assert reduced.denominator_q_step == 1


def test_convergent_common_factor_reduction_finds_exact_hero_reduction():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )
    b0, a_terms, b_terms = _template_reciprocal_coeffs(target, q=t, depth=12)
    reduction = convergent_common_factor_reduction(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
    )

    assert reduction.gcd_factors[0] == 1
    for n in range(1, 9):
        assert sp.simplify(reduction.gcd_factors[n] - (1 + t**n)) == 0

    reduced = reduction.reduced_coeffs
    assert sp.simplify(reduced.b0 - 1) == 0
    assert sp.simplify(reduced.a_terms[1] - t) == 0
    assert sp.simplify(reduced.b_terms[1] - 1) == 0
    assert sp.simplify(reduced.a_terms[2] - t**2) == 0
    assert sp.simplify(reduced.b_terms[2] - (1 + t)) == 0
    assert sp.simplify(reduced.a_terms[3] - t**3 * (1 + t)) == 0
    assert sp.simplify(reduced.b_terms[3] - (1 + t**2)) == 0
    assert sp.simplify(reduced.a_terms[4] - t**4 * (1 + t**2)) == 0
    assert sp.simplify(reduced.b_terms[4] - (1 + t**3)) == 0


def test_direct_bauer_muir_obstruction_rules_out_rr_and_cubic_sources():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )

    rr_check = direct_bauer_muir_obstruction(
        source_label="RR reciprocal",
        source_template=RR_TEMPLATE,
        target_template=target,
        q=t,
    )
    assert rr_check.obstruction_stage == 1
    assert sp.simplify(rr_check.transformed_a1 - t) == 0
    assert sp.simplify(rr_check.target_a1 - (t + t**2)) == 0

    cubic_check = direct_bauer_muir_obstruction(
        source_label="cubic reciprocal",
        source_template=CUBIC_TEMPLATE,
        target_template=target,
        q=t,
    )
    assert cubic_check.obstruction_stage == 2
    assert sp.simplify(cubic_check.forced_w2 - t**2) == 0
    assert sp.simplify(cubic_check.target_a2 - (t**2 + t**4)) == 0
    assert sp.simplify(cubic_check.transformed_a2 - (-t + t**2 - t**3 + t**4)) == 0


def test_bauer_muir_pattern_search_finds_identity_zero_chain():
    t = sp.Symbol("t")

    one_step_hits = bauer_muir_pattern_search(
        source_label="RR reciprocal",
        source_template=RR_TEMPLATE,
        target_template=RR_TEMPLATE,
        q=t,
        depth=3,
        steps=1,
    )
    assert any(hit["pattern_chain"] == ["w_n = 0"] for hit in one_step_hits)

    two_step_hits = bauer_muir_pattern_search(
        source_label="RR reciprocal",
        source_template=RR_TEMPLATE,
        target_template=RR_TEMPLATE,
        q=t,
        depth=3,
        steps=2,
    )
    assert any(hit["pattern_chain"] == ["w_n = 0", "w_n = 0"] for hit in two_step_hits)

    three_step_hits = bauer_muir_pattern_search(
        source_label="RR reciprocal",
        source_template=RR_TEMPLATE,
        target_template=RR_TEMPLATE,
        q=t,
        depth=3,
        steps=3,
    )
    assert any(hit["pattern_chain"] == ["w_n = 0", "w_n = 0", "w_n = 0"] for hit in three_step_hits)


def test_bauer_muir_pattern_search_has_no_small_chain_hit_for_hero_template():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )

    assert (
        bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=RR_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=1,
        )
        == []
    )
    assert (
        bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=RR_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=2,
        )
        == []
    )
    assert (
        bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=CUBIC_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=1,
        )
        == []
    )
    assert (
        bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=CUBIC_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=2,
        )
        == []
    )
    assert (
        bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=RR_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=3,
        )
        == []
    )
    assert (
        bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=CUBIC_TEMPLATE,
            target_template=target,
            q=t,
            depth=4,
            steps=3,
        )
        == []
    )


def test_parity_contraction_reproduces_cubic_odd_even_convergents():
    t = sp.Symbol("t")
    b0, a_terms, b_terms = _template_reciprocal_coeffs(CUBIC_TEMPLATE, q=t, depth=6)
    source_convergents = continued_fraction_convergents(b0=b0, a_terms=a_terms, b_terms=b_terms)

    odd_part = parity_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        parity="odd",
    )
    even_part = parity_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        parity="even",
    )

    odd_values = [
        sp.simplify(numerator / denominator)
        for numerator, denominator in continued_fraction_convergents(
            b0=odd_part.b0,
            a_terms=odd_part.a_terms,
            b_terms=odd_part.b_terms,
        )
    ]
    even_values = [
        sp.simplify(numerator / denominator)
        for numerator, denominator in continued_fraction_convergents(
            b0=even_part.b0,
            a_terms=even_part.a_terms,
            b_terms=even_part.b_terms,
        )
    ]

    assert odd_values == [
        sp.simplify(numerator / denominator) for numerator, denominator in source_convergents[1::2]
    ]
    assert even_values == [
        sp.simplify(numerator / denominator) for numerator, denominator in source_convergents[0::2]
    ]


def test_arithmetic_subsequence_contraction_specializes_to_parity_cases():
    t = sp.Symbol("t")
    b0, a_terms, b_terms = _template_reciprocal_coeffs(CUBIC_TEMPLATE, q=t, depth=6)
    even_part = parity_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        parity="even",
    )
    odd_part = parity_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        parity="odd",
    )

    even_general = arithmetic_subsequence_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        stride=2,
        offset=0,
    )
    odd_general = arithmetic_subsequence_contraction_coeffs(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
        stride=2,
        offset=1,
    )

    assert sp.simplify(even_general.b0 - even_part.b0) == 0
    assert even_general.a_terms[:3] == even_part.a_terms[:3]
    assert even_general.b_terms[:3] == even_part.b_terms[:3]
    assert sp.simplify(odd_general.b0 - odd_part.b0) == 0
    assert odd_general.a_terms[:3] == odd_part.a_terms[:3]
    assert odd_general.b_terms[:3] == odd_part.b_terms[:3]


def test_cubic_parity_contractions_expose_hero_denominator_mismatch():
    t = sp.Symbol("t")
    source_b0, source_a_terms, source_b_terms = _template_reciprocal_coeffs(CUBIC_TEMPLATE, q=t, depth=6)
    odd_part = parity_contraction_coeffs(
        b0=source_b0,
        a_terms=source_a_terms,
        b_terms=source_b_terms,
        parity="odd",
    )
    even_part = parity_contraction_coeffs(
        b0=source_b0,
        a_terms=source_a_terms,
        b_terms=source_b_terms,
        parity="even",
    )

    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )
    target_b0, target_a_terms, target_b_terms = _template_reciprocal_coeffs(target, q=t, depth=4)

    assert sp.simplify(odd_part.b0 - (1 + t + t**2)) == 0
    assert sp.simplify(target_b0 - 1) == 0
    assert sp.simplify(even_part.b0 - target_b0) == 0
    assert sp.simplify(even_part.a_terms[1] - target_a_terms[1]) == 0
    assert sp.simplify(even_part.b_terms[1] - (1 + t**2 + t**4)) == 0
    assert sp.simplify(target_b_terms[1] - (1 + t)) == 0
    assert sp.expand(sp.simplify(even_part.b_terms[1] - target_b_terms[1])) == t**4 + t**2 - t


def test_page43_monomial_search_finds_simple_f2_specialization():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )
    hits = page43_monomial_parameter_search(
        family="f2",
        target_template=target,
        q=t,
        max_shift=2,
        stages=3,
    )
    assert any(
        hit.a_shift == 0
        and hit.b_shift == 0
        and hit.lambda_shift == 0
        and sp.simplify(hit.a_coeff) == 0
        and sp.simplify(hit.b_coeff - 1) == 0
        and sp.simplify(hit.lambda_coeff - 1) == 0
        for hit in hits
    )


def test_page43_monomial_search_has_no_hit_for_hero_template():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )
    assert (
        page43_monomial_parameter_search(
            family="f2",
            target_template=target,
            q=t,
            max_shift=3,
            stages=3,
        )
        == []
    )
    assert (
        page43_monomial_parameter_search(
            family="f4",
            target_template=target,
            q=t,
            max_shift=3,
            stages=3,
        )
        == []
    )


def test_arithmetic_subsequence_contraction_search_has_no_hit_for_hero_template():
    t = sp.Symbol("t")
    target = QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=1,
        numerator_q_step=1,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2,
        numerator_extra_q_step=2,
        denominator_constant=1,
        denominator_scale=1,
        denominator_q_shift=1,
        denominator_q_step=1,
    )
    assert (
        arithmetic_subsequence_contraction_search(
            source_label="RR reciprocal",
            source_template=RR_TEMPLATE,
            target_template=target,
            q=t,
            max_stride=4,
            stages=3,
        )
        == []
    )
    assert (
        arithmetic_subsequence_contraction_search(
            source_label="cubic reciprocal",
            source_template=CUBIC_TEMPLATE,
            target_template=target,
            q=t,
            max_stride=4,
            stages=3,
        )
        == []
    )


def test_cli_research_writes_note(tmp_path: Path):
    verified = tmp_path / "verified.jsonl"
    output_path = tmp_path / "note.md"

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
                "research",
                "--in",
                str(verified),
                "--candidate-id",
                "hero",
                "--depth",
                "12",
                "--series-order",
                "61",
                "--out",
                str(output_path),
            ]
        )
        == 0
    )
    text = output_path.read_text(encoding="utf-8")
    assert "Research Note: `hero`" in text
    assert "Euler Product Exponents" in text
    assert "Direct 1-Step Bauer-Muir Obstruction" in text
    assert "Page-43 Monomial Substitution Check" in text
    assert "Cubic Odd/Even Contraction Check" in text
    assert "Arithmetic Subsequence Contraction Scan" in text
    assert "Constrained Bauer-Muir Search" in text
    assert "Heine hcf2 Specialization Check (c=bz=-1)" in text
