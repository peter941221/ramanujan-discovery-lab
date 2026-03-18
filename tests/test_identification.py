from pathlib import Path

import pytest
import sympy as sp

from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.cli import main
from ramanujan_discovery.identification import (
    benchmark_power_substitution_series,
    signed_argument_substitution_series,
    scan_explicit_source_family_eta_correction_templates,
    scan_quotient_core_source_family_eta_corrections,
    scan_source_family_eta_corrections,
    scan_gg_modular_equation_box,
    scan_self_fractional_linear_uniqueness_relations,
    scan_self_polynomial_uniqueness_relations,
    scan_source_correction_self_fractional_linear_uniqueness_relations,
    scan_source_correction_self_polynomial_uniqueness_relations,
    scan_two_quotient_core_source_family_self_fractional_linear_relations,
    scan_two_quotient_core_source_family_self_eta_corrections,
    scan_two_quotient_core_source_family_self_polynomial_relations,
    scan_two_quotient_core_source_family_self_quotient_products,
    scan_two_quotient_core_source_family_eta_corrections,
    scan_two_core_source_family_eta_corrections,
    scan_named_fractional_linear_prefixes,
    scan_benchmark_power_relation_prefixes,
    scan_named_multiplicative_prefixes,
    scan_named_polynomial_prefixes,
    scan_named_prefix_boxes,
    scan_explicit_source_family_transform_templates,
    scan_parameterized_source_family_power_boxes,
    scan_named_two_layer_fractional_linear_prefixes,
    scan_ratio_benchmark_fractional_linear_prefixes,
    scan_ratio_benchmark_multiplicative_prefixes,
    scan_ratio_benchmark_power_relation_prefixes,
    scan_ratio_benchmark_two_layer_fractional_linear_prefixes,
    scan_ratio_eta_quotient_relations,
    scan_ratio_self_quotient_product_relations,
    search_eta_quotient_relation,
    search_fractional_linear_relation,
    search_multiplicative_relation,
    search_polynomial_relation,
    search_self_polynomial_uniqueness_relation,
    search_self_quotient_product_relation,
    search_self_t_polynomial_fractional_linear_relation,
    search_two_layer_fractional_linear_relation,
)
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.series import (
    continued_fraction_series_coeffs,
    series_div,
    series_invert,
    series_mul,
    series_pow,
)
from ramanujan_discovery.storage import write_candidates


def _eval_relation_series(relation, x_series, y_series, order: int):
    x_series = x_series[:order]
    y_series = y_series[:order]
    max_degree = relation.max_total_degree

    x_pows = [series_pow(x_series, i) for i in range(max_degree + 1)]
    y_pows = [series_pow(y_series, j) for j in range(max_degree + 1)]

    residual = [sp.Integer(0) for _ in range(order)]
    for (i, j), coeff in relation.coefficients.items():
        term = series_mul(x_pows[i], y_pows[j])
        for n in range(order):
            if term[n] == 0:
                continue
            residual[n] = sp.simplify(residual[n] + coeff * term[n])
    return residual


def _eval_three_variable_relation_series(relation, t_series, f_series, g_series, order: int):
    t_series = t_series[:order]
    f_series = f_series[:order]
    g_series = g_series[:order]
    max_t_degree = max(exponents[0] for exponents in relation.coefficients)
    max_f_degree = max(exponents[1] for exponents in relation.coefficients)
    max_g_degree = max(exponents[2] for exponents in relation.coefficients)

    t_pows = [series_pow(t_series, i) for i in range(max_t_degree + 1)]
    f_pows = [series_pow(f_series, i) for i in range(max_f_degree + 1)]
    g_pows = [series_pow(g_series, i) for i in range(max_g_degree + 1)]

    residual = [sp.Integer(0) for _ in range(order)]
    for (t_degree, f_degree, g_degree), coeff in relation.coefficients.items():
        term = series_mul(series_mul(t_pows[t_degree], f_pows[f_degree]), g_pows[g_degree])
        for n in range(order):
            if term[n] == 0:
                continue
            residual[n] = sp.simplify(residual[n] + coeff * term[n])
    return residual


def _one_minus_power_series(power: int, order: int):
    factor = [sp.Integer(0) for _ in range(order)]
    factor[0] = sp.Integer(1)
    if power < order:
        factor[power] = sp.Integer(-1)
    return factor


def _build_self_quotient_product_target(*, modulus: int, exponents_by_residue: dict[int, int], order: int):
    product_series = [sp.Integer(0) for _ in range(order)]
    product_series[0] = sp.Integer(1)
    for residue, exponent in sorted(exponents_by_residue.items()):
        factor = _one_minus_power_series(residue, order)
        if exponent >= 0:
            product_series = series_mul(product_series, series_pow(factor, exponent))
        else:
            product_series = series_mul(product_series, series_invert(series_pow(factor, -exponent)))

    target = [sp.Integer(0) for _ in range(order)]
    target[0] = sp.Integer(1)
    current = product_series
    while any(sp.simplify(value) != 0 for value in current[1:]):
        target = series_mul(target, current)
        current = benchmark_power_substitution_series(current, power=modulus, order=order)
    return target


def _build_self_eta_target(*, modulus: int, eta_series, order: int):
    target = [sp.Integer(0) for _ in range(order)]
    target[0] = sp.Integer(1)
    current = eta_series
    while any(sp.simplify(value) != 0 for value in current[1:]):
        target = series_mul(target, current)
        current = benchmark_power_substitution_series(current, power=modulus, order=order)
    return target


def _build_self_fractional_linear_eta_target(
    *,
    modulus: int,
    self_coeff: sp.Expr,
    eta_coeff: sp.Expr,
    eta_series,
    order: int,
):
    target = [sp.Integer(0) for _ in range(order)]
    target[0] = sp.Integer(1)
    eta_shifted = [sp.simplify(value) for value in eta_series]
    eta_shifted[0] = sp.Integer(0)
    for n in range(1, order):
        value = sp.simplify(eta_coeff * eta_shifted[n])
        if n % modulus == 0:
            value = sp.simplify(value + self_coeff * target[n // modulus])
        target[n] = value
    return target


def _eta_pochhammer_series(divisor: int, order: int):
    series = [sp.Integer(0) for _ in range(order)]
    series[0] = sp.Integer(1)
    for power in range(divisor, order, divisor):
        series = series_mul(series, _one_minus_power_series(power, order))
    return series


def test_search_bivariate_relation_finds_linear_identity():
    # x = 1 + q, y = 1 - q => x + y - 2 == 0.
    order = 8
    x = [sp.Integer(0) for _ in range(order)]
    y = [sp.Integer(0) for _ in range(order)]
    x[0] = 1
    x[1] = 1
    y[0] = 1
    y[1] = -1

    relation = search_polynomial_relation(series_by_variable={"x": x, "y": y}, max_total_degree=1, order=order)
    assert relation is not None
    # Map back to our helper which assumes variable order [x,y].
    # relation.variables preserves insertion order from the dict above.
    assert relation.variables == ("x", "y")
    residual = _eval_relation_series(relation, x, y, order)
    assert all(sp.simplify(value) == 0 for value in residual)


def test_search_bivariate_relation_returns_none_when_no_relation_in_box():
    # Choose sparse low-order series so the columns {1, x, y, x*y} are independent.
    order = 8
    x = [sp.Integer(0) for _ in range(order)]
    y = [sp.Integer(0) for _ in range(order)]
    x[0] = 1
    x[1] = 1
    x[3] = 1
    y[0] = 1
    y[2] = 1
    y[4] = 1

    relation = search_polynomial_relation(series_by_variable={"x": x, "y": y}, max_total_degree=1, order=order)
    assert relation is None


def test_search_relation_can_require_variable_to_appear():
    # y and z are identical, so a relation exists (y - z == 0) that does not involve x.
    order = 8
    x = [sp.Integer(0) for _ in range(order)]
    y = [sp.Integer(0) for _ in range(order)]
    z = [sp.Integer(0) for _ in range(order)]
    x[0] = 1
    x[1] = 1
    y[0] = 1
    y[2] = 1
    z[0] = 1
    z[2] = 1

    relation = search_polynomial_relation(series_by_variable={"x": x, "y": y, "z": z}, max_total_degree=1, order=order)
    assert relation is not None

    relation_x = search_polynomial_relation(
        series_by_variable={"x": x, "y": y, "z": z},
        max_total_degree=1,
        order=order,
        required_variable="x",
    )
    assert relation_x is None


def test_search_relation_raises_when_underdetermined():
    order = 10
    x = [sp.Integer(0) for _ in range(order)]
    y = [sp.Integer(0) for _ in range(order)]
    z = [sp.Integer(0) for _ in range(order)]
    w = [sp.Integer(0) for _ in range(order)]
    x[0] = 1
    y[0] = 1
    z[0] = 1
    w[0] = 1

    with pytest.raises(ValueError, match="underdetermined"):
        search_polynomial_relation(
            series_by_variable={"x": x, "y": y, "z": z, "w": w},
            max_total_degree=3,
            order=order,
        )


def test_scan_benchmark_power_relation_prefixes_finds_identity_relation():
    order = 10
    candidate = [sp.Integer(0) for _ in range(order)]
    benchmark = [sp.Integer(0) for _ in range(order)]
    candidate[0] = 1
    candidate[1] = 1
    benchmark[0] = 1
    benchmark[1] = 1

    scans = scan_benchmark_power_relation_prefixes(
        candidate_recip=candidate,
        benchmark_recip=benchmark,
        powers=(2, 3, 4),
        order=order,
        degree_values=(1, 2),
        required_variable="C",
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.powers == (2,)
    assert first_hit.max_total_degree == 1
    assert first_hit.relation is not None
    assert first_hit.relation.variables == ("C", "B1", "B2")


def test_scan_ratio_benchmark_power_relation_prefixes_finds_identity_relation():
    order = 10
    ratio = [sp.Integer(0) for _ in range(order)]
    benchmark = [sp.Integer(0) for _ in range(order)]
    ratio[0] = 1
    ratio[1] = 1
    benchmark[0] = 1
    benchmark[1] = 1

    scans = scan_ratio_benchmark_power_relation_prefixes(
        ratio_series=ratio,
        benchmark_series=benchmark,
        powers=(2, 3, 4),
        order=order,
        degree_values=(1, 2),
        required_variable="F",
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.powers == (2,)
    assert first_hit.max_total_degree == 1
    assert first_hit.relation is not None
    assert first_hit.relation.variables == ("F", "B1", "B2")


def test_search_fractional_linear_relation_finds_structured_ratio():
    order = 10
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1
    benchmark_q2 = benchmark_power_substitution_series(benchmark, power=2, order=order)

    numerator = [sp.Integer(0) for _ in range(order)]
    denominator = [sp.Integer(0) for _ in range(order)]
    numerator[0] = 1
    denominator[0] = 1
    numerator[1] = 2
    denominator[2] = -1
    ratio = series_mul(numerator, series_invert(denominator))

    relation = search_fractional_linear_relation(
        target_series=ratio,
        basis_series_by_variable={"B1": benchmark, "B2": benchmark_q2},
        order=order,
    )
    assert relation is not None
    assert relation.basis_variables == ("B1", "B2")
    assert relation.numerator_coefficients == {"B1": sp.Integer(2)}
    assert relation.denominator_coefficients == {"B2": sp.Integer(-1)}


def test_search_self_polynomial_uniqueness_relation_finds_tfg_relation():
    order = 16
    target = [sp.Integer(0) for _ in range(order)]
    target[0] = 1
    for exponent in (1, 2, 4, 8):
        target[exponent] = 1

    relation = search_self_polynomial_uniqueness_relation(
        target_series=target,
        modulus=2,
        order=order,
        max_fg_total_degree=1,
        max_t_degree=1,
    )
    assert relation is not None
    assert relation.variables == ("T", "F", "G2")

    t_series = [sp.Integer(0) for _ in range(order)]
    t_series[1] = 1
    residual = _eval_three_variable_relation_series(
        relation,
        t_series,
        target,
        benchmark_power_substitution_series(target, power=2, order=order),
        order,
    )
    assert all(sp.simplify(value) == 0 for value in residual)


def test_search_self_t_polynomial_fractional_linear_relation_finds_rhs_equation():
    order = 20
    eta = [sp.Integer(0) for _ in range(order)]
    eta[0] = 1
    eta[1] = 1
    target = _build_self_fractional_linear_eta_target(
        modulus=2,
        self_coeff=sp.Integer(2),
        eta_coeff=sp.Integer(1),
        eta_series=eta,
        order=order,
    )

    relation = search_self_t_polynomial_fractional_linear_relation(
        target_series=target,
        modulus=2,
        order=order,
        max_t_degree=1,
    )
    assert relation is not None
    assert relation.numerator_t_coefficients == (sp.Integer(1), sp.Integer(1))
    assert relation.numerator_self_coefficients == (sp.Integer(2), sp.Integer(0))
    assert relation.denominator_t_coefficients == (sp.Integer(1), sp.Integer(0))
    assert relation.denominator_self_coefficients == (sp.Integer(0), sp.Integer(0))


def test_rhs_uniqueness_searches_skip_degenerate_non_self_relations():
    order = 12
    target = [sp.Integer(0) for _ in range(order)]
    target[0] = 1
    target[2] = 1

    assert (
        search_self_polynomial_uniqueness_relation(
            target_series=target,
            modulus=2,
            order=order,
            max_fg_total_degree=1,
            max_t_degree=1,
        )
        is None
    )
    assert (
        search_self_t_polynomial_fractional_linear_relation(
            target_series=target,
            modulus=2,
            order=order,
            max_t_degree=1,
        )
        is None
    )

    polynomial_scan = scan_self_polynomial_uniqueness_relations(
        target_series=target,
        moduli=(2, 3, 4),
        order=order,
        fg_degree_values=(1,),
        t_degree_values=(1,),
    )
    fractional_scan = scan_self_fractional_linear_uniqueness_relations(
        target_series=target,
        moduli=(2, 3, 4),
        order=order,
        t_degree_values=(1,),
    )
    assert polynomial_scan.hits == ()
    assert fractional_scan.hits == ()


def test_scan_source_correction_self_polynomial_uniqueness_relations_finds_one_core_hit():
    order = 16
    rr = [sp.Integer(0) for _ in range(order)]
    rr[0] = 1
    rr[1] = 1
    cubic = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[2] = 1

    correction = [sp.Integer(0) for _ in range(order)]
    correction[0] = 1
    for exponent in (1, 2, 4, 8):
        if exponent < order:
            correction[exponent] = 1
    target = series_mul(rr, correction)

    scan = scan_source_correction_self_polynomial_uniqueness_relations(
        target_series=target,
        ordered_base_families=(
            ("RR", "rogers_ramanujan_normalized", rr),
            ("cubic", "ramanujan_cubic_normalized", cubic),
        ),
        correction_size=1,
        moduli=(2, 3),
        order=order,
        fg_degree_values=(1,),
        t_degree_values=(1,),
    )
    assert scan.total_corrections_checked == 2
    assert scan.hits
    hit = next(item for item in scan.hits if item.basis_labels == ("RR",) and item.modulus == 2)
    assert hit.relation.variables == ("T", "F", "G2")


def test_scan_source_correction_self_fractional_linear_uniqueness_relations_finds_two_core_hit():
    order = 20
    rr = [sp.Integer(0) for _ in range(order)]
    rr[0] = 1
    rr[1] = 1
    cubic = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[2] = 1
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[3] = 1

    eta = [sp.Integer(0) for _ in range(order)]
    eta[0] = 1
    eta[1] = 1
    correction = _build_self_fractional_linear_eta_target(
        modulus=2,
        self_coeff=sp.Integer(2),
        eta_coeff=sp.Integer(1),
        eta_series=eta,
        order=order,
    )
    target = series_mul(series_mul(rr, cubic), correction)

    scan = scan_source_correction_self_fractional_linear_uniqueness_relations(
        target_series=target,
        ordered_base_families=(
            ("RR", "rogers_ramanujan_normalized", rr),
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        correction_size=2,
        moduli=(2, 3),
        order=order,
        t_degree_values=(1,),
    )
    assert scan.total_corrections_checked == 3
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.basis_labels == ("RR", "cubic") and item.modulus == 2
    )
    assert hit.relation.numerator_t_coefficients == (sp.Integer(1), sp.Integer(1))
    assert hit.relation.numerator_self_coefficients == (sp.Integer(2), sp.Integer(0))


def test_scan_ratio_benchmark_fractional_linear_prefixes_finds_identity_relation():
    order = 10
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1

    numerator = [sp.Integer(0) for _ in range(order)]
    denominator = [sp.Integer(0) for _ in range(order)]
    numerator[0] = 1
    denominator[0] = 1
    numerator[1] = 2
    denominator[2] = -1
    ratio = series_mul(numerator, series_invert(denominator))

    scans = scan_ratio_benchmark_fractional_linear_prefixes(
        ratio_series=ratio,
        benchmark_series=benchmark,
        powers=(2, 3, 4),
        order=order,
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.powers == (2,)
    assert first_hit.relation is not None
    assert first_hit.relation.numerator_coefficients == {"B1": sp.Integer(2)}
    assert first_hit.relation.denominator_coefficients == {"B2": sp.Integer(-1)}


def test_search_multiplicative_relation_finds_integer_exponents():
    order = 12
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1
    benchmark_q2 = benchmark_power_substitution_series(benchmark, power=2, order=order)

    ratio = series_mul(series_pow(benchmark, 2), series_invert(benchmark_q2))

    relation = search_multiplicative_relation(
        target_series=ratio,
        basis_series_by_variable={"B1": benchmark, "B2": benchmark_q2},
        order=order,
    )
    assert relation is not None
    assert relation.basis_variables == ("B1", "B2")
    assert relation.exponents == {"B1": 2, "B2": -1}


def test_scan_ratio_benchmark_multiplicative_prefixes_finds_hit():
    order = 12
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1

    benchmark_q2 = benchmark_power_substitution_series(benchmark, power=2, order=order)
    ratio = series_mul(series_pow(benchmark, 2), series_invert(benchmark_q2))

    scans = scan_ratio_benchmark_multiplicative_prefixes(
        ratio_series=ratio,
        benchmark_series=benchmark,
        powers=(2, 3, 4),
        order=order,
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.powers == (2,)
    assert first_hit.relation is not None
    assert first_hit.relation.exponents == {"B1": 2, "B2": -1}


def test_search_self_quotient_product_relation_finds_periodic_product_identity():
    order = 20
    expected = {1: 1, 2: -1, 3: -1, 4: 1}
    ratio = _build_self_quotient_product_target(modulus=5, exponents_by_residue=expected, order=order)

    relation = search_self_quotient_product_relation(
        target_series=ratio,
        modulus=5,
        order=order,
    )
    assert relation is not None
    assert relation.modulus == 5
    assert relation.exponents_by_residue == expected


def test_scan_ratio_self_quotient_product_relations_hits_matching_modulus():
    order = 20
    expected = {1: 1, 2: -1, 3: -1, 4: 1}
    ratio = _build_self_quotient_product_target(modulus=5, exponents_by_residue=expected, order=order)

    scans = scan_ratio_self_quotient_product_relations(
        ratio_series=ratio,
        moduli=(2, 5, 7),
        order=order,
    )
    assert scans
    hit = next(scan for scan in scans if scan.modulus == 5)
    assert hit.relation is not None
    assert hit.relation.exponents_by_residue == expected


def test_search_eta_quotient_relation_finds_small_level_hit():
    order = 18
    e1 = _eta_pochhammer_series(1, order)
    e2 = _eta_pochhammer_series(2, order)
    ratio = series_mul(series_pow(e1, 2), series_invert(e2))

    relation = search_eta_quotient_relation(
        target_series=ratio,
        level=2,
        order=order,
    )
    assert relation is not None
    assert relation.basis_variables == ("E1", "E2")
    assert relation.exponents == {"E1": 2, "E2": -1}


def test_scan_ratio_eta_quotient_relations_hits_matching_level():
    order = 18
    e1 = _eta_pochhammer_series(1, order)
    e2 = _eta_pochhammer_series(2, order)
    ratio = series_mul(series_pow(e1, 2), series_invert(e2))

    scans = scan_ratio_eta_quotient_relations(
        ratio_series=ratio,
        levels=(2, 3, 4),
        order=order,
    )
    assert scans
    hit = next(scan for scan in scans if scan.level == 2)
    assert hit.relation is not None
    assert hit.relation.exponents == {"E1": 2, "E2": -1}


def test_scan_source_family_eta_corrections_finds_power_times_eta_hit():
    order = 18
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    e1 = _eta_pochhammer_series(1, order)
    e2 = _eta_pochhammer_series(2, order)
    eta_ratio = series_mul(series_pow(e1, 2), series_invert(e2))
    target = series_mul(gg2, eta_ratio)

    scans = scan_source_family_eta_corrections(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        eta_levels=(2,),
        order=order,
    )
    assert len(scans) == 1
    family_scan = scans[0]
    gg2_scan = next(scan for scan in family_scan.direct_basis_scans if scan.basis_label == "GG2")
    level_2_hit = next(scan for scan in gg2_scan.eta_scans if scan.level == 2)
    assert level_2_hit.relation is not None
    assert level_2_hit.relation.exponents == {"E1": 2, "E2": -1}


def test_scan_two_core_source_family_eta_corrections_finds_cross_family_hit():
    order = 18
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    eta_ratio = series_mul(
        series_pow(_eta_pochhammer_series(1, order), 2),
        series_invert(_eta_pochhammer_series(2, order)),
    )
    target = series_mul(series_mul(cubic2, series_invert(gg2)), eta_ratio)

    scan = scan_two_core_source_family_eta_corrections(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        eta_levels=(2,),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = scan.hits[0]
    assert hit.basis_labels == ("cubic2", "GG2")
    assert hit.level == 2
    assert hit.relation.exponents == {"cubic2": 1, "GG2": -1, "E1": 2, "E2": -1}


def test_scan_quotient_core_source_family_eta_corrections_finds_cross_family_hit():
    order = 18
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    gg_q2 = series_div(gg2, gg)
    eta_ratio = series_mul(
        series_pow(_eta_pochhammer_series(1, order), 2),
        series_invert(_eta_pochhammer_series(2, order)),
    )
    target = series_mul(series_mul(gg_q2, cubic), eta_ratio)

    scan = scan_quotient_core_source_family_eta_corrections(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        eta_levels=(2,),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_label == "GG_Q2" and item.raw_label == "cubic"
    )
    assert hit.quotient_expression == "GG2 / GG"
    assert hit.level == 2
    assert hit.relation.exponents == {"GG_Q2": 1, "cubic": 1, "E1": 2, "E2": -1}


def test_scan_two_quotient_core_source_family_eta_corrections_finds_cross_family_hit():
    order = 18
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    cubic_q2 = series_div(cubic2, cubic)
    gg_q2 = series_div(gg2, gg)
    eta_ratio = series_mul(
        series_pow(_eta_pochhammer_series(1, order), 2),
        series_invert(_eta_pochhammer_series(2, order)),
    )
    target = series_mul(series_mul(cubic_q2, gg_q2), eta_ratio)

    scan = scan_two_quotient_core_source_family_eta_corrections(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        eta_levels=(2,),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_labels == ("cubic_Q2", "GG_Q2")
    )
    assert hit.quotient_expressions == ("cubic2 / cubic", "GG2 / GG")
    assert hit.level == 2
    assert hit.relation.exponents == {"cubic_Q2": 1, "GG_Q2": 1, "E1": 2, "E2": -1}


def test_scan_two_quotient_core_source_family_self_quotient_products_finds_cross_family_hit():
    order = 20
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    cubic_q2 = series_div(cubic2, cubic)
    gg_q2 = series_div(gg2, gg)
    correction = _build_self_quotient_product_target(
        modulus=2,
        exponents_by_residue={1: 2},
        order=order,
    )
    target = series_mul(series_mul(cubic_q2, gg_q2), correction)

    scan = scan_two_quotient_core_source_family_self_quotient_products(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        moduli=(2, 3, 4),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_labels == ("cubic_Q2", "GG_Q2")
    )
    assert hit.modulus == 2
    assert hit.relation.exponents_by_residue == {1: 2}


def test_scan_two_quotient_core_source_family_self_polynomial_relations_finds_cross_family_hit():
    order = 16
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    cubic_q2 = series_div(cubic2, cubic)
    gg_q2 = series_div(gg2, gg)
    correction = [sp.Integer(0) for _ in range(order)]
    correction[0] = 1
    correction[1] = 1
    target = series_mul(series_mul(cubic_q2, gg_q2), correction)

    scan = scan_two_quotient_core_source_family_self_polynomial_relations(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        moduli=(2, 3),
        order=order,
        degree_values=(1, 2),
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_labels == ("cubic_Q2", "GG_Q2") and item.modulus == 2
    )
    assert hit.max_total_degree == 2
    assert hit.relation.variables == ("G", "G2")
    residual = _eval_relation_series(
        hit.relation,
        correction,
        benchmark_power_substitution_series(correction, power=2, order=order),
        order,
    )
    assert all(sp.simplify(value) == 0 for value in residual)


def test_scan_two_quotient_core_source_family_self_eta_corrections_finds_cross_family_hit():
    order = 20
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    cubic_q2 = series_div(cubic2, cubic)
    gg_q2 = series_div(gg2, gg)
    eta_ratio = series_mul(series_pow(_eta_pochhammer_series(1, order), 2), series_invert(_eta_pochhammer_series(2, order)))
    correction = _build_self_eta_target(modulus=2, eta_series=eta_ratio, order=order)
    target = series_mul(series_mul(cubic_q2, gg_q2), correction)

    scan = scan_two_quotient_core_source_family_self_eta_corrections(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        moduli=(2, 3),
        eta_levels=(2, 3),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_labels == ("cubic_Q2", "GG_Q2") and item.modulus == 2 and item.level == 2
    )
    assert hit.relation.exponents == {"G2": 1, "E1": 2, "E2": -1}


def test_scan_two_quotient_core_source_family_self_fractional_linear_relations_finds_cross_family_hit():
    order = 20
    cubic = [sp.Integer(0) for _ in range(order)]
    gg = [sp.Integer(0) for _ in range(order)]
    cubic[0] = 1
    cubic[1] = 1
    gg[0] = 1
    gg[2] = 1

    cubic2 = benchmark_power_substitution_series(cubic, power=2, order=order)
    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    cubic_q2 = series_div(cubic2, cubic)
    gg_q2 = series_div(gg2, gg)
    e1 = _eta_pochhammer_series(1, order)
    correction = _build_self_fractional_linear_eta_target(
        modulus=2,
        self_coeff=sp.Integer(2),
        eta_coeff=sp.Integer(3),
        eta_series=e1,
        order=order,
    )
    target = series_mul(series_mul(cubic_q2, gg_q2), correction)

    scan = scan_two_quotient_core_source_family_self_fractional_linear_relations(
        target_series=target,
        ordered_base_families=(
            ("cubic", "ramanujan_cubic_normalized", cubic),
            ("GG", "gollnitz_gordon_normalized", gg),
        ),
        powers=(2, 3),
        moduli=(2, 3),
        eta_levels=(1, 2),
        order=order,
    )
    assert scan.total_basis_pairs_checked > 0
    assert scan.hits
    hit = next(
        item
        for item in scan.hits
        if item.quotient_labels == ("cubic_Q2", "GG_Q2") and item.modulus == 2 and item.level == 1
    )
    assert hit.relation.numerator_coefficients == {"G2": sp.Integer(2), "E1": sp.Integer(3)}
    assert hit.relation.denominator_coefficients == {}


def test_scan_named_multiplicative_prefixes_finds_hit():
    order = 12
    rr = [sp.Integer(0) for _ in range(order)]
    cubic = [sp.Integer(0) for _ in range(order)]
    rr[0] = 1
    rr[1] = 1
    cubic[0] = 1
    cubic[2] = 1

    ratio = series_mul(series_pow(rr, 2), series_invert(cubic))

    scans = scan_named_multiplicative_prefixes(
        target_series=ratio,
        ordered_basis_series=(("RR", rr), ("cubic", cubic)),
        order=order,
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.basis_labels == ("RR", "cubic")
    assert first_hit.relation is not None
    assert first_hit.relation.exponents == {"RR": 2, "cubic": -1}


def test_scan_named_fractional_linear_prefixes_finds_hit():
    order = 12
    rr = [sp.Integer(0) for _ in range(order)]
    cubic = [sp.Integer(0) for _ in range(order)]
    rr[0] = 1
    rr[1] = 1
    cubic[0] = 1
    cubic[2] = 1

    numerator = [sp.Integer(0) for _ in range(order)]
    denominator = [sp.Integer(0) for _ in range(order)]
    numerator[0] = 1
    denominator[0] = 1
    numerator[1] = 2
    denominator[2] = -1
    ratio = series_mul(numerator, series_invert(denominator))

    scans = scan_named_fractional_linear_prefixes(
        target_series=ratio,
        ordered_basis_series=(("RR", rr), ("cubic", cubic)),
        order=order,
    )
    assert scans
    assert any(scan.relation is not None for scan in scans)
    first_hit = next(scan for scan in scans if scan.relation is not None)
    assert first_hit.basis_labels == ("RR", "cubic")
    assert first_hit.relation is not None
    assert first_hit.relation.numerator_coefficients == {"RR": sp.Integer(2)}
    assert first_hit.relation.denominator_coefficients == {"cubic": sp.Integer(-1)}


def test_search_two_layer_fractional_linear_relation_finds_structured_ratio():
    order = 14
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1
    benchmark_q2 = benchmark_power_substitution_series(benchmark, power=2, order=order)

    u1 = [sp.Integer(0) for _ in range(order)]
    u2 = [sp.Integer(0) for _ in range(order)]
    u1[1] = 1
    u2[2] = 1

    num_factor_1 = [sp.Integer(0) for _ in range(order)]
    den_factor_1 = [sp.Integer(0) for _ in range(order)]
    num_factor_2 = [sp.Integer(0) for _ in range(order)]
    den_factor_2 = [sp.Integer(0) for _ in range(order)]
    num_factor_1[0] = 1
    den_factor_1[0] = 1
    num_factor_2[0] = 1
    den_factor_2[0] = 1
    num_factor_1[1] = 2
    den_factor_1[2] = -1
    num_factor_2[2] = 3
    den_factor_2[1] = 4

    ratio = series_mul(
        series_mul(num_factor_1, num_factor_2),
        series_invert(series_mul(den_factor_1, den_factor_2)),
    )

    relation = search_two_layer_fractional_linear_relation(
        target_series=ratio,
        basis_series_by_variable={"B1": benchmark, "B2": benchmark_q2},
        numerator_variables=("B1", "B2"),
        denominator_variables=("B2", "B1"),
        order=order,
        solve_order=10,
    )
    assert relation is not None
    assert relation.numerator_variables == ("B1", "B2")
    assert relation.denominator_variables == ("B2", "B1")
    assert relation.numerator_coefficients == (sp.Integer(2), sp.Integer(3))
    assert relation.denominator_coefficients == (sp.Integer(-1), sp.Integer(4))


def test_scan_ratio_benchmark_two_layer_fractional_linear_prefixes_finds_hit():
    order = 14
    benchmark = [sp.Integer(0) for _ in range(order)]
    benchmark[0] = 1
    benchmark[1] = 1

    num_factor_1 = [sp.Integer(0) for _ in range(order)]
    den_factor_1 = [sp.Integer(0) for _ in range(order)]
    num_factor_2 = [sp.Integer(0) for _ in range(order)]
    den_factor_2 = [sp.Integer(0) for _ in range(order)]
    num_factor_1[0] = 1
    den_factor_1[0] = 1
    num_factor_2[0] = 1
    den_factor_2[0] = 1
    num_factor_1[1] = 2
    den_factor_1[2] = -1
    num_factor_2[2] = 3
    den_factor_2[1] = 4

    ratio = series_mul(
        series_mul(num_factor_1, num_factor_2),
        series_invert(series_mul(den_factor_1, den_factor_2)),
    )

    scans = scan_ratio_benchmark_two_layer_fractional_linear_prefixes(
        ratio_series=ratio,
        benchmark_series=benchmark,
        powers=(2, 3, 4),
        order=order,
        solve_order=10,
    )
    assert scans
    assert any(scan.total_hits > 0 for scan in scans)
    first_hit = next(scan for scan in scans if scan.total_hits > 0)
    assert first_hit.powers == (2,)
    assert first_hit.relations
    assert first_hit.total_hits >= 1
    assert all(
        variable in {"B1", "B2"}
        for relation in first_hit.relations
        for variable in relation.numerator_variables + relation.denominator_variables
    )


def test_scan_named_two_layer_fractional_linear_prefixes_finds_hit():
    order = 14
    rr = [sp.Integer(0) for _ in range(order)]
    cubic = [sp.Integer(0) for _ in range(order)]
    rr[0] = 1
    rr[1] = 1
    cubic[0] = 1
    cubic[2] = 1

    num_factor_1 = [sp.Integer(0) for _ in range(order)]
    den_factor_1 = [sp.Integer(0) for _ in range(order)]
    num_factor_2 = [sp.Integer(0) for _ in range(order)]
    den_factor_2 = [sp.Integer(0) for _ in range(order)]
    num_factor_1[0] = 1
    den_factor_1[0] = 1
    num_factor_2[0] = 1
    den_factor_2[0] = 1
    num_factor_1[1] = 2
    den_factor_1[2] = -1
    num_factor_2[2] = 3
    den_factor_2[1] = 4

    ratio = series_mul(
        series_mul(num_factor_1, num_factor_2),
        series_invert(series_mul(den_factor_1, den_factor_2)),
    )

    scans = scan_named_two_layer_fractional_linear_prefixes(
        target_series=ratio,
        ordered_basis_series=(("RR", rr), ("cubic", cubic)),
        order=order,
        solve_order=10,
    )
    assert scans
    assert any(scan.total_hits > 0 for scan in scans)
    first_hit = next(scan for scan in scans if scan.total_hits > 0)
    assert first_hit.basis_labels == ("RR", "cubic")
    assert first_hit.relations
    assert first_hit.total_hits >= 1
    assert all(
        variable in {"RR", "cubic"}
        for relation in first_hit.relations
        for variable in relation.numerator_variables + relation.denominator_variables
    )


def test_scan_named_prefix_boxes_matches_individual_scans():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    target = benchmark_power_substitution_series(gg, power=2, order=order)
    ordered_basis_series = (
        ("GG", gg),
        ("GG2", benchmark_power_substitution_series(gg, power=2, order=order)),
        ("GG3", benchmark_power_substitution_series(gg, power=3, order=order)),
    )

    bundle = scan_named_prefix_boxes(
        target_series=target,
        ordered_basis_series=ordered_basis_series,
        order=order,
        degree_values=(1, 2),
        max_abs_exponent=4,
        solve_order=10,
    )

    assert bundle.polynomial_scans == tuple(
        scan_named_polynomial_prefixes(
            target_series=target,
            ordered_basis_series=ordered_basis_series,
            order=order,
            degree_values=(1, 2),
        )
    )
    assert bundle.multiplicative_scans == tuple(
        scan_named_multiplicative_prefixes(
            target_series=target,
            ordered_basis_series=ordered_basis_series,
            order=order,
            max_abs_exponent=4,
        )
    )
    assert bundle.fractional_linear_scans == tuple(
        scan_named_fractional_linear_prefixes(
            target_series=target,
            ordered_basis_series=ordered_basis_series,
            order=order,
        )
    )
    assert bundle.two_layer_fractional_linear_scans == tuple(
        scan_named_two_layer_fractional_linear_prefixes(
            target_series=target,
            ordered_basis_series=ordered_basis_series,
            order=order,
            solve_order=10,
        )
    )


def test_scan_parameterized_source_family_power_boxes_finds_power_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    target = benchmark_power_substitution_series(gg, power=2, order=order)
    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert family_scan.family_label == "GG"
    assert family_scan.ordered_basis_series[1][0] == "GG2"
    assert any(scan.relation is not None for scan in family_scan.polynomial_scans)
    polynomial_hit = next(
        scan
        for scan in family_scan.polynomial_scans
        if scan.relation is not None and scan.basis_labels[-1] == "GG2" and scan.max_total_degree == 1
    )
    assert polynomial_hit.basis_labels == ("GG", "GG2")
    assert polynomial_hit.max_total_degree == 1
    assert any(scan.relation is not None for scan in family_scan.multiplicative_scans)
    multiplicative_hit = next(scan for scan in family_scan.multiplicative_scans if scan.relation is not None)
    assert multiplicative_hit.basis_labels == ("GG", "GG2")
    assert multiplicative_hit.relation.exponents == {"GG2": 1}
    assert any(scan.relation is not None for scan in family_scan.fractional_linear_scans)
    assert family_scan.two_layer_fractional_linear_scans


def test_scan_parameterized_source_family_power_boxes_finds_two_layer_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    num_factor_1 = [sp.Integer(0) for _ in range(order)]
    den_factor_1 = [sp.Integer(0) for _ in range(order)]
    num_factor_2 = [sp.Integer(0) for _ in range(order)]
    den_factor_2 = [sp.Integer(0) for _ in range(order)]
    num_factor_1[0] = 1
    den_factor_1[0] = 1
    num_factor_2[0] = 1
    den_factor_2[0] = 1
    num_factor_1[1] = 2
    den_factor_1[2] = -1
    num_factor_2[2] = 3
    den_factor_2[1] = 4

    target = series_mul(
        series_mul(num_factor_1, num_factor_2),
        series_invert(series_mul(den_factor_1, den_factor_2)),
    )

    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
        solve_order=10,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert any(scan.total_hits > 0 for scan in family_scan.two_layer_fractional_linear_scans)
    first_hit = next(
        scan for scan in family_scan.two_layer_fractional_linear_scans if scan.total_hits > 0
    )
    assert first_hit.basis_labels == ("GG", "GG2")
    assert first_hit.relations
    assert first_hit.total_hits >= 1
    assert all(
        variable in {"GG", "GG2"}
        for relation in first_hit.relations
        for variable in relation.numerator_variables + relation.denominator_variables
    )


def test_scan_parameterized_source_family_power_boxes_finds_quotient_ladder_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    target = benchmark_power_substitution_series(gg, power=2, order=order)
    target = series_div(target, gg)
    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert tuple(label for label, _, _ in family_scan.quotient_basis_series) == ("Q2", "Q3")
    assert family_scan.quotient_basis_series[0][1] == "GG2 / GG"
    assert any(scan.relation is not None for scan in family_scan.quotient_polynomial_scans)
    quotient_polynomial_hit = next(
        scan
        for scan in family_scan.quotient_polynomial_scans
        if scan.relation is not None and scan.basis_labels[-1] == "Q2" and scan.max_total_degree == 1
    )
    assert quotient_polynomial_hit.basis_labels == ("Q2",)
    assert any(scan.relation is not None for scan in family_scan.quotient_multiplicative_scans)
    quotient_multiplicative_hit = next(
        scan for scan in family_scan.quotient_multiplicative_scans if scan.relation is not None
    )
    assert quotient_multiplicative_hit.basis_labels == ("Q2",)
    assert quotient_multiplicative_hit.relation.exponents == {"Q2": 1}
    assert any(scan.relation is not None for scan in family_scan.quotient_fractional_linear_scans)
    assert family_scan.quotient_two_layer_fractional_linear_scans


def test_scan_parameterized_source_family_power_boxes_finds_quotient_ladder_two_layer_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = series_div(benchmark_power_substitution_series(gg, power=2, order=order), gg)
    gg3 = series_div(benchmark_power_substitution_series(gg, power=3, order=order), gg)

    num_factor_1 = [sp.Integer(1)] + [sp.simplify(2 * value) for value in gg2[1:]]
    den_factor_1 = [sp.Integer(1)] + [sp.simplify(-value) for value in gg2[1:]]
    num_factor_2 = [sp.Integer(1)] + [sp.simplify(3 * value) for value in gg3[1:]]
    den_factor_2 = [sp.Integer(1)] + [sp.simplify(4 * value) for value in gg3[1:]]

    target = series_mul(
        series_mul(num_factor_1, num_factor_2),
        series_invert(series_mul(den_factor_1, den_factor_2)),
    )

    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
        solve_order=10,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert any(scan.total_hits > 0 for scan in family_scan.quotient_two_layer_fractional_linear_scans)
    first_hit = next(
        scan
        for scan in family_scan.quotient_two_layer_fractional_linear_scans
        if scan.total_hits > 0
    )
    assert first_hit.basis_labels == ("Q2", "Q3")
    assert first_hit.relations
    assert first_hit.total_hits >= 1
    assert all(
        variable in {"Q2", "Q3"}
        for relation in first_hit.relations
        for variable in relation.numerator_variables + relation.denominator_variables
    )


def test_scan_parameterized_source_family_power_boxes_finds_mixed_quotient_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    target = series_div(benchmark_power_substitution_series(gg, power=2, order=order), gg)
    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert tuple(label for label, _, _ in family_scan.mixed_quotient_basis_series) == ("GG", "Q2", "Q3")
    assert any(scan.relation is not None for scan in family_scan.mixed_quotient_polynomial_scans)
    mixed_polynomial_hit = next(
        scan
        for scan in family_scan.mixed_quotient_polynomial_scans
        if scan.relation is not None and scan.basis_labels[-1] == "Q2" and scan.max_total_degree == 1
    )
    assert mixed_polynomial_hit.basis_labels == ("GG", "Q2")
    assert any(scan.relation is not None for scan in family_scan.mixed_quotient_multiplicative_scans)
    mixed_multiplicative_hit = next(
        scan for scan in family_scan.mixed_quotient_multiplicative_scans if scan.relation is not None
    )
    assert mixed_multiplicative_hit.basis_labels == ("GG", "Q2")
    assert mixed_multiplicative_hit.relation.exponents == {"Q2": 1}
    assert any(scan.relation is not None for scan in family_scan.mixed_quotient_fractional_linear_scans)


def test_scan_parameterized_source_family_power_boxes_finds_mixed_quotient_two_layer_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = series_div(benchmark_power_substitution_series(gg, power=2, order=order), gg)
    gg3 = series_div(benchmark_power_substitution_series(gg, power=3, order=order), gg)

    factor_num_1 = [sp.Integer(1)] + [sp.simplify(2 * value) for value in gg[1:]]
    factor_den_1 = [sp.Integer(1)] + [sp.simplify(-value) for value in gg2[1:]]
    factor_num_2 = [sp.Integer(1)] + [sp.simplify(3 * value) for value in gg3[1:]]
    factor_den_2 = [sp.Integer(1)] + [sp.simplify(4 * value) for value in gg[1:]]

    target = series_mul(
        series_mul(factor_num_1, factor_num_2),
        series_invert(series_mul(factor_den_1, factor_den_2)),
    )

    scans = scan_parameterized_source_family_power_boxes(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
        solve_order=10,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert any(scan.total_hits > 0 for scan in family_scan.mixed_quotient_two_layer_fractional_linear_scans)
    first_hit = next(
        scan
        for scan in family_scan.mixed_quotient_two_layer_fractional_linear_scans
        if scan.total_hits > 0
    )
    assert first_hit.basis_labels == ("GG", "Q2", "Q3")
    assert first_hit.relations
    assert first_hit.total_hits >= 1
    assert all(
        variable in {"GG", "Q2", "Q3"}
        for relation in first_hit.relations
        for variable in relation.numerator_variables + relation.denominator_variables
    )


def test_scan_explicit_source_family_transform_templates_finds_quotient_hit():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    target = series_div(gg2, gg)
    scans = scan_explicit_source_family_transform_templates(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert family_scan.family_label == "GG"
    assert "GG2 / GG" in family_scan.checked_templates
    assert family_scan.hit_templates == ("GG2 / GG",)


def test_scan_explicit_source_family_eta_correction_templates_finds_quotient_times_eta_hit():
    order = 18
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    gg2 = benchmark_power_substitution_series(gg, power=2, order=order)
    eta_ratio = series_mul(series_pow(_eta_pochhammer_series(1, order), 2), series_invert(_eta_pochhammer_series(2, order)))
    target = series_mul(series_div(gg2, gg), eta_ratio)

    scans = scan_explicit_source_family_eta_correction_templates(
        target_series=target,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        eta_levels=(2,),
        order=order,
    )

    assert len(scans) == 1
    family_scan = scans[0]
    assert "GG2 / GG" in family_scan.checked_templates
    assert family_scan.hits
    first_hit = family_scan.hits[0]
    assert first_hit.template_label == "GG2 / GG"
    assert first_hit.level == 2
    assert first_hit.relation.exponents == {"E1": 2, "E2": -1}


def test_scan_gg_modular_equation_box_finds_signed_quotient_hit():
    order = 18
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1
    gg[2] = 2

    ggneg = signed_argument_substitution_series(gg, order=order)
    gg3 = benchmark_power_substitution_series(gg, power=3, order=order)
    target = series_div(gg3, ggneg)

    scan = scan_gg_modular_equation_box(
        target_series=target,
        benchmark_name="gollnitz_gordon_normalized",
        gg_series=gg,
        order=order,
        max_abs_exponent=4,
        solve_order=14,
    )

    assert scan.benchmark_name == "gollnitz_gordon_normalized"
    assert tuple(label for label, _, _ in scan.ordered_basis_series) == ("GG", "GGneg", "GG2", "GG3", "GG4")
    assert tuple(label for label, _, _ in scan.quotient_basis_series) == ("Q_neg", "Q_2", "Q_3", "Q_4")
    assert tuple(label for label, _, _ in scan.mixed_quotient_basis_series) == (
        "GG",
        "Q_neg",
        "Q_2",
        "Q_3",
        "Q_4",
    )
    assert "GG3 / GGneg" in scan.checked_templates
    assert scan.hit_templates == ("GG3 / GGneg",)
    mixed_hit = next(
        relation_scan
        for relation_scan in scan.mixed_quotient_multiplicative_scans
        if relation_scan.relation is not None
    )
    assert mixed_hit.basis_labels == ("GG", "Q_neg", "Q_2", "Q_3")
    assert mixed_hit.relation.exponents == {"Q_neg": -1, "Q_3": 1}


def test_scan_gg_modular_equation_box_supports_odd_prime_descendants():
    order = 8
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    scan = scan_gg_modular_equation_box(
        target_series=gg,
        benchmark_name="gollnitz_gordon_normalized",
        gg_series=gg,
        order=order,
        degree_values=(1,),
        max_abs_exponent=2,
        solve_order=6,
        supplemental_powers=(5, 7, 11),
    )

    assert tuple(label for label, _, _ in scan.ordered_basis_series) == (
        "GG",
        "GGneg",
        "GG2",
        "GG3",
        "GG4",
        "GG5",
        "GG7",
        "GG11",
    )
    assert tuple(label for label, _, _ in scan.quotient_basis_series) == (
        "Q_neg",
        "Q_2",
        "Q_3",
        "Q_4",
        "Q_5",
        "Q_7",
        "Q_11",
    )


def test_scan_gg_modular_equation_box_finds_chan_huang_exact_templates_on_true_gg():
    order = 24
    gg_template = get_benchmark("gollnitz_gordon_normalized").canonical_template.normalized()
    gg = continued_fraction_series_coeffs(gg_template, depth=24, order=order)

    scan = scan_gg_modular_equation_box(
        target_series=gg,
        benchmark_name="gollnitz_gordon_normalized",
        gg_series=gg,
        order=order,
        degree_values=(1,),
        max_abs_exponent=4,
        solve_order=10,
    )

    assert scan.exact_polynomial_template_labels == (
        "Chan--Huang Cor. 3.2(i) on (F, GG3)",
        "Chan--Huang Cor. 3.2(ii) on (F, GG4)",
    )
    assert scan.exact_polynomial_template_hits == scan.exact_polynomial_template_labels
    assert scan.quotient_exact_polynomial_template_labels == (
        "Chan--Huang Cor. 3.2(i) on (F, Q_3)",
        "Chan--Huang Cor. 3.2(ii) on (F, Q_4)",
    )
    assert scan.quotient_exact_polynomial_template_hits == scan.quotient_exact_polynomial_template_labels


def test_scan_parameterized_source_family_power_boxes_supports_family_specific_powers():
    order = 12
    gg = [sp.Integer(0) for _ in range(order)]
    gg[0] = 1
    gg[1] = 1

    scans = scan_parameterized_source_family_power_boxes(
        target_series=gg,
        ordered_base_families=(("GG", "gollnitz_gordon_normalized", gg),),
        powers=(2, 3),
        order=order,
        supplemental_powers_by_family={"GG": (5, 7, 11)},
    )

    assert len(scans) == 1
    labels = tuple(label for label, _ in scans[0].ordered_basis_series)
    assert labels == ("GG", "GG2", "GG3", "GG5", "GG7", "GG11")


def test_cli_identify_writes_power_tower_scan(tmp_path: Path):
    verified = tmp_path / "verified.jsonl"
    output_path = tmp_path / "identify.md"

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
                "identify",
                "--in",
                str(verified),
                "--candidate-id",
                "hero",
                "--benchmark-powers",
                "2,3,4",
                "--smoke",
                "--out",
                str(output_path),
            ]
        )
        == 0
    )

    text = output_path.read_text(encoding="utf-8")
    assert "Identification Note: `hero`" in text
    assert "## Build Timing" in text
    assert "- `rhs-uniqueness-search`:" in text
    assert "- `source-family-scans`:" in text
    assert "- `final-render`:" in text
    assert "Extra Multivariate Search" in text
    assert "## RHS Uniqueness Search" in text
    assert "Polynomial Functional Box" in text
    assert "Fractional-Linear Functional Box" in text
    assert "One-Source-Core Correction Objects" in text
    assert "Two-Source-Core Correction Objects" in text
    assert "Rational-Equivalence Reduced Object" in text
    assert ("F_red = B1 / R" in text) or ("Reduced-object bridge construction failed" in text)
    assert "Benchmark Power-Tower Prefix Scan" in text
    assert "Ratio-Object Source-Family Multiplicative Scan" in text
    assert "Ratio-Object Source-Family Fractional-Linear Scan" in text
    assert "Ratio-Object Source-Family Two-Layer Fractional-Linear Scan" in text
    assert "Ratio-Object Parameterized Source-Family Power Scan" in text
    assert "Ratio-Object Source-Family Eta-Correction Scan" in text
    assert "Ratio-Object Two-Core Source-Family Eta Scan" in text
    assert "Ratio-Object Quotient-Core Source-Family Eta Scan" in text
    assert "Ratio-Object Two-Quotient-Core Source-Family Eta Scan" in text
    assert "Ratio-Object Two-Quotient-Core Self-Quotient Finite-Product Scan" in text
    assert "Ratio-Object Two-Quotient-Core Self-Eta Functional Scan" in text
    assert "Ratio-Object Two-Quotient-Core Self-Fractional-Linear Scan" in text
    assert "Ratio-Object Two-Quotient-Core Self-Polynomial Functional Scan" in text
    assert "Ratio-Object Explicit GG/S Template Eta-Correction Scan" in text
    assert "Ratio-Object GG Modular-Equation Template Scan" in text
    assert "`GGneg = GG(-t)`" in text
    assert "Quotient basis: `Q_neg = GG(-t) / GG(t)`" in text
    assert "Quotient-coordinate polynomial" in text
    assert "Two-layer fractional-linear scan" in text
    assert "Quotient ladder" in text
    assert "Quotient-ladder polynomial" in text
    assert "Quotient-ladder multiplicative scan" in text
    assert "Quotient-ladder fractional-linear scan" in text
    assert "Quotient-ladder two-layer fractional-linear scan" in text
    assert "Mixed quotient basis" in text
    assert "Mixed-quotient polynomial" in text
    assert "Mixed-quotient multiplicative scan" in text
    assert "Mixed-quotient fractional-linear scan" in text
    assert "Mixed-quotient two-layer fractional-linear scan" in text
    assert "Polynomial `total degree <= 1`" in text
    assert "Ratio-Object Explicit GG/S Transform Template Scan" in text
    assert "### `GG` Explicit Transform Box" in text
    assert "### `GG` Family" in text
    assert "Ratio-Object RR-Tower Prefix Scan" in text
    assert "Ratio-Object Self-Quotient Finite-Product Scan" in text
    assert "Ratio-Object Eta-Quotient Scan" in text
    assert "Ratio-Object Multiplicative RR-Tower Scan" in text
    assert "Ratio-Object Fractional-Linear RR-Tower Scan" in text
    assert "Ratio-Object Two-Layer Fractional-Linear RR-Tower Scan" in text
    assert "`F = candidate / rogers_ramanujan_q3_normalized`" in text
    assert "`B2 = B1(t^2)`" in text
    assert "`B4 = B1(t^4)`" in text
    assert "total degree <= 1" in text
