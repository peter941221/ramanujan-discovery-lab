from pathlib import Path

import pytest
import sympy as sp

from ramanujan_discovery.cli import main
from ramanujan_discovery.identification import (
    benchmark_power_substitution_series,
    scan_named_fractional_linear_prefixes,
    scan_benchmark_power_relation_prefixes,
    scan_named_multiplicative_prefixes,
    scan_named_two_layer_fractional_linear_prefixes,
    scan_ratio_benchmark_fractional_linear_prefixes,
    scan_ratio_benchmark_multiplicative_prefixes,
    scan_ratio_benchmark_power_relation_prefixes,
    scan_ratio_benchmark_two_layer_fractional_linear_prefixes,
    search_fractional_linear_relation,
    search_multiplicative_relation,
    search_polynomial_relation,
    search_two_layer_fractional_linear_relation,
)
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.series import series_invert, series_mul, series_pow
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
    assert "Extra Multivariate Search" in text
    assert "Benchmark Power-Tower Prefix Scan" in text
    assert "Ratio-Object Source-Family Multiplicative Scan" in text
    assert "Ratio-Object Source-Family Fractional-Linear Scan" in text
    assert "Ratio-Object Source-Family Two-Layer Fractional-Linear Scan" in text
    assert "Ratio-Object RR-Tower Prefix Scan" in text
    assert "Ratio-Object Multiplicative RR-Tower Scan" in text
    assert "Ratio-Object Fractional-Linear RR-Tower Scan" in text
    assert "Ratio-Object Two-Layer Fractional-Linear RR-Tower Scan" in text
    assert "`F = candidate / rogers_ramanujan_q3_normalized`" in text
    assert "`B2 = B1(t^2)`" in text
    assert "`B4 = B1(t^4)`" in text
    assert "total degree <= 1" in text
