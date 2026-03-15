from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from pathlib import Path

import sympy as sp

from ramanujan_discovery.analysis import _format_expr
from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.research import reduce_template_by_step
from ramanujan_discovery.series import (
    Series,
    continued_fraction_series_coeffs,
    series_invert,
    series_mul,
)
from ramanujan_discovery.storage import read_candidates


@dataclass(frozen=True)
class PolynomialRelation:
    order_checked: int
    max_deg_x: int
    max_deg_y: int
    coefficients: dict[tuple[int, int], sp.Integer]  # (i,j) -> c_ij for X^i Y^j

    def as_sympy(self, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
        expr = sp.Integer(0)
        for (i, j), coeff in sorted(self.coefficients.items()):
            expr += coeff * x**i * y**j
        return sp.expand(expr)


def _lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def _series_active_exponents(template: QCFTemplate) -> list[int]:
    parts = [abs(template.numerator_q_shift), abs(template.numerator_q_step)]
    if template.numerator_extra_scale != 0:
        parts.extend([abs(template.numerator_extra_q_shift), abs(template.numerator_extra_q_step)])
    if template.denominator_scale != 0:
        parts.extend([abs(template.denominator_q_shift), abs(template.denominator_q_step)])
    return [value for value in parts if value != 0]


def guess_bivariate_polynomial_relation(
    *,
    x_series: Series,
    y_series: Series,
    max_deg_x: int,
    max_deg_y: int,
    order: int,
) -> PolynomialRelation | None:
    """Find integer coefficients c_ij such that sum c_ij X^i Y^j == 0 modulo q^order."""
    if order < 2:
        raise ValueError("order must be at least 2")
    if max_deg_x < 0 or max_deg_y < 0:
        raise ValueError("max degrees must be non-negative")
    if len(x_series) < order or len(y_series) < order:
        raise ValueError("series are shorter than requested order")

    x_series = x_series[:order]
    y_series = y_series[:order]

    x_pows = [None] * (max_deg_x + 1)
    x_pows[0] = [sp.Integer(0) for _ in range(order)]
    x_pows[0][0] = sp.Integer(1)
    for i in range(1, max_deg_x + 1):
        x_pows[i] = series_mul(x_pows[i - 1], x_series)  # type: ignore[arg-type]

    y_pows = [None] * (max_deg_y + 1)
    y_pows[0] = [sp.Integer(0) for _ in range(order)]
    y_pows[0][0] = sp.Integer(1)
    for j in range(1, max_deg_y + 1):
        y_pows[j] = series_mul(y_pows[j - 1], y_series)  # type: ignore[arg-type]

    monomials: list[tuple[int, int]] = []
    columns: list[Series] = []
    for i in range(max_deg_x + 1):
        for j in range(max_deg_y + 1):
            monomials.append((i, j))
            columns.append(series_mul(x_pows[i], y_pows[j]))  # type: ignore[arg-type]

    matrix = sp.Matrix([[columns[col][row] for col in range(len(columns))] for row in range(order)])
    nullspace = matrix.nullspace()
    if not nullspace:
        return None

    # Pick a small-ish basis vector and scale it to integer coefficients.
    basis_vec = min(nullspace, key=lambda v: sum(1 for item in v if item != 0))
    den_lcm = 1
    nums: list[sp.Integer] = []
    dens: list[int] = []
    for entry in basis_vec:
        num, den = sp.together(entry).as_numer_denom()
        if not (num.is_Integer and den.is_Integer):
            return None
        num_i = int(num)
        den_i = int(den)
        nums.append(sp.Integer(num_i))
        dens.append(abs(den_i))
        den_lcm = _lcm(den_lcm, abs(den_i))

    scaled = [sp.Integer(den_lcm) * num // sp.Integer(den) for num, den in zip(nums, dens)]
    if all(value == 0 for value in scaled):
        return None

    # Normalize by gcd and sign.
    int_scaled = [int(value) for value in scaled]
    overall_gcd = 0
    for value in int_scaled:
        overall_gcd = gcd(overall_gcd, abs(value))
    if overall_gcd > 1:
        scaled = [sp.Integer(int(value) // overall_gcd) for value in scaled]

    for value in scaled:
        if value != 0:
            if value < 0:
                scaled = [-v for v in scaled]
            break

    coeff_map: dict[tuple[int, int], sp.Integer] = {}
    for (i, j), coeff in zip(monomials, scaled):
        coeff_s = sp.simplify(coeff)
        if coeff_s != 0:
            coeff_map[(i, j)] = sp.Integer(int(coeff_s))

    if not coeff_map:
        return None

    return PolynomialRelation(
        order_checked=order,
        max_deg_x=max_deg_x,
        max_deg_y=max_deg_y,
        coefficients=coeff_map,
    )


def search_bivariate_relation(
    *,
    x_series: Series,
    y_series: Series,
    max_degree: int,
    order: int,
) -> PolynomialRelation | None:
    if max_degree < 1:
        raise ValueError("max_degree must be at least 1")
    for degree in range(1, max_degree + 1):
        relation = guess_bivariate_polynomial_relation(
            x_series=x_series,
            y_series=y_series,
            max_deg_x=degree,
            max_deg_y=degree,
            order=order,
        )
        if relation is not None:
            return relation
    return None


def _relation_residual_series(
    relation: PolynomialRelation,
    *,
    x_series: Series,
    y_series: Series,
    order: int,
) -> Series:
    x_pows = {0: [sp.Integer(0) for _ in range(order)]}
    x_pows[0][0] = sp.Integer(1)
    for i in range(1, relation.max_deg_x + 1):
        x_pows[i] = series_mul(x_pows[i - 1], x_series[:order])

    y_pows = {0: [sp.Integer(0) for _ in range(order)]}
    y_pows[0][0] = sp.Integer(1)
    for j in range(1, relation.max_deg_y + 1):
        y_pows[j] = series_mul(y_pows[j - 1], y_series[:order])

    residual: Series = [sp.Integer(0) for _ in range(order)]
    for (i, j), coeff in relation.coefficients.items():
        term = series_mul(x_pows[i], y_pows[j])
        for n in range(order):
            if term[n] == 0:
                continue
            residual[n] = sp.simplify(residual[n] + coeff * term[n])
    return residual


def build_candidate_identification_note(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    depth: int = 40,
    series_order: int = 90,
    max_degree: int = 4,
    smoke: bool = False,
) -> None:
    records = read_candidates(input_path)
    record: CandidateRecord | None = None
    for item in records:
        if item.id == candidate_id:
            record = item
            break
    if record is None:
        raise KeyError(f"unknown candidate id: {candidate_id}")

    benchmark = get_benchmark(record.closest_benchmark)

    profile_order = series_order
    profile_degree = max_degree
    profile_depth = depth
    if smoke:
        profile_order = min(profile_order, 60)
        profile_degree = min(profile_degree, 3)
        profile_depth = min(profile_depth, 28)

    active = _series_active_exponents(record.template) + _series_active_exponents(benchmark.canonical_template)
    step = 0
    for value in active:
        step = gcd(step, value)
    if step <= 0:
        step = 1
    reduced_candidate = record.template
    reduced_benchmark = benchmark.canonical_template
    variable_label = "q"
    series_symbol = "q"
    if step > 1:
        maybe_candidate = reduce_template_by_step(record.template, step)
        maybe_benchmark = reduce_template_by_step(benchmark.canonical_template, step)
        if maybe_candidate is not None and maybe_benchmark is not None:
            reduced_candidate = maybe_candidate
            reduced_benchmark = maybe_benchmark
            variable_label = f"t = q^{step}"
            series_symbol = "t"
        else:
            step = 1

    candidate_series = continued_fraction_series_coeffs(reduced_candidate, depth=profile_depth, order=profile_order)
    benchmark_series = continued_fraction_series_coeffs(reduced_benchmark, depth=profile_depth, order=profile_order)
    if candidate_series[0] == 0 or benchmark_series[0] == 0:
        raise ValueError("series constant term was zero; cannot build reciprocals for identification")

    candidate_recip = series_invert(candidate_series)
    benchmark_recip = series_invert(benchmark_series)

    relation = search_bivariate_relation(
        x_series=candidate_recip,
        y_series=benchmark_recip,
        max_degree=profile_degree,
        order=profile_order,
    )

    lines: list[str] = [
        f"# Identification Note: `{record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{record.id}`",
        f"- Closest benchmark: `{record.closest_benchmark}`",
        f"- Variable view: `{variable_label}`",
        f"- Depth: `{profile_depth}`",
        f"- Series order: `{profile_order}`",
        f"- Algebraic relation search: degrees `<= {profile_degree}` (both variables)",
        "",
        "## Objects",
        "",
        "- Candidate template:",
        f"  - `{reduced_candidate.signature()}`",
        "- Benchmark template:",
        f"  - `{reduced_benchmark.signature()}`",
        "",
        "We run the relation search on the **reciprocal** continued fractions (the `1 + ...` objects):",
        "",
        "- `C = 1 / candidate`",
        f"- `B = 1 / {record.closest_benchmark}`",
        "",
        "## Result",
        "",
    ]

    if relation is None:
        lines.extend(
            [
                "No nontrivial bivariate polynomial relation",
                "",
                f"```text\nP(C, B) = 0\n```",
                "",
                f"was found in the search box `deg_C, deg_B <= {profile_degree}` when checked modulo `{series_symbol}^{profile_order}`.",
            ]
        )
    else:
        X, Y = sp.symbols("C B")
        poly = relation.as_sympy(X, Y)
        residual = _relation_residual_series(
            relation,
            x_series=candidate_recip,
            y_series=benchmark_recip,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                "Found a candidate bivariate polynomial relation:",
                "",
                "```text",
                _format_expr(poly),
                "```",
                "",
                f"- Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
            ]
        )

    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
