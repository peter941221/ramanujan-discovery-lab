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
    """Multivariate polynomial relation P(X_0, ..., X_{k-1}) == 0 mod q^N."""

    order_checked: int
    variables: tuple[str, ...]
    max_total_degree: int
    coefficients: dict[tuple[int, ...], sp.Integer]  # exponent tuple -> coefficient

    def as_sympy(self, symbols: tuple[sp.Symbol, ...]) -> sp.Expr:
        if len(symbols) != len(self.variables):
            raise ValueError("symbol count must match variable count")
        expr = sp.Integer(0)
        for exponents, coeff in sorted(self.coefficients.items()):
            term = sp.Integer(1)
            for sym, exp in zip(symbols, exponents):
                if exp:
                    term *= sym**exp
            expr += coeff * term
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


def guess_polynomial_relation(
    *,
    series_by_variable: dict[str, Series],
    order: int,
    max_total_degree: int,
) -> PolynomialRelation | None:
    """Find integer coefficients for a small polynomial relation among series variables.

    The search space is monomials with total degree <= max_total_degree.
    """
    if order < 2:
        raise ValueError("order must be at least 2")
    if max_total_degree < 1:
        raise ValueError("max_total_degree must be at least 1")
    if len(series_by_variable) < 2:
        raise ValueError("need at least two variables for a relation search")
    if any(len(series) < order for series in series_by_variable.values()):
        raise ValueError("series are shorter than requested order")

    variables = tuple(series_by_variable.keys())
    series_list = [series_by_variable[name][:order] for name in variables]

    # Precompute powers up to max_total_degree for each variable.
    one = [sp.Integer(0) for _ in range(order)]
    one[0] = sp.Integer(1)
    powers: list[list[Series]] = []
    for series in series_list:
        var_pows: list[Series] = [one]
        for _ in range(max_total_degree):
            var_pows.append(series_mul(var_pows[-1], series))
        powers.append(var_pows)

    # Enumerate exponent tuples with total degree <= max_total_degree.
    exponent_tuples: list[tuple[int, ...]] = []
    columns: list[Series] = []

    def _recurse(idx: int, remaining: int, current: list[int]) -> None:
        if idx == len(variables):
            exponent_tuples.append(tuple(current))
            # Multiply the precomputed powers for this monomial.
            term = one
            for var_idx, exp in enumerate(current):
                if exp == 0:
                    continue
                term = series_mul(term, powers[var_idx][exp])
            columns.append(term)
            return
        for exp in range(remaining + 1):
            current.append(exp)
            _recurse(idx + 1, remaining - exp, current)
            current.pop()

    _recurse(0, max_total_degree, [])

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

    coeff_map: dict[tuple[int, ...], sp.Integer] = {}
    for exponents, coeff in zip(exponent_tuples, scaled):
        coeff_s = sp.simplify(coeff)
        if coeff_s == 0:
            continue
        coeff_map[exponents] = sp.Integer(int(coeff_s))

    if not coeff_map:
        return None

    return PolynomialRelation(
        order_checked=order,
        variables=variables,
        max_total_degree=max_total_degree,
        coefficients=coeff_map,
    )


def search_polynomial_relation(
    *,
    series_by_variable: dict[str, Series],
    order: int,
    max_total_degree: int,
) -> PolynomialRelation | None:
    if max_total_degree < 1:
        raise ValueError("max_total_degree must be at least 1")
    return guess_polynomial_relation(
        series_by_variable=series_by_variable,
        order=order,
        max_total_degree=max_total_degree,
    )


def _relation_residual_series(
    relation: PolynomialRelation,
    *,
    series_by_variable: dict[str, Series],
    order: int,
) -> Series:
    variables = relation.variables
    series_list = [series_by_variable[name][:order] for name in variables]

    one = [sp.Integer(0) for _ in range(order)]
    one[0] = sp.Integer(1)
    powers: list[list[Series]] = []
    for series in series_list:
        var_pows: list[Series] = [one]
        for _ in range(relation.max_total_degree):
            var_pows.append(series_mul(var_pows[-1], series))
        powers.append(var_pows)

    residual: Series = [sp.Integer(0) for _ in range(order)]
    for exponents, coeff in relation.coefficients.items():
        term = one
        for var_idx, exp in enumerate(exponents):
            if exp == 0:
                continue
            term = series_mul(term, powers[var_idx][exp])
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
    benchmark_powers: tuple[int, ...] = (),
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

    relation = search_polynomial_relation(
        series_by_variable={"C": candidate_recip, "B1": benchmark_recip},
        order=profile_order,
        max_total_degree=profile_degree,
    )

    extra_relation: PolynomialRelation | None = None
    benchmark_power_series: dict[int, Series] = {}
    if benchmark_powers:
        for power in sorted(set(benchmark_powers)):
            if power < 2:
                continue
            powered: Series = [sp.Integer(0) for _ in range(profile_order)]
            for idx, coeff in enumerate(benchmark_recip):
                j = power * idx
                if j >= profile_order:
                    break
                powered[j] = sp.simplify(coeff)
            benchmark_power_series[power] = powered

        variables: dict[str, Series] = {"C": candidate_recip, "B1": benchmark_recip}
        for power, series in benchmark_power_series.items():
            variables[f"B{power}"] = series

        extra_relation = search_polynomial_relation(
            series_by_variable=variables,
            order=profile_order,
            max_total_degree=min(profile_degree, 3 if smoke else profile_degree),
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
        f"- Polynomial relation search: total degree `<= {profile_degree}`",
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
        f"- `B1 = 1 / {record.closest_benchmark}`",
        "",
        "## Result",
        "",
    ]

    if relation is None:
        lines.extend(
            [
                "No nontrivial polynomial relation",
                "",
                "```text",
                "P(C, B1) = 0",
                "```",
                "",
                f"was found in the search box `total degree <= {profile_degree}` when checked modulo `{series_symbol}^{profile_order}`.",
            ]
        )
    else:
        sym_map = {name: sp.Symbol(name) for name in relation.variables}
        symbols = tuple(sym_map[name] for name in relation.variables)
        poly = relation.as_sympy(symbols)
        residual = _relation_residual_series(
            relation,
            series_by_variable={"C": candidate_recip, "B1": benchmark_recip},
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                "Found a candidate polynomial relation:",
                "",
                "```text",
                _format_expr(poly),
                "```",
                "",
                f"- Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
            ]
        )

    if benchmark_power_series:
        lines.extend(
            [
                "",
                "## Extra Multivariate Search",
                "",
                "We also tried a small multivariate search that includes benchmark power substitutions:",
                "",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.append("")

        if extra_relation is None:
            lines.extend(
                [
                    "No nontrivial multivariate polynomial relation was found",
                    "",
                    f"under `total degree <= {min(profile_degree, 3 if smoke else profile_degree)}` when checked modulo `{series_symbol}^{profile_order}`.",
                ]
            )
        else:
            sym_map = {name: sp.Symbol(name) for name in extra_relation.variables}
            symbols = tuple(sym_map[name] for name in extra_relation.variables)
            poly = extra_relation.as_sympy(symbols)
            residual = _relation_residual_series(
                extra_relation,
                series_by_variable={
                    "C": candidate_recip,
                    "B1": benchmark_recip,
                    **{f"B{p}": series for p, series in benchmark_power_series.items()},
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    "Found a candidate multivariate polynomial relation:",
                    "",
                    "```text",
                    _format_expr(poly),
                    "```",
                    "",
                    f"- Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                ]
            )

    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
