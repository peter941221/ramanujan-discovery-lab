from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb, gcd
from pathlib import Path

import sympy as sp

from ramanujan_discovery.analysis import _format_expr
from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.research import reduce_template_by_step
from ramanujan_discovery.series import (
    Series,
    continued_fraction_series_coeffs,
    series_div,
    series_invert,
    series_mul,
    series_pow,
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


@dataclass(frozen=True)
class BenchmarkPowerRelationScan:
    """Outcome for one prefix scan against benchmark power substitutions."""

    powers: tuple[int, ...]
    max_total_degree: int
    relation: PolynomialRelation | None
    error: str | None = None


@dataclass(frozen=True)
class FractionalLinearRelation:
    """A structured relation F = (1 + sum a_i U_i) / (1 + sum b_i U_i)."""

    order_checked: int
    basis_variables: tuple[str, ...]
    numerator_coefficients: dict[str, sp.Expr]
    denominator_coefficients: dict[str, sp.Expr]


@dataclass(frozen=True)
class FractionalLinearRelationScan:
    """Outcome for one prefix scan of fractional-linear templates."""

    powers: tuple[int, ...]
    relation: FractionalLinearRelation | None
    error: str | None = None


@dataclass(frozen=True)
class MultiplicativeRelation:
    """A structured relation F = prod_i B_i^e_i with small integer exponents."""

    order_checked: int
    basis_variables: tuple[str, ...]
    exponents: dict[str, int]


@dataclass(frozen=True)
class MultiplicativeRelationScan:
    """Outcome for one prefix scan of multiplicative benchmark-tower templates."""

    powers: tuple[int, ...]
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class TwoLayerFractionalLinearRelation:
    """A product of two low-complexity fractional-linear factors."""

    order_checked: int
    numerator_variables: tuple[str, str]
    denominator_variables: tuple[str, str]
    numerator_coefficients: tuple[sp.Expr, sp.Expr]
    denominator_coefficients: tuple[sp.Expr, sp.Expr]


@dataclass(frozen=True)
class TwoLayerFractionalLinearRelationScan:
    """Outcome for one prefix scan of two-layer fractional-linear templates."""

    powers: tuple[int, ...]
    relations: tuple[TwoLayerFractionalLinearRelation, ...]
    total_hits: int
    tuples_checked: int
    error: str | None = None


def _lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def _monomial_count(num_variables: int, max_total_degree: int) -> int:
    if num_variables < 1:
        raise ValueError("num_variables must be at least 1")
    if max_total_degree < 0:
        raise ValueError("max_total_degree must be non-negative")
    return comb(num_variables + max_total_degree, num_variables)


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
    required_variable: str | None = None,
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

    if required_variable is not None and required_variable not in series_by_variable:
        raise ValueError("required_variable must be one of the series variable names")

    num_monomials = _monomial_count(len(series_by_variable), max_total_degree)
    if num_monomials > order:
        raise ValueError(
            "underdetermined polynomial relation search: "
            f"{num_monomials} monomials > {order} constraints "
            "(increase order, lower max_total_degree, or reduce variables)"
        )

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
    candidate_vecs = nullspace
    if required_variable is not None:
        required_index = variables.index(required_variable)
        required_columns = [idx for idx, exps in enumerate(exponent_tuples) if exps[required_index] > 0]
        candidate_vecs = [
            vec
            for vec in nullspace
            if any(vec[col] != 0 for col in required_columns)
        ]
        if not candidate_vecs:
            return None

    basis_vec = min(candidate_vecs, key=lambda v: sum(1 for item in v if item != 0))
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
    required_variable: str | None = None,
) -> PolynomialRelation | None:
    if max_total_degree < 1:
        raise ValueError("max_total_degree must be at least 1")
    return guess_polynomial_relation(
        series_by_variable=series_by_variable,
        order=order,
        max_total_degree=max_total_degree,
        required_variable=required_variable,
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


def benchmark_power_substitution_series(base_series: Series, *, power: int, order: int) -> Series:
    if power < 2:
        raise ValueError("power must be at least 2")
    if len(base_series) < order:
        raise ValueError("base_series is shorter than requested order")

    powered: Series = [sp.Integer(0) for _ in range(order)]
    for idx, coeff in enumerate(base_series[:order]):
        mapped_idx = power * idx
        if mapped_idx >= order:
            break
        powered[mapped_idx] = sp.simplify(coeff)
    return powered


def _series_subtract_one(series: Series) -> Series:
    shifted = [sp.simplify(value) for value in series]
    shifted[0] = sp.simplify(shifted[0] - 1)
    return shifted


def _series_log_coeffs(series: Series) -> Series:
    """Return coefficients of log(F) for a series F with constant term 1."""
    if not series:
        raise ValueError("series must be non-empty")
    if sp.simplify(series[0] - 1) != 0:
        raise ValueError("series constant term must be 1 for a logarithm expansion")

    coeffs: Series = [sp.Integer(0) for _ in range(len(series))]
    for n in range(1, len(series)):
        rhs = sp.Integer(0)
        for j in range(1, n):
            rhs += sp.Integer(j) * coeffs[j] * series[n - j]
        coeffs[n] = sp.simplify(series[n] - rhs / sp.Integer(n))
    return coeffs


def _signed_series_pow(base: Series, exponent: int) -> Series:
    if exponent >= 0:
        return series_pow(base, exponent)
    return series_invert(series_pow(base, -exponent))


def search_multiplicative_relation(
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = prod_i B_i^e_i with bounded integer exponents."""
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if not basis_series_by_variable:
        raise ValueError("need at least one basis series for a multiplicative search")
    if any(len(series) < order for series in basis_series_by_variable.values()):
        raise ValueError("series are shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a multiplicative search")
    if any(sp.simplify(series[0] - 1) != 0 for series in basis_series_by_variable.values()):
        raise ValueError("basis series must have constant term 1 for a multiplicative search")
    if max_abs_exponent < 1:
        raise ValueError("max_abs_exponent must be at least 1")

    basis_variables = tuple(basis_series_by_variable.keys())
    num_unknowns = len(basis_variables)
    num_constraints = order - 1
    if num_unknowns > num_constraints:
        raise ValueError(
            "underdetermined multiplicative relation search: "
            f"{num_unknowns} exponents > {num_constraints} constraints "
            "(increase order or reduce variables)"
        )

    target_log = _series_log_coeffs(target_series[:order])
    basis_logs = {
        name: _series_log_coeffs(series[:order])
        for name, series in basis_series_by_variable.items()
    }

    matrix = sp.Matrix(
        [
            [basis_logs[name][n] for name in basis_variables]
            for n in range(1, order)
        ]
    )
    rhs_vector = sp.Matrix([target_log[n] for n in range(1, order)])

    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs_vector).rank()
    if augmented_rank > rank:
        return None
    if rank < num_unknowns:
        raise ValueError(
            "underdetermined multiplicative relation search: "
            f"rank {rank} < {num_unknowns} exponents "
            "(increase order or reduce variables)"
        )

    unknowns = sp.symbols(f"e0:{len(basis_variables)}")
    solution_set = sp.linsolve((matrix, rhs_vector), unknowns)
    if not solution_set:
        return None
    solution = next(iter(solution_set))
    if any(value.free_symbols for value in solution):
        raise ValueError("multiplicative relation search returned a parametric solution")

    exponent_map: dict[str, int] = {}
    for name, value in zip(basis_variables, solution):
        value_simplified = sp.simplify(value)
        if not value_simplified.is_integer:
            return None
        exponent = int(value_simplified)
        if abs(exponent) > max_abs_exponent:
            return None
        if exponent != 0:
            exponent_map[name] = exponent

    if not exponent_map:
        return None

    relation = MultiplicativeRelation(
        order_checked=order,
        basis_variables=basis_variables,
        exponents=exponent_map,
    )
    residual = _multiplicative_relation_residual_series(
        relation,
        target_series=target_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return relation


def _format_multiplicative_relation(
    relation: MultiplicativeRelation,
    *,
    target_variable: str,
) -> str:
    terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        terms.append(f"{name}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable} = {rhs}"


def _multiplicative_relation_residual_series(
    relation: MultiplicativeRelation,
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    order: int,
) -> Series:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if any(name not in basis_series_by_variable for name in relation.basis_variables):
        raise ValueError("basis series are missing variables from the relation")
    if any(len(basis_series_by_variable[name]) < order for name in relation.basis_variables):
        raise ValueError("basis series are shorter than requested order")

    product_series: Series = [sp.Integer(0) for _ in range(order)]
    product_series[0] = sp.Integer(1)
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name, 0)
        if exponent == 0:
            continue
        factor = _signed_series_pow(basis_series_by_variable[name][:order], exponent)
        product_series = series_mul(product_series, factor)
    return [sp.simplify(target_series[n] - product_series[n]) for n in range(order)]


def search_fractional_linear_relation(
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    order: int,
) -> FractionalLinearRelation | None:
    """Search F = (1 + sum a_i U_i) / (1 + sum b_i U_i) with U_i = B_i - 1."""
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if not basis_series_by_variable:
        raise ValueError("need at least one basis series for a fractional-linear search")
    if any(len(series) < order for series in basis_series_by_variable.values()):
        raise ValueError("series are shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a fractional-linear search")

    basis_variables = tuple(basis_series_by_variable.keys())
    num_unknowns = 2 * len(basis_variables)
    num_constraints = order - 1
    if num_unknowns > num_constraints:
        raise ValueError(
            "underdetermined fractional-linear relation search: "
            f"{num_unknowns} coefficients > {num_constraints} constraints "
            "(increase order or reduce variables)"
        )

    shifted_basis = {
        name: _series_subtract_one(series[:order])
        for name, series in basis_series_by_variable.items()
    }
    rhs = [sp.simplify(-value) for value in _series_subtract_one(target_series[:order])[1:]]

    columns: list[list[sp.Expr]] = []
    for name in basis_variables:
        basis = shifted_basis[name]
        columns.append([sp.simplify(-basis[n]) for n in range(1, order)])
    for name in basis_variables:
        basis = shifted_basis[name]
        product = series_mul(target_series[:order], basis)
        columns.append([sp.simplify(product[n]) for n in range(1, order)])

    matrix = sp.Matrix(
        [
            [column[row_idx] for column in columns]
            for row_idx in range(num_constraints)
        ]
    )
    rhs_vector = sp.Matrix(rhs)

    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs_vector).rank()
    if augmented_rank > rank:
        return None
    if rank < num_unknowns:
        raise ValueError(
            "underdetermined fractional-linear relation search: "
            f"rank {rank} < {num_unknowns} coefficients "
            "(increase order or reduce variables)"
        )

    unknowns = sp.symbols(f"a0:{len(basis_variables)} b0:{len(basis_variables)}")
    solution_set = sp.linsolve((matrix, rhs_vector), unknowns)
    if not solution_set:
        return None
    solution = next(iter(solution_set))
    if any(value.free_symbols for value in solution):
        raise ValueError("fractional-linear relation search returned a parametric solution")

    numerator_coefficients: dict[str, sp.Expr] = {}
    denominator_coefficients: dict[str, sp.Expr] = {}
    for idx, name in enumerate(basis_variables):
        value = sp.simplify(solution[idx])
        if value != 0:
            numerator_coefficients[name] = value
    for idx, name in enumerate(basis_variables, start=len(basis_variables)):
        value = sp.simplify(solution[idx])
        if value != 0:
            denominator_coefficients[name] = value

    if not numerator_coefficients and not denominator_coefficients:
        return None

    return FractionalLinearRelation(
        order_checked=order,
        basis_variables=basis_variables,
        numerator_coefficients=numerator_coefficients,
        denominator_coefficients=denominator_coefficients,
    )


def _format_fractional_linear_relation(relation: FractionalLinearRelation, *, target_variable: str) -> str:
    def _side(coefficients: dict[str, sp.Expr]) -> str:
        terms = ["1"]
        for name in relation.basis_variables:
            coeff = coefficients.get(name)
            if coeff is None:
                continue
            coeff_str = _format_expr(coeff)
            terms.append(f"{coeff_str}*({name} - 1)")
        return " + ".join(terms)

    return f"{target_variable} = ({_side(relation.numerator_coefficients)}) / ({_side(relation.denominator_coefficients)})"


def _fractional_linear_relation_residual_series(
    relation: FractionalLinearRelation,
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    order: int,
) -> Series:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if any(name not in basis_series_by_variable for name in relation.basis_variables):
        raise ValueError("basis series are missing variables from the relation")
    if any(len(basis_series_by_variable[name]) < order for name in relation.basis_variables):
        raise ValueError("basis series are shorter than requested order")

    numerator: Series = [sp.Integer(0) for _ in range(order)]
    denominator: Series = [sp.Integer(0) for _ in range(order)]
    numerator[0] = sp.Integer(1)
    denominator[0] = sp.Integer(1)

    for name in relation.basis_variables:
        shifted = _series_subtract_one(basis_series_by_variable[name][:order])
        numer_coeff = relation.numerator_coefficients.get(name)
        if numer_coeff is not None:
            for idx, value in enumerate(shifted):
                if value == 0:
                    continue
                numerator[idx] = sp.simplify(numerator[idx] + numer_coeff * value)
        denom_coeff = relation.denominator_coefficients.get(name)
        if denom_coeff is not None:
            for idx, value in enumerate(shifted):
                if value == 0:
                    continue
                denominator[idx] = sp.simplify(denominator[idx] + denom_coeff * value)

    lhs = series_mul(target_series[:order], denominator)
    return [sp.simplify(lhs[idx] - numerator[idx]) for idx in range(order)]


def search_two_layer_fractional_linear_relation(
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    numerator_variables: tuple[str, str],
    denominator_variables: tuple[str, str],
    order: int,
    solve_order: int | None = None,
) -> TwoLayerFractionalLinearRelation | None:
    """Search F = prod_j (1 + a_j U_xj) / (1 + b_j U_yj) with U = B - 1."""
    if order < 7:
        raise ValueError("order must be at least 7 for a two-layer fractional-linear search")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if not basis_series_by_variable:
        raise ValueError("need at least one basis series for a two-layer fractional-linear search")
    if any(len(series) < order for series in basis_series_by_variable.values()):
        raise ValueError("series are shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a two-layer fractional-linear search")

    basis_variables = tuple(basis_series_by_variable.keys())
    if any(name not in basis_variables for name in numerator_variables + denominator_variables):
        raise ValueError("all numerator and denominator variables must be present in the basis")

    solve_bound = solve_order if solve_order is not None else min(order, 14)
    if solve_bound > order:
        raise ValueError("solve_order cannot exceed order")
    if solve_bound < 7:
        raise ValueError("solve_order must be at least 7")

    solve_target = target_series[:solve_bound]
    shifted_basis = {
        name: _series_subtract_one(series[:solve_bound])
        for name, series in basis_series_by_variable.items()
    }
    rhs = [sp.simplify(-value) for value in _series_subtract_one(solve_target)[1:]]

    ux1 = shifted_basis[numerator_variables[0]]
    ux2 = shifted_basis[numerator_variables[1]]
    uy1 = shifted_basis[denominator_variables[0]]
    uy2 = shifted_basis[denominator_variables[1]]

    product_num = series_mul(ux1, ux2)
    product_den = series_mul(uy1, uy2)
    target_times_den_1 = series_mul(solve_target, uy1)
    target_times_den_2 = series_mul(solve_target, uy2)
    target_times_product_den = series_mul(solve_target, product_den)

    columns = [
        [sp.simplify(-ux1[n]) for n in range(1, solve_bound)],
        [sp.simplify(target_times_den_1[n]) for n in range(1, solve_bound)],
        [sp.simplify(-ux2[n]) for n in range(1, solve_bound)],
        [sp.simplify(target_times_den_2[n]) for n in range(1, solve_bound)],
        [sp.simplify(-product_num[n]) for n in range(1, solve_bound)],
        [sp.simplify(target_times_product_den[n]) for n in range(1, solve_bound)],
    ]

    matrix = sp.Matrix(
        [
            [column[row_idx] for column in columns]
            for row_idx in range(solve_bound - 1)
        ]
    )
    rhs_vector = sp.Matrix(rhs)

    num_unknowns = 6
    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs_vector).rank()
    if augmented_rank > rank:
        return None
    if rank < num_unknowns:
        raise ValueError(
            "underdetermined two-layer fractional-linear relation search: "
            f"rank {rank} < {num_unknowns} coefficients "
            "(increase solve_order or reduce the template family)"
        )

    unknowns = sp.symbols("a0 b0 a1 b1 p q")
    solution_set = sp.linsolve((matrix, rhs_vector), unknowns)
    if not solution_set:
        return None
    solution = next(iter(solution_set))
    if any(value.free_symbols for value in solution):
        raise ValueError("two-layer fractional-linear relation search returned a parametric solution")

    a0, b0, a1, b1, product_num_coeff, product_den_coeff = map(sp.simplify, solution)
    if sp.simplify(product_num_coeff - a0 * a1) != 0:
        return None
    if sp.simplify(product_den_coeff - b0 * b1) != 0:
        return None
    if all(value == 0 for value in (a0, b0, a1, b1)):
        return None

    relation = TwoLayerFractionalLinearRelation(
        order_checked=order,
        numerator_variables=numerator_variables,
        denominator_variables=denominator_variables,
        numerator_coefficients=(a0, a1),
        denominator_coefficients=(b0, b1),
    )

    residual = _two_layer_fractional_linear_relation_residual_series(
        relation,
        target_series=target_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return relation


def _format_two_layer_fractional_linear_relation(
    relation: TwoLayerFractionalLinearRelation,
    *,
    target_variable: str,
) -> str:
    factors: list[str] = []
    for idx in range(2):
        num_var = relation.numerator_variables[idx]
        den_var = relation.denominator_variables[idx]
        num_coeff = _format_expr(relation.numerator_coefficients[idx])
        den_coeff = _format_expr(relation.denominator_coefficients[idx])
        factors.append(
            f"((1 + {num_coeff}*({num_var} - 1)) / (1 + {den_coeff}*({den_var} - 1)))"
        )
    return f"{target_variable} = {' * '.join(factors)}"


def _two_layer_fractional_linear_relation_residual_series(
    relation: TwoLayerFractionalLinearRelation,
    *,
    target_series: Series,
    basis_series_by_variable: dict[str, Series],
    order: int,
) -> Series:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if any(name not in basis_series_by_variable for name in relation.numerator_variables + relation.denominator_variables):
        raise ValueError("basis series are missing variables from the relation")

    numerator: Series = [sp.Integer(0) for _ in range(order)]
    denominator: Series = [sp.Integer(0) for _ in range(order)]
    numerator[0] = sp.Integer(1)
    denominator[0] = sp.Integer(1)

    for idx in range(2):
        num_factor: Series = [sp.Integer(0) for _ in range(order)]
        den_factor: Series = [sp.Integer(0) for _ in range(order)]
        num_factor[0] = sp.Integer(1)
        den_factor[0] = sp.Integer(1)

        num_shifted = _series_subtract_one(basis_series_by_variable[relation.numerator_variables[idx]][:order])
        den_shifted = _series_subtract_one(basis_series_by_variable[relation.denominator_variables[idx]][:order])

        for n, value in enumerate(num_shifted):
            if value == 0:
                continue
            num_factor[n] = sp.simplify(num_factor[n] + relation.numerator_coefficients[idx] * value)
        for n, value in enumerate(den_shifted):
            if value == 0:
                continue
            den_factor[n] = sp.simplify(den_factor[n] + relation.denominator_coefficients[idx] * value)

        numerator = series_mul(numerator, num_factor)
        denominator = series_mul(denominator, den_factor)

    lhs = series_mul(target_series[:order], denominator)
    return [sp.simplify(lhs[idx] - numerator[idx]) for idx in range(order)]


def scan_benchmark_power_relation_prefixes(
    *,
    candidate_recip: Series,
    benchmark_recip: Series,
    powers: tuple[int, ...],
    order: int,
    degree_values: tuple[int, ...],
    required_variable: str | None = "C",
) -> list[BenchmarkPowerRelationScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not unique_powers:
        return []

    benchmark_power_series = {
        power: benchmark_power_substitution_series(benchmark_recip, power=power, order=order)
        for power in unique_powers
    }

    scans: list[BenchmarkPowerRelationScan] = []
    prefix: list[int] = []
    for power in unique_powers:
        prefix.append(power)
        variables = {"C": candidate_recip, "B1": benchmark_recip}
        for prefix_power in prefix:
            variables[f"B{prefix_power}"] = benchmark_power_series[prefix_power]
        for degree in tuple(sorted({value for value in degree_values if value >= 1})):
            try:
                relation = search_polynomial_relation(
                    series_by_variable=variables,
                    order=order,
                    max_total_degree=degree,
                    required_variable=required_variable,
                )
                scans.append(
                    BenchmarkPowerRelationScan(
                        powers=tuple(prefix),
                        max_total_degree=degree,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                scans.append(
                    BenchmarkPowerRelationScan(
                        powers=tuple(prefix),
                        max_total_degree=degree,
                        relation=None,
                        error=str(exc),
                    )
                )
    return scans


def scan_ratio_benchmark_power_relation_prefixes(
    *,
    ratio_series: Series,
    benchmark_series: Series,
    powers: tuple[int, ...],
    order: int,
    degree_values: tuple[int, ...],
    required_variable: str | None = "F",
) -> list[BenchmarkPowerRelationScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not unique_powers:
        return []

    benchmark_power_series = {
        power: benchmark_power_substitution_series(benchmark_series, power=power, order=order)
        for power in unique_powers
    }

    scans: list[BenchmarkPowerRelationScan] = []
    prefix: list[int] = []
    for power in unique_powers:
        prefix.append(power)
        variables = {"F": ratio_series, "B1": benchmark_series}
        for prefix_power in prefix:
            variables[f"B{prefix_power}"] = benchmark_power_series[prefix_power]
        for degree in tuple(sorted({value for value in degree_values if value >= 1})):
            try:
                relation = search_polynomial_relation(
                    series_by_variable=variables,
                    order=order,
                    max_total_degree=degree,
                    required_variable=required_variable,
                )
                scans.append(
                    BenchmarkPowerRelationScan(
                        powers=tuple(prefix),
                        max_total_degree=degree,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                scans.append(
                    BenchmarkPowerRelationScan(
                        powers=tuple(prefix),
                        max_total_degree=degree,
                        relation=None,
                        error=str(exc),
                    )
                )
    return scans


def scan_ratio_benchmark_fractional_linear_prefixes(
    *,
    ratio_series: Series,
    benchmark_series: Series,
    powers: tuple[int, ...],
    order: int,
) -> list[FractionalLinearRelationScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not unique_powers:
        return []

    benchmark_power_series = {
        power: benchmark_power_substitution_series(benchmark_series, power=power, order=order)
        for power in unique_powers
    }

    scans: list[FractionalLinearRelationScan] = []
    prefix: list[int] = []
    for power in unique_powers:
        prefix.append(power)
        variables = {"B1": benchmark_series}
        for prefix_power in prefix:
            variables[f"B{prefix_power}"] = benchmark_power_series[prefix_power]
        try:
            relation = search_fractional_linear_relation(
                target_series=ratio_series,
                basis_series_by_variable=variables,
                order=order,
            )
            scans.append(
                FractionalLinearRelationScan(
                    powers=tuple(prefix),
                    relation=relation,
                )
            )
        except ValueError as exc:
            scans.append(
                FractionalLinearRelationScan(
                    powers=tuple(prefix),
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def scan_ratio_benchmark_multiplicative_prefixes(
    *,
    ratio_series: Series,
    benchmark_series: Series,
    powers: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[MultiplicativeRelationScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not unique_powers:
        return []

    benchmark_power_series = {
        power: benchmark_power_substitution_series(benchmark_series, power=power, order=order)
        for power in unique_powers
    }

    scans: list[MultiplicativeRelationScan] = []
    prefix: list[int] = []
    for power in unique_powers:
        prefix.append(power)
        variables = {"B1": benchmark_series}
        for prefix_power in prefix:
            variables[f"B{prefix_power}"] = benchmark_power_series[prefix_power]
        try:
            relation = search_multiplicative_relation(
                target_series=ratio_series,
                basis_series_by_variable=variables,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(
                MultiplicativeRelationScan(
                    powers=tuple(prefix),
                    relation=relation,
                )
            )
        except ValueError as exc:
            scans.append(
                MultiplicativeRelationScan(
                    powers=tuple(prefix),
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def scan_ratio_benchmark_two_layer_fractional_linear_prefixes(
    *,
    ratio_series: Series,
    benchmark_series: Series,
    powers: tuple[int, ...],
    order: int,
    solve_order: int | None = None,
    max_reported_hits: int = 3,
) -> list[TwoLayerFractionalLinearRelationScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not unique_powers:
        return []

    if max_reported_hits < 1:
        raise ValueError("max_reported_hits must be at least 1")

    benchmark_power_series = {
        power: benchmark_power_substitution_series(benchmark_series, power=power, order=order)
        for power in unique_powers
    }

    scans: list[TwoLayerFractionalLinearRelationScan] = []
    prefix: list[int] = []
    for power in unique_powers:
        prefix.append(power)
        variables = {"B1": benchmark_series}
        for prefix_power in prefix:
            variables[f"B{prefix_power}"] = benchmark_power_series[prefix_power]

        basis_names = tuple(variables.keys())
        hits: list[TwoLayerFractionalLinearRelation] = []
        seen_signatures: set[str] = set()
        total_hits = 0
        tuples_checked = 0
        try:
            for numerator_variables in product(basis_names, repeat=2):
                for denominator_variables in product(basis_names, repeat=2):
                    factor_1 = (numerator_variables[0], denominator_variables[0])
                    factor_2 = (numerator_variables[1], denominator_variables[1])
                    if factor_1 > factor_2:
                        continue
                    tuples_checked += 1
                    try:
                        relation = search_two_layer_fractional_linear_relation(
                            target_series=ratio_series,
                            basis_series_by_variable=variables,
                            numerator_variables=numerator_variables,
                            denominator_variables=denominator_variables,
                            order=order,
                            solve_order=solve_order,
                        )
                    except ValueError:
                        continue
                    if relation is None:
                        continue
                    signature = _format_two_layer_fractional_linear_relation(
                        relation,
                        target_variable="F",
                    )
                    if signature in seen_signatures:
                        continue
                    seen_signatures.add(signature)
                    total_hits += 1
                    if len(hits) < max_reported_hits:
                        hits.append(relation)

            scans.append(
                TwoLayerFractionalLinearRelationScan(
                    powers=tuple(prefix),
                    relations=tuple(hits),
                    total_hits=total_hits,
                    tuples_checked=tuples_checked,
                )
            )
        except ValueError as exc:
            scans.append(
                TwoLayerFractionalLinearRelationScan(
                    powers=tuple(prefix),
                    relations=(),
                    total_hits=0,
                    tuples_checked=tuples_checked,
                    error=str(exc),
                )
            )
    return scans


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

    ratio_series = series_div(candidate_series, benchmark_series)
    candidate_recip = series_invert(candidate_series)
    benchmark_recip = series_invert(benchmark_series)

    relation: PolynomialRelation | None = None
    relation_error: str | None = None
    try:
        relation = search_polynomial_relation(
            series_by_variable={"C": candidate_recip, "B1": benchmark_recip},
            order=profile_order,
            max_total_degree=profile_degree,
            required_variable="C",
        )
    except ValueError as exc:
        relation_error = str(exc)

    extra_relation: PolynomialRelation | None = None
    extra_relation_error: str | None = None
    benchmark_power_series: dict[int, Series] = {}
    extra_search_degree = min(profile_degree, 3 if smoke else profile_degree)
    power_tower_scans: list[BenchmarkPowerRelationScan] = []
    ratio_power_tower_scans: list[BenchmarkPowerRelationScan] = []
    ratio_multiplicative_scans: list[MultiplicativeRelationScan] = []
    ratio_fractional_linear_scans: list[FractionalLinearRelationScan] = []
    ratio_two_layer_fractional_linear_scans: list[TwoLayerFractionalLinearRelationScan] = []
    if benchmark_powers:
        for power in sorted(set(benchmark_powers)):
            if power < 2:
                continue
            benchmark_power_series[power] = benchmark_power_substitution_series(
                benchmark_recip,
                power=power,
                order=profile_order,
            )

        variables: dict[str, Series] = {"C": candidate_recip, "B1": benchmark_recip}
        for power, series in benchmark_power_series.items():
            variables[f"B{power}"] = series

        try:
            extra_relation = search_polynomial_relation(
                series_by_variable=variables,
                order=profile_order,
                max_total_degree=extra_search_degree,
                required_variable="C",
            )
        except ValueError as exc:
            extra_relation_error = str(exc)

        power_tower_scans = scan_benchmark_power_relation_prefixes(
            candidate_recip=candidate_recip,
            benchmark_recip=benchmark_recip,
            powers=tuple(benchmark_power_series),
            order=profile_order,
            degree_values=tuple(value for value in (1, min(profile_degree, 2)) if value >= 1),
            required_variable="C",
        )
        ratio_power_tower_scans = scan_ratio_benchmark_power_relation_prefixes(
            ratio_series=ratio_series,
            benchmark_series=benchmark_series,
            powers=tuple(benchmark_power_series),
            order=profile_order,
            degree_values=tuple(value for value in (1, min(profile_degree, 2)) if value >= 1),
            required_variable="F",
        )
        ratio_multiplicative_scans = scan_ratio_benchmark_multiplicative_prefixes(
            ratio_series=ratio_series,
            benchmark_series=benchmark_series,
            powers=tuple(benchmark_power_series),
            order=profile_order,
            max_abs_exponent=6 if smoke else 8,
        )
        ratio_fractional_linear_scans = scan_ratio_benchmark_fractional_linear_prefixes(
            ratio_series=ratio_series,
            benchmark_series=benchmark_series,
            powers=tuple(benchmark_power_series),
            order=profile_order,
        )
        ratio_two_layer_fractional_linear_scans = scan_ratio_benchmark_two_layer_fractional_linear_prefixes(
            ratio_series=ratio_series,
            benchmark_series=benchmark_series,
            powers=tuple(benchmark_power_series),
            order=profile_order,
            solve_order=min(profile_order, 14 if smoke else 18),
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

    if relation_error is not None:
        lines.extend(
            [
                "Skipped polynomial relation search:",
                "",
                "```text",
                relation_error,
                "```",
            ]
        )
    elif relation is None:
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

        if extra_relation_error is not None:
            lines.extend(
                [
                    "Skipped multivariate relation search:",
                    "",
                    "```text",
                    extra_relation_error,
                    "```",
                ]
            )
        elif extra_relation is None:
            lines.extend(
                [
                    "No candidate-dependent multivariate polynomial relation was found",
                    "",
                    f"under `total degree <= {extra_search_degree}` when checked modulo `{series_symbol}^{profile_order}`.",
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

    if power_tower_scans:
        lines.extend(
            [
                "",
                "## Benchmark Power-Tower Prefix Scan",
                "",
                "We also ran a structured low-degree scan against prefixes of the benchmark power tower:",
                "",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(C, B1, B2)`, then `(C, B1, B2, B3)`, and so on through the final listed power.",
                "- Degrees scanned here are intentionally low (`1` and `2`) to stay in candidate-dependent, well-determined boxes.",
                "",
            ]
        )

        grouped_scans: dict[int, list[BenchmarkPowerRelationScan]] = {}
        for scan in power_tower_scans:
            grouped_scans.setdefault(scan.max_total_degree, []).append(scan)

        any_hit = any(scan.relation is not None for scan in power_tower_scans)
        if not any_hit:
            lines.append("No candidate-dependent relation was found in any scanned prefix box.")
            lines.append("")

        for degree in sorted(grouped_scans):
            scans = grouped_scans[degree]
            prefix_labels = [f"`B{scan.powers[-1]}`" for scan in scans if scan.error is None]
            if not any_hit and prefix_labels:
                lines.append(
                    f"- `total degree <= {degree}`: no hit for prefixes ending at {', '.join(prefix_labels)}."
                )
            for scan in scans:
                if scan.error is not None:
                    lines.append(
                        f"- `total degree <= {degree}` prefix ending at `B{scan.powers[-1]}` skipped: {scan.error}"
                    )
                elif scan.relation is not None:
                    sym_map = {name: sp.Symbol(name) for name in scan.relation.variables}
                    poly = scan.relation.as_sympy(tuple(sym_map[name] for name in scan.relation.variables))
                    lines.extend(
                        [
                            f"- `total degree <= {degree}` prefix ending at `B{scan.powers[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_expr(poly),
                            "```",
                            "",
                        ]
                    )

    if ratio_power_tower_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object RR-Tower Prefix Scan",
                "",
                "We also scanned the multiplicative correction object against prefixes of the benchmark tower:",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                f"- `B1 = {record.closest_benchmark}`",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(F, B1, B2)`, then `(F, B1, B2, B3)`, and so on through the final listed power.",
                "- Degrees scanned here are intentionally low (`1` and `2`) to stay in candidate-dependent, well-determined boxes.",
                "",
            ]
        )

        grouped_ratio_scans: dict[int, list[BenchmarkPowerRelationScan]] = {}
        for scan in ratio_power_tower_scans:
            grouped_ratio_scans.setdefault(scan.max_total_degree, []).append(scan)

        any_ratio_hit = any(scan.relation is not None for scan in ratio_power_tower_scans)
        if not any_ratio_hit:
            lines.append("No candidate-dependent relation was found for the ratio object in any scanned prefix box.")
            lines.append("")

        for degree in sorted(grouped_ratio_scans):
            scans = grouped_ratio_scans[degree]
            prefix_labels = [f"`B{scan.powers[-1]}`" for scan in scans if scan.error is None]
            if not any_ratio_hit and prefix_labels:
                lines.append(
                    f"- `total degree <= {degree}`: no hit for ratio-object prefixes ending at {', '.join(prefix_labels)}."
                )
            for scan in scans:
                if scan.error is not None:
                    lines.append(
                        f"- `total degree <= {degree}` ratio-object prefix ending at `B{scan.powers[-1]}` skipped: {scan.error}"
                    )
                elif scan.relation is not None:
                    sym_map = {name: sp.Symbol(name) for name in scan.relation.variables}
                    poly = scan.relation.as_sympy(tuple(sym_map[name] for name in scan.relation.variables))
                    residual = _relation_residual_series(
                        scan.relation,
                        series_by_variable={
                            "F": ratio_series,
                            "B1": benchmark_series,
                            **{f"B{p}": benchmark_power_substitution_series(benchmark_series, power=p, order=profile_order) for p in scan.powers},
                        },
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- `total degree <= {degree}` ratio-object prefix ending at `B{scan.powers[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_expr(poly),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

    if ratio_fractional_linear_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Multiplicative RR-Tower Scan",
                "",
                "We also searched for exact multiplicative corrections built from the benchmark tower:",
                "",
                "```text",
                "F = prod_i B_i^e_i",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                f"- `B1 = {record.closest_benchmark}`",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.",
                "- Exponents are solved exactly from the log-series constraints, then verified by exact series re-expansion.",
                "",
            ]
        )

        any_multiplicative_hit = any(scan.relation is not None for scan in ratio_multiplicative_scans)
        if not any_multiplicative_hit:
            lines.append("No multiplicative ratio-object relation was found in any scanned prefix box.")
            lines.append("")

        multiplicative_no_hit_labels = [f"`B{scan.powers[-1]}`" for scan in ratio_multiplicative_scans if scan.error is None]
        if not any_multiplicative_hit and multiplicative_no_hit_labels:
            lines.append(
                f"- No hit for multiplicative prefixes ending at {', '.join(multiplicative_no_hit_labels)}."
            )

        for scan in ratio_multiplicative_scans:
            if scan.error is not None:
                lines.append(
                    f"- Multiplicative prefix ending at `B{scan.powers[-1]}` skipped: {scan.error}"
                )
            elif scan.relation is not None:
                residual = _multiplicative_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        "B1": benchmark_series,
                        **{
                            f"B{p}": benchmark_power_substitution_series(
                                benchmark_series,
                                power=p,
                                order=profile_order,
                            )
                            for p in scan.powers
                        },
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Multiplicative prefix ending at `B{scan.powers[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_multiplicative_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

    if ratio_fractional_linear_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Fractional-Linear RR-Tower Scan",
                "",
                "We also searched for low-complexity fractional-linear corrections built from the benchmark tower:",
                "",
                "```text",
                "F = (1 + sum a_i*(B_i - 1)) / (1 + sum b_i*(B_i - 1))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                f"- `B1 = {record.closest_benchmark}`",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.",
                "- Each prefix solves an exact linear system for the numerator and denominator correction coefficients.",
                "",
            ]
        )

        any_fractional_hit = any(scan.relation is not None for scan in ratio_fractional_linear_scans)
        if not any_fractional_hit:
            lines.append("No fractional-linear ratio-object relation was found in any scanned prefix box.")
            lines.append("")

        prefix_labels = [f"`B{scan.powers[-1]}`" for scan in ratio_fractional_linear_scans if scan.error is None]
        if not any_fractional_hit and prefix_labels:
            lines.append(
                f"- No hit for fractional-linear prefixes ending at {', '.join(prefix_labels)}."
            )

        for scan in ratio_fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Fractional-linear prefix ending at `B{scan.powers[-1]}` skipped: {scan.error}"
                )
            elif scan.relation is not None:
                residual = _fractional_linear_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        "B1": benchmark_series,
                        **{
                            f"B{p}": benchmark_power_substitution_series(
                                benchmark_series,
                                power=p,
                                order=profile_order,
                            )
                            for p in scan.powers
                        },
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Fractional-linear prefix ending at `B{scan.powers[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_fractional_linear_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

    if ratio_two_layer_fractional_linear_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Layer Fractional-Linear RR-Tower Scan",
                "",
                "We then expanded to a second-ring nonlinear box built from two single-basis fractional-linear factors:",
                "",
                "```text",
                "F = ((1 + a0*(X1 - 1)) / (1 + b0*(Y1 - 1))) * ((1 + a1*(X2 - 1)) / (1 + b1*(Y2 - 1)))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                f"- `B1 = {record.closest_benchmark}`",
            ]
        )
        for power in sorted(benchmark_power_series):
            lines.append(f"- `B{power} = B1({series_symbol}^{power})`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(B1, B2)`, then `(B1, B2, B3)`, and so on through the final listed power.",
                "- Each prefix scans low-complexity two-factor templates and verifies any candidate hit by exact series re-expansion.",
                "",
            ]
        )

        any_two_layer_hit = any(scan.total_hits > 0 for scan in ratio_two_layer_fractional_linear_scans)
        if not any_two_layer_hit:
            lines.append("No two-layer fractional-linear ratio-object relation was found in any scanned prefix box.")
            lines.append("")

        no_hit_labels = [f"`B{scan.powers[-1]}`" for scan in ratio_two_layer_fractional_linear_scans if scan.error is None and scan.total_hits == 0]
        if not any_two_layer_hit and no_hit_labels:
            lines.append(
                f"- No hit for two-layer fractional-linear prefixes ending at {', '.join(no_hit_labels)}."
            )

        for scan in ratio_two_layer_fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Two-layer fractional-linear prefix ending at `B{scan.powers[-1]}` skipped: {scan.error}"
                )
                continue
            if scan.total_hits == 0:
                continue
            lines.append(
                f"- Two-layer fractional-linear prefix ending at `B{scan.powers[-1]}` found `{scan.total_hits}` exact hit(s) after checking `{scan.tuples_checked}` template(s):"
            )
            lines.append("")
            for relation in scan.relations:
                residual = _two_layer_fractional_linear_relation_residual_series(
                    relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        "B1": benchmark_series,
                        **{
                            f"B{p}": benchmark_power_substitution_series(
                                benchmark_series,
                                power=p,
                                order=profile_order,
                            )
                            for p in scan.powers
                        },
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        "```text",
                        _format_two_layer_fractional_linear_relation(relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
