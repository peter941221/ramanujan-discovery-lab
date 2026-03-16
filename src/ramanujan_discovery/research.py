from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Callable

import sympy as sp

from ramanujan_discovery.analysis import _series_expr
from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.storage import read_candidates


def reduce_template_by_step(template: QCFTemplate, step: int) -> QCFTemplate | None:
    """Return a new template expressed in t=q^step if all exponents are divisible by step."""
    if step <= 0:
        raise ValueError("step must be positive")

    def _maybe_div(value: int) -> int | None:
        if value % step != 0:
            return None
        return value // step

    numerator_q_shift = _maybe_div(template.numerator_q_shift)
    numerator_q_step = _maybe_div(template.numerator_q_step)
    if numerator_q_shift is None or numerator_q_step is None:
        return None

    numerator_extra_q_shift = template.numerator_extra_q_shift
    numerator_extra_q_step = template.numerator_extra_q_step
    if template.numerator_extra_scale != 0:
        numerator_extra_q_shift = _maybe_div(template.numerator_extra_q_shift)
        numerator_extra_q_step = _maybe_div(template.numerator_extra_q_step)
        if numerator_extra_q_shift is None or numerator_extra_q_step is None:
            return None
    else:
        numerator_extra_q_shift = 0
        numerator_extra_q_step = 0

    denominator_q_shift = template.denominator_q_shift
    denominator_q_step = template.denominator_q_step
    if template.denominator_scale != 0:
        denominator_q_shift = _maybe_div(template.denominator_q_shift)
        denominator_q_step = _maybe_div(template.denominator_q_step)
        if denominator_q_shift is None or denominator_q_step is None:
            return None
    else:
        denominator_q_shift = 0
        denominator_q_step = 0

    return QCFTemplate(
        top_constant=template.top_constant,
        base_denominator=template.base_denominator,
        numerator_scale=template.numerator_scale,
        numerator_q_shift=numerator_q_shift,
        numerator_q_step=numerator_q_step,
        numerator_extra_scale=template.numerator_extra_scale,
        numerator_extra_q_shift=numerator_extra_q_shift,
        numerator_extra_q_step=numerator_extra_q_step,
        denominator_constant=template.denominator_constant,
        denominator_scale=template.denominator_scale,
        denominator_q_shift=denominator_q_shift,
        denominator_q_step=denominator_q_step,
    ).normalized()


def _series_coeffs(expr, q: sp.Symbol, order: int) -> list[sp.Expr]:
    coeffs: list[sp.Expr] = []
    expanded = sp.expand(expr)
    for n in range(order):
        coeffs.append(sp.simplify(expanded.coeff(q, n)))
    return coeffs


def _monomial_coeff_map(expr, q: sp.Symbol) -> dict[int, sp.Expr]:
    coefficients: dict[int, sp.Expr] = {}
    expanded = sp.expand(expr)
    if expanded == 0:
        return coefficients
    for term in expanded.as_ordered_terms():
        coefficient, exponent = term.as_coeff_exponent(q)
        if coefficient == 0:
            continue
        if not exponent.is_integer:
            raise ValueError(f"non-integer exponent encountered: {exponent}")
        exponent_int = int(exponent)
        coefficients[exponent_int] = sp.simplify(coefficients.get(exponent_int, 0) + coefficient)
    return {exponent: coeff for exponent, coeff in coefficients.items() if sp.simplify(coeff) != 0}


def _coefficient_equations(left, right, q: sp.Symbol) -> list[sp.Equality]:
    left_map = _monomial_coeff_map(left, q=q)
    right_map = _monomial_coeff_map(right, q=q)
    equations: list[sp.Equality] = []
    for exponent in sorted(set(left_map) | set(right_map)):
        equations.append(sp.Eq(left_map.get(exponent, 0), right_map.get(exponent, 0)))
    return equations


def _exact_zero_equations(expr, q: sp.Symbol) -> list[sp.Equality]:
    numerator = sp.expand(sp.together(expr).as_numer_denom()[0])
    numerator_map = _monomial_coeff_map(numerator, q=q)
    return [sp.Eq(value, 0) for _, value in sorted(numerator_map.items())]


def _exact_match_equations(left, right, q: sp.Symbol) -> list[sp.Equality]:
    return _exact_zero_equations(left - right, q=q)


def _log_series_coeffs(series_coeffs: list[sp.Expr]) -> list[sp.Expr]:
    """Return coefficients g_n for log(F) where F = sum f_n q^n and f_0 == 1."""
    if not series_coeffs:
        raise ValueError("series_coeffs must be non-empty")
    if sp.simplify(series_coeffs[0] - 1) != 0:
        raise ValueError("series constant term must be 1 to expand log(F)")

    # g_0 is unused; return list indexed by n with g[0]=0.
    g: list[sp.Expr] = [sp.Integer(0)]
    for n in range(1, len(series_coeffs)):
        rhs = sp.Integer(0)
        for j in range(1, n):
            rhs += sp.Integer(j) * g[j] * series_coeffs[n - j]
        g_n = sp.simplify(series_coeffs[n] - rhs / sp.Integer(n))
        g.append(g_n)
    return g


def euler_product_exponents(expr, q: sp.Symbol, order: int) -> list[sp.Expr]:
    """Compute c_1..c_{order-1} such that F(q)=prod_{n>=1} (1-q^n)^{c_n} modulo q^order.

    This is the Euler transform of the q-series with constant term 1.
    """
    if order < 2:
        raise ValueError("order must be at least 2")

    coeffs = _series_coeffs(expr, q=q, order=order)
    log_coeffs = _log_series_coeffs(coeffs)

    # s_m := -m * [q^m] log(F) = sum_{d|m} d c_d
    c: list[sp.Expr] = [sp.Integer(0)] * order
    for m in range(1, order):
        s_m = sp.simplify(-sp.Integer(m) * log_coeffs[m])
        sub = sp.Integer(0)
        for d in range(1, m):
            if m % d == 0:
                sub += sp.Integer(d) * c[d]
        c[m] = sp.simplify((s_m - sub) / sp.Integer(m))
    return c[1:]


def _periodic_exponents(exponents: list[sp.Expr], period: int) -> list[sp.Expr] | None:
    if period <= 0:
        raise ValueError("period must be positive")
    if len(exponents) < period:
        return None
    for i in range(len(exponents) - period):
        if sp.simplify(exponents[i] - exponents[i + period]) != 0:
            return None
    # Return one period, 1-indexed residues shown as 1..period.
    return exponents[:period]


def try_fit_periodic_pochhammer(
    exponents: list[sp.Expr],
    max_period: int = 12,
    max_abs: int = 8,
) -> tuple[int, list[int]] | None:
    """Attempt to fit Euler exponents as a pure (q^r; q^m)_inf product with bounded integer powers.

    If F(q)=prod_{r=1..m} (q^r; q^m)_inf^{e_r}, then Euler exponents are periodic with period m.
    """
    for period in range(1, max_period + 1):
        period_values = _periodic_exponents(exponents, period=period)
        if period_values is None:
            continue
        integers: list[int] = []
        ok = True
        for value in period_values:
            value_simplified = sp.simplify(value)
            if not value_simplified.is_integer:
                ok = False
                break
            value_int = int(value_simplified)
            if abs(value_int) > max_abs:
                ok = False
                break
            integers.append(value_int)
        if ok:
            return period, integers
    return None


def try_fit_eta_quotient(
    exponents: list[sp.Expr],
    max_level: int = 12,
    max_abs: int = 8,
) -> tuple[int, dict[int, int]] | None:
    """Attempt to fit Euler exponents as an eta-quotient over a small level.

    For F(q)=prod_{d|N} (q^d; q^d)_inf^{e_d}, the Euler exponent c_n satisfies c_n=sum_{d|n, d|N} e_d.
    """
    if not exponents:
        return None

    # Use the earliest constraints; if exponents are huge or non-integers, it will fail quickly.
    for level in range(1, max_level + 1):
        divisors = [d for d in range(1, level + 1) if level % d == 0]
        unknowns = sp.symbols(" ".join(f"e{d}" for d in divisors), integer=True, seq=True)
        equations = []
        for n, c_n in enumerate(exponents, start=1):
            if n > len(exponents):
                break
            rhs = sp.Integer(0)
            for d, unknown in zip(divisors, unknowns):
                if n % d == 0:
                    rhs += unknown
            equations.append(sp.Eq(rhs, sp.simplify(c_n)))
            # Keep it small to avoid overconstraining with noisy high-order behavior.
            if n >= max_level:
                break

        solution = sp.linsolve(equations, unknowns)
        if not solution:
            continue
        candidate = next(iter(solution), None)
        if candidate is None:
            continue
        # linsolve may return parameters; reject unless fully determined integers.
        params = list(solution.free_symbols)
        if params:
            continue
        e_map: dict[int, int] = {}
        ok = True
        for d, value in zip(divisors, candidate):
            value_simplified = sp.simplify(value)
            if not value_simplified.is_integer:
                ok = False
                break
            value_int = int(value_simplified)
            if abs(value_int) > max_abs:
                ok = False
                break
            if value_int != 0:
                e_map[d] = value_int
        if ok:
            return level, e_map
    return None


def _two_modulus_residue(n: int, modulus: int) -> int:
    """Return the 1-indexed residue class used by the q-Pochhammer display helpers."""
    if n <= 0:
        raise ValueError("n must be positive")
    if modulus <= 0:
        raise ValueError("modulus must be positive")
    return ((n - 1) % modulus) + 1


def _choose_balancing_shift(*, left_values: list[int], right_values: list[int], max_abs: int) -> int | None:
    """Choose a gauge shift that keeps both residue lists inside the bounded box."""
    if max_abs < 0:
        raise ValueError("max_abs must be non-negative")
    lower = -max_abs
    upper = max_abs
    for value in left_values:
        lower = max(lower, -max_abs - value)
        upper = min(upper, max_abs - value)
    for value in right_values:
        lower = max(lower, value - max_abs)
        upper = min(upper, value + max_abs)
    if lower > upper:
        return None

    best_shift = lower
    best_key: tuple[int, int, int] | None = None
    for shift in range(lower, upper + 1):
        shifted_left = [value + shift for value in left_values]
        shifted_right = [value - shift for value in right_values]
        key = (
            max(abs(value) for value in shifted_left + shifted_right),
            sum(abs(value) for value in shifted_left + shifted_right),
            abs(shift),
        )
        if best_key is None or key < best_key:
            best_key = key
            best_shift = shift
    return best_shift


def _solve_two_modulus_pochhammer_fit(
    *,
    exponents: list[int],
    first_modulus: int,
    second_modulus: int,
    max_abs: int,
) -> TwoModulusPochhammerFit | None:
    """Solve c_n = e^(m1)_{n mod m1} + e^(m2)_{n mod m2} over one full lcm period."""
    period = math.lcm(first_modulus, second_modulus)
    if len(exponents) < period:
        return None

    left_edges: dict[int, list[tuple[int, int]]] = {residue: [] for residue in range(1, first_modulus + 1)}
    right_edges: dict[int, list[tuple[int, int]]] = {residue: [] for residue in range(1, second_modulus + 1)}
    for n in range(1, period + 1):
        left_residue = _two_modulus_residue(n, first_modulus)
        right_residue = _two_modulus_residue(n, second_modulus)
        weight = exponents[n - 1]
        left_edges[left_residue].append((right_residue, weight))
        right_edges[right_residue].append((left_residue, weight))

    left_base: list[int | None] = [None] * (first_modulus + 1)
    right_base: list[int | None] = [None] * (second_modulus + 1)
    left_values = [0] * (first_modulus + 1)
    right_values = [0] * (second_modulus + 1)

    for start_residue in range(1, first_modulus + 1):
        if left_base[start_residue] is not None:
            continue

        left_component: list[int] = []
        right_component: list[int] = []
        queue: list[tuple[str, int]] = [("left", start_residue)]
        left_base[start_residue] = 0
        left_component.append(start_residue)
        index = 0
        while index < len(queue):
            side, residue = queue[index]
            index += 1
            if side == "left":
                current = left_base[residue]
                assert current is not None
                for right_residue, weight in left_edges[residue]:
                    candidate = weight - current
                    existing = right_base[right_residue]
                    if existing is None:
                        right_base[right_residue] = candidate
                        right_component.append(right_residue)
                        queue.append(("right", right_residue))
                    elif existing != candidate:
                        return None
            else:
                current = right_base[residue]
                assert current is not None
                for left_residue, weight in right_edges[residue]:
                    candidate = weight - current
                    existing = left_base[left_residue]
                    if existing is None:
                        left_base[left_residue] = candidate
                        left_component.append(left_residue)
                        queue.append(("left", left_residue))
                    elif existing != candidate:
                        return None

        shift = _choose_balancing_shift(
            left_values=[left_base[residue] for residue in left_component if left_base[residue] is not None],
            right_values=[right_base[residue] for residue in right_component if right_base[residue] is not None],
            max_abs=max_abs,
        )
        if shift is None:
            return None
        for residue in left_component:
            base_value = left_base[residue]
            assert base_value is not None
            left_values[residue] = base_value + shift
        for residue in right_component:
            base_value = right_base[residue]
            assert base_value is not None
            right_values[residue] = base_value - shift

    return TwoModulusPochhammerFit(
        first_modulus=first_modulus,
        first_exponents=left_values[1:],
        second_modulus=second_modulus,
        second_exponents=right_values[1:],
    )


def try_fit_two_modulus_pochhammer(
    exponents: list[sp.Expr],
    max_modulus: int = 12,
    max_abs: int = 8,
) -> TwoModulusPochhammerFit | None:
    """Attempt to fit Euler exponents as a product of two shifted q-Pochhammer blocks.

    This solves

        c_n = e^(m1)_{n mod m1} + e^(m2)_{n mod m2}

    for bounded integer residue exponents at moduli `m1 < m2 <= max_modulus`.
    It is aimed at mixed-modulus products whose combined period can exceed the
    single-period search box, e.g. moduli `5` and `6` with lcm `30`.
    """
    if not exponents:
        return None

    integer_exponents: list[int] = []
    for value in exponents:
        value_simplified = sp.simplify(value)
        if not value_simplified.is_integer:
            return None
        integer_exponents.append(int(value_simplified))

    for first_modulus in range(1, max_modulus + 1):
        for second_modulus in range(first_modulus + 1, max_modulus + 1):
            period = math.lcm(first_modulus, second_modulus)
            if period <= max(first_modulus, second_modulus):
                continue
            fit = _solve_two_modulus_pochhammer_fit(
                exponents=integer_exponents,
                first_modulus=first_modulus,
                second_modulus=second_modulus,
                max_abs=max_abs,
            )
            if fit is None:
                continue
            mismatches = _verify_two_modulus_pochhammer_fit(
                euler_exponents=exponents,
                fit=fit,
                check_count=len(exponents),
            )
            if not mismatches:
                return fit
    return None


def _pochhammer_inf_symbol(*, variable: str, shift: int, step: int) -> str:
    """Return a compact ASCII representation of (variable^shift; variable^step)_inf."""
    if shift <= 0 or step <= 0:
        raise ValueError("shift and step must be positive")
    base = variable if shift == 1 else f"{variable}^{shift}"
    qstep = variable if step == 1 else f"{variable}^{step}"
    return f"({base}; {qstep})_inf"


def _format_periodic_pochhammer_closed_form(*, period: int, exponents: list[int], variable: str) -> str:
    """Format Π_{r=1..period} (variable^r; variable^period)_inf^{e_r}, skipping zero powers."""
    if period <= 0:
        raise ValueError("period must be positive")
    if len(exponents) != period:
        raise ValueError("exponents length must equal period")
    terms: list[str] = []
    for residue, exponent in enumerate(exponents, start=1):
        if exponent == 0:
            continue
        symbol = _pochhammer_inf_symbol(variable=variable, shift=residue, step=period)
        if exponent == 1:
            terms.append(symbol)
        else:
            terms.append(f"{symbol}^{exponent}")
    return " * ".join(terms) if terms else "1"


def _format_eta_quotient_closed_form(*, level: int, exponents_by_divisor: dict[int, int], variable: str) -> str:
    """Format Π_{d|level} (variable^d; variable^d)_inf^{e_d}, skipping zero powers."""
    if level <= 0:
        raise ValueError("level must be positive")
    terms: list[str] = []
    for divisor in sorted(exponents_by_divisor):
        exponent = exponents_by_divisor[divisor]
        if exponent == 0:
            continue
        symbol = _pochhammer_inf_symbol(variable=variable, shift=divisor, step=divisor)
        if exponent == 1:
            terms.append(symbol)
        else:
            terms.append(f"{symbol}^{exponent}")
    return " * ".join(terms) if terms else "1"


def _format_two_modulus_pochhammer_closed_form(*, fit: TwoModulusPochhammerFit, variable: str) -> str:
    """Format a product of two residue-class q-Pochhammer blocks."""
    terms = [
        _format_periodic_pochhammer_closed_form(
            period=fit.first_modulus,
            exponents=fit.first_exponents,
            variable=variable,
        ),
        _format_periodic_pochhammer_closed_form(
            period=fit.second_modulus,
            exponents=fit.second_exponents,
            variable=variable,
        ),
    ]
    nontrivial_terms = [term for term in terms if term != "1"]
    return " * ".join(nontrivial_terms) if nontrivial_terms else "1"


def _verify_periodic_pochhammer_fit(
    *,
    euler_exponents: list[sp.Expr],
    period: int,
    residue_exponents: list[int],
    check_count: int,
) -> list[int]:
    """Return the first mismatch indices n where c_n disagrees with the periodic rule."""
    if check_count < 1:
        raise ValueError("check_count must be positive")
    check_count = min(check_count, len(euler_exponents))
    mismatches: list[int] = []
    for n in range(1, check_count + 1):
        predicted = sp.Integer(residue_exponents[(n - 1) % period])
        actual = sp.simplify(euler_exponents[n - 1])
        if sp.simplify(actual - predicted) != 0:
            mismatches.append(n)
            if len(mismatches) >= 8:
                break
    return mismatches


def _verify_eta_quotient_fit(
    *,
    euler_exponents: list[sp.Expr],
    level: int,
    exponents_by_divisor: dict[int, int],
    check_count: int,
) -> list[int]:
    """Return the first mismatch indices n where c_n disagrees with the divisor-sum rule."""
    if check_count < 1:
        raise ValueError("check_count must be positive")
    check_count = min(check_count, len(euler_exponents))
    items = list(exponents_by_divisor.items())
    mismatches: list[int] = []
    for n in range(1, check_count + 1):
        predicted = sp.Integer(0)
        for divisor, exponent in items:
            if n % divisor == 0:
                predicted += sp.Integer(exponent)
        actual = sp.simplify(euler_exponents[n - 1])
        if sp.simplify(actual - predicted) != 0:
            mismatches.append(n)
            if len(mismatches) >= 8:
                break
    return mismatches


def _verify_two_modulus_pochhammer_fit(
    *,
    euler_exponents: list[sp.Expr],
    fit: TwoModulusPochhammerFit,
    check_count: int,
) -> list[int]:
    """Return the first mismatch indices n where c_n disagrees with a two-modulus rule."""
    if check_count < 1:
        raise ValueError("check_count must be positive")
    check_count = min(check_count, len(euler_exponents))
    mismatches: list[int] = []
    for n in range(1, check_count + 1):
        predicted = sp.Integer(fit.first_exponents[(n - 1) % fit.first_modulus])
        predicted += sp.Integer(fit.second_exponents[(n - 1) % fit.second_modulus])
        actual = sp.simplify(euler_exponents[n - 1])
        if sp.simplify(actual - predicted) != 0:
            mismatches.append(n)
            if len(mismatches) >= 8:
                break
    return mismatches


def _benchmark_product_closed_form_in_reduced_variable(benchmark_name: str, variable: str) -> str | None:
    """Return a closed-form product string for a known benchmark in the step-reduced variable.

    This is intentionally a small, high-confidence mapping aligned with `benchmarks.py`:
    the product is stated in the reduced variable (so q -> t) and avoids numeric evaluation.
    """
    if benchmark_name.startswith("rogers_ramanujan"):
        # (t; t^5)_inf (t^4; t^5)_inf / ((t^2; t^5)_inf (t^3; t^5)_inf)
        return (
            f"{_pochhammer_inf_symbol(variable=variable, shift=1, step=5)}"
            f" * {_pochhammer_inf_symbol(variable=variable, shift=4, step=5)}"
            f" / ({_pochhammer_inf_symbol(variable=variable, shift=2, step=5)}"
            f" * {_pochhammer_inf_symbol(variable=variable, shift=3, step=5)})"
        )
    if benchmark_name.startswith("ramanujan_cubic"):
        # (t; t^6)_inf (t^5; t^6)_inf / (t^3; t^6)_inf^2
        base = (
            f"{_pochhammer_inf_symbol(variable=variable, shift=1, step=6)}"
            f" * {_pochhammer_inf_symbol(variable=variable, shift=5, step=6)}"
        )
        denom = _pochhammer_inf_symbol(variable=variable, shift=3, step=6)
        return f"{base} / ({denom}^2)"
    return None


@dataclass(frozen=True)
class HeineHCF2Coeffs:
    b0: sp.Expr
    a_terms: list[sp.Expr]  # 1-indexed, a_terms[0] unused
    b_terms: list[sp.Expr]  # 1-indexed, b_terms[0] unused


@dataclass(frozen=True)
class ContinuedFractionCoeffs:
    b0: sp.Expr
    a_terms: list[sp.Expr]  # 1-indexed, a_terms[0] unused
    b_terms: list[sp.Expr]  # 1-indexed, b_terms[0] unused


@dataclass(frozen=True)
class ConvergentCommonFactorReduction:
    gcd_factors: list[sp.Expr]
    reduced_coeffs: ContinuedFractionCoeffs


@dataclass(frozen=True)
class EquivalenceTransformedCoeffs:
    b0: sp.Expr
    a_terms: list[sp.Expr]  # 1-indexed, a_terms[0] unused
    b_terms: list[sp.Expr]  # 1-indexed, b_terms[0] unused


@dataclass(frozen=True)
class ConvergentFactorEquivalenceWitness:
    reduction: ConvergentCommonFactorReduction
    scale_terms: list[sp.Expr]  # 1-indexed, scale_terms[0] unused
    retransformed_coeffs: EquivalenceTransformedCoeffs


@dataclass(frozen=True)
class ResearchBuildProfile:
    label: str
    raw_series_cap: int
    reduced_order_cap: int
    euler_order_cap: int
    fit_order_cap: int
    factor_depth_cap: int
    page43_max_shift: int
    page43_stages: int
    max_subsequence_stride: int
    subsequence_stages: int
    max_bauer_muir_steps: int
    bauer_muir_depth: int


@dataclass(frozen=True)
class Page43ParameterHit:
    family: str
    a_shift: int
    b_shift: int
    lambda_shift: int
    a_coeff: sp.Expr
    b_coeff: sp.Expr
    lambda_coeff: sp.Expr
    a_profile: str = "1"
    b_profile: str = "1"
    lambda_profile: str = "1"


Page43MonomialHit = Page43ParameterHit


@dataclass(frozen=True)
class SubsequenceContractionHit:
    source_label: str
    stride: int
    offset: int


@dataclass(frozen=True)
class HeineCor2CFContractionObstruction:
    source_coeffs: ContinuedFractionCoeffs
    odd_part: ContinuedFractionCoeffs
    even_part: ContinuedFractionCoeffs
    even_odd_part: ContinuedFractionCoeffs
    even_even_part: ContinuedFractionCoeffs


@dataclass(frozen=True)
class Page43F2EquivalenceObstruction:
    residual_polynomial: sp.Expr
    m_coefficients: dict[int, sp.Expr]
    forced_ab_solution: dict[sp.Symbol, sp.Expr]
    reduced_m1_coefficient: sp.Expr
    forced_lambda_solution: dict[sp.Symbol, sp.Expr]
    final_m2_coefficient: sp.Expr


@dataclass(frozen=True)
class Page43F4EquivalenceObstruction:
    residual_polynomial: sp.Expr
    m_coefficients: dict[int, sp.Expr]
    forced_a_solution: dict[sp.Symbol, sp.Expr]
    forced_b_solution: dict[sp.Symbol, sp.Expr]
    reduced_m1_coefficient: sp.Expr
    forced_lambda_solution: dict[sp.Symbol, sp.Expr]
    final_m2_coefficient: sp.Expr


@dataclass(frozen=True)
class TwoModulusPochhammerFit:
    first_modulus: int
    first_exponents: list[int]
    second_modulus: int
    second_exponents: list[int]


def heine_hcf2_standardized_coeffs(
    *,
    a: sp.Expr,
    b: sp.Expr,
    c: sp.Expr,
    z: sp.Expr,
    q: sp.Symbol,
    depth: int,
) -> HeineHCF2Coeffs:
    """Coefficient extraction for the theorem-level hcf2 continued fraction (standardized).

    This follows the displayed hcf2 in Lee–Mc Laughlin–Sohn 2019 (arXiv:1906.11991),
    with an equivalence normalization that removes the extra (1-c) factor in b1.
    """
    if depth < 1:
        raise ValueError("depth must be at least 1")

    b0 = sp.simplify((1 - b * z) / (1 - c))
    a_terms = [sp.Integer(0)]
    b_terms = [sp.Integer(0)]

    # Normalization: choose r1 = 1/(1-c), r_n = 1 for n>=2 so that b1 becomes 1 - b z q.
    a1 = sp.simplify((c - a * b * z) * (z - 1) / (1 - c))
    b1 = sp.simplify(1 - b * z * q)
    a_terms.append(a1)
    b_terms.append(b1)

    for n in range(2, depth + 1):
        if n % 2 == 0:
            k = n // 2
            a_n = sp.simplify((1 - b * q**k) * (c * q**k - a) * z * q ** (k - 1))
        else:
            k = (n - 1) // 2
            a_n = sp.simplify((c - a * b * z * q**k) * (z * q**k - 1) * q**k)
        b_n = sp.simplify(1 - b * z * q**n)
        a_terms.append(a_n)
        b_terms.append(b_n)

    return HeineHCF2Coeffs(b0=b0, a_terms=a_terms, b_terms=b_terms)


def heine_cor2cf_a_zero_specialized_coeffs(
    *,
    b: sp.Expr,
    lam: sp.Expr,
    q: sp.Symbol,
    depth: int,
) -> ContinuedFractionCoeffs:
    """Return the `cor2cf` lane that remains after forcing the initial term to be `1`.

    In the notation of the transform audit, this is the `a = 0` specialization

        1 + (lam*q)/1 + (b*q + lam*q^2)/1 + (lam*q^3)/1 + (b*q^2 + lam*q^4)/1 + ...

    which is the only nearby `cor2cf` branch compatible with the target initial term.
    """
    if depth < 1:
        raise ValueError("depth must be at least 1")

    a_terms = [sp.Integer(0)]
    b_terms = [sp.Integer(0)]
    for n in range(1, depth + 1):
        if n % 2 == 1:
            stage = (n + 1) // 2
            a_n = lam * q ** (2 * stage - 1)
        else:
            stage = n // 2
            a_n = b * q**stage + lam * q ** (2 * stage)
        a_terms.append(sp.simplify(a_n))
        b_terms.append(sp.Integer(1))

    return ContinuedFractionCoeffs(
        b0=sp.Integer(1),
        a_terms=a_terms,
        b_terms=b_terms,
    )


def heine_cor2cf_a_zero_contraction_obstruction(
    *,
    b: sp.Expr,
    lam: sp.Expr,
    q: sp.Symbol,
    depth: int = 12,
) -> HeineCor2CFContractionObstruction:
    """Compute exact odd/even contraction data for the relevant `cor2cf` specialization."""
    if depth < 8:
        raise ValueError("depth must be at least 8 to recover the two-step contraction branches")

    source = heine_cor2cf_a_zero_specialized_coeffs(b=b, lam=lam, q=q, depth=depth)
    odd = parity_contraction_coeffs(
        b0=source.b0,
        a_terms=source.a_terms,
        b_terms=source.b_terms,
        parity="odd",
    )
    even = parity_contraction_coeffs(
        b0=source.b0,
        a_terms=source.a_terms,
        b_terms=source.b_terms,
        parity="even",
    )
    even_odd = parity_contraction_coeffs(
        b0=even.b0,
        a_terms=even.a_terms,
        b_terms=even.b_terms,
        parity="odd",
    )
    even_even = parity_contraction_coeffs(
        b0=even.b0,
        a_terms=even.a_terms,
        b_terms=even.b_terms,
        parity="even",
    )

    return HeineCor2CFContractionObstruction(
        source_coeffs=source,
        odd_part=odd,
        even_part=even,
        even_odd_part=even_odd,
        even_even_part=even_even,
    )


def page43_f2_zero_shift_equivalence_obstruction(
    *,
    q: sp.Symbol,
) -> Page43F2EquivalenceObstruction:
    """Return the exact zero-shift `f2` / `gcf3` equivalence obstruction.

    This is the source-family-specific lane where the page-43 `f2` family is tested
    against the hero-case reciprocal under an arbitrary `n`-dependent equivalence
    transformation, with the low-complexity zero-shift specialization

        a = alpha,  b = beta,  lambda = gamma.

    Writing `m = q^(n-1)`, the necessary identity becomes a polynomial in `m`
    and `q`. Vanishing identically would force the source-family parameters to
    satisfy impossible exact coefficient constraints.
    """
    a, b, lam, m = sp.symbols("a b lambda m")

    alpha_n = sp.simplify(lam * m * q - a * b * m**2 * q**2)
    beta_prev = sp.simplify(1 + b * m + a * m * q)
    beta_curr = sp.simplify(1 + b * m * q + a * m * q**2)
    residual = sp.expand(alpha_n * (1 + m) - m * q * beta_prev * beta_curr)

    residual_poly = sp.Poly(residual, m)
    m_coefficients = {
        degree: sp.expand(residual_poly.nth(degree))
        for degree in range(1, residual_poly.degree() + 1)
        if sp.simplify(residual_poly.nth(degree)) != 0
    }

    leading_coefficient = m_coefficients[max(m_coefficients)]
    forced_ab_solutions = sp.solve(_exact_zero_equations(leading_coefficient, q=q), (a, b), dict=True)
    if forced_ab_solutions != [{a: sp.Integer(0), b: sp.Integer(0)}]:
        raise ValueError("unexpected exact-solution set for the zero-shift `f2` equivalence leading coefficient")
    forced_ab_solution = forced_ab_solutions[0]

    reduced_m1_coefficient = sp.simplify(m_coefficients[1].subs(forced_ab_solution))
    forced_lambda_solutions = sp.solve(_exact_zero_equations(reduced_m1_coefficient, q=q), (lam,), dict=True)
    if forced_lambda_solutions != [{lam: sp.Integer(1)}]:
        raise ValueError("unexpected exact-solution set for the zero-shift `f2` equivalence linear coefficient")
    forced_lambda_solution = forced_lambda_solutions[0]

    final_m2_coefficient = sp.simplify(
        m_coefficients[2].subs(forced_ab_solution).subs(forced_lambda_solution)
    )
    if _exact_zero_equations(final_m2_coefficient, q=q) == []:
        raise ValueError("unexpected vanishing final obstruction in the zero-shift `f2` equivalence lane")

    return Page43F2EquivalenceObstruction(
        residual_polynomial=residual,
        m_coefficients=m_coefficients,
        forced_ab_solution=forced_ab_solution,
        reduced_m1_coefficient=reduced_m1_coefficient,
        forced_lambda_solution=forced_lambda_solution,
        final_m2_coefficient=final_m2_coefficient,
    )


def page43_f4_zero_shift_equivalence_obstruction(
    *,
    q: sp.Symbol,
) -> Page43F4EquivalenceObstruction:
    """Return the exact zero-shift `f4` / `gcf2` equivalence obstruction."""
    a, b, lam, m = sp.symbols("a b lambda m")

    alpha_n = sp.simplify(a * q + lam * m * q)
    beta_prev = sp.simplify(1 - a * q + b * m)
    beta_curr = sp.simplify(1 - a * q + b * m * q)
    residual = sp.expand(alpha_n * (1 + m) - m * q * beta_prev * beta_curr)

    residual_poly = sp.Poly(residual, m)
    m_coefficients = {
        degree: sp.expand(residual_poly.nth(degree))
        for degree in range(0, residual_poly.degree() + 1)
        if sp.simplify(residual_poly.nth(degree)) != 0
    }

    forced_a_solutions = sp.solve(_exact_zero_equations(m_coefficients[0], q=q), (a,), dict=True)
    if forced_a_solutions != [{a: sp.Integer(0)}]:
        raise ValueError("unexpected exact-solution set for the zero-shift `f4` equivalence constant coefficient")
    forced_a_solution = forced_a_solutions[0]

    reduced_m3_coefficient = sp.simplify(m_coefficients[3].subs(forced_a_solution))
    forced_b_solutions = sp.solve(_exact_zero_equations(reduced_m3_coefficient, q=q), (b,), dict=True)
    if forced_b_solutions != [{b: sp.Integer(0)}]:
        raise ValueError("unexpected exact-solution set for the zero-shift `f4` equivalence cubic coefficient")
    forced_b_solution = forced_b_solutions[0]

    reduced_m1_coefficient = sp.simplify(
        m_coefficients[1].subs(forced_a_solution).subs(forced_b_solution)
    )
    forced_lambda_solutions = sp.solve(_exact_zero_equations(reduced_m1_coefficient, q=q), (lam,), dict=True)
    if forced_lambda_solutions != [{lam: sp.Integer(1)}]:
        raise ValueError("unexpected exact-solution set for the zero-shift `f4` equivalence linear coefficient")
    forced_lambda_solution = forced_lambda_solutions[0]

    final_m2_coefficient = sp.simplify(
        m_coefficients[2]
        .subs(forced_a_solution)
        .subs(forced_b_solution)
        .subs(forced_lambda_solution)
    )
    if _exact_zero_equations(final_m2_coefficient, q=q) == []:
        raise ValueError("unexpected vanishing final obstruction in the zero-shift `f4` equivalence lane")

    return Page43F4EquivalenceObstruction(
        residual_polynomial=residual,
        m_coefficients=m_coefficients,
        forced_a_solution=forced_a_solution,
        forced_b_solution=forced_b_solution,
        reduced_m1_coefficient=reduced_m1_coefficient,
        forced_lambda_solution=forced_lambda_solution,
        final_m2_coefficient=final_m2_coefficient,
    )


def continued_fraction_convergents(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
) -> list[tuple[sp.Expr, sp.Expr]]:
    """Return convergent numerator/denominator pairs for b0 + K a_n / b_n."""
    if len(a_terms) != len(b_terms):
        raise ValueError("a_terms and b_terms must have the same length")
    if len(a_terms) < 2:
        raise ValueError("a_terms and b_terms must be 1-indexed with at least one term")

    convergents: list[tuple[sp.Expr, sp.Expr]] = [(sp.simplify(b0), sp.Integer(1))]
    a_prev2 = sp.Integer(1)
    a_prev1 = sp.simplify(b0)
    b_prev2 = sp.Integer(0)
    b_prev1 = sp.Integer(1)

    for n in range(1, len(a_terms)):
        a_n = sp.expand(b_terms[n] * a_prev1 + a_terms[n] * a_prev2)
        b_n = sp.expand(b_terms[n] * b_prev1 + a_terms[n] * b_prev2)
        convergents.append((a_n, b_n))
        a_prev2, a_prev1 = a_prev1, a_n
        b_prev2, b_prev1 = b_prev1, b_n

    return convergents


def recover_cf_from_convergents(convergents: list[tuple[sp.Expr, sp.Expr]]) -> ContinuedFractionCoeffs:
    """Recover one generalized continued fraction whose convergents match the supplied sequence."""
    if not convergents:
        raise ValueError("convergents must be non-empty")

    b0 = sp.simplify(convergents[0][0] / convergents[0][1])
    a_terms = [sp.Integer(0)]
    b_terms = [sp.Integer(0)]

    a_prev2 = sp.Integer(1)
    b_prev2 = sp.Integer(0)
    a_prev1 = convergents[0][0]
    b_prev1 = convergents[0][1]

    for numerator, denominator in convergents[1:]:
        determinant = sp.simplify(a_prev1 * b_prev2 - b_prev1 * a_prev2)
        if determinant == 0:
            raise ValueError("degenerate convergent sequence")

        b_n = sp.simplify((numerator * b_prev2 - denominator * a_prev2) / determinant)
        a_n = sp.simplify((denominator * a_prev1 - numerator * b_prev1) / determinant)
        a_terms.append(sp.expand(a_n))
        b_terms.append(sp.expand(b_n))

        next_a = sp.expand(b_n * a_prev1 + a_n * a_prev2)
        next_b = sp.expand(b_n * b_prev1 + a_n * b_prev2)
        a_prev2, a_prev1 = a_prev1, next_a
        b_prev2, b_prev1 = b_prev1, next_b

    return ContinuedFractionCoeffs(b0=b0, a_terms=a_terms, b_terms=b_terms)


def convergent_common_factor_reduction(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
) -> ConvergentCommonFactorReduction:
    """Cancel the exact gcd of each convergent pair and recover the induced continued fraction."""
    convergents = continued_fraction_convergents(b0=b0, a_terms=a_terms, b_terms=b_terms)
    reduced_convergents: list[tuple[sp.Expr, sp.Expr]] = []
    gcd_factors: list[sp.Expr] = []

    for numerator, denominator in convergents:
        common_factor = sp.factor(sp.gcd(sp.expand(numerator), sp.expand(denominator)))
        gcd_factors.append(common_factor)
        reduced_convergents.append(
            (
                sp.expand(sp.cancel(numerator / common_factor)),
                sp.expand(sp.cancel(denominator / common_factor)),
            )
        )

    return ConvergentCommonFactorReduction(
        gcd_factors=gcd_factors,
        reduced_coeffs=recover_cf_from_convergents(reduced_convergents),
    )


def apply_equivalence_transform(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
    scale_terms: list[sp.Expr],
) -> EquivalenceTransformedCoeffs:
    """Apply an equivalence transformation with stage scales r_n to b0 + K a_n / b_n.

    The input coefficient lists are 1-indexed. `scale_terms[n]` is the multiplier r_n
    for the n-th stage, with `scale_terms[0]` unused.
    """
    if len(a_terms) != len(b_terms):
        raise ValueError("a_terms and b_terms must have the same length")
    if len(a_terms) < 2:
        raise ValueError("a_terms and b_terms must be 1-indexed with at least one term")
    if len(scale_terms) != len(a_terms):
        raise ValueError("scale_terms must be 1-indexed and match a_terms length")

    transformed_a = [sp.Integer(0)]
    transformed_b = [sp.Integer(0)]
    for n in range(1, len(a_terms)):
        transformed_b.append(sp.expand(scale_terms[n] * b_terms[n]))
        if n == 1:
            transformed_a.append(sp.expand(scale_terms[n] * a_terms[n]))
        else:
            transformed_a.append(sp.expand(scale_terms[n - 1] * scale_terms[n] * a_terms[n]))

    return EquivalenceTransformedCoeffs(
        b0=sp.simplify(b0),
        a_terms=transformed_a,
        b_terms=transformed_b,
    )


def convergent_factor_equivalence_witness(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
) -> ConvergentFactorEquivalenceWitness:
    """Recover the convergent-factor reduction plus the reverse equivalence scales."""
    reduction = convergent_common_factor_reduction(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
    )

    scale_terms = [sp.Integer(0)]
    if len(reduction.gcd_factors) > 1:
        scale_terms.append(sp.simplify(reduction.gcd_factors[1]))
        for n in range(2, len(reduction.gcd_factors)):
            previous_factor = sp.simplify(reduction.gcd_factors[n - 1])
            if previous_factor == 0:
                raise ValueError("gcd factors must stay nonzero to recover equivalence scales")
            scale_terms.append(sp.simplify(reduction.gcd_factors[n] / previous_factor))

    retransformed = apply_equivalence_transform(
        b0=reduction.reduced_coeffs.b0,
        a_terms=reduction.reduced_coeffs.a_terms,
        b_terms=reduction.reduced_coeffs.b_terms,
        scale_terms=scale_terms,
    )
    return ConvergentFactorEquivalenceWitness(
        reduction=reduction,
        scale_terms=scale_terms,
        retransformed_coeffs=retransformed,
    )


def parity_contraction_coeffs(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
    parity: str,
) -> ContinuedFractionCoeffs:
    """Recover the odd or even canonical contraction from a truncated source continued fraction."""
    if parity == "odd":
        return arithmetic_subsequence_contraction_coeffs(
            b0=b0,
            a_terms=a_terms,
            b_terms=b_terms,
            stride=2,
            offset=1,
        )
    if parity == "even":
        return arithmetic_subsequence_contraction_coeffs(
            b0=b0,
            a_terms=a_terms,
            b_terms=b_terms,
            stride=2,
            offset=0,
        )
    raise ValueError("parity must be 'odd' or 'even'")


def arithmetic_subsequence_contraction_coeffs(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
    stride: int,
    offset: int,
) -> ContinuedFractionCoeffs:
    """Recover the contraction determined by every `stride`-th convergent starting at `offset`."""
    if stride < 1:
        raise ValueError("stride must be positive")
    if offset < 0 or offset >= stride:
        raise ValueError("offset must satisfy 0 <= offset < stride")

    convergents = continued_fraction_convergents(b0=b0, a_terms=a_terms, b_terms=b_terms)
    selected = [convergents[index] for index in range(offset, len(convergents), stride)]
    if not selected:
        raise ValueError(f"source continued fraction has no subsequence for stride={stride}, offset={offset}")

    return recover_cf_from_convergents(selected)


def arithmetic_subsequence_contraction_search(
    *,
    source_label: str,
    source_template: QCFTemplate,
    target_template: QCFTemplate,
    q: sp.Symbol,
    max_stride: int = 4,
    stages: int = 3,
) -> list[SubsequenceContractionHit]:
    """Search simple arithmetic convergent subsequences of a nearby source reciprocal."""
    if max_stride < 2:
        raise ValueError("max_stride must be at least 2")
    if stages < 1:
        raise ValueError("stages must be at least 1")

    source_depth = stages * max_stride + (max_stride - 1)
    source_b0, source_a_terms, source_b_terms = _template_reciprocal_coeffs(
        source_template.normalized(),
        q=q,
        depth=source_depth,
    )
    target_b0, target_a_terms, target_b_terms = _template_reciprocal_coeffs(
        target_template.normalized(),
        q=q,
        depth=stages,
    )

    hits: list[SubsequenceContractionHit] = []
    for stride in range(2, max_stride + 1):
        for offset in range(stride):
            contraction = arithmetic_subsequence_contraction_coeffs(
                b0=source_b0,
                a_terms=source_a_terms,
                b_terms=source_b_terms,
                stride=stride,
                offset=offset,
            )
            if sp.simplify(contraction.b0 - target_b0) != 0:
                continue
            exact = True
            for n in range(1, stages + 1):
                if sp.simplify(contraction.a_terms[n] - target_a_terms[n]) != 0:
                    exact = False
                    break
                if sp.simplify(contraction.b_terms[n] - target_b_terms[n]) != 0:
                    exact = False
                    break
            if exact:
                hits.append(SubsequenceContractionHit(source_label=source_label, stride=stride, offset=offset))

    return hits


def page43_monomial_parameter_search(
    *,
    family: str,
    target_template: QCFTemplate,
    q: sp.Symbol,
    max_shift: int = 3,
    stages: int = 3,
) -> list[Page43MonomialHit]:
    """Search bounded monomial parameter substitutions in the page-43 gcf2/gcf3 families."""
    if family not in {"f2", "f4"}:
        raise ValueError("family must be 'f2' or 'f4'")
    if stages < 1:
        raise ValueError("stages must be at least 1")

    _, target_a_terms, target_b_terms = _template_reciprocal_coeffs(target_template.normalized(), q=q, depth=stages)
    alpha, beta, gamma = sp.symbols("alpha beta gamma")
    hits: list[Page43MonomialHit] = []

    for a_shift in range(-max_shift, max_shift + 1):
        for b_shift in range(-max_shift, max_shift + 1):
            for lambda_shift in range(-max_shift, max_shift + 1):
                equations: list[sp.Equality] = []
                for n in range(1, stages + 1):
                    source_a, source_b = _page43_family_terms(
                        family=family,
                        alpha=alpha,
                        beta=beta,
                        gamma=gamma,
                        a_shift=a_shift,
                        b_shift=b_shift,
                        lambda_shift=lambda_shift,
                        q=q,
                        n=n,
                    )
                    equations.extend(_coefficient_equations(source_a, target_a_terms[n], q=q))
                    equations.extend(_coefficient_equations(source_b, target_b_terms[n], q=q))

                solutions = sp.solve(equations, (alpha, beta, gamma), dict=True)
                if not solutions:
                    continue

                for solution in solutions:
                    if any(value.free_symbols for value in solution.values()):
                        continue

                    exact_match = True
                    for n in range(1, stages + 1):
                        source_a, source_b = _page43_family_terms(
                            family=family,
                            alpha=alpha,
                            beta=beta,
                            gamma=gamma,
                            a_shift=a_shift,
                            b_shift=b_shift,
                            lambda_shift=lambda_shift,
                            q=q,
                            n=n,
                        )
                        if sp.simplify(source_a.subs(solution) - target_a_terms[n]) != 0:
                            exact_match = False
                            break
                        if sp.simplify(source_b.subs(solution) - target_b_terms[n]) != 0:
                            exact_match = False
                            break

                    if exact_match:
                        hits.append(
                            Page43MonomialHit(
                                family=family,
                                a_shift=a_shift,
                                b_shift=b_shift,
                                lambda_shift=lambda_shift,
                                a_coeff=sp.simplify(solution[alpha]),
                                b_coeff=sp.simplify(solution[beta]),
                                lambda_coeff=sp.simplify(solution[gamma]),
                            )
                        )

    return hits


def _page43_plus_profile_catalog(q: sp.Symbol) -> list[tuple[str, sp.Expr]]:
    variable = str(q)
    return [
        ("1", sp.Integer(1)),
        (f"1 + {variable}", 1 + q),
        (f"1 / (1 + {variable})", sp.simplify(1 / (1 + q))),
    ]


def _page43_family_terms(
    *,
    family: str,
    alpha: sp.Expr,
    beta: sp.Expr,
    gamma: sp.Expr,
    a_shift: int,
    b_shift: int,
    lambda_shift: int,
    q: sp.Symbol,
    n: int,
) -> tuple[sp.Expr, sp.Expr]:
    if family not in {"f2", "f4"}:
        raise ValueError("family must be 'f2' or 'f4'")
    if family == "f2":
        source_a = gamma * q ** (lambda_shift + n) - alpha * beta * q ** (a_shift + b_shift + 2 * n)
        source_b = 1 + beta * q ** (b_shift + n) + alpha * q ** (a_shift + n + 1)
    else:
        source_a = alpha * q ** (a_shift + 1) + gamma * q ** (lambda_shift + n)
        source_b = 1 - alpha * q ** (a_shift + 1) + beta * q ** (b_shift + n)
    return sp.simplify(source_a), sp.simplify(source_b)


def page43_rational_parameter_search(
    *,
    family: str,
    target_template: QCFTemplate,
    q: sp.Symbol,
    max_shift: int = 3,
    stages: int = 3,
    profile_catalog: list[tuple[str, sp.Expr]] | None = None,
    max_nontrivial_profiles: int = 1,
) -> list[Page43ParameterHit]:
    """Search bounded low-complexity rational substitutions in the page-43 gcf2/gcf3 families."""
    if family not in {"f2", "f4"}:
        raise ValueError("family must be 'f2' or 'f4'")
    if stages < 1:
        raise ValueError("stages must be at least 1")
    if max_nontrivial_profiles < 0:
        raise ValueError("max_nontrivial_profiles must be non-negative")

    _, target_a_terms, target_b_terms = _template_reciprocal_coeffs(target_template.normalized(), q=q, depth=stages)
    alpha, beta, gamma = sp.symbols("alpha beta gamma")
    profiles = profile_catalog or _page43_plus_profile_catalog(q)
    hits: list[Page43ParameterHit] = []
    seen: set[tuple[str, int, int, int, str, str, str, str, str, str]] = set()

    for a_profile_label, a_profile in profiles:
        for b_profile_label, b_profile in profiles:
            for lambda_profile_label, lambda_profile in profiles:
                if (
                    sum(
                        label != "1"
                        for label in (a_profile_label, b_profile_label, lambda_profile_label)
                    )
                    > max_nontrivial_profiles
                ):
                    continue
                for a_shift in range(-max_shift, max_shift + 1):
                    for b_shift in range(-max_shift, max_shift + 1):
                        for lambda_shift in range(-max_shift, max_shift + 1):
                            equations: list[sp.Equality] = []
                            for n in range(1, stages + 1):
                                source_a, source_b = _page43_family_terms(
                                    family=family,
                                    alpha=alpha * a_profile,
                                    beta=beta * b_profile,
                                    gamma=gamma * lambda_profile,
                                    a_shift=a_shift,
                                    b_shift=b_shift,
                                    lambda_shift=lambda_shift,
                                    q=q,
                                    n=n,
                                )
                                equations.extend(_exact_match_equations(source_a, target_a_terms[n], q=q))
                                equations.extend(_exact_match_equations(source_b, target_b_terms[n], q=q))

                            solutions = sp.solve(equations, (alpha, beta, gamma), dict=True)
                            if not solutions:
                                continue

                            for solution in solutions:
                                if any(value.free_symbols for value in solution.values()):
                                    continue

                                exact_match = True
                                for n in range(1, stages + 1):
                                    source_a, source_b = _page43_family_terms(
                                        family=family,
                                        alpha=alpha * a_profile,
                                        beta=beta * b_profile,
                                        gamma=gamma * lambda_profile,
                                        a_shift=a_shift,
                                        b_shift=b_shift,
                                        lambda_shift=lambda_shift,
                                        q=q,
                                        n=n,
                                    )
                                    if sp.simplify(source_a.subs(solution) - target_a_terms[n]) != 0:
                                        exact_match = False
                                        break
                                    if sp.simplify(source_b.subs(solution) - target_b_terms[n]) != 0:
                                        exact_match = False
                                        break

                                if exact_match:
                                    hit = Page43ParameterHit(
                                        family=family,
                                        a_shift=a_shift,
                                        b_shift=b_shift,
                                        lambda_shift=lambda_shift,
                                        a_coeff=sp.simplify(solution[alpha]),
                                        b_coeff=sp.simplify(solution[beta]),
                                        lambda_coeff=sp.simplify(solution[gamma]),
                                        a_profile=a_profile_label,
                                        b_profile=b_profile_label,
                                        lambda_profile=lambda_profile_label,
                                    )
                                    key = (
                                        hit.family,
                                        hit.a_shift,
                                        hit.b_shift,
                                        hit.lambda_shift,
                                        sp.srepr(hit.a_coeff),
                                        sp.srepr(hit.b_coeff),
                                        sp.srepr(hit.lambda_coeff),
                                        hit.a_profile if sp.simplify(hit.a_coeff) != 0 else "0",
                                        hit.b_profile if sp.simplify(hit.b_coeff) != 0 else "0",
                                        hit.lambda_profile if sp.simplify(hit.lambda_coeff) != 0 else "0",
                                    )
                                    if key not in seen:
                                        seen.add(key)
                                        hits.append(hit)

    return hits


@dataclass(frozen=True)
class BauerMuirTransformResult:
    b0: sp.Expr
    a_terms: list[sp.Expr]  # 1-indexed
    b_terms: list[sp.Expr]  # 1-indexed
    lambdas: list[sp.Expr]  # 1-indexed


@dataclass(frozen=True)
class BauerMuirDirectObstruction:
    source_label: str
    forced_w0: sp.Expr
    forced_w1: sp.Expr
    forced_w2: sp.Expr | None
    transformed_a1: sp.Expr
    target_a1: sp.Expr
    transformed_a2: sp.Expr | None
    target_a2: sp.Expr | None
    obstruction_stage: int | None


def bauer_muir_transform_trunc(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
    w_terms: list[sp.Expr],
) -> BauerMuirTransformResult:
    """Truncated Bauer–Muir transform for a generalized continued fraction.

    Input sequences are 1-indexed for a_n,b_n and 0-indexed for w_n.
    """
    if len(a_terms) != len(b_terms):
        raise ValueError("a_terms and b_terms must have the same length")
    if len(a_terms) < 2:
        raise ValueError("a_terms and b_terms must be 1-indexed with at least one term")
    depth = len(a_terms) - 1
    if len(w_terms) < depth + 1:
        raise ValueError("w_terms must include indices 0..depth")

    lam = [sp.Integer(0)]
    for n in range(1, depth + 1):
        lam.append(sp.simplify(a_terms[n] - w_terms[n - 1] * (b_terms[n] + w_terms[n])))

    b0_prime = sp.simplify(b0 + w_terms[0])
    a_prime = [sp.Integer(0)]
    b_prime = [sp.Integer(0)]
    # n=1
    a_prime.append(lam[1])
    b_prime.append(sp.simplify(b_terms[1] + w_terms[1]))
    for n in range(2, depth + 1):
        a_prime.append(sp.simplify(a_terms[n - 1] * lam[n] / lam[n - 1]))
        b_prime.append(sp.simplify(b_terms[n] + w_terms[n] - w_terms[n - 2] * lam[n] / lam[n - 1]))

    return BauerMuirTransformResult(b0=b0_prime, a_terms=a_prime, b_terms=b_prime, lambdas=lam)


def direct_bauer_muir_obstruction(
    *,
    source_label: str,
    source_template: QCFTemplate,
    target_template: QCFTemplate,
    q: sp.Symbol,
) -> BauerMuirDirectObstruction:
    """Return the earliest forced obstruction for a direct 1-step Bauer-Muir match when visible."""
    source_b0, source_a, source_b = _template_reciprocal_coeffs(source_template.normalized(), q=q, depth=3)
    target_b0, target_a, target_b = _template_reciprocal_coeffs(target_template.normalized(), q=q, depth=3)

    forced_w0 = sp.simplify(target_b0 - source_b0)
    forced_w1 = sp.simplify(target_b[1] - source_b[1])
    transformed_a1 = sp.simplify(source_a[1] - forced_w0 * (source_b[1] + forced_w1))

    forced_w2 = None
    transformed_a2 = None
    obstruction_stage = None

    if sp.simplify(transformed_a1 - target_a[1]) != 0:
        obstruction_stage = 1
    elif sp.simplify(forced_w0) == 0:
        forced_w2 = sp.simplify(target_b[2] - source_b[2])
        lambda2 = sp.simplify(source_a[2] - forced_w1 * (source_b[2] + forced_w2))
        if sp.simplify(transformed_a1) != 0:
            transformed_a2 = sp.simplify(source_a[1] * lambda2 / transformed_a1)
            if sp.simplify(transformed_a2 - target_a[2]) != 0:
                obstruction_stage = 2

    return BauerMuirDirectObstruction(
        source_label=source_label,
        forced_w0=forced_w0,
        forced_w1=forced_w1,
        forced_w2=forced_w2,
        transformed_a1=transformed_a1,
        target_a1=target_a[1],
        transformed_a2=transformed_a2,
        target_a2=target_a[2],
        obstruction_stage=obstruction_stage,
    )


def _template_reciprocal_coeffs(template: QCFTemplate, q: sp.Symbol, depth: int) -> tuple[sp.Expr, list[sp.Expr], list[sp.Expr]]:
    """Return b0 and coefficient lists for 1/value(template) = b0 + K a_n / b_n (truncated)."""
    b0 = sp.Rational(template.base_denominator, template.top_constant)
    a_terms = [sp.Integer(0)]
    b_terms = [sp.Integer(0)]
    for n in range(1, depth + 1):
        position = n - 1
        numerator_exponent = template.numerator_q_shift + template.numerator_q_step * position
        a_n = template.numerator_scale * q**numerator_exponent
        if template.numerator_extra_scale != 0:
            extra_exponent = template.numerator_extra_q_shift + template.numerator_extra_q_step * position
            a_n += template.numerator_extra_scale * q**extra_exponent
        b_n = sp.Integer(template.denominator_constant)
        if template.denominator_scale != 0:
            denominator_exponent = template.denominator_q_shift + template.denominator_q_step * position
            b_n += template.denominator_scale * q**denominator_exponent
        a_terms.append(sp.simplify(a_n))
        b_terms.append(sp.simplify(b_n))
    return b0, a_terms, b_terms


def _format_expr(expr) -> str:
    return str(sp.expand(expr)).replace("**", "^")


def _format_fraction_expr(expr) -> str:
    return str(sp.cancel(expr)).replace("**", "^")


def _research_build_profile(*, smoke: bool) -> ResearchBuildProfile:
    if smoke:
        return ResearchBuildProfile(
            label="smoke",
            raw_series_cap=31,
            reduced_order_cap=24,
            euler_order_cap=36,
            fit_order_cap=24,
            factor_depth_cap=4,
            page43_max_shift=1,
            page43_stages=2,
            max_subsequence_stride=2,
            subsequence_stages=2,
            max_bauer_muir_steps=1,
            bauer_muir_depth=3,
        )
    return ResearchBuildProfile(
        label="full",
        raw_series_cap=61,
        reduced_order_cap=80,
        euler_order_cap=80,
        fit_order_cap=60,
        factor_depth_cap=8,
        page43_max_shift=3,
        page43_stages=3,
        max_subsequence_stride=4,
        subsequence_stages=3,
        max_bauer_muir_steps=3,
        bauer_muir_depth=4,
    )


def _find_candidate(records: list[CandidateRecord], candidate_id: str) -> CandidateRecord:
    for record in records:
        if record.id == candidate_id:
            return record
    raise KeyError(f"unknown candidate id: {candidate_id}")


def _residue_support(expr, q: sp.Symbol, modulus: int) -> set[int]:
    residues: set[int] = set()
    expanded = sp.expand(expr)
    if expanded == 0:
        return residues
    for term in expanded.as_ordered_terms():
        coefficient, exponent = term.as_coeff_exponent(q)
        if coefficient == 0:
            continue
        if not exponent.is_integer:
            continue
        residues.add(int(exponent) % modulus)
    return residues


def _bauer_muir_patterns() -> list[tuple[str, Callable[[int, sp.Expr], sp.Expr]]]:
    patterns: list[tuple[str, Callable[[int, sp.Expr], sp.Expr]]] = [("w_n = 0", lambda n, q: sp.Integer(0))]
    for scale in (-2, -1, 1, 2):
        patterns.append((f"w_n = {scale}*(q^n - 1)", lambda n, q, s=scale: s * (q**n - 1)))
        patterns.append((f"w_n = {scale}*(q^(2n) - 1)", lambda n, q, s=scale: s * (q ** (2 * n) - 1)))
        patterns.append((f"w_n = {scale}*(q^n - q^(2n))", lambda n, q, s=scale: s * (q**n - q ** (2 * n))))
    return patterns


def _template_reciprocal_coeffs_at_value(
    template: QCFTemplate,
    *,
    q_value: sp.Expr,
    depth: int,
) -> tuple[sp.Expr, list[sp.Expr], list[sp.Expr]]:
    b0 = sp.Rational(template.base_denominator, template.top_constant)
    a_terms = [sp.Integer(0)]
    b_terms = [sp.Integer(0)]
    for n in range(1, depth + 1):
        position = n - 1
        numerator_exponent = template.numerator_q_shift + template.numerator_q_step * position
        a_n = sp.Integer(template.numerator_scale) * q_value**numerator_exponent
        if template.numerator_extra_scale != 0:
            extra_exponent = template.numerator_extra_q_shift + template.numerator_extra_q_step * position
            a_n += sp.Integer(template.numerator_extra_scale) * q_value**extra_exponent

        b_n = sp.Integer(template.denominator_constant)
        if template.denominator_scale != 0:
            denominator_exponent = template.denominator_q_shift + template.denominator_q_step * position
            b_n += sp.Integer(template.denominator_scale) * q_value**denominator_exponent
        a_terms.append(sp.simplify(a_n))
        b_terms.append(sp.simplify(b_n))
    return b0, a_terms, b_terms


def _bauer_muir_matches_numeric(
    transformed: BauerMuirTransformResult,
    target_b0: sp.Expr,
    target_a_terms: list[sp.Expr],
    target_b_terms: list[sp.Expr],
    depth: int,
) -> bool:
    if sp.simplify(transformed.b0 - target_b0) != 0:
        return False
    for n in range(1, depth + 1):
        if sp.simplify(transformed.a_terms[n] - target_a_terms[n]) != 0:
            return False
        if sp.simplify(transformed.b_terms[n] - target_b_terms[n]) != 0:
            return False
    return True


def bauer_muir_pattern_search(
    *,
    source_label: str,
    source_template: QCFTemplate,
    target_template: QCFTemplate,
    q: sp.Symbol,
    depth: int,
    steps: int,
) -> list[dict[str, object]]:
    """Search a tiny fixed Bauer-Muir pattern family against exact rational sample points."""
    if steps < 1 or steps > 3:
        raise ValueError("steps must be between 1 and 3")

    patterns = _bauer_muir_patterns()
    sample_points = (sp.Rational(1, 10), sp.Rational(1, 7), sp.Rational(1, 5))
    normalized_source = source_template.normalized()
    normalized_target = target_template.normalized()
    sample_data = {
        sample: (
            _template_reciprocal_coeffs_at_value(normalized_source, q_value=sample, depth=depth),
            _template_reciprocal_coeffs_at_value(normalized_target, q_value=sample, depth=depth),
        )
        for sample in sample_points
    }

    hits: list[dict[str, object]] = []

    def _start_transforms() -> dict[sp.Expr, BauerMuirTransformResult]:
        initial: dict[sp.Expr, BauerMuirTransformResult] = {}
        for sample in sample_points:
            source_coeffs, _ = sample_data[sample]
            source_b0, source_a_terms, source_b_terms = source_coeffs
            initial[sample] = BauerMuirTransformResult(
                b0=source_b0,
                a_terms=source_a_terms,
                b_terms=source_b_terms,
                lambdas=[sp.Integer(0)],
            )
        return initial

    def _search(
        *,
        transforms_by_sample: dict[sp.Expr, BauerMuirTransformResult],
        chain: list[str],
        remaining_steps: int,
    ) -> None:
        for label, func in patterns:
            next_transforms: dict[sp.Expr, BauerMuirTransformResult] = {}
            chain_matches = True
            for sample in sample_points:
                _, (target_b0, target_a_terms, target_b_terms) = sample_data[sample]
                current_transform = transforms_by_sample[sample]
                next_terms = [func(n, sample) for n in range(0, depth + 1)]
                next_transform = bauer_muir_transform_trunc(
                    b0=current_transform.b0,
                    a_terms=current_transform.a_terms,
                    b_terms=current_transform.b_terms,
                    w_terms=next_terms,
                )
                next_transforms[sample] = next_transform
                if not _bauer_muir_matches_numeric(next_transform, target_b0, target_a_terms, target_b_terms, depth):
                    chain_matches = False
            next_chain = chain + [label]
            if chain_matches and len(next_chain) == steps:
                hits.append(
                    {
                        "source": source_label,
                        "steps": len(next_chain),
                        "pattern_chain": next_chain,
                        "samples": [str(sample) for sample in sample_points],
                    }
                )
            if remaining_steps > 1:
                _search(
                    transforms_by_sample=next_transforms,
                    chain=next_chain,
                    remaining_steps=remaining_steps - 1,
                )

    _search(
        transforms_by_sample=_start_transforms(),
        chain=[],
        remaining_steps=steps,
    )

    return hits


def build_candidate_research_note(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    depth: int = 40,
    series_order: int = 151,
    smoke: bool = False,
) -> None:
    records = read_candidates(input_path)
    record = _find_candidate(records, candidate_id)
    benchmark = get_benchmark(record.closest_benchmark)
    profile = _research_build_profile(smoke=smoke)

    q = sp.Symbol("q")
    step = benchmark.canonical_template.numerator_q_step

    # Computing very high-order q-series directly is expensive in Sympy. For the hero-case
    # RR(q^3) neighborhood, it's better to step-reduce to t=q^3 and work there.
    raw_order = min(series_order, profile.raw_series_cap)
    candidate_series = _series_expr(record.template, depth=depth, order=raw_order, q_symbol=q)
    benchmark_series = _series_expr(benchmark.canonical_template, depth=depth, order=raw_order, q_symbol=q)
    ratio_series = sp.expand(sp.series(candidate_series / benchmark_series, q, 0, raw_order).removeO())

    ratio_residues = _residue_support(
        ratio_series - ratio_series.subs(q, 0),
        q=q,
        modulus=step or 1,
    )

    lines: list[str] = [
        f"# Research Note: `{record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{record.id}`",
        f"- Closest benchmark: `{record.closest_benchmark}` ({record.closest_benchmark_digits} shared digits)",
        f"- Candidate template: `{record.template.signature()}`",
        f"- Benchmark template: `{benchmark.canonical_template.signature()}`",
        f"- Depth: `{depth}`",
        f"- Series order request: `{series_order}`",
        f"- Raw q-series computed order: `{raw_order}`",
        f"- Build profile: `{profile.label}`",
        "",
        "## Ratio Series",
        "",
        f"- Computed `candidate / {record.closest_benchmark}` as a truncated q-series.",
        f"- Residue support mod step `{step}` (excluding constant term): `{sorted(ratio_residues)}`",
        "",
        "```text",
        _format_expr(sp.series(ratio_series, q, 0, min(series_order, 40)).removeO()),
        "```",
    ]

    # Step-reduced t=q^step view (only when it is cleanly divisible).
    reduced_candidate = reduce_template_by_step(record.template.normalized(), step=step)
    reduced_benchmark = reduce_template_by_step(benchmark.canonical_template.normalized(), step=step)
    if reduced_candidate is not None and reduced_benchmark is not None:
        t = sp.Symbol("t")
        reduced_order = min((series_order // step) + 3, profile.reduced_order_cap)
        candidate_t = _series_expr(
            reduced_candidate,
            depth=depth,
            order=reduced_order,
            q_symbol=t,
        )
        benchmark_t = _series_expr(
            reduced_benchmark,
            depth=depth,
            order=reduced_order,
            q_symbol=t,
        )
        ratio_t = sp.expand(sp.series(candidate_t / benchmark_t, t, 0, reduced_order).removeO())
        euler_exponents = euler_product_exponents(ratio_t, q=t, order=min(reduced_order, profile.euler_order_cap))

        max_abs_exponent = max(int(abs(sp.N(value))) for value in euler_exponents[:30] if value.is_number) if euler_exponents else 0
        nonzero_count = sum(1 for value in euler_exponents[:30] if sp.simplify(value) != 0)

        lines.extend(
            [
                "",
                "## Step-Reduced View (t-Series)",
                "",
                f"- Reduced variable: `t = q^{step}`",
                f"- Reduced candidate template: `{reduced_candidate.signature()}`",
                f"- Reduced benchmark template: `{reduced_benchmark.signature()}`",
                "",
                "### Ratio(t)",
                "",
                "```text",
                _format_expr(sp.series(ratio_t, t, 0, 40).removeO()),
                "```",
                "",
                "### Euler Product Exponents",
                "",
                "- Using the truncated Euler transform: `Ratio(t) = Π_{n>=1} (1 - t^n)^{c_n}` modulo the visible order.",
                f"- Nonzero count in first 30 exponents: `{nonzero_count}`",
                f"- Max |c_n| among first 30 exponents (rough): `{max_abs_exponent}`",
                "",
                "```text",
                "c_1..c_30:",
                ", ".join(str(sp.simplify(value)) for value in euler_exponents[:30]),
                "```",
            ]
        )

        fit_window = euler_exponents[: profile.fit_order_cap]
        periodic_fit = try_fit_periodic_pochhammer(fit_window, max_period=12, max_abs=8)
        two_modulus_fit = try_fit_two_modulus_pochhammer(fit_window, max_modulus=12, max_abs=8)
        eta_fit = try_fit_eta_quotient(fit_window, max_level=12, max_abs=8)
        two_modulus_summary = (
            "no fit"
            if two_modulus_fit is None
            else (
                f"(m1={two_modulus_fit.first_modulus}, e1={two_modulus_fit.first_exponents}, "
                f"m2={two_modulus_fit.second_modulus}, e2={two_modulus_fit.second_exponents})"
            )
        )
        lines.extend(
            [
                "",
                "### Structured Fit Attempts",
                "",
                (
                    f"- Periodic Pochhammer fit (bounded |e_r|<=8, period<=12): "
                    f"`{periodic_fit if periodic_fit is not None else 'no fit'}`"
                ),
                (
                    "- Two-modulus Pochhammer fit (bounded |e_r|<=8, moduli<=12): "
                    f"`{two_modulus_summary}`"
                ),
                (
                    f"- Eta-quotient fit (bounded |e_d|<=8, level<=12): "
                    f"`{eta_fit if eta_fit is not None else 'no fit'}`"
                ),
            ]
        )

        benchmark_closed_form = _benchmark_product_closed_form_in_reduced_variable(
            benchmark_name=record.closest_benchmark,
            variable="t",
        )

        # Upgrade fits into "closed-form drafts" with explicit, checkable certificates.
        # We *solve* using the early segment, but we *certify* against the full computed range.
        fit_check_count = len(euler_exponents)
        lines.extend(
            [
                "",
                "### Closed-Form Drafts (Product Guesses)",
                "",
                "- These are *hypotheses* derived from the Euler-exponent signature of `Ratio(t)`.",
                f"- Certificate rule: verify the predicted Euler exponents `c_1..c_{fit_check_count}` exactly.",
                (
                    "- Closest-benchmark product (in reduced variable) is known: "
                    f"`{benchmark_closed_form}`"
                    if benchmark_closed_form is not None
                    else "- Closest-benchmark product (in reduced variable) is not in the small built-in mapping."
                ),
            ]
        )

        if periodic_fit is None and two_modulus_fit is None and eta_fit is None:
            lines.extend(
                [
                    "- Current outcome: no periodic-Pochhammer, two-modulus Pochhammer, or eta-quotient fit was found in the bounded search box.",
                    "",
                ]
            )
        else:
            lines.append("")

        if periodic_fit is not None:
            period, residue_exponents = periodic_fit
            poch_symbolic = _format_periodic_pochhammer_closed_form(
                period=period,
                exponents=residue_exponents,
                variable="t",
            )
            mismatches = _verify_periodic_pochhammer_fit(
                euler_exponents=euler_exponents,
                period=period,
                residue_exponents=residue_exponents,
                check_count=fit_check_count,
            )
            if not mismatches and benchmark_closed_form is not None:
                lines.extend(
                    [
                        "- Candidate closed-form draft (if the fit persists beyond the checked order):",
                        "  Candidate(t) = Benchmark(t) * Ratio(t)",
                        f"  Benchmark(t) = {benchmark_closed_form}",
                        f"  Ratio(t) = {poch_symbolic}",
                    ]
                )
            periodic_cert = {
                "object": "Ratio(t)",
                "fit_type": "periodic_pochhammer",
                "period": period,
                "residue_exponents": {str(r): residue_exponents[r - 1] for r in range(1, period + 1)},
                "checked_c_n": fit_check_count,
                "mismatch_indices": mismatches,
            }
            lines.extend(
                [
                    f"- Periodic Pochhammer closed form: `{poch_symbolic}`",
                    f"- Certificate: mismatches in c_1..c_{fit_check_count}: `{len(mismatches)}`",
                    "",
                    "```json",
                    json.dumps(periodic_cert, indent=2, ensure_ascii=True),
                    "```",
                ]
            )

        if two_modulus_fit is not None:
            two_modulus_symbolic = _format_two_modulus_pochhammer_closed_form(
                fit=two_modulus_fit,
                variable="t",
            )
            mismatches = _verify_two_modulus_pochhammer_fit(
                euler_exponents=euler_exponents,
                fit=two_modulus_fit,
                check_count=fit_check_count,
            )
            if not mismatches and benchmark_closed_form is not None:
                lines.extend(
                    [
                        "- Candidate closed-form draft (if the fit persists beyond the checked order):",
                        "  Candidate(t) = Benchmark(t) * Ratio(t)",
                        f"  Benchmark(t) = {benchmark_closed_form}",
                        f"  Ratio(t) = {two_modulus_symbolic}",
                    ]
                )
            two_modulus_cert = {
                "object": "Ratio(t)",
                "fit_type": "two_modulus_pochhammer",
                "first_modulus": two_modulus_fit.first_modulus,
                "first_residue_exponents": {
                    str(r): two_modulus_fit.first_exponents[r - 1]
                    for r in range(1, two_modulus_fit.first_modulus + 1)
                },
                "second_modulus": two_modulus_fit.second_modulus,
                "second_residue_exponents": {
                    str(r): two_modulus_fit.second_exponents[r - 1]
                    for r in range(1, two_modulus_fit.second_modulus + 1)
                },
                "checked_c_n": fit_check_count,
                "mismatch_indices": mismatches,
            }
            lines.extend(
                [
                    f"- Two-modulus Pochhammer closed form: `{two_modulus_symbolic}`",
                    f"- Certificate: mismatches in c_1..c_{fit_check_count}: `{len(mismatches)}`",
                    "",
                    "```json",
                    json.dumps(two_modulus_cert, indent=2, ensure_ascii=True),
                    "```",
                ]
            )

        if eta_fit is not None:
            level, exponents_by_divisor = eta_fit
            eta_symbolic = _format_eta_quotient_closed_form(
                level=level,
                exponents_by_divisor=exponents_by_divisor,
                variable="t",
            )
            mismatches = _verify_eta_quotient_fit(
                euler_exponents=euler_exponents,
                level=level,
                exponents_by_divisor=exponents_by_divisor,
                check_count=fit_check_count,
            )
            if not mismatches and benchmark_closed_form is not None:
                lines.extend(
                    [
                        "- Candidate closed-form draft (if the fit persists beyond the checked order):",
                        "  Candidate(t) = Benchmark(t) * Ratio(t)",
                        f"  Benchmark(t) = {benchmark_closed_form}",
                        f"  Ratio(t) = {eta_symbolic}",
                    ]
                )
            eta_cert = {
                "object": "Ratio(t)",
                "fit_type": "eta_quotient",
                "level": level,
                "exponents_by_divisor": {str(d): e for d, e in sorted(exponents_by_divisor.items())},
                "checked_c_n": fit_check_count,
                "mismatch_indices": mismatches,
            }
            lines.extend(
                [
                    f"- Eta-quotient closed form: `{eta_symbolic}`",
                    f"- Certificate: mismatches in c_1..c_{fit_check_count}: `{len(mismatches)}`",
                    "",
                    "```json",
                    json.dumps(eta_cert, indent=2, ensure_ascii=True),
                    "```",
                ]
            )

        # Heine hcf2 c=bz=-1 coefficient check when the reduced reciprocal matches 1+K (t^n+t^(2n))/(1+t^n).
        b0_target, a_target, b_target = _template_reciprocal_coeffs(reduced_candidate, q=t, depth=4)
        factor_depth = min(depth, profile.factor_depth_cap)
        factor_b0, factor_a_terms, factor_b_terms = _template_reciprocal_coeffs(
            reduced_candidate,
            q=t,
            depth=factor_depth,
        )
        factor_witness = convergent_factor_equivalence_witness(
            b0=factor_b0,
            a_terms=factor_a_terms,
            b_terms=factor_b_terms,
        )
        has_nontrivial_factor = any(
            sp.simplify(common_factor - 1) != 0 for common_factor in factor_witness.reduction.gcd_factors[1:]
        )
        if has_nontrivial_factor:
            display_stages = min(4, len(factor_witness.reduction.reduced_coeffs.a_terms) - 1)
            factor_lines = []
            for n in range(1, min(5, len(factor_witness.reduction.gcd_factors))):
                factor_lines.append(f"g{n} = {_format_expr(factor_witness.reduction.gcd_factors[n])}")

            reduced_lines = [f"b0_red = {_format_expr(factor_witness.reduction.reduced_coeffs.b0)}"]
            for n in range(1, display_stages + 1):
                reduced_lines.append(
                    f"a{n}_red = {_format_expr(factor_witness.reduction.reduced_coeffs.a_terms[n])}"
                )
                reduced_lines.append(
                    f"b{n}_red = {_format_expr(factor_witness.reduction.reduced_coeffs.b_terms[n])}"
                )

            scale_lines = []
            for n in range(1, min(5, len(factor_witness.scale_terms))):
                scale_lines.append(f"r{n} = {_format_fraction_expr(factor_witness.scale_terms[n])}")

            retransformed_matches = all(
                sp.simplify(factor_witness.retransformed_coeffs.a_terms[n] - factor_a_terms[n]) == 0
                and sp.simplify(factor_witness.retransformed_coeffs.b_terms[n] - factor_b_terms[n]) == 0
                for n in range(1, len(factor_a_terms))
            )
            lines.extend(
                [
                    "",
                    "## Exact Convergent-Factor Reduction",
                    "",
                    (
                        f"- Checked exact convergent gcd factors through stage `{factor_depth}`. "
                        "The first visible factors are:"
                    ),
                    "",
                    "```text",
                    *factor_lines,
                    "```",
                    "",
                    "- After cancellation, the induced reduced continued fraction begins:",
                    "",
                    "```text",
                    *reduced_lines,
                    "```",
                    "",
                    "- The original reduced target is recovered by the reverse equivalence transform with stage scales:",
                    "",
                    "```text",
                    *scale_lines,
                    "```",
                    "",
                    (
                        f"- Applying those scales back to the cancelled fraction reproduces the target "
                        f"coefficients exactly through the checked depth: `{retransformed_matches}`."
                    ),
                    "- These reverse scales are rational functions in `t`, so they point toward a future fraction-field formalization layer rather than a purely polynomial one.",
                ]
            )

        rr_direct = direct_bauer_muir_obstruction(
            source_label="RR reciprocal",
            source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
        )
        cubic_direct = direct_bauer_muir_obstruction(
            source_label="cubic reciprocal",
            source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
        )
        lines.extend(
            [
                "",
                "## Direct 1-Step Bauer-Muir Obstruction",
                "",
                (
                    f"- `{rr_direct.source_label}`: matching `b0` forces `w0 = {_format_expr(rr_direct.forced_w0)}` "
                    f"and matching `b1` forces `w1 = {_format_expr(rr_direct.forced_w1)}`. "
                    f"That leaves the first transformed numerator at "
                    f"`{_format_expr(rr_direct.transformed_a1)}` instead of the target "
                    f"`{_format_expr(rr_direct.target_a1)}`, so no direct 1-step transform exists."
                ),
                (
                    f"- `{cubic_direct.source_label}`: the forced first-stage data "
                    f"`w0 = {_format_expr(cubic_direct.forced_w0)}`, "
                    f"`w1 = {_format_expr(cubic_direct.forced_w1)}` does preserve the first numerator, "
                    f"but matching `b2` also forces `w2 = {_format_expr(cubic_direct.forced_w2)}`. "
                    f"Then the second transformed numerator becomes "
                    f"`{_format_expr(cubic_direct.transformed_a2)}` instead of "
                    f"`{_format_expr(cubic_direct.target_a2)}`, so this direct path is also ruled out."
                ),
            ]
        )

        looks_like_hybrid = (
            sp.simplify(b0_target - 1) == 0
            and sp.simplify(b_target[1] - (1 + t)) == 0
            and sp.simplify(a_target[1] - (t + t**2)) == 0
        )
        if looks_like_hybrid:
            a_sym, b_sym = sp.symbols("a b")
            f2_equivalence_obstruction = page43_f2_zero_shift_equivalence_obstruction(q=t)
            f4_equivalence_obstruction = page43_f4_zero_shift_equivalence_obstruction(q=t)
            coeffs = heine_hcf2_standardized_coeffs(
                a=a_sym,
                b=b_sym,
                c=sp.Integer(-1),
                z=-1 / b_sym,
                q=t,
                depth=4,
            )
            lines.extend(
                [
                    "",
                    "## Heine hcf2 Specialization Check (c=bz=-1)",
                    "",
                    "- In hcf2, setting `c = b z = -1` forces denominators `b_n = 1 + t^n` and `b0 = 1`.",
                    "- Under the standardized coefficient extraction used here, the first numerator term becomes a constant in `t`:",
                    "",
                    "```text",
                    f"a1_hcf2 = {_format_expr(sp.simplify(coeffs.a_terms[1]))}",
                    "```",
                    "",
                    "- But the target reciprocal for this candidate has `a1 = t + t^2`, so this specialization cannot match at the coefficient level.",
                ]
            )

            b_cor, lam_cor = sp.symbols("b lambda")
            cor2cf_obstruction = heine_cor2cf_a_zero_contraction_obstruction(
                b=b_cor,
                lam=lam_cor,
                q=t,
                depth=12,
            )
            lines.extend(
                [
                    "",
                    "## Heine `cor2cf` Contraction Check (`a = 0` lane)",
                    "",
                    "- For the nearby `cor2cf` family, matching the target initial term `1` forces the relevant branch to the `a = 0` specialization:",
                    "",
                    "```text",
                    f"b0 = {_format_expr(cor2cf_obstruction.source_coeffs.b0)}",
                    f"a1 = {_format_expr(cor2cf_obstruction.source_coeffs.a_terms[1])}",
                    f"b1 = {_format_expr(cor2cf_obstruction.source_coeffs.b_terms[1])}",
                    f"a2 = {_format_expr(cor2cf_obstruction.source_coeffs.a_terms[2])}",
                    f"b2 = {_format_expr(cor2cf_obstruction.source_coeffs.b_terms[2])}",
                    f"a3 = {_format_expr(cor2cf_obstruction.source_coeffs.a_terms[3])}",
                    f"b3 = {_format_expr(cor2cf_obstruction.source_coeffs.b_terms[3])}",
                    f"a4 = {_format_expr(cor2cf_obstruction.source_coeffs.a_terms[4])}",
                    f"b4 = {_format_expr(cor2cf_obstruction.source_coeffs.b_terms[4])}",
                    "```",
                    "",
                    (
                        f"- Odd part: the initial term is `d0 = {_format_expr(cor2cf_obstruction.odd_part.b0)}` "
                        f"instead of `{_format_expr(b0_target)}`."
                    ),
                    (
                        f"- Even part: the initial term stays `{_format_expr(cor2cf_obstruction.even_part.b0)}`, "
                        f"but the first numerator is `{_format_expr(cor2cf_obstruction.even_part.a_terms[1])}` "
                        f"instead of `{_format_expr(a_target[1])}`."
                    ),
                    (
                        "- Odd-of-even branch: the new initial term is "
                        f"`{_format_fraction_expr(cor2cf_obstruction.even_odd_part.b0)}`, "
                        f"so this two-step branch also fails at stage 0."
                    ),
                    (
                        "- Even-of-even branch: the initial term stays `1`, but the first numerator is "
                        f"`{_format_expr(cor2cf_obstruction.even_even_part.a_terms[1])}`; "
                        "its `t^2` coefficient is `0`, so it cannot equal the target `t + t^2`."
                    ),
                    "- So the relevant 1-step and 2-step odd/even contraction branches around `cor2cf` are ruled out exactly at low stage.",
                ]
            )

            f2_monomial_hits = page43_monomial_parameter_search(
                family="f2",
                target_template=reduced_candidate,
                q=t,
                max_shift=profile.page43_max_shift,
                stages=profile.page43_stages,
            )
            f4_monomial_hits = page43_monomial_parameter_search(
                family="f4",
                target_template=reduced_candidate,
                q=t,
                max_shift=profile.page43_max_shift,
                stages=profile.page43_stages,
            )
            lines.extend(
                [
                    "",
                    "## Page-43 Monomial Substitution Check",
                    "",
                    (
                        "- Search shape: `a = alpha*t^A`, `b = beta*t^B`, `lambda = gamma*t^L` "
                        f"with integer shifts `A,B,L in [-{profile.page43_max_shift},{profile.page43_max_shift}]`."
                    ),
                    (
                        "- Matching rule: solve exactly for `alpha, beta, gamma` so the first "
                        f"`{profile.page43_stages}` reciprocal stages match the reduced target."
                    ),
                    f"- `f2` / `gcf3` hits in this box: `{len(f2_monomial_hits)}`",
                    f"- `f4` / `gcf2` hits in this box: `{len(f4_monomial_hits)}`",
                ]
            )
            monomial_hits = f2_monomial_hits + f4_monomial_hits
            if monomial_hits:
                lines.extend(["", "```text"])
                for hit in monomial_hits:
                    lines.append(
                        f"{hit.family}: A={hit.a_shift}, B={hit.b_shift}, L={hit.lambda_shift}, "
                        f"alpha={_format_expr(hit.a_coeff)}, beta={_format_expr(hit.b_coeff)}, "
                        f"lambda={_format_expr(hit.lambda_coeff)}"
                    )
                lines.extend(["```"])

            rational_page43_shift = 0
            rational_page43_stages = min(profile.page43_stages, 2)
            f2_rational_hits = page43_rational_parameter_search(
                family="f2",
                target_template=reduced_candidate,
                q=t,
                max_shift=rational_page43_shift,
                stages=rational_page43_stages,
            )
            f4_rational_hits = page43_rational_parameter_search(
                family="f4",
                target_template=reduced_candidate,
                q=t,
                max_shift=rational_page43_shift,
                stages=rational_page43_stages,
            )
            lines.extend(
                [
                    "",
                    "## Page-43 Low-Complexity Rational Prefactor Check",
                    "",
                    (
                        "- Search shape: `a = alpha*phi_a(t)*t^A`, `b = beta*phi_b(t)*t^B`, "
                        "`lambda = gamma*phi_lambda(t)*t^L`."
                    ),
                    (
                        "- Prefactor box: `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active "
                        f"and integer shifts fixed to `A=B=L=0`."
                    ),
                    (
                        "- Matching rule: solve exactly for scalar `alpha, beta, gamma` so the first "
                        f"`{rational_page43_stages}` reciprocal stages match the reduced target."
                    ),
                    f"- `f2` / `gcf3` hits in this prefactor box: `{len(f2_rational_hits)}`",
                    f"- `f4` / `gcf2` hits in this prefactor box: `{len(f4_rational_hits)}`",
                ]
            )
            rational_hits = f2_rational_hits + f4_rational_hits
            if rational_hits:
                lines.extend(["", "```text"])
                for hit in rational_hits:
                    lines.append(
                        f"{hit.family}: A={hit.a_shift}, B={hit.b_shift}, L={hit.lambda_shift}, "
                        f"phi_a={hit.a_profile}, phi_b={hit.b_profile}, phi_lambda={hit.lambda_profile}, "
                        f"alpha={_format_expr(hit.a_coeff)}, beta={_format_expr(hit.b_coeff)}, "
                        f"lambda={_format_expr(hit.lambda_coeff)}"
                    )
                lines.extend(["```"])

            lines.extend(
                [
                    "",
                    "## Exact `f2` / `gcf3` n-Dependent Equivalence Check",
                    "",
                    "- Prioritized source-family-specific lane: the zero-shift `f2` / `gcf3` family under an arbitrary `n`-dependent equivalence transformation.",
                    "- Write `m = t^(n-1)` and enforce the necessary identity for matching the hero-case reciprocal:",
                    "",
                    "```text",
                    "alpha_n * (1 + t^(n-1)) = t^n * beta_(n-1) * beta_n",
                    "```",
                    "",
                    "- In the zero-shift `f2` lane this becomes:",
                    "",
                    "```text",
                    "alpha_n = lambda*m*t - a*b*m^2*t^2",
                    "beta_(n-1) = 1 + b*m + a*m*t",
                    "beta_n = 1 + b*m*t + a*m*t^2",
                    "```",
                    "",
                    "- Residual polynomial:",
                    "",
                    "```text",
                    _format_expr(f2_equivalence_obstruction.residual_polynomial),
                    "```",
                    "",
                    (
                        "- `m^3` coefficient: "
                        f"`{_format_expr(f2_equivalence_obstruction.m_coefficients[3])}`; "
                        "exact vanishing forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient is "
                        f"`{_format_expr(f2_equivalence_obstruction.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes "
                        f"`{_format_expr(f2_equivalence_obstruction.final_m2_coefficient)}`, "
                        "still nonzero."
                    ),
                    "- So the current hero case is not in this zero-shift `f2` / `gcf3` equivalence lane.",
                ]
            )
            lines.extend(
                [
                    "",
                    "## Exact `f4` / `gcf2` n-Dependent Equivalence Check",
                    "",
                    "- Next source-family-specific lane: the zero-shift `f4` / `gcf2` family under an arbitrary `n`-dependent equivalence transformation.",
                    "- Write `m = t^(n-1)` and enforce the same necessary identity:",
                    "",
                    "```text",
                    "alpha_n * (1 + t^(n-1)) = t^n * beta_(n-1) * beta_n",
                    "```",
                    "",
                    "- In the zero-shift `f4` lane this becomes:",
                    "",
                    "```text",
                    "alpha_n = a*t + lambda*m*t",
                    "beta_(n-1) = 1 - a*t + b*m",
                    "beta_n = 1 - a*t + b*m*t",
                    "```",
                    "",
                    "- Residual polynomial:",
                    "",
                    "```text",
                    _format_expr(f4_equivalence_obstruction.residual_polynomial),
                    "```",
                    "",
                    (
                        "- `m^0` coefficient: "
                        f"`{_format_expr(f4_equivalence_obstruction.m_coefficients[0])}`; "
                        "exact vanishing forces `a = 0`."
                    ),
                    (
                        "- After `a = 0`, the `m^3` coefficient is "
                        f"`{_format_expr(sp.simplify(f4_equivalence_obstruction.m_coefficients[3].subs(f4_equivalence_obstruction.forced_a_solution)))}`; "
                        "exact vanishing forces `b = 0`."
                    ),
                    (
                        "- After that, the `m^1` coefficient is "
                        f"`{_format_expr(f4_equivalence_obstruction.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes "
                        f"`{_format_expr(f4_equivalence_obstruction.final_m2_coefficient)}`, "
                        "still nonzero."
                    ),
                    "- So the current hero case is not in this zero-shift `f4` / `gcf2` equivalence lane either.",
                ]
            )

            cubic_b0, cubic_a_terms, cubic_b_terms = _template_reciprocal_coeffs(
                get_benchmark("ramanujan_cubic_normalized").canonical_template,
                q=t,
                depth=6,
            )
            cubic_odd = parity_contraction_coeffs(
                b0=cubic_b0,
                a_terms=cubic_a_terms,
                b_terms=cubic_b_terms,
                parity="odd",
            )
            cubic_even = parity_contraction_coeffs(
                b0=cubic_b0,
                a_terms=cubic_a_terms,
                b_terms=cubic_b_terms,
                parity="even",
            )
            lines.extend(
                [
                    "",
                    "## Cubic Odd/Even Contraction Check",
                    "",
                    "- Treat the reduced cubic benchmark reciprocal as `1 + K (t^n + t^(2n)) / 1`.",
                    (
                        f"- Odd part: the initial term is `d0 = {_format_expr(cubic_odd.b0)}`, "
                        f"already incompatible with the target initial term `{_format_expr(b0_target)}`."
                    ),
                    (
                        f"- Even part: the initial term `{_format_expr(cubic_even.b0)}` and first numerator "
                        f"`{_format_expr(cubic_even.a_terms[1])}` do match the target, "
                        f"but the first denominator is `{_format_expr(cubic_even.b_terms[1])}` "
                        f"instead of `{_format_expr(b_target[1])}`."
                    ),
                    "- So the reduced candidate is not a simple odd/even canonical contraction of the cubic reciprocal.",
                ]
            )

            rr_subsequence_hits = arithmetic_subsequence_contraction_search(
                source_label="RR reciprocal",
                source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
                target_template=reduced_candidate,
                q=t,
                max_stride=profile.max_subsequence_stride,
                stages=profile.subsequence_stages,
            )
            cubic_subsequence_hits = arithmetic_subsequence_contraction_search(
                source_label="cubic reciprocal",
                source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
                target_template=reduced_candidate,
                q=t,
                max_stride=profile.max_subsequence_stride,
                stages=profile.subsequence_stages,
            )
            stride_values = ", ".join(str(value) for value in range(2, profile.max_subsequence_stride + 1))
            lines.extend(
                [
                    "",
                    "## Arithmetic Subsequence Contraction Scan",
                    "",
                    (
                        "- Search shape: every `stride`-th convergent subsequence with "
                        f"`stride in {{{stride_values}}}` and all offsets."
                    ),
                    (
                        "- Matching rule: recover the induced contracted fraction and compare "
                        f"`b0, a1..a{profile.subsequence_stages}, b1..b{profile.subsequence_stages}` exactly against the reduced target."
                    ),
                    f"- RR source hits in this box: `{len(rr_subsequence_hits)}`",
                    f"- Cubic source hits in this box: `{len(cubic_subsequence_hits)}`",
                ]
            )
            subsequence_hits = rr_subsequence_hits + cubic_subsequence_hits
            if subsequence_hits:
                lines.extend(["", "```text"])
                for hit in subsequence_hits:
                    lines.append(f"{hit.source_label}: stride={hit.stride}, offset={hit.offset}")
                lines.extend(["```"])

        bm_pattern_count = len(_bauer_muir_patterns())
        rr_hits_by_step: dict[int, list[dict[str, object]]] = {}
        cubic_hits_by_step: dict[int, list[dict[str, object]]] = {}
        for step_count in range(1, profile.max_bauer_muir_steps + 1):
            rr_hits_by_step[step_count] = bauer_muir_pattern_search(
                source_label="RR reciprocal",
                source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
                target_template=reduced_candidate,
                q=t,
                depth=profile.bauer_muir_depth,
                steps=step_count,
            )
            cubic_hits_by_step[step_count] = bauer_muir_pattern_search(
                source_label="cubic reciprocal",
                source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
                target_template=reduced_candidate,
                q=t,
                depth=profile.bauer_muir_depth,
                steps=step_count,
            )

        def _bm_line(source_label: str, hits_by_step: dict[int, list[dict[str, object]]], step_count: int) -> str:
            if step_count not in hits_by_step:
                return (
                    f"- {source_label}, {step_count}-step search space: "
                    f"`skipped in {profile.label} profile`"
                )
            return (
                f"- {source_label}, {step_count}-step search space: "
                f"`{bm_pattern_count**step_count}`; hits: `{len(hits_by_step[step_count])}`"
            )

        lines.extend(
            [
                "",
                "## Constrained Bauer-Muir Search",
                "",
                (
                    "- Match target: reciprocal coefficients through depth "
                    f"`{profile.bauer_muir_depth}`, checked at exact rational sample points `t = 1/10, 1/7, 1/5`."
                ),
                (
                    f"- Pattern family per step: `{bm_pattern_count}` low-complexity modifiers "
                    "(`0`, `±(t^n-1)`, `±(t^(2n)-1)`, `±(t^n-t^(2n))`, with scales `1` or `2`)."
                ),
                _bm_line("RR source", rr_hits_by_step, 1),
                _bm_line("RR source", rr_hits_by_step, 2),
                _bm_line("RR source", rr_hits_by_step, 3),
                _bm_line("Cubic source", cubic_hits_by_step, 1),
                _bm_line("Cubic source", cubic_hits_by_step, 2),
                _bm_line("Cubic source", cubic_hits_by_step, 3),
            ]
        )
        bm_hits: list[dict[str, object]] = []
        for hits_by_step in (rr_hits_by_step, cubic_hits_by_step):
            for step_count in sorted(hits_by_step):
                bm_hits.extend(hits_by_step[step_count])
        if bm_hits:
            lines.extend(["", "```text"])
            for hit in bm_hits:
                lines.append(str(hit))
            lines.extend(["```"])

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
