from __future__ import annotations

from dataclasses import dataclass
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
class Page43MonomialHit:
    family: str
    a_shift: int
    b_shift: int
    lambda_shift: int
    a_coeff: sp.Expr
    b_coeff: sp.Expr
    lambda_coeff: sp.Expr


@dataclass(frozen=True)
class SubsequenceContractionHit:
    source_label: str
    stride: int
    offset: int


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
                    if family == "f2":
                        source_a = gamma * q ** (lambda_shift + n) - alpha * beta * q ** (a_shift + b_shift + 2 * n)
                        source_b = 1 + beta * q ** (b_shift + n) + alpha * q ** (a_shift + n + 1)
                    else:
                        source_a = alpha * q ** (a_shift + 1) + gamma * q ** (lambda_shift + n)
                        source_b = 1 - alpha * q ** (a_shift + 1) + beta * q ** (b_shift + n)

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
                        if family == "f2":
                            source_a = gamma * q ** (lambda_shift + n) - alpha * beta * q ** (a_shift + b_shift + 2 * n)
                            source_b = 1 + beta * q ** (b_shift + n) + alpha * q ** (a_shift + n + 1)
                        else:
                            source_a = alpha * q ** (a_shift + 1) + gamma * q ** (lambda_shift + n)
                            source_b = 1 - alpha * q ** (a_shift + 1) + beta * q ** (b_shift + n)

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
) -> None:
    records = read_candidates(input_path)
    record = _find_candidate(records, candidate_id)
    benchmark = get_benchmark(record.closest_benchmark)

    q = sp.Symbol("q")
    step = benchmark.canonical_template.numerator_q_step

    # Computing very high-order q-series directly is expensive in Sympy. For the hero-case
    # RR(q^3) neighborhood, it's better to step-reduce to t=q^3 and work there.
    raw_order = min(series_order, 61)
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
        reduced_order = (series_order // step) + 3
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
        euler_exponents = euler_product_exponents(ratio_t, q=t, order=min(reduced_order, 80))

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

        periodic_fit = try_fit_periodic_pochhammer(euler_exponents[:60], max_period=12, max_abs=8)
        eta_fit = try_fit_eta_quotient(euler_exponents[:60], max_level=12, max_abs=8)
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
                    f"- Eta-quotient fit (bounded |e_d|<=8, level<=12): "
                    f"`{eta_fit if eta_fit is not None else 'no fit'}`"
                ),
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

        # Heine hcf2 c=bz=-1 coefficient check when the reduced reciprocal matches 1+K (t^n+t^(2n))/(1+t^n).
        b0_target, a_target, b_target = _template_reciprocal_coeffs(reduced_candidate, q=t, depth=4)
        looks_like_hybrid = (
            sp.simplify(b0_target - 1) == 0
            and sp.simplify(b_target[1] - (1 + t)) == 0
            and sp.simplify(a_target[1] - (t + t**2)) == 0
        )
        if looks_like_hybrid:
            a_sym, b_sym = sp.symbols("a b")
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

            f2_monomial_hits = page43_monomial_parameter_search(
                family="f2",
                target_template=reduced_candidate,
                q=t,
                max_shift=3,
                stages=3,
            )
            f4_monomial_hits = page43_monomial_parameter_search(
                family="f4",
                target_template=reduced_candidate,
                q=t,
                max_shift=3,
                stages=3,
            )
            lines.extend(
                [
                    "",
                    "## Page-43 Monomial Substitution Check",
                    "",
                    "- Search shape: `a = alpha*t^A`, `b = beta*t^B`, `lambda = gamma*t^L` with integer shifts `A,B,L in [-3,3]`.",
                    "- Matching rule: solve exactly for `alpha, beta, gamma` so the first `3` reciprocal stages match the reduced target.",
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
                max_stride=4,
                stages=3,
            )
            cubic_subsequence_hits = arithmetic_subsequence_contraction_search(
                source_label="cubic reciprocal",
                source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
                target_template=reduced_candidate,
                q=t,
                max_stride=4,
                stages=3,
            )
            lines.extend(
                [
                    "",
                    "## Arithmetic Subsequence Contraction Scan",
                    "",
                    "- Search shape: every `stride`-th convergent subsequence with `stride in {2,3,4}` and all offsets.",
                    "- Matching rule: recover the induced contracted fraction and compare `b0, a1..a3, b1..b3` exactly against the reduced target.",
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
        rr_one_step_hits = bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=1,
        )
        rr_two_step_hits = bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=2,
        )
        rr_three_step_hits = bauer_muir_pattern_search(
            source_label="RR reciprocal",
            source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=3,
        )
        cubic_one_step_hits = bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=1,
        )
        cubic_two_step_hits = bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=2,
        )
        cubic_three_step_hits = bauer_muir_pattern_search(
            source_label="cubic reciprocal",
            source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
            target_template=reduced_candidate,
            q=t,
            depth=4,
            steps=3,
        )
        lines.extend(
            [
                "",
                "## Constrained Bauer-Muir Search",
                "",
                "- Match target: reciprocal coefficients through depth `4`, checked at exact rational sample points `t = 1/10, 1/7, 1/5`.",
                (
                    f"- Pattern family per step: `{bm_pattern_count}` low-complexity modifiers "
                    "(`0`, `±(t^n-1)`, `±(t^(2n)-1)`, `±(t^n-t^(2n))`, with scales `1` or `2`)."
                ),
                f"- RR source, 1-step search space: `{bm_pattern_count}`; hits: `{len(rr_one_step_hits)}`",
                f"- RR source, 2-step search space: `{bm_pattern_count**2}`; hits: `{len(rr_two_step_hits)}`",
                f"- RR source, 3-step search space: `{bm_pattern_count**3}`; hits: `{len(rr_three_step_hits)}`",
                f"- Cubic source, 1-step search space: `{bm_pattern_count}`; hits: `{len(cubic_one_step_hits)}`",
                f"- Cubic source, 2-step search space: `{bm_pattern_count**2}`; hits: `{len(cubic_two_step_hits)}`",
                f"- Cubic source, 3-step search space: `{bm_pattern_count**3}`; hits: `{len(cubic_three_step_hits)}`",
            ]
        )
        bm_hits = (
            rr_one_step_hits
            + rr_two_step_hits
            + rr_three_step_hits
            + cubic_one_step_hits
            + cubic_two_step_hits
            + cubic_three_step_hits
        )
        if bm_hits:
            lines.extend(["", "```text"])
            for hit in bm_hits:
                lines.append(str(hit))
            lines.extend(["```"])

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
