from __future__ import annotations

from typing import Iterable

import sympy as sp

from ramanujan_discovery.models import QCFTemplate


Series = list[sp.Expr]  # coefficients [c_0, c_1, ..., c_{order-1}] for sum c_n q^n


def _zero_series(order: int) -> Series:
    return [sp.Integer(0) for _ in range(order)]


def series_add(a: Series, b: Series) -> Series:
    if len(a) != len(b):
        raise ValueError("series must have the same length")
    return [sp.simplify(x + y) for x, y in zip(a, b)]


def series_add_shifted_scaled(target: Series, src: Series, *, shift: int, scale: sp.Expr) -> None:
    """Add scale * q^shift * src into target (in place), truncated to target length."""
    if shift < 0:
        raise ValueError("shift must be non-negative")
    order = len(target)
    if len(src) != order:
        raise ValueError("series must have the same length")
    if shift >= order or scale == 0:
        return
    for idx in range(order - shift):
        if src[idx] == 0:
            continue
        target[idx + shift] = sp.simplify(target[idx + shift] + scale * src[idx])


def series_mul(a: Series, b: Series) -> Series:
    if len(a) != len(b):
        raise ValueError("series must have the same length")
    order = len(a)
    out = _zero_series(order)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        max_j = order - i
        for j in range(max_j):
            bj = b[j]
            if bj == 0:
                continue
            out[i + j] = sp.simplify(out[i + j] + ai * bj)
    return out


def series_pow(base: Series, exponent: int) -> Series:
    if exponent < 0:
        raise ValueError("exponent must be non-negative")
    order = len(base)
    out = _zero_series(order)
    out[0] = sp.Integer(1)
    if exponent == 0:
        return out
    acc = base
    exp = exponent
    while exp > 0:
        if exp & 1:
            out = series_mul(out, acc)
        exp >>= 1
        if exp:
            acc = series_mul(acc, acc)
    return out


def series_invert(den: Series) -> Series:
    """Return inv so that den * inv = 1 modulo q^order, assuming den[0] != 0."""
    order = len(den)
    if order < 1:
        raise ValueError("series must be non-empty")
    d0 = sp.simplify(den[0])
    if d0 == 0:
        raise ValueError("series constant term must be nonzero for inversion")

    inv = _zero_series(order)
    inv0 = sp.simplify(sp.Integer(1) / d0)
    inv[0] = inv0
    for n in range(1, order):
        acc = sp.Integer(0)
        for k in range(1, n + 1):
            dk = den[k]
            if dk == 0:
                continue
            acc = sp.simplify(acc + dk * inv[n - k])
        inv[n] = sp.simplify(-inv0 * acc)
    return inv


def series_div(num: Series, den: Series) -> Series:
    if len(num) != len(den):
        raise ValueError("series must have the same length")
    return series_mul(num, series_invert(den))


def series_to_sympy(series: Series, q: sp.Symbol) -> sp.Expr:
    return sp.expand(sum(coeff * q**idx for idx, coeff in enumerate(series) if coeff != 0))


def _template_exponent_parts(template: QCFTemplate) -> Iterable[int]:
    yield template.numerator_q_shift
    yield template.numerator_q_step
    if template.numerator_extra_scale != 0:
        yield template.numerator_extra_q_shift
        yield template.numerator_extra_q_step
    if template.denominator_scale != 0:
        yield template.denominator_q_shift
        yield template.denominator_q_step


def continued_fraction_series_coeffs(template: QCFTemplate, *, depth: int, order: int) -> Series:
    """Compute q-series coefficients for the truncated template continued fraction.

    This mirrors ramanujan_discovery.analysis._series_expr, but avoids Sympy's series()
    inside the depth loop for speed. Coefficients are exact Sympy Integers/Rationals.
    """
    if depth < 1:
        raise ValueError("depth must be at least 1")
    if order < 2:
        raise ValueError("order must be at least 2")
    if any(value < 0 for value in _template_exponent_parts(template)):
        raise ValueError("negative exponents are not supported in series evaluation")

    tail = _zero_series(order)
    for index in range(depth, 0, -1):
        position = index - 1
        numerator_exponent = template.numerator_q_shift + template.numerator_q_step * position
        extra_exponent = None
        if template.numerator_extra_scale != 0:
            extra_exponent = template.numerator_extra_q_shift + template.numerator_extra_q_step * position

        den_series = _zero_series(order)
        den_series[0] = sp.Integer(template.denominator_constant)
        if template.denominator_scale != 0:
            denom_exponent = template.denominator_q_shift + template.denominator_q_step * position
            if 0 <= denom_exponent < order:
                den_series[denom_exponent] = sp.simplify(den_series[denom_exponent] + template.denominator_scale)
        # den_series := den + tail
        for idx, coeff in enumerate(tail):
            if coeff == 0:
                continue
            den_series[idx] = sp.simplify(den_series[idx] + coeff)

        inv = series_invert(den_series)
        next_tail = _zero_series(order)
        series_add_shifted_scaled(
            next_tail,
            inv,
            shift=numerator_exponent,
            scale=sp.Integer(template.numerator_scale),
        )
        if extra_exponent is not None:
            series_add_shifted_scaled(
                next_tail,
                inv,
                shift=extra_exponent,
                scale=sp.Integer(template.numerator_extra_scale),
            )
        tail = next_tail

    final_den = _zero_series(order)
    final_den[0] = sp.Integer(template.base_denominator)
    for idx, coeff in enumerate(tail):
        if coeff == 0:
            continue
        final_den[idx] = sp.simplify(final_den[idx] + coeff)

    inv_final = series_invert(final_den)
    return [sp.simplify(sp.Integer(template.top_constant) * coeff) for coeff in inv_final]

