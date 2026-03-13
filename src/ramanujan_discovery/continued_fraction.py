from __future__ import annotations

from mpmath import mp

from ramanujan_discovery.models import QCFTemplate


def _to_mpf(value):
    if isinstance(value, mp.mpf):
        return value
    if isinstance(value, float):
        return mp.mpf(repr(value))
    return mp.mpf(value)


def evaluate_qcf(template: QCFTemplate, q: float, depth: int, precision: int):
    """Evaluate a truncated q-continued fraction using bottom-up recursion."""
    if depth < 1:
        raise ValueError("depth must be at least 1")

    with mp.workdps(precision):
        q_value = _to_mpf(q)
        numerators = []
        denominators = []

        for index in range(depth):
            exponent_numerator = template.numerator_q_shift + template.numerator_q_step * index
            numerator = template.numerator_scale * (q_value ** exponent_numerator)
            if template.numerator_extra_scale != 0:
                exponent_extra = template.numerator_extra_q_shift + template.numerator_extra_q_step * index
                numerator += template.numerator_extra_scale * (q_value ** exponent_extra)
            exponent_denominator = template.denominator_q_shift + template.denominator_q_step * index
            numerators.append(numerator)
            denominators.append(
                template.denominator_constant
                + template.denominator_scale * (q_value ** exponent_denominator)
            )

        tail = numerators[-1] / denominators[-1]
        for index in range(depth - 2, -1, -1):
            tail = numerators[index] / (denominators[index] + tail)

        return template.top_constant / (template.base_denominator + tail)


def q_pochhammer(a, q, precision: int, max_terms: int | None = None):
    """Compute (a; q)_inf using a truncated product that is stable for |q| < 1."""
    if max_terms is None:
        max_terms = max(precision * 2, 80)

    with mp.workdps(precision):
        a_value = _to_mpf(a)
        q_value = _to_mpf(q)
        threshold = mp.power(10, -(precision - 8))
        term = a_value
        product = mp.mpf(1)

        for _ in range(max_terms):
            product *= 1 - term
            term *= q_value
            if abs(term) < threshold:
                break

        return product


def agreement_digits(lhs, rhs, cap: int | None = None) -> int:
    with mp.workdps(50):
        error = abs(lhs - rhs)
        if error == 0:
            return cap if cap is not None else 999

        digits = int(mp.floor(-mp.log10(error)))
        digits = max(0, digits)
        if cap is None:
            return digits
        return min(digits, cap)


def format_mpf(value, digits: int = 24) -> str:
    return mp.nstr(value, n=digits)
