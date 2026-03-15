import sympy as sp

from ramanujan_discovery.analysis import _series_expr
from ramanujan_discovery.benchmarks import CUBIC_TEMPLATE, RR_TEMPLATE
from ramanujan_discovery.series import continued_fraction_series_coeffs, series_to_sympy


def test_fast_series_matches_sympy_series_for_rr():
    q = sp.Symbol("q")
    order = 24
    depth = 14
    fast = continued_fraction_series_coeffs(RR_TEMPLATE, depth=depth, order=order)
    fast_expr = series_to_sympy(fast, q=q)
    sympy_expr = _series_expr(RR_TEMPLATE, depth=depth, order=order, q_symbol=q)
    assert sp.simplify(fast_expr - sympy_expr) == 0


def test_fast_series_matches_sympy_series_for_cubic():
    q = sp.Symbol("q")
    order = 24
    depth = 14
    fast = continued_fraction_series_coeffs(CUBIC_TEMPLATE, depth=depth, order=order)
    fast_expr = series_to_sympy(fast, q=q)
    sympy_expr = _series_expr(CUBIC_TEMPLATE, depth=depth, order=order, q_symbol=q)
    assert sp.simplify(fast_expr - sympy_expr) == 0

