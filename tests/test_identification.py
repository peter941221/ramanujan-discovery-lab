import sympy as sp

from ramanujan_discovery.identification import search_polynomial_relation
from ramanujan_discovery.series import series_mul, series_pow


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
