from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import comb, gcd
from pathlib import Path
from time import perf_counter

import sympy as sp

from ramanujan_discovery.analysis import _format_expr
from ramanujan_discovery.benchmarks import get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.research import (
    ContinuedFractionCoeffs,
    _template_reciprocal_coeffs,
    convergent_common_factor_reduction,
    convergent_factor_equivalence_witness,
    reduce_template_by_step,
)
from ramanujan_discovery.series import (
    Series,
    continued_fraction_series_coeffs,
    series_add,
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
class SelfPolynomialUniquenessHit:
    """A theorem-facing self-polynomial candidate P(T, F, G) = 0."""

    modulus: int
    max_fg_total_degree: int
    max_t_degree: int
    relation: PolynomialRelation


@dataclass(frozen=True)
class SelfPolynomialUniquenessScan:
    """Summary of low-degree self-polynomial uniqueness scans."""

    moduli_checked: tuple[int, ...]
    fg_degree_values: tuple[int, ...]
    t_degree_values: tuple[int, ...]
    hits: tuple[SelfPolynomialUniquenessHit, ...]


@dataclass(frozen=True)
class SelfTPolynomialFractionalLinearRelation:
    """A relation F = (A(T) + B(T)*(G-1)) / (C(T) + D(T)*(G-1))."""

    order_checked: int
    max_t_degree: int
    numerator_t_coefficients: tuple[sp.Expr, ...]
    numerator_self_coefficients: tuple[sp.Expr, ...]
    denominator_t_coefficients: tuple[sp.Expr, ...]
    denominator_self_coefficients: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class SelfFractionalLinearUniquenessHit:
    """A theorem-facing low-degree self-fractional-linear candidate."""

    modulus: int
    max_t_degree: int
    relation: SelfTPolynomialFractionalLinearRelation


@dataclass(frozen=True)
class SelfFractionalLinearUniquenessScan:
    """Summary of low-degree self-fractional-linear uniqueness scans."""

    moduli_checked: tuple[int, ...]
    t_degree_values: tuple[int, ...]
    hits: tuple[SelfFractionalLinearUniquenessHit, ...]


@dataclass(frozen=True)
class SelfMahlerLinearHit:
    """A bounded linear Mahler-style relation using deeper power substitutions."""

    modulus: int
    levels: int
    max_t_degree: int
    relation: PolynomialRelation


@dataclass(frozen=True)
class SelfMahlerLinearScan:
    """Summary of bounded linear Mahler-style scans."""

    moduli_checked: tuple[int, ...]
    levels_checked: tuple[int, ...]
    t_degree_values: tuple[int, ...]
    hits: tuple[SelfMahlerLinearHit, ...]


@dataclass(frozen=True)
class SignedSelfQuotientProductRelationScan:
    """Outcome for one modulus scan of a mixed (1-t^r)/(1+t^r) self-quotient box."""

    modulus: int
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class SelfSignedEtaRelationScan:
    """Outcome for one modulus/level scan of a mixed signed modular-unit self-functional box."""

    modulus: int
    level: int
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class ReducedTailTransferEquation:
    """Exact stationary tail law extracted from reduced continued-fraction coefficients."""

    start_stage: int
    stages_checked: int
    state_variable: str
    denominator_expr: sp.Expr
    numerator_expr: sp.Expr
    next_state_expr: sp.Expr


@dataclass(frozen=True)
class SourceCorrectionSelfPolynomialHit:
    """A self-polynomial uniqueness hit after factoring source cores."""

    basis_labels: tuple[str, ...]
    basis_expressions: tuple[str, ...]
    modulus: int
    max_fg_total_degree: int
    max_t_degree: int
    relation: PolynomialRelation


@dataclass(frozen=True)
class SourceCorrectionSelfPolynomialScan:
    """Summary of source-core-corrected self-polynomial scans."""

    correction_size: int
    moduli_checked: tuple[int, ...]
    fg_degree_values: tuple[int, ...]
    t_degree_values: tuple[int, ...]
    total_corrections_checked: int
    hits: tuple[SourceCorrectionSelfPolynomialHit, ...]


@dataclass(frozen=True)
class SourceCorrectionSelfFractionalLinearHit:
    """A self-fractional-linear uniqueness hit after factoring source cores."""

    basis_labels: tuple[str, ...]
    basis_expressions: tuple[str, ...]
    modulus: int
    max_t_degree: int
    relation: SelfTPolynomialFractionalLinearRelation


@dataclass(frozen=True)
class SourceCorrectionSelfFractionalLinearScan:
    """Summary of source-core-corrected self-fractional-linear scans."""

    correction_size: int
    moduli_checked: tuple[int, ...]
    t_degree_values: tuple[int, ...]
    total_corrections_checked: int
    hits: tuple[SourceCorrectionSelfFractionalLinearHit, ...]


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
class SelfQuotientProductRelation:
    """A periodic-product-style quotient relation F(t) / F(t^m) = prod_r (1 - t^r)^e_r."""

    order_checked: int
    modulus: int
    exponents_by_residue: dict[int, int]


@dataclass(frozen=True)
class SelfQuotientProductRelationScan:
    """Outcome for one modulus scan of the self-quotient finite-product box."""

    modulus: int
    relation: SelfQuotientProductRelation | None
    error: str | None = None


@dataclass(frozen=True)
class EtaQuotientRelationScan:
    """Outcome for one eta-quotient level scan."""

    level: int
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class EtaCorrectionBasisScan:
    """Eta-quotient correction scans for one fixed source-family basis choice."""

    basis_label: str
    basis_expression: str
    basis_series: Series
    eta_scans: tuple[EtaQuotientRelationScan, ...]


@dataclass(frozen=True)
class SourceFamilyEtaCorrectionScan:
    """Per-family eta-correction scans over raw and quotient basis choices."""

    family_label: str
    benchmark_name: str
    direct_basis_scans: tuple[EtaCorrectionBasisScan, ...]
    quotient_basis_scans: tuple[EtaCorrectionBasisScan, ...]


@dataclass(frozen=True)
class TwoCoreSourceFamilyEtaCorrectionHit:
    """A low-complexity cross-family two-core eta-correction hit."""

    basis_labels: tuple[str, str]
    basis_expressions: tuple[str, str]
    level: int
    relation: MultiplicativeRelation


@dataclass(frozen=True)
class TwoCoreSourceFamilyEtaCorrectionScan:
    """Cross-family raw-basis two-core eta-correction summary."""

    levels_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoCoreSourceFamilyEtaCorrectionHit, ...]


@dataclass(frozen=True)
class QuotientCoreSourceFamilyEtaCorrectionHit:
    """A low-complexity cross-family quotient-core eta-correction hit."""

    quotient_label: str
    quotient_expression: str
    raw_label: str
    raw_expression: str
    level: int
    relation: MultiplicativeRelation


@dataclass(frozen=True)
class QuotientCoreSourceFamilyEtaCorrectionScan:
    """Cross-family quotient-core eta-correction summary."""

    levels_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[QuotientCoreSourceFamilyEtaCorrectionHit, ...]


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilyEtaCorrectionHit:
    """A low-complexity cross-family quotient-pair eta-correction hit."""

    quotient_labels: tuple[str, str]
    quotient_expressions: tuple[str, str]
    level: int
    relation: MultiplicativeRelation


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilyEtaCorrectionScan:
    """Cross-family quotient-pair eta-correction summary."""

    levels_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoQuotientCoreSourceFamilyEtaCorrectionHit, ...]


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfQuotientProductHit:
    """A cross-family quotient-pair finite-product self-quotient hit."""

    quotient_labels: tuple[str, str]
    quotient_expressions: tuple[str, str]
    modulus: int
    relation: SelfQuotientProductRelation


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfQuotientProductScan:
    """Cross-family quotient-pair finite-product self-quotient summary."""

    moduli_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoQuotientCoreSourceFamilySelfQuotientProductHit, ...]


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfPolynomialHit:
    """A cross-family quotient-pair low-degree self-polynomial hit."""

    quotient_labels: tuple[str, str]
    quotient_expressions: tuple[str, str]
    modulus: int
    max_total_degree: int
    relation: PolynomialRelation


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfPolynomialScan:
    """Cross-family quotient-pair low-degree self-polynomial summary."""

    moduli_checked: tuple[int, ...]
    degree_values: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoQuotientCoreSourceFamilySelfPolynomialHit, ...]


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfEtaCorrectionHit:
    """A cross-family quotient-pair self-eta functional hit."""

    quotient_labels: tuple[str, str]
    quotient_expressions: tuple[str, str]
    modulus: int
    level: int
    relation: MultiplicativeRelation


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfEtaCorrectionScan:
    """Cross-family quotient-pair self-eta functional summary."""

    moduli_checked: tuple[int, ...]
    levels_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoQuotientCoreSourceFamilySelfEtaCorrectionHit, ...]


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfFractionalLinearHit:
    """A cross-family quotient-pair self-fractional-linear eta hit."""

    quotient_labels: tuple[str, str]
    quotient_expressions: tuple[str, str]
    modulus: int
    level: int
    relation: FractionalLinearRelation


@dataclass(frozen=True)
class TwoQuotientCoreSourceFamilySelfFractionalLinearScan:
    """Cross-family quotient-pair self-fractional-linear eta summary."""

    moduli_checked: tuple[int, ...]
    levels_checked: tuple[int, ...]
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    total_basis_pairs_checked: int
    total_boxes_checked: int
    hits: tuple[TwoQuotientCoreSourceFamilySelfFractionalLinearHit, ...]


@dataclass(frozen=True)
class NamedMultiplicativeRelationScan:
    """Outcome for one prefix scan of named source-family multiplicative templates."""

    basis_labels: tuple[str, ...]
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class NamedFractionalLinearRelationScan:
    """Outcome for one prefix scan of named source-family fractional-linear templates."""

    basis_labels: tuple[str, ...]
    relation: FractionalLinearRelation | None
    error: str | None = None


@dataclass(frozen=True)
class NamedPolynomialRelationScan:
    """Outcome for one prefix scan of named source-family polynomial templates."""

    basis_labels: tuple[str, ...]
    max_total_degree: int
    relation: PolynomialRelation | None
    error: str | None = None


@dataclass(frozen=True)
class NamedTwoLayerFractionalLinearRelationScan:
    """Outcome for one prefix scan of named source-family two-layer templates."""

    basis_labels: tuple[str, ...]
    relations: tuple[TwoLayerFractionalLinearRelation, ...]
    total_hits: int
    tuples_checked: int
    error: str | None = None


@dataclass(frozen=True)
class NamedPrefixBoxScans:
    """Bundled prefix scans over one ordered named basis ladder."""

    polynomial_scans: tuple[NamedPolynomialRelationScan, ...] = ()
    multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...] = ()
    fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...] = ()
    two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...] = ()


@dataclass(frozen=True)
class ParameterizedSourceFamilyScan:
    """Per-family powered-basis scans that keep the source-family label explicit."""

    family_label: str
    benchmark_name: str
    ordered_basis_series: tuple[tuple[str, Series], ...]
    polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    quotient_basis_series: tuple[tuple[str, str, Series], ...]
    quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    mixed_quotient_basis_series: tuple[tuple[str, str, Series], ...]
    mixed_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    mixed_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    mixed_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    mixed_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]


@dataclass(frozen=True)
class ExplicitSourceFamilyTransformScan:
    """Exact direct/reciprocal/quotient checks that keep source-family meaning explicit."""

    family_label: str
    benchmark_name: str
    ordered_basis_series: tuple[tuple[str, Series], ...]
    checked_templates: tuple[str, ...]
    hit_templates: tuple[str, ...]


@dataclass(frozen=True)
class ExplicitTransformEtaCorrectionHit:
    """One explicit-transform template whose residual factor is an eta quotient."""

    template_label: str
    level: int
    relation: MultiplicativeRelation


@dataclass(frozen=True)
class ExplicitSourceFamilyEtaCorrectionScan:
    """Eta-correction scans over the explicit direct/reciprocal/quotient template orbit."""

    family_label: str
    benchmark_name: str
    ordered_basis_series: tuple[tuple[str, Series], ...]
    checked_templates: tuple[str, ...]
    hits: tuple[ExplicitTransformEtaCorrectionHit, ...]


@dataclass(frozen=True)
class GGModularEquationScan:
    """A literature-driven GG box built from sign and low-power substitutions."""

    benchmark_name: str
    ordered_basis_series: tuple[tuple[str, str, Series], ...]
    checked_templates: tuple[str, ...]
    hit_templates: tuple[str, ...]
    exact_polynomial_template_labels: tuple[str, ...]
    exact_polynomial_template_hits: tuple[str, ...]
    polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    quotient_basis_series: tuple[tuple[str, str, Series], ...]
    quotient_exact_polynomial_template_labels: tuple[str, ...]
    quotient_exact_polynomial_template_hits: tuple[str, ...]
    quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    mixed_quotient_basis_series: tuple[tuple[str, str, Series], ...]
    mixed_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    mixed_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    mixed_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    mixed_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]


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


def _guess_polynomial_relation_from_exponent_tuples(
    *,
    series_by_variable: dict[str, Series],
    order: int,
    exponent_tuples: tuple[tuple[int, ...], ...],
    required_variables: tuple[str, ...] = (),
) -> PolynomialRelation | None:
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(series_by_variable) < 2:
        raise ValueError("need at least two variables for a relation search")
    if any(len(series) < order for series in series_by_variable.values()):
        raise ValueError("series are shorter than requested order")
    if not exponent_tuples:
        raise ValueError("need at least one monomial for a relation search")

    variables = tuple(series_by_variable.keys())
    if any(required not in series_by_variable for required in required_variables):
        raise ValueError("required_variables must be series variable names")
    if any(len(exponents) != len(variables) for exponents in exponent_tuples):
        raise ValueError("exponent tuple arity must match the number of variables")
    if any(any(exp < 0 for exp in exponents) for exponents in exponent_tuples):
        raise ValueError("exponents must be non-negative")

    num_monomials = len(exponent_tuples)
    if num_monomials > order:
        raise ValueError(
            "underdetermined polynomial relation search: "
            f"{num_monomials} monomials > {order} constraints "
            "(increase order, lower the search box, or reduce variables)"
        )

    series_list = [series_by_variable[name][:order] for name in variables]
    max_exponents_by_variable = [0 for _ in variables]
    for exponents in exponent_tuples:
        for index, exponent in enumerate(exponents):
            max_exponents_by_variable[index] = max(max_exponents_by_variable[index], exponent)

    one = [sp.Integer(0) for _ in range(order)]
    one[0] = sp.Integer(1)
    powers: list[list[Series]] = []
    for series, max_exponent in zip(series_list, max_exponents_by_variable):
        variable_powers: list[Series] = [one]
        for _ in range(max_exponent):
            variable_powers.append(series_mul(variable_powers[-1], series))
        powers.append(variable_powers)

    columns: list[Series] = []
    for exponents in exponent_tuples:
        term = one
        for variable_index, exponent in enumerate(exponents):
            if exponent == 0:
                continue
            term = series_mul(term, powers[variable_index][exponent])
        columns.append(term)

    matrix = sp.Matrix([[columns[col][row] for col in range(len(columns))] for row in range(order)])
    nullspace = matrix.nullspace()
    if not nullspace:
        return None

    candidate_vecs = nullspace
    for required_variable in required_variables:
        required_index = variables.index(required_variable)
        required_columns = [
            index
            for index, exponents in enumerate(exponent_tuples)
            if exponents[required_index] > 0
        ]
        candidate_vecs = [
            vec
            for vec in candidate_vecs
            if any(vec[column_index] != 0 for column_index in required_columns)
        ]
        if not candidate_vecs:
            return None

    basis_vec = min(candidate_vecs, key=lambda vec: sum(1 for item in vec if item != 0))
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

    int_scaled = [int(value) for value in scaled]
    overall_gcd = 0
    for value in int_scaled:
        overall_gcd = gcd(overall_gcd, abs(value))
    if overall_gcd > 1:
        scaled = [sp.Integer(int(value) // overall_gcd) for value in scaled]

    for value in scaled:
        if value != 0:
            if value < 0:
                scaled = [-value for value in scaled]
            break

    coeff_map: dict[tuple[int, ...], sp.Integer] = {}
    max_total_degree = 0
    for exponents, coeff in zip(exponent_tuples, scaled):
        coeff_s = sp.simplify(coeff)
        if coeff_s == 0:
            continue
        normalized_exponents = tuple(int(exponent) for exponent in exponents)
        coeff_map[normalized_exponents] = sp.Integer(int(coeff_s))
        max_total_degree = max(max_total_degree, sum(normalized_exponents))

    if not coeff_map:
        return None

    return PolynomialRelation(
        order_checked=order,
        variables=variables,
        max_total_degree=max_total_degree,
        coefficients=coeff_map,
    )


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

    variables = tuple(series_by_variable.keys())
    exponent_tuples: list[tuple[int, ...]] = []

    def _recurse(idx: int, remaining: int, current: list[int]) -> None:
        if idx == len(variables):
            exponent_tuples.append(tuple(current))
            return
        for exp in range(remaining + 1):
            current.append(exp)
            _recurse(idx + 1, remaining - exp, current)
            current.pop()

    _recurse(0, max_total_degree, [])
    return _guess_polynomial_relation_from_exponent_tuples(
        series_by_variable=series_by_variable,
        order=order,
        exponent_tuples=tuple(exponent_tuples),
        required_variables=() if required_variable is None else (required_variable,),
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


def _polynomial_relation_from_sympy_expr(
    *,
    expr: sp.Expr,
    variables: tuple[str, ...],
    order: int,
) -> PolynomialRelation:
    symbols = tuple(sp.Symbol(name) for name in variables)
    poly = sp.Poly(sp.expand(expr), *symbols)
    coeff_map: dict[tuple[int, ...], sp.Integer] = {}
    max_total_degree = 0
    for exponents, coeff in poly.terms():
        coeff_s = sp.simplify(coeff)
        if coeff_s == 0:
            continue
        normalized_exponents = tuple(int(value) for value in exponents)
        coeff_map[normalized_exponents] = sp.Integer(coeff_s)
        max_total_degree = max(max_total_degree, sum(normalized_exponents))
    if not coeff_map:
        raise ValueError("polynomial expression was identically zero")
    return PolynomialRelation(
        order_checked=order,
        variables=variables,
        max_total_degree=max_total_degree,
        coefficients=coeff_map,
    )


@lru_cache(maxsize=None)
def _two_layer_factor_index_pairs(
    basis_size: int,
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    if basis_size < 1:
        return ()
    factor_indices = tuple(product(range(basis_size), repeat=2))
    pairs: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for index, factor_1 in enumerate(factor_indices):
        for factor_2 in factor_indices[index:]:
            pairs.append((factor_1, factor_2))
    return tuple(pairs)


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


def signed_argument_substitution_series(base_series: Series, *, order: int) -> Series:
    if len(base_series) < order:
        raise ValueError("base_series is shorter than requested order")

    signed: Series = [sp.Integer(0) for _ in range(order)]
    for idx, coeff in enumerate(base_series[:order]):
        signed[idx] = sp.simplify(coeff if idx % 2 == 0 else -coeff)
    return signed


def _series_subtract_one(series: Series) -> Series:
    shifted = [sp.simplify(value) for value in series]
    shifted[0] = sp.simplify(shifted[0] - 1)
    return shifted


def _expr_to_series(expr: sp.Expr, *, symbol: sp.Symbol, order: int) -> Series:
    expanded = sp.expand(expr)
    return [sp.simplify(expanded.coeff(symbol, index)) for index in range(order)]


def _generalized_continued_fraction_series_from_coeffs(
    *,
    b0: sp.Expr,
    a_terms: list[sp.Expr],
    b_terms: list[sp.Expr],
    symbol: sp.Symbol,
    order: int,
) -> Series:
    if len(a_terms) != len(b_terms):
        raise ValueError("a_terms and b_terms must have the same length")
    if len(a_terms) < 2:
        raise ValueError("a_terms and b_terms must be 1-indexed with at least one term")
    if order < 2:
        raise ValueError("order must be at least 2")

    tail: Series = [sp.Integer(0) for _ in range(order)]
    for stage in range(len(a_terms) - 1, 0, -1):
        denominator_series = series_add(_expr_to_series(b_terms[stage], symbol=symbol, order=order), tail)
        numerator_series = _expr_to_series(a_terms[stage], symbol=symbol, order=order)
        tail = series_mul(numerator_series, series_invert(denominator_series))

    return series_add(_expr_to_series(b0, symbol=symbol, order=order), tail)


def _reduced_reciprocal_bridge(
    *,
    template: QCFTemplate,
    symbol: sp.Symbol,
    depth: int,
    order: int,
):
    b0, a_terms, b_terms = _template_reciprocal_coeffs(template.normalized(), q=symbol, depth=depth)
    witness = convergent_factor_equivalence_witness(
        b0=b0,
        a_terms=a_terms,
        b_terms=b_terms,
    )
    reduced_coeffs = witness.reduction.reduced_coeffs
    reduced_series = _generalized_continued_fraction_series_from_coeffs(
        b0=reduced_coeffs.b0,
        a_terms=reduced_coeffs.a_terms,
        b_terms=reduced_coeffs.b_terms,
        symbol=symbol,
        order=order,
    )
    return witness, reduced_series


def detect_reduced_tail_transfer_equation(
    *,
    reduced_coeffs: ContinuedFractionCoeffs,
    symbol: sp.Symbol,
    start_stage: int = 3,
    state_variable: str = "x",
) -> ReducedTailTransferEquation | None:
    """Detect an exact stationary tail law of the form b = 1 + x, a = x*(t + x), x' = t*x."""
    if start_stage < 1:
        raise ValueError("start_stage must be positive")
    if len(reduced_coeffs.a_terms) != len(reduced_coeffs.b_terms):
        raise ValueError("reduced coefficient lists must have the same length")
    if len(reduced_coeffs.a_terms) <= start_stage:
        return None

    state = sp.Symbol(state_variable)
    stages_checked = 0
    for stage in range(start_stage, len(reduced_coeffs.a_terms)):
        state_value = sp.simplify(reduced_coeffs.b_terms[stage] - 1)
        if sp.simplify(reduced_coeffs.b_terms[stage] - (1 + state_value)) != 0:
            return None
        if sp.simplify(reduced_coeffs.a_terms[stage] - state_value * (symbol + state_value)) != 0:
            return None
        stages_checked += 1
        if stage + 1 < len(reduced_coeffs.b_terms):
            next_state_value = sp.simplify(reduced_coeffs.b_terms[stage + 1] - 1)
            if sp.simplify(next_state_value - symbol * state_value) != 0:
                return None

    if stages_checked < 1:
        return None

    return ReducedTailTransferEquation(
        start_stage=start_stage,
        stages_checked=stages_checked,
        state_variable=state_variable,
        denominator_expr=1 + state,
        numerator_expr=sp.expand(state * (symbol + state)),
        next_state_expr=sp.expand(symbol * state),
    )


@lru_cache(maxsize=None)
def _series_log_coeffs_cached(series: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    """Cached coefficients of log(F) for a series F with constant term 1."""
    if not series:
        raise ValueError("series must be non-empty")
    if sp.simplify(series[0] - 1) != 0:
        raise ValueError("series constant term must be 1 for a logarithm expansion")

    coeffs: list[sp.Expr] = [sp.Integer(0) for _ in range(len(series))]
    for n in range(1, len(series)):
        rhs = sp.Integer(0)
        for j in range(1, n):
            rhs += sp.Integer(j) * coeffs[j] * series[n - j]
        coeffs[n] = sp.simplify(series[n] - rhs / sp.Integer(n))
    return tuple(coeffs)


def _series_log_coeffs(series: Series) -> Series:
    """Return coefficients of log(F) for a series F with constant term 1."""
    return list(_series_log_coeffs_cached(tuple(series)))


def _signed_series_pow(base: Series, exponent: int) -> Series:
    if exponent >= 0:
        return series_pow(base, exponent)
    return series_invert(series_pow(base, -exponent))


@lru_cache(maxsize=None)
def _one_minus_power_series_tuple(*, power: int, order: int) -> tuple[sp.Expr, ...]:
    if power < 1:
        raise ValueError("power must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    series: Series = [sp.Integer(0) for _ in range(order)]
    series[0] = sp.Integer(1)
    if power < order:
        series[power] = sp.Integer(-1)
    return tuple(series)


def _one_minus_power_series(*, power: int, order: int) -> Series:
    return list(_one_minus_power_series_tuple(power=power, order=order))


@lru_cache(maxsize=None)
def _one_plus_power_series_tuple(*, power: int, order: int) -> tuple[sp.Expr, ...]:
    if power < 1:
        raise ValueError("power must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    series: Series = [sp.Integer(0) for _ in range(order)]
    series[0] = sp.Integer(1)
    if power < order:
        series[power] = sp.Integer(1)
    return tuple(series)


def _one_plus_power_series(*, power: int, order: int) -> Series:
    return list(_one_plus_power_series_tuple(power=power, order=order))


@lru_cache(maxsize=None)
def _eta_pochhammer_series_tuple(*, divisor: int, order: int) -> tuple[sp.Expr, ...]:
    if divisor < 1:
        raise ValueError("divisor must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")

    series: Series = [sp.Integer(0) for _ in range(order)]
    series[0] = sp.Integer(1)
    for power in range(divisor, order, divisor):
        series = series_mul(series, list(_one_minus_power_series_tuple(power=power, order=order)))
    return tuple(series)


def _eta_pochhammer_series(*, divisor: int, order: int) -> Series:
    return list(_eta_pochhammer_series_tuple(divisor=divisor, order=order))


@lru_cache(maxsize=None)
def _eta_quotient_basis_series_tuple(*, level: int, order: int) -> tuple[tuple[str, tuple[sp.Expr, ...]], ...]:
    if level < 1:
        raise ValueError("level must be positive")
    divisors = tuple(divisor for divisor in range(1, level + 1) if level % divisor == 0)
    return tuple(
        (
            f"E{divisor}",
            _eta_pochhammer_series_tuple(divisor=divisor, order=order),
        )
        for divisor in divisors
    )


def _eta_quotient_basis_series(*, level: int, order: int) -> dict[str, Series]:
    return {
        name: list(series)
        for name, series in _eta_quotient_basis_series_tuple(level=level, order=order)
    }


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


def search_self_quotient_product_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_abs_exponent: int = 8,
) -> SelfQuotientProductRelation | None:
    """Search F(t) / F(t^m) = prod_{r=1}^{m-1} (1 - t^r)^e_r with bounded integer exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a self-quotient product search")

    target_power_series = benchmark_power_substitution_series(target_series, power=modulus, order=order)
    quotient_series = series_div(target_series[:order], target_power_series)
    basis_series_by_variable = {
        f"U{residue}": _one_minus_power_series(power=residue, order=order)
        for residue in range(1, modulus)
    }

    relation = search_multiplicative_relation(
        target_series=quotient_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )
    if relation is None:
        return None

    exponents_by_residue = {
        int(name.removeprefix("U")): exponent
        for name, exponent in relation.exponents.items()
    }
    if not exponents_by_residue:
        return None

    self_relation = SelfQuotientProductRelation(
        order_checked=order,
        modulus=modulus,
        exponents_by_residue=exponents_by_residue,
    )
    residual = _self_quotient_product_relation_residual_series(
        self_relation,
        target_series=target_series,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return self_relation


def scan_ratio_self_quotient_product_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SelfQuotientProductRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    if not unique_moduli:
        return []

    scans: list[SelfQuotientProductRelationScan] = []
    for modulus in unique_moduli:
        try:
            relation = search_self_quotient_product_relation(
                target_series=ratio_series,
                modulus=modulus,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(SelfQuotientProductRelationScan(modulus=modulus, relation=relation))
        except ValueError as exc:
            scans.append(
                SelfQuotientProductRelationScan(
                    modulus=modulus,
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def search_self_quotient_plus_product_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_abs_exponent: int = 8,
) -> SelfQuotientProductRelation | None:
    """Search F(t) / F(t^m) = prod_{r=1}^{m-1} (1 + t^r)^e_r with bounded integer exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a self-quotient plus-product search")

    target_power_series = benchmark_power_substitution_series(target_series, power=modulus, order=order)
    quotient_series = series_div(target_series[:order], target_power_series)
    basis_series_by_variable = {
        f"U{residue}": _one_plus_power_series(power=residue, order=order)
        for residue in range(1, modulus)
    }

    relation = search_multiplicative_relation(
        target_series=quotient_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )
    if relation is None:
        return None

    exponents_by_residue = {
        int(name.removeprefix("U")): exponent
        for name, exponent in relation.exponents.items()
    }
    if not exponents_by_residue:
        return None

    plus_relation = SelfQuotientProductRelation(
        order_checked=order,
        modulus=modulus,
        exponents_by_residue=exponents_by_residue,
    )
    product_series: Series = [sp.Integer(0) for _ in range(order)]
    product_series[0] = sp.Integer(1)
    for residue, exponent in sorted(exponents_by_residue.items()):
        factor = _signed_series_pow(_one_plus_power_series(power=residue, order=order), exponent)
        product_series = series_mul(product_series, factor)
    residual = [sp.simplify(quotient_series[index] - product_series[index]) for index in range(order)]
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return plus_relation


def scan_ratio_self_plus_product_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SelfQuotientProductRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    if not unique_moduli:
        return []

    scans: list[SelfQuotientProductRelationScan] = []
    for modulus in unique_moduli:
        try:
            relation = search_self_quotient_plus_product_relation(
                target_series=ratio_series,
                modulus=modulus,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(SelfQuotientProductRelationScan(modulus=modulus, relation=relation))
        except ValueError as exc:
            scans.append(
                SelfQuotientProductRelationScan(
                    modulus=modulus,
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def search_self_quotient_signed_product_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F(t)/F(t^m) = prod_r (1-t^r)^a_r (1+t^r)^b_r with bounded exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a signed self-quotient product search")

    target_power_series = benchmark_power_substitution_series(target_series, power=modulus, order=order)
    quotient_series = series_div(target_series[:order], target_power_series)
    basis_series_by_variable = {
        **{f"M{residue}": _one_minus_power_series(power=residue, order=order) for residue in range(1, modulus)},
        **{f"P{residue}": _one_plus_power_series(power=residue, order=order) for residue in range(1, modulus)},
    }
    return search_multiplicative_relation(
        target_series=quotient_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )


def scan_ratio_self_signed_product_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SignedSelfQuotientProductRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    if not unique_moduli:
        return []

    scans: list[SignedSelfQuotientProductRelationScan] = []
    for modulus in unique_moduli:
        try:
            relation = search_self_quotient_signed_product_relation(
                target_series=ratio_series,
                modulus=modulus,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(SignedSelfQuotientProductRelationScan(modulus=modulus, relation=relation))
        except ValueError as exc:
            scans.append(
                SignedSelfQuotientProductRelationScan(
                    modulus=modulus,
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def search_self_signed_eta_relation(
    *,
    target_series: Series,
    modulus: int,
    level: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = F(t^m)^a * prod (1-t^r)^u_r (1+t^r)^v_r * eta(t)^w with bounded exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if level < 1:
        raise ValueError("level must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a self signed-eta search")

    g_power_label = f"G{modulus}"
    basis_series_by_variable = {
        g_power_label: benchmark_power_substitution_series(target_series, power=modulus, order=order),
        **{f"M{residue}": _one_minus_power_series(power=residue, order=order) for residue in range(1, modulus)},
        **{f"P{residue}": _one_plus_power_series(power=residue, order=order) for residue in range(1, modulus)},
        **_eta_quotient_basis_series(level=level, order=order),
    }
    relation = search_multiplicative_relation(
        target_series=target_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )
    if relation is None:
        return None
    g_power_exponent = relation.exponents.get(g_power_label)
    if g_power_exponent not in {-1, 1}:
        return None
    return relation


def scan_ratio_self_signed_eta_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SelfSignedEtaRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    unique_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not unique_moduli or not unique_levels:
        return []

    scans: list[SelfSignedEtaRelationScan] = []
    for modulus in unique_moduli:
        for level in unique_levels:
            try:
                relation = search_self_signed_eta_relation(
                    target_series=ratio_series,
                    modulus=modulus,
                    level=level,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
                scans.append(SelfSignedEtaRelationScan(modulus=modulus, level=level, relation=relation))
            except ValueError as exc:
                scans.append(
                    SelfSignedEtaRelationScan(
                        modulus=modulus,
                        level=level,
                        relation=None,
                        error=str(exc),
                    )
                )
    return scans


def _t_series(*, order: int) -> Series:
    series = [sp.Integer(0) for _ in range(order)]
    if order > 1:
        series[1] = sp.Integer(1)
    return series


def search_self_polynomial_uniqueness_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_fg_total_degree: int,
    max_t_degree: int,
) -> PolynomialRelation | None:
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if max_fg_total_degree < 1:
        raise ValueError("max_fg_total_degree must be at least 1")
    if max_t_degree < 0:
        raise ValueError("max_t_degree must be non-negative")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    g_label = f"G{modulus}"
    series_by_variable = {
        "T": _t_series(order=order),
        "F": target_series[:order],
        g_label: benchmark_power_substitution_series(target_series, power=modulus, order=order),
    }
    exponent_tuples: list[tuple[int, int, int]] = []
    for t_degree in range(max_t_degree + 1):
        for f_degree in range(max_fg_total_degree + 1):
            for g_degree in range(max_fg_total_degree - f_degree + 1):
                exponent_tuples.append((t_degree, f_degree, g_degree))

    relation = _guess_polynomial_relation_from_exponent_tuples(
        series_by_variable=series_by_variable,
        order=order,
        exponent_tuples=tuple(exponent_tuples),
        required_variables=("F", g_label),
    )
    if relation is None:
        return None

    residual = _relation_residual_series(
        relation,
        series_by_variable=series_by_variable,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return relation


def scan_self_polynomial_uniqueness_relations(
    *,
    target_series: Series,
    moduli: tuple[int, ...],
    order: int,
    fg_degree_values: tuple[int, ...] = (1, 2),
    t_degree_values: tuple[int, ...] = (1, 2),
) -> SelfPolynomialUniquenessScan:
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_fg_degrees = tuple(sorted({degree for degree in fg_degree_values if degree >= 1}))
    normalized_t_degrees = tuple(sorted({degree for degree in t_degree_values if degree >= 0}))
    if not normalized_moduli or not normalized_fg_degrees or not normalized_t_degrees:
        return SelfPolynomialUniquenessScan(
            moduli_checked=normalized_moduli,
            fg_degree_values=normalized_fg_degrees,
            t_degree_values=normalized_t_degrees,
            hits=(),
        )

    hits: list[SelfPolynomialUniquenessHit] = []
    for modulus in normalized_moduli:
        for fg_degree in normalized_fg_degrees:
            for t_degree in normalized_t_degrees:
                try:
                    relation = search_self_polynomial_uniqueness_relation(
                        target_series=target_series,
                        modulus=modulus,
                        order=order,
                        max_fg_total_degree=fg_degree,
                        max_t_degree=t_degree,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                hits.append(
                    SelfPolynomialUniquenessHit(
                        modulus=modulus,
                        max_fg_total_degree=fg_degree,
                        max_t_degree=t_degree,
                        relation=relation,
                    )
                )

    return SelfPolynomialUniquenessScan(
        moduli_checked=normalized_moduli,
        fg_degree_values=normalized_fg_degrees,
        t_degree_values=normalized_t_degrees,
        hits=tuple(hits),
    )


def search_eta_quotient_relation(
    *,
    target_series: Series,
    level: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = prod_{d|N} (t^d; t^d)_inf^e_d with bounded integer exponents."""
    if level < 1:
        raise ValueError("level must be positive")
    basis_series_by_variable = _eta_quotient_basis_series(level=level, order=order)
    return search_multiplicative_relation(
        target_series=target_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )


def scan_ratio_eta_quotient_relations(
    *,
    ratio_series: Series,
    levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[EtaQuotientRelationScan]:
    unique_levels = tuple(sorted({level for level in levels if level >= 1}))
    if not unique_levels:
        return []

    scans: list[EtaQuotientRelationScan] = []
    for level in unique_levels:
        try:
            relation = search_eta_quotient_relation(
                target_series=ratio_series,
                level=level,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(EtaQuotientRelationScan(level=level, relation=relation))
        except ValueError as exc:
            scans.append(
                EtaQuotientRelationScan(
                    level=level,
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def scan_source_family_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> list[SourceFamilyEtaCorrectionScan]:
    if not ordered_base_families:
        return []

    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not normalized_eta_levels:
        return []

    scans: list[SourceFamilyEtaCorrectionScan] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        family_powers = set(unique_powers)
        if supplemental_powers_by_family is not None:
            family_powers.update(
                power
                for power in supplemental_powers_by_family.get(family_label, ())
                if power >= 2
            )

        direct_basis_entries: list[tuple[str, str, Series]] = [(family_label, family_label, base_series)]
        for power in tuple(sorted(family_powers)):
            direct_basis_entries.append(
                (
                    f"{family_label}{power}",
                    f"{family_label}{power}",
                    benchmark_power_substitution_series(base_series, power=power, order=order),
                )
            )

        quotient_basis_entries: list[tuple[str, str, Series]] = []
        for label, _, basis_series in direct_basis_entries[1:]:
            quotient_basis_entries.append(
                (
                    f"Q{int(label.removeprefix(family_label))}",
                    f"{label} / {family_label}",
                    series_div(basis_series, base_series),
                )
            )

        direct_basis_scans = tuple(
            EtaCorrectionBasisScan(
                basis_label=label,
                basis_expression=expr,
                basis_series=basis_series,
                eta_scans=tuple(
                    scan_ratio_eta_quotient_relations(
                        ratio_series=series_div(target_series, basis_series),
                        levels=normalized_eta_levels,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                ),
            )
            for label, expr, basis_series in direct_basis_entries
        )
        quotient_basis_scans = tuple(
            EtaCorrectionBasisScan(
                basis_label=label,
                basis_expression=expr,
                basis_series=basis_series,
                eta_scans=tuple(
                    scan_ratio_eta_quotient_relations(
                        ratio_series=series_div(target_series, basis_series),
                        levels=normalized_eta_levels,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                ),
            )
            for label, expr, basis_series in quotient_basis_entries
        )

        scans.append(
            SourceFamilyEtaCorrectionScan(
                family_label=family_label,
                benchmark_name=benchmark_name,
                direct_basis_scans=direct_basis_scans,
                quotient_basis_scans=quotient_basis_scans,
            )
        )
    return scans


def _source_family_raw_basis_entries(
    *,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    order: int,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> tuple[tuple[str, str, str, Series], ...]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    entries: list[tuple[str, str, str, Series]] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        family_powers = set(unique_powers)
        if supplemental_powers_by_family is not None:
            family_powers.update(
                power
                for power in supplemental_powers_by_family.get(family_label, ())
                if power >= 2
            )
        entries.append((family_label, benchmark_name, family_label, base_series))
        for power in tuple(sorted(family_powers)):
            entries.append(
                (
                    family_label,
                    benchmark_name,
                    f"{family_label}{power}",
                    benchmark_power_substitution_series(base_series, power=power, order=order),
                )
            )
    return tuple(entries)


def _source_family_quotient_basis_entries(
    *,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    order: int,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> tuple[tuple[str, str, str, str, Series], ...]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    entries: list[tuple[str, str, str, str, Series]] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        family_powers = set(unique_powers)
        if supplemental_powers_by_family is not None:
            family_powers.update(
                power
                for power in supplemental_powers_by_family.get(family_label, ())
                if power >= 2
            )
        for power in tuple(sorted(family_powers)):
            powered_label = f"{family_label}{power}"
            entries.append(
                (
                    family_label,
                    benchmark_name,
                    f"{family_label}_Q{power}",
                    f"{powered_label} / {family_label}",
                    series_div(
                        benchmark_power_substitution_series(base_series, power=power, order=order),
                        base_series,
                    ),
                )
            )
    return tuple(entries)


def _gg_modular_equation_quotient_basis_series(
    ordered_basis_entries: tuple[tuple[str, str, Series], ...],
) -> tuple[tuple[str, str, Series], ...]:
    if not ordered_basis_entries:
        return ()

    base_label, _, base_series = ordered_basis_entries[0]
    entries: list[tuple[str, str, Series]] = []
    for label, expression, basis_series in ordered_basis_entries[1:]:
        entries.append(
            (
                f"Q_{label.removeprefix(base_label)}",
                f"{expression} / {base_label}(t)",
                series_div(basis_series, base_series),
            )
        )
    return tuple(entries)


def scan_two_core_source_family_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    raw_basis_entries: tuple[tuple[str, str, str, Series], ...] | None = None,
) -> TwoCoreSourceFamilyEtaCorrectionScan:
    if raw_basis_entries is None:
        raw_basis_entries = _source_family_raw_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not raw_basis_entries or not normalized_eta_levels:
        return TwoCoreSourceFamilyEtaCorrectionScan(
            levels_checked=normalized_eta_levels,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    eta_basis_by_level = {
        level: _eta_quotient_basis_series(level=level, order=order)
        for level in normalized_eta_levels
    }
    family_pair_counts: dict[str, int] = {}
    hits: list[TwoCoreSourceFamilyEtaCorrectionHit] = []
    total_basis_pairs_checked = 0
    for left_index, left_entry in enumerate(raw_basis_entries):
        left_family, _, left_label, left_series = left_entry
        for right_entry in raw_basis_entries[left_index + 1 :]:
            right_family, _, right_label, right_series = right_entry
            if left_family == right_family:
                continue
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
            total_basis_pairs_checked += 1

            for level, eta_basis in eta_basis_by_level.items():
                basis_series_by_variable = {
                    left_label: left_series,
                    right_label: right_series,
                    **eta_basis,
                }
                try:
                    relation = search_multiplicative_relation(
                        target_series=target_series,
                        basis_series_by_variable=basis_series_by_variable,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                left_exponent = relation.exponents.get(left_label)
                right_exponent = relation.exponents.get(right_label)
                if left_exponent not in {-1, 1} or right_exponent not in {-1, 1}:
                    continue
                hits.append(
                    TwoCoreSourceFamilyEtaCorrectionHit(
                        basis_labels=(left_label, right_label),
                        basis_expressions=(left_label, right_label),
                        level=level,
                        relation=relation,
                    )
                )

    return TwoCoreSourceFamilyEtaCorrectionScan(
        levels_checked=normalized_eta_levels,
        family_pair_basis_counts=tuple(sorted(family_pair_counts.items())),
        total_basis_pairs_checked=total_basis_pairs_checked,
        total_boxes_checked=total_basis_pairs_checked * len(normalized_eta_levels),
        hits=tuple(hits),
    )


def scan_quotient_core_source_family_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
) -> QuotientCoreSourceFamilyEtaCorrectionScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    raw_basis_entries = _source_family_raw_basis_entries(
        ordered_base_families=ordered_base_families,
        powers=powers,
        order=order,
        supplemental_powers_by_family=supplemental_powers_by_family,
    )
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not quotient_basis_entries or not raw_basis_entries or not normalized_eta_levels:
        return QuotientCoreSourceFamilyEtaCorrectionScan(
            levels_checked=normalized_eta_levels,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    eta_basis_by_level = {
        level: _eta_quotient_basis_series(level=level, order=order)
        for level in normalized_eta_levels
    }
    family_pair_counts: dict[str, int] = {}
    hits: list[QuotientCoreSourceFamilyEtaCorrectionHit] = []
    total_basis_pairs_checked = 0
    for quotient_entry in quotient_basis_entries:
        quotient_family, _, quotient_label, quotient_expression, quotient_series = quotient_entry
        for raw_entry in raw_basis_entries:
            raw_family, _, raw_label, raw_series = raw_entry
            if quotient_family == raw_family:
                continue
            pair_label = f"{quotient_family}->{raw_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
            total_basis_pairs_checked += 1

            for level, eta_basis in eta_basis_by_level.items():
                basis_series_by_variable = {
                    quotient_label: quotient_series,
                    raw_label: raw_series,
                    **eta_basis,
                }
                try:
                    relation = search_multiplicative_relation(
                        target_series=target_series,
                        basis_series_by_variable=basis_series_by_variable,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                quotient_exponent = relation.exponents.get(quotient_label)
                raw_exponent = relation.exponents.get(raw_label)
                if quotient_exponent not in {-1, 1} or raw_exponent not in {-1, 1}:
                    continue
                hits.append(
                    QuotientCoreSourceFamilyEtaCorrectionHit(
                        quotient_label=quotient_label,
                        quotient_expression=quotient_expression,
                        raw_label=raw_label,
                        raw_expression=raw_label,
                        level=level,
                        relation=relation,
                    )
                )

    return QuotientCoreSourceFamilyEtaCorrectionScan(
        levels_checked=normalized_eta_levels,
        family_pair_basis_counts=tuple(sorted(family_pair_counts.items())),
        total_basis_pairs_checked=total_basis_pairs_checked,
        total_boxes_checked=total_basis_pairs_checked * len(normalized_eta_levels),
        hits=tuple(hits),
    )


def scan_two_quotient_core_source_family_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
) -> TwoQuotientCoreSourceFamilyEtaCorrectionScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not quotient_basis_entries or not normalized_eta_levels:
        return TwoQuotientCoreSourceFamilyEtaCorrectionScan(
            levels_checked=normalized_eta_levels,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    eta_basis_by_level = {
        level: _eta_quotient_basis_series(level=level, order=order)
        for level in normalized_eta_levels
    }
    family_pair_counts: dict[str, int] = {}
    hits: list[TwoQuotientCoreSourceFamilyEtaCorrectionHit] = []
    total_basis_pairs_checked = 0
    for left_index, left_entry in enumerate(quotient_basis_entries):
        left_family, _, left_label, left_expression, left_series = left_entry
        for right_entry in quotient_basis_entries[left_index + 1 :]:
            right_family, _, right_label, right_expression, right_series = right_entry
            if left_family == right_family:
                continue
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
            total_basis_pairs_checked += 1

            for level, eta_basis in eta_basis_by_level.items():
                basis_series_by_variable = {
                    left_label: left_series,
                    right_label: right_series,
                    **eta_basis,
                }
                try:
                    relation = search_multiplicative_relation(
                        target_series=target_series,
                        basis_series_by_variable=basis_series_by_variable,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                left_exponent = relation.exponents.get(left_label)
                right_exponent = relation.exponents.get(right_label)
                if left_exponent not in {-1, 1} or right_exponent not in {-1, 1}:
                    continue
                hits.append(
                    TwoQuotientCoreSourceFamilyEtaCorrectionHit(
                        quotient_labels=(left_label, right_label),
                        quotient_expressions=(left_expression, right_expression),
                        level=level,
                        relation=relation,
                    )
                )

    return TwoQuotientCoreSourceFamilyEtaCorrectionScan(
        levels_checked=normalized_eta_levels,
        family_pair_basis_counts=tuple(sorted(family_pair_counts.items())),
        total_basis_pairs_checked=total_basis_pairs_checked,
        total_boxes_checked=total_basis_pairs_checked * len(normalized_eta_levels),
        hits=tuple(hits),
    )


def _two_quotient_core_correction_entries(
    *,
    target_series: Series,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...],
) -> tuple[
    tuple[tuple[str, int], ...],
    tuple[tuple[str, str, str, str, str, str, Series], ...],
]:
    family_pair_counts: dict[str, int] = {}
    entries: list[tuple[str, str, str, str, str, str, Series]] = []
    for left_index, left_entry in enumerate(quotient_basis_entries):
        left_family, _, left_label, left_expression, left_series = left_entry
        for right_entry in quotient_basis_entries[left_index + 1 :]:
            right_family, _, right_label, right_expression, right_series = right_entry
            if left_family == right_family:
                continue
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
            correction_series = series_div(target_series, series_mul(left_series, right_series))
            entries.append(
                (
                    left_family,
                    right_family,
                    left_label,
                    right_label,
                    left_expression,
                    right_expression,
                    correction_series,
                )
            )
    return tuple(sorted(family_pair_counts.items())), tuple(entries)


def scan_two_quotient_core_source_family_self_quotient_products(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    moduli: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
    correction_entries: tuple[tuple[str, str, str, str, str, str, Series], ...] | None = None,
) -> TwoQuotientCoreSourceFamilySelfQuotientProductScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    if not quotient_basis_entries or not normalized_moduli:
        return TwoQuotientCoreSourceFamilySelfQuotientProductScan(
            moduli_checked=normalized_moduli,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    family_pair_basis_counts: tuple[tuple[str, int], ...]
    if correction_entries is None:
        family_pair_basis_counts, correction_entries = _two_quotient_core_correction_entries(
            target_series=target_series,
            quotient_basis_entries=quotient_basis_entries,
        )
    else:
        family_pair_counts: dict[str, int] = {}
        for left_family, right_family, *_ in correction_entries:
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
        family_pair_basis_counts = tuple(sorted(family_pair_counts.items()))
    hits: list[TwoQuotientCoreSourceFamilySelfQuotientProductHit] = []
    for _, _, left_label, right_label, left_expression, right_expression, correction_series in correction_entries:
        for modulus in normalized_moduli:
            try:
                relation = search_self_quotient_product_relation(
                    target_series=correction_series,
                    modulus=modulus,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            except ValueError:
                continue
            if relation is None:
                continue
            hits.append(
                TwoQuotientCoreSourceFamilySelfQuotientProductHit(
                    quotient_labels=(left_label, right_label),
                    quotient_expressions=(left_expression, right_expression),
                    modulus=modulus,
                    relation=relation,
                )
            )

    return TwoQuotientCoreSourceFamilySelfQuotientProductScan(
        moduli_checked=normalized_moduli,
        family_pair_basis_counts=family_pair_basis_counts,
        total_basis_pairs_checked=len(correction_entries),
        total_boxes_checked=len(correction_entries) * len(normalized_moduli),
        hits=tuple(hits),
    )


def scan_two_quotient_core_source_family_self_polynomial_relations(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    moduli: tuple[int, ...],
    order: int,
    degree_values: tuple[int, ...] = (1, 2),
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
    correction_entries: tuple[tuple[str, str, str, str, str, str, Series], ...] | None = None,
) -> TwoQuotientCoreSourceFamilySelfPolynomialScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_degrees = tuple(sorted({degree for degree in degree_values if degree >= 1}))
    if not quotient_basis_entries or not normalized_moduli or not normalized_degrees:
        return TwoQuotientCoreSourceFamilySelfPolynomialScan(
            moduli_checked=normalized_moduli,
            degree_values=normalized_degrees,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    family_pair_basis_counts: tuple[tuple[str, int], ...]
    if correction_entries is None:
        family_pair_basis_counts, correction_entries = _two_quotient_core_correction_entries(
            target_series=target_series,
            quotient_basis_entries=quotient_basis_entries,
        )
    else:
        family_pair_counts: dict[str, int] = {}
        for left_family, right_family, *_ in correction_entries:
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
        family_pair_basis_counts = tuple(sorted(family_pair_counts.items()))
    hits: list[TwoQuotientCoreSourceFamilySelfPolynomialHit] = []
    for _, _, left_label, right_label, left_expression, right_expression, correction_series in correction_entries:
        for modulus in normalized_moduli:
            powered_correction = benchmark_power_substitution_series(
                correction_series,
                power=modulus,
                order=order,
            )
            series_by_variable = {
                "G": correction_series,
                f"G{modulus}": powered_correction,
            }
            for degree in normalized_degrees:
                try:
                    relation = search_polynomial_relation(
                        series_by_variable=series_by_variable,
                        order=order,
                        max_total_degree=degree,
                        required_variable="G",
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                hits.append(
                    TwoQuotientCoreSourceFamilySelfPolynomialHit(
                        quotient_labels=(left_label, right_label),
                        quotient_expressions=(left_expression, right_expression),
                        modulus=modulus,
                        max_total_degree=degree,
                        relation=relation,
                    )
                )

    return TwoQuotientCoreSourceFamilySelfPolynomialScan(
        moduli_checked=normalized_moduli,
        degree_values=normalized_degrees,
        family_pair_basis_counts=family_pair_basis_counts,
        total_basis_pairs_checked=len(correction_entries),
        total_boxes_checked=len(correction_entries) * len(normalized_moduli) * len(normalized_degrees),
        hits=tuple(hits),
    )


def scan_two_quotient_core_source_family_self_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
    correction_entries: tuple[tuple[str, str, str, str, str, str, Series], ...] | None = None,
) -> TwoQuotientCoreSourceFamilySelfEtaCorrectionScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not quotient_basis_entries or not normalized_moduli or not normalized_levels:
        return TwoQuotientCoreSourceFamilySelfEtaCorrectionScan(
            moduli_checked=normalized_moduli,
            levels_checked=normalized_levels,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    eta_basis_by_level = {
        level: _eta_quotient_basis_series(level=level, order=order)
        for level in normalized_levels
    }
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    if correction_entries is None:
        family_pair_basis_counts, correction_entries = _two_quotient_core_correction_entries(
            target_series=target_series,
            quotient_basis_entries=quotient_basis_entries,
        )
    else:
        family_pair_counts: dict[str, int] = {}
        for left_family, right_family, *_ in correction_entries:
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
        family_pair_basis_counts = tuple(sorted(family_pair_counts.items()))
    hits: list[TwoQuotientCoreSourceFamilySelfEtaCorrectionHit] = []
    for _, _, left_label, right_label, left_expression, right_expression, correction_series in correction_entries:
        for modulus in normalized_moduli:
            powered_correction = benchmark_power_substitution_series(
                correction_series,
                power=modulus,
                order=order,
            )
            for level, eta_basis in eta_basis_by_level.items():
                g_power_label = f"G{modulus}"
                basis_series_by_variable = {
                    g_power_label: powered_correction,
                    **eta_basis,
                }
                try:
                    relation = search_multiplicative_relation(
                        target_series=correction_series,
                        basis_series_by_variable=basis_series_by_variable,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                g_power_exponent = relation.exponents.get(g_power_label)
                if g_power_exponent not in {-1, 1}:
                    continue
                hits.append(
                    TwoQuotientCoreSourceFamilySelfEtaCorrectionHit(
                        quotient_labels=(left_label, right_label),
                        quotient_expressions=(left_expression, right_expression),
                        modulus=modulus,
                        level=level,
                        relation=relation,
                    )
                )

    return TwoQuotientCoreSourceFamilySelfEtaCorrectionScan(
        moduli_checked=normalized_moduli,
        levels_checked=normalized_levels,
        family_pair_basis_counts=family_pair_basis_counts,
        total_basis_pairs_checked=len(correction_entries),
        total_boxes_checked=len(correction_entries) * len(normalized_moduli) * len(normalized_levels),
        hits=tuple(hits),
    )


def scan_two_quotient_core_source_family_self_fractional_linear_relations(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    quotient_basis_entries: tuple[tuple[str, str, str, str, Series], ...] | None = None,
    correction_entries: tuple[tuple[str, str, str, str, str, str, Series], ...] | None = None,
) -> TwoQuotientCoreSourceFamilySelfFractionalLinearScan:
    if quotient_basis_entries is None:
        quotient_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=ordered_base_families,
            powers=powers,
            order=order,
            supplemental_powers_by_family=supplemental_powers_by_family,
        )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not quotient_basis_entries or not normalized_moduli or not normalized_levels:
        return TwoQuotientCoreSourceFamilySelfFractionalLinearScan(
            moduli_checked=normalized_moduli,
            levels_checked=normalized_levels,
            family_pair_basis_counts=(),
            total_basis_pairs_checked=0,
            total_boxes_checked=0,
            hits=(),
        )

    eta_basis_by_level = {
        level: _eta_quotient_basis_series(level=level, order=order)
        for level in normalized_levels
    }
    family_pair_basis_counts: tuple[tuple[str, int], ...]
    if correction_entries is None:
        family_pair_basis_counts, correction_entries = _two_quotient_core_correction_entries(
            target_series=target_series,
            quotient_basis_entries=quotient_basis_entries,
        )
    else:
        family_pair_counts: dict[str, int] = {}
        for left_family, right_family, *_ in correction_entries:
            pair_label = f"{left_family}×{right_family}"
            family_pair_counts[pair_label] = family_pair_counts.get(pair_label, 0) + 1
        family_pair_basis_counts = tuple(sorted(family_pair_counts.items()))
    hits: list[TwoQuotientCoreSourceFamilySelfFractionalLinearHit] = []
    for _, _, left_label, right_label, left_expression, right_expression, correction_series in correction_entries:
        for modulus in normalized_moduli:
            powered_correction = benchmark_power_substitution_series(
                correction_series,
                power=modulus,
                order=order,
            )
            g_power_label = f"G{modulus}"
            for level, eta_basis in eta_basis_by_level.items():
                basis_series_by_variable = {
                    g_power_label: powered_correction,
                    **eta_basis,
                }
                try:
                    relation = search_fractional_linear_relation(
                        target_series=correction_series,
                        basis_series_by_variable=basis_series_by_variable,
                        order=order,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                uses_g_power = (
                    relation.numerator_coefficients.get(g_power_label) is not None
                    or relation.denominator_coefficients.get(g_power_label) is not None
                )
                if not uses_g_power:
                    continue
                hits.append(
                    TwoQuotientCoreSourceFamilySelfFractionalLinearHit(
                        quotient_labels=(left_label, right_label),
                        quotient_expressions=(left_expression, right_expression),
                        modulus=modulus,
                        level=level,
                        relation=relation,
                    )
                )

    return TwoQuotientCoreSourceFamilySelfFractionalLinearScan(
        moduli_checked=normalized_moduli,
        levels_checked=normalized_levels,
        family_pair_basis_counts=family_pair_basis_counts,
        total_basis_pairs_checked=len(correction_entries),
        total_boxes_checked=len(correction_entries) * len(normalized_moduli) * len(normalized_levels),
        hits=tuple(hits),
    )


def _format_self_quotient_product_relation(
    relation: SelfQuotientProductRelation,
    *,
    target_variable: str,
    series_symbol: str,
) -> str:
    terms: list[str] = []
    for residue in sorted(relation.exponents_by_residue):
        exponent = relation.exponents_by_residue[residue]
        if exponent == 0:
            continue
        base = f"(1 - {series_symbol})" if residue == 1 else f"(1 - {series_symbol}^{residue})"
        if exponent == 1:
            terms.append(base)
        else:
            terms.append(f"{base}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable}({series_symbol}) / {target_variable}({series_symbol}^{relation.modulus}) = {rhs}"


def _format_self_plus_product_relation(
    relation: SelfQuotientProductRelation,
    *,
    target_variable: str,
    series_symbol: str,
) -> str:
    terms: list[str] = []
    for residue in sorted(relation.exponents_by_residue):
        exponent = relation.exponents_by_residue[residue]
        if exponent == 0:
            continue
        base = f"(1 + {series_symbol})" if residue == 1 else f"(1 + {series_symbol}^{residue})"
        if exponent == 1:
            terms.append(base)
        else:
            terms.append(f"{base}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable}({series_symbol}) / {target_variable}({series_symbol}^{relation.modulus}) = {rhs}"


def _format_self_signed_product_relation(
    relation: MultiplicativeRelation,
    *,
    modulus: int,
    target_variable: str,
    series_symbol: str,
) -> str:
    terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name.startswith("M"):
            residue = int(name.removeprefix("M"))
            base = f"(1 - {series_symbol})" if residue == 1 else f"(1 - {series_symbol}^{residue})"
        elif name.startswith("P"):
            residue = int(name.removeprefix("P"))
            base = f"(1 + {series_symbol})" if residue == 1 else f"(1 + {series_symbol}^{residue})"
        else:
            base = name
        terms.append(base if exponent == 1 else f"{base}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable}({series_symbol}) / {target_variable}({series_symbol}^{modulus}) = {rhs}"


def _format_self_signed_eta_relation(
    relation: MultiplicativeRelation,
    *,
    modulus: int,
    target_variable: str,
    series_symbol: str,
) -> str:
    g_power_label = f"G{modulus}"
    terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name == g_power_label:
            base = f"{target_variable}({series_symbol}^{modulus})"
            terms.append(base if exponent == 1 else f"({base})^{exponent}")
            continue
        if name.startswith("M"):
            residue = int(name.removeprefix("M"))
            base = f"(1 - {series_symbol})" if residue == 1 else f"(1 - {series_symbol}^{residue})"
            terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        if name.startswith("P"):
            residue = int(name.removeprefix("P"))
            base = f"(1 + {series_symbol})" if residue == 1 else f"(1 + {series_symbol}^{residue})"
            terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        if name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
            terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        terms.append(name if exponent == 1 else f"{name}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable} = {rhs}"


def _format_reduced_tail_transfer_equation(
    relation: ReducedTailTransferEquation,
) -> tuple[str, str, str, str]:
    state_variable = relation.state_variable
    state_symbol = sp.Symbol(state_variable)
    base_symbol = _format_expr(sp.simplify(relation.next_state_expr / state_symbol))
    next_state = _format_expr(relation.next_state_expr)
    state_n = f"{state_variable}_n"
    next_state_n = f"{base_symbol}*{state_n}"
    return (
        f"For n >= {relation.start_stage}, let {state_variable}_n = b_n_red - 1.",
        f"Then b_n_red = 1 + {state_n}, a_n_red = {state_n}*({base_symbol} + {state_n}), {state_variable}_{{n+1}} = {next_state_n}.",
        f"T({state_variable}) = 1 + {state_variable} + ({state_variable}*({base_symbol} + {state_variable}))/T({next_state})",
        f"T({state_variable})*T({next_state}) - (1 + {state_variable})*T({next_state}) - {state_variable}*({base_symbol} + {state_variable}) = 0",
    )


def _format_eta_quotient_relation(
    relation: MultiplicativeRelation,
    *,
    target_variable: str,
    series_symbol: str,
) -> str:
    terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        divisor = int(name.removeprefix("E"))
        base = (
            f"({series_symbol}; {series_symbol})_inf"
            if divisor == 1
            else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
        )
        if exponent == 1:
            terms.append(base)
        else:
            terms.append(f"{base}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable} = {rhs}"


def _format_source_family_eta_correction(
    *,
    basis_expression: str,
    relation: MultiplicativeRelation,
    target_variable: str,
    series_symbol: str,
) -> str:
    eta_product = _format_eta_quotient_relation(
        relation,
        target_variable="G",
        series_symbol=series_symbol,
    ).split(" = ", 1)[1]
    basis = basis_expression if "/" not in basis_expression else f"({basis_expression})"
    if eta_product == "1":
        return f"{target_variable} = {basis}"
    return f"{target_variable} = {basis} * {eta_product}"


def _format_two_core_source_family_eta_correction(
    *,
    relation: MultiplicativeRelation,
    target_variable: str,
    series_symbol: str,
) -> str:
    source_terms: list[str] = []
    eta_terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
            eta_terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        source_terms.append(name if exponent == 1 else f"{name}^{exponent}")
    rhs_terms = source_terms + eta_terms
    rhs = " * ".join(rhs_terms) if rhs_terms else "1"
    return f"{target_variable} = {rhs}"


def _format_quotient_core_source_family_eta_correction(
    *,
    quotient_label: str,
    quotient_expression: str,
    raw_label: str,
    raw_expression: str,
    relation: MultiplicativeRelation,
    target_variable: str,
    series_symbol: str,
) -> str:
    source_terms: list[str] = []
    eta_terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
            eta_terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        if name == quotient_label:
            base = f"({quotient_expression})"
        elif name == raw_label:
            base = raw_expression
        else:
            base = name
        source_terms.append(base if exponent == 1 else f"{base}^{exponent}")
    rhs_terms = source_terms + eta_terms
    rhs = " * ".join(rhs_terms) if rhs_terms else "1"
    return f"{target_variable} = {rhs}"


def _format_two_quotient_core_source_family_eta_correction(
    *,
    quotient_labels: tuple[str, str],
    quotient_expressions: tuple[str, str],
    relation: MultiplicativeRelation,
    target_variable: str,
    series_symbol: str,
) -> str:
    expression_by_label = dict(zip(quotient_labels, quotient_expressions))
    source_terms: list[str] = []
    eta_terms: list[str] = []
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
            eta_terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        base = f"({expression_by_label.get(name, name)})"
        source_terms.append(base if exponent == 1 else f"{base}^{exponent}")
    rhs_terms = source_terms + eta_terms
    rhs = " * ".join(rhs_terms) if rhs_terms else "1"
    return f"{target_variable} = {rhs}"


def _format_self_eta_correction(
    *,
    relation: MultiplicativeRelation,
    modulus: int,
    target_variable: str,
    series_symbol: str,
) -> str:
    source_terms: list[str] = []
    eta_terms: list[str] = []
    g_power_label = f"G{modulus}"
    for name in relation.basis_variables:
        exponent = relation.exponents.get(name)
        if exponent is None:
            continue
        if name == g_power_label:
            base = f"{target_variable}({series_symbol}^{modulus})"
            source_terms.append(base if exponent == 1 else f"({base})^{exponent}")
            continue
        if name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
            eta_terms.append(base if exponent == 1 else f"{base}^{exponent}")
    rhs_terms = source_terms + eta_terms
    rhs = " * ".join(rhs_terms) if rhs_terms else "1"
    return f"{target_variable} = {rhs}"


def _self_quotient_product_relation_residual_series(
    relation: SelfQuotientProductRelation,
    *,
    target_series: Series,
    order: int,
) -> Series:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    lhs = series_div(
        target_series[:order],
        benchmark_power_substitution_series(target_series, power=relation.modulus, order=order),
    )
    rhs: Series = [sp.Integer(0) for _ in range(order)]
    rhs[0] = sp.Integer(1)
    for residue, exponent in sorted(relation.exponents_by_residue.items()):
        factor = _signed_series_pow(_one_minus_power_series(power=residue, order=order), exponent)
        rhs = series_mul(rhs, factor)
    return [sp.simplify(lhs[idx] - rhs[idx]) for idx in range(order)]


def _series_match(lhs: Series, rhs: Series, *, order: int) -> bool:
    return all(sp.simplify(lhs[index] - rhs[index]) == 0 for index in range(order))


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


def search_self_t_polynomial_fractional_linear_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_t_degree: int,
) -> SelfTPolynomialFractionalLinearRelation | None:
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if order < 2:
        raise ValueError("order must be at least 2")
    if max_t_degree < 0:
        raise ValueError("max_t_degree must be non-negative")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1")

    g_series = benchmark_power_substitution_series(target_series, power=modulus, order=order)
    g_shifted = _series_subtract_one(g_series)
    t_series = _t_series(order=order)
    t_powers = [series_pow(t_series, degree) for degree in range(max_t_degree + 1)]

    rhs_vector = sp.Matrix([sp.simplify(-target_series[n]) for n in range(1, order)])
    columns: list[list[sp.Expr]] = []

    for degree in range(1, max_t_degree + 1):
        columns.append([sp.simplify(-t_powers[degree][n]) for n in range(1, order)])
    for degree in range(max_t_degree + 1):
        numer_self_term = series_mul(t_powers[degree], g_shifted)
        columns.append([sp.simplify(-numer_self_term[n]) for n in range(1, order)])
    for degree in range(1, max_t_degree + 1):
        denom_t_term = series_mul(target_series[:order], t_powers[degree])
        columns.append([sp.simplify(denom_t_term[n]) for n in range(1, order)])
    for degree in range(max_t_degree + 1):
        denom_self_term = series_mul(series_mul(target_series[:order], t_powers[degree]), g_shifted)
        columns.append([sp.simplify(denom_self_term[n]) for n in range(1, order)])

    num_unknowns = len(columns)
    num_constraints = order - 1
    if num_unknowns > num_constraints:
        raise ValueError(
            "underdetermined self-fractional-linear relation search: "
            f"{num_unknowns} coefficients > {num_constraints} constraints "
            "(increase order or reduce max_t_degree)"
        )

    matrix = sp.Matrix(
        [[column[row_index] for column in columns] for row_index in range(num_constraints)]
    )
    rank = matrix.rank()
    augmented_rank = matrix.row_join(rhs_vector).rank()
    if augmented_rank > rank:
        return None
    if rank < num_unknowns:
        raise ValueError(
            "underdetermined self-fractional-linear relation search: "
            f"rank {rank} < {num_unknowns} coefficients "
            "(increase order or reduce max_t_degree)"
        )

    unknowns = sp.symbols(f"a1:{max_t_degree + 1} b0:{max_t_degree + 1} c1:{max_t_degree + 1} d0:{max_t_degree + 1}")
    solution_set = sp.linsolve((matrix, rhs_vector), unknowns)
    if not solution_set:
        return None
    solution = next(iter(solution_set))
    if any(value.free_symbols for value in solution):
        raise ValueError("self-fractional-linear relation search returned a parametric solution")

    index = 0
    numerator_t_coefficients = [sp.Integer(1)]
    for _ in range(1, max_t_degree + 1):
        numerator_t_coefficients.append(sp.simplify(solution[index]))
        index += 1
    numerator_self_coefficients: list[sp.Expr] = []
    for _ in range(max_t_degree + 1):
        numerator_self_coefficients.append(sp.simplify(solution[index]))
        index += 1
    denominator_t_coefficients = [sp.Integer(1)]
    for _ in range(1, max_t_degree + 1):
        denominator_t_coefficients.append(sp.simplify(solution[index]))
        index += 1
    denominator_self_coefficients: list[sp.Expr] = []
    for _ in range(max_t_degree + 1):
        denominator_self_coefficients.append(sp.simplify(solution[index]))
        index += 1

    if all(value == 0 for value in numerator_self_coefficients + denominator_self_coefficients):
        return None

    relation = SelfTPolynomialFractionalLinearRelation(
        order_checked=order,
        max_t_degree=max_t_degree,
        numerator_t_coefficients=tuple(numerator_t_coefficients),
        numerator_self_coefficients=tuple(numerator_self_coefficients),
        denominator_t_coefficients=tuple(denominator_t_coefficients),
        denominator_self_coefficients=tuple(denominator_self_coefficients),
    )

    residual = _self_t_polynomial_fractional_linear_relation_residual_series(
        relation,
        target_series=target_series,
        modulus=modulus,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return relation


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


def _format_self_t_polynomial_fractional_linear_relation(
    relation: SelfTPolynomialFractionalLinearRelation,
    *,
    modulus: int,
    target_variable: str,
    series_symbol: str,
) -> str:
    def _polynomial(coefficients: tuple[sp.Expr, ...]) -> str:
        terms: list[str] = []
        for degree, coefficient in enumerate(coefficients):
            coefficient_s = sp.simplify(coefficient)
            if coefficient_s == 0:
                continue
            coefficient_str = _format_expr(coefficient_s)
            if degree == 0:
                terms.append(coefficient_str)
                continue
            power = series_symbol if degree == 1 else f"{series_symbol}^{degree}"
            if coefficient_s == 1:
                terms.append(power)
            elif coefficient_s == -1:
                terms.append(f"-{power}")
            else:
                terms.append(f"{coefficient_str}*{power}")
        return " + ".join(terms) if terms else "0"

    g_term = f"{target_variable}({series_symbol}^{modulus}) - 1"
    numerator = _polynomial(relation.numerator_t_coefficients)
    numerator_self = _polynomial(relation.numerator_self_coefficients)
    denominator = _polynomial(relation.denominator_t_coefficients)
    denominator_self = _polynomial(relation.denominator_self_coefficients)
    return (
        f"{target_variable} = "
        f"(({numerator}) + ({numerator_self})*({g_term})) / "
        f"(({denominator}) + ({denominator_self})*({g_term}))"
    )


def _format_source_correction_expression(
    *,
    basis_labels: tuple[str, ...],
    basis_expressions: tuple[str, ...],
    target_variable: str,
) -> str:
    factors: list[str] = []
    for label, expression in zip(basis_labels, basis_expressions):
        base = label if expression == label else f"{label} = {expression}"
        factors.append(base)
    source_product = " * ".join(label for label in basis_labels)
    if len(factors) == 1:
        return f"G = {target_variable} / {source_product}    with    {factors[0]}"
    return (
        f"G = {target_variable} / ({source_product})    with    "
        + ", ".join(factors)
    )


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


def _self_t_polynomial_fractional_linear_relation_residual_series(
    relation: SelfTPolynomialFractionalLinearRelation,
    *,
    target_series: Series,
    modulus: int,
    order: int,
) -> Series:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    g_series = benchmark_power_substitution_series(target_series, power=modulus, order=order)
    g_shifted = _series_subtract_one(g_series)
    t_series = _t_series(order=order)
    t_powers = [series_pow(t_series, degree) for degree in range(relation.max_t_degree + 1)]

    numerator: Series = [sp.Integer(0) for _ in range(order)]
    denominator: Series = [sp.Integer(0) for _ in range(order)]
    for degree, coefficient in enumerate(relation.numerator_t_coefficients):
        term = t_powers[degree]
        for index, value in enumerate(term):
            if value == 0:
                continue
            numerator[index] = sp.simplify(numerator[index] + coefficient * value)
    for degree, coefficient in enumerate(relation.denominator_t_coefficients):
        term = t_powers[degree]
        for index, value in enumerate(term):
            if value == 0:
                continue
            denominator[index] = sp.simplify(denominator[index] + coefficient * value)
    for degree, coefficient in enumerate(relation.numerator_self_coefficients):
        term = series_mul(t_powers[degree], g_shifted)
        for index, value in enumerate(term):
            if value == 0:
                continue
            numerator[index] = sp.simplify(numerator[index] + coefficient * value)
    for degree, coefficient in enumerate(relation.denominator_self_coefficients):
        term = series_mul(t_powers[degree], g_shifted)
        for index, value in enumerate(term):
            if value == 0:
                continue
            denominator[index] = sp.simplify(denominator[index] + coefficient * value)

    lhs = series_mul(target_series[:order], denominator)
    return [sp.simplify(lhs[index] - numerator[index]) for index in range(order)]


def scan_self_fractional_linear_uniqueness_relations(
    *,
    target_series: Series,
    moduli: tuple[int, ...],
    order: int,
    t_degree_values: tuple[int, ...] = (1, 2),
) -> SelfFractionalLinearUniquenessScan:
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_t_degrees = tuple(sorted({degree for degree in t_degree_values if degree >= 0}))
    if not normalized_moduli or not normalized_t_degrees:
        return SelfFractionalLinearUniquenessScan(
            moduli_checked=normalized_moduli,
            t_degree_values=normalized_t_degrees,
            hits=(),
        )

    hits: list[SelfFractionalLinearUniquenessHit] = []
    for modulus in normalized_moduli:
        for t_degree in normalized_t_degrees:
            try:
                relation = search_self_t_polynomial_fractional_linear_relation(
                    target_series=target_series,
                    modulus=modulus,
                    order=order,
                    max_t_degree=t_degree,
                )
            except ValueError:
                continue
            if relation is None:
                continue
            hits.append(
                SelfFractionalLinearUniquenessHit(
                    modulus=modulus,
                    max_t_degree=t_degree,
                    relation=relation,
                )
            )

    return SelfFractionalLinearUniquenessScan(
        moduli_checked=normalized_moduli,
        t_degree_values=normalized_t_degrees,
        hits=tuple(hits),
    )


def search_self_mahler_linear_relation(
    *,
    target_series: Series,
    modulus: int,
    levels: int,
    order: int,
    max_t_degree: int,
) -> PolynomialRelation | None:
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if levels < 2:
        raise ValueError("levels must be at least 2")
    if max_t_degree < 0:
        raise ValueError("max_t_degree must be non-negative")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    variable_names = ["T", "F"]
    series_by_variable: dict[str, Series] = {"T": _t_series(order=order), "F": target_series[:order]}
    current_power = modulus
    for _ in range(levels):
        label = f"G{current_power}"
        variable_names.append(label)
        series_by_variable[label] = benchmark_power_substitution_series(
            target_series,
            power=current_power,
            order=order,
        )
        current_power *= modulus

    exponent_tuples: list[tuple[int, ...]] = []
    num_variables = len(variable_names)
    for t_degree in range(max_t_degree + 1):
        exponent_tuples.append((t_degree,) + (0,) * (num_variables - 1))
        for variable_index in range(1, num_variables):
            exponents = [0] * num_variables
            exponents[0] = t_degree
            exponents[variable_index] = 1
            exponent_tuples.append(tuple(exponents))

    relation = _guess_polynomial_relation_from_exponent_tuples(
        series_by_variable=series_by_variable,
        order=order,
        exponent_tuples=tuple(exponent_tuples),
        required_variables=("F", variable_names[-1]),
    )
    if relation is None:
        return None

    residual = _relation_residual_series(
        relation,
        series_by_variable=series_by_variable,
        order=order,
    )
    if any(sp.simplify(value) != 0 for value in residual):
        return None
    return relation


def scan_self_mahler_linear_relations(
    *,
    target_series: Series,
    moduli: tuple[int, ...],
    levels_checked: tuple[int, ...] = (2,),
    order: int,
    t_degree_values: tuple[int, ...] = (1, 2),
) -> SelfMahlerLinearScan:
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_levels = tuple(sorted({level for level in levels_checked if level >= 2}))
    normalized_t_degrees = tuple(sorted({degree for degree in t_degree_values if degree >= 0}))
    if not normalized_moduli or not normalized_levels or not normalized_t_degrees:
        return SelfMahlerLinearScan(
            moduli_checked=normalized_moduli,
            levels_checked=normalized_levels,
            t_degree_values=normalized_t_degrees,
            hits=(),
        )

    hits: list[SelfMahlerLinearHit] = []
    for modulus in normalized_moduli:
        for levels in normalized_levels:
            for t_degree in normalized_t_degrees:
                try:
                    relation = search_self_mahler_linear_relation(
                        target_series=target_series,
                        modulus=modulus,
                        levels=levels,
                        order=order,
                        max_t_degree=t_degree,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                hits.append(
                    SelfMahlerLinearHit(
                        modulus=modulus,
                        levels=levels,
                        max_t_degree=t_degree,
                        relation=relation,
                    )
                )

    return SelfMahlerLinearScan(
        moduli_checked=normalized_moduli,
        levels_checked=normalized_levels,
        t_degree_values=normalized_t_degrees,
        hits=tuple(hits),
    )


def _source_core_correction_entries(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    correction_size: int,
) -> tuple[tuple[tuple[str, ...], tuple[str, ...], Series], ...]:
    if correction_size < 1:
        raise ValueError("correction_size must be at least 1")
    if correction_size > len(ordered_base_families):
        return ()

    entries: list[tuple[tuple[str, ...], tuple[str, ...], Series]] = []
    if correction_size == 1:
        for label, benchmark_name, basis_series in ordered_base_families:
            entries.append(((label,), (benchmark_name,), series_div(target_series, basis_series)))
        return tuple(entries)

    if correction_size == 2:
        for left_index, left_entry in enumerate(ordered_base_families):
            left_label, left_name, left_series = left_entry
            for right_entry in ordered_base_families[left_index + 1 :]:
                right_label, right_name, right_series = right_entry
                entries.append(
                    (
                        (left_label, right_label),
                        (left_name, right_name),
                        series_div(target_series, series_mul(left_series, right_series)),
                    )
                )
        return tuple(entries)

    raise ValueError("supported correction sizes are 1 and 2")


def scan_source_correction_self_polynomial_uniqueness_relations(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    correction_size: int,
    moduli: tuple[int, ...],
    order: int,
    fg_degree_values: tuple[int, ...] = (1, 2),
    t_degree_values: tuple[int, ...] = (1, 2),
) -> SourceCorrectionSelfPolynomialScan:
    correction_entries = _source_core_correction_entries(
        target_series=target_series,
        ordered_base_families=ordered_base_families,
        correction_size=correction_size,
    )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_fg_degrees = tuple(sorted({degree for degree in fg_degree_values if degree >= 1}))
    normalized_t_degrees = tuple(sorted({degree for degree in t_degree_values if degree >= 0}))
    if not correction_entries or not normalized_moduli or not normalized_fg_degrees or not normalized_t_degrees:
        return SourceCorrectionSelfPolynomialScan(
            correction_size=correction_size,
            moduli_checked=normalized_moduli,
            fg_degree_values=normalized_fg_degrees,
            t_degree_values=normalized_t_degrees,
            total_corrections_checked=len(correction_entries),
            hits=(),
        )

    hits: list[SourceCorrectionSelfPolynomialHit] = []
    for basis_labels, basis_expressions, correction_series in correction_entries:
        for modulus in normalized_moduli:
            for fg_degree in normalized_fg_degrees:
                for t_degree in normalized_t_degrees:
                    try:
                        relation = search_self_polynomial_uniqueness_relation(
                            target_series=correction_series,
                            modulus=modulus,
                            order=order,
                            max_fg_total_degree=fg_degree,
                            max_t_degree=t_degree,
                        )
                    except ValueError:
                        continue
                    if relation is None:
                        continue
                    hits.append(
                        SourceCorrectionSelfPolynomialHit(
                            basis_labels=basis_labels,
                            basis_expressions=basis_expressions,
                            modulus=modulus,
                            max_fg_total_degree=fg_degree,
                            max_t_degree=t_degree,
                            relation=relation,
                        )
                    )

    return SourceCorrectionSelfPolynomialScan(
        correction_size=correction_size,
        moduli_checked=normalized_moduli,
        fg_degree_values=normalized_fg_degrees,
        t_degree_values=normalized_t_degrees,
        total_corrections_checked=len(correction_entries),
        hits=tuple(hits),
    )


def scan_source_correction_self_fractional_linear_uniqueness_relations(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    correction_size: int,
    moduli: tuple[int, ...],
    order: int,
    t_degree_values: tuple[int, ...] = (1, 2),
) -> SourceCorrectionSelfFractionalLinearScan:
    correction_entries = _source_core_correction_entries(
        target_series=target_series,
        ordered_base_families=ordered_base_families,
        correction_size=correction_size,
    )
    normalized_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_t_degrees = tuple(sorted({degree for degree in t_degree_values if degree >= 0}))
    if not correction_entries or not normalized_moduli or not normalized_t_degrees:
        return SourceCorrectionSelfFractionalLinearScan(
            correction_size=correction_size,
            moduli_checked=normalized_moduli,
            t_degree_values=normalized_t_degrees,
            total_corrections_checked=len(correction_entries),
            hits=(),
        )

    hits: list[SourceCorrectionSelfFractionalLinearHit] = []
    for basis_labels, basis_expressions, correction_series in correction_entries:
        for modulus in normalized_moduli:
            for t_degree in normalized_t_degrees:
                try:
                    relation = search_self_t_polynomial_fractional_linear_relation(
                        target_series=correction_series,
                        modulus=modulus,
                        order=order,
                        max_t_degree=t_degree,
                    )
                except ValueError:
                    continue
                if relation is None:
                    continue
                hits.append(
                    SourceCorrectionSelfFractionalLinearHit(
                        basis_labels=basis_labels,
                        basis_expressions=basis_expressions,
                        modulus=modulus,
                        max_t_degree=t_degree,
                        relation=relation,
                    )
                )

    return SourceCorrectionSelfFractionalLinearScan(
        correction_size=correction_size,
        moduli_checked=normalized_moduli,
        t_degree_values=normalized_t_degrees,
        total_corrections_checked=len(correction_entries),
        hits=tuple(hits),
    )


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


def scan_named_multiplicative_prefixes(
    *,
    target_series: Series,
    ordered_basis_series: tuple[tuple[str, Series], ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[NamedMultiplicativeRelationScan]:
    if not ordered_basis_series:
        return []

    scans: list[NamedMultiplicativeRelationScan] = []
    prefix: list[tuple[str, Series]] = []
    for label, series in ordered_basis_series:
        prefix.append((label, series))
        basis_series_by_variable = {name: basis for name, basis in prefix}
        try:
            relation = search_multiplicative_relation(
                target_series=target_series,
                basis_series_by_variable=basis_series_by_variable,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(
                NamedMultiplicativeRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relation=relation,
                )
            )
        except ValueError as exc:
            scans.append(
                NamedMultiplicativeRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def scan_named_fractional_linear_prefixes(
    *,
    target_series: Series,
    ordered_basis_series: tuple[tuple[str, Series], ...],
    order: int,
) -> list[NamedFractionalLinearRelationScan]:
    if not ordered_basis_series:
        return []

    scans: list[NamedFractionalLinearRelationScan] = []
    prefix: list[tuple[str, Series]] = []
    for label, series in ordered_basis_series:
        prefix.append((label, series))
        basis_series_by_variable = {name: basis for name, basis in prefix}
        try:
            relation = search_fractional_linear_relation(
                target_series=target_series,
                basis_series_by_variable=basis_series_by_variable,
                order=order,
            )
            scans.append(
                NamedFractionalLinearRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relation=relation,
                )
            )
        except ValueError as exc:
            scans.append(
                NamedFractionalLinearRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def scan_named_two_layer_fractional_linear_prefixes(
    *,
    target_series: Series,
    ordered_basis_series: tuple[tuple[str, Series], ...],
    order: int,
    solve_order: int | None = None,
    max_reported_hits: int = 3,
) -> list[NamedTwoLayerFractionalLinearRelationScan]:
    if len(ordered_basis_series) < 2:
        return []

    if max_reported_hits < 1:
        raise ValueError("max_reported_hits must be at least 1")

    scans: list[NamedTwoLayerFractionalLinearRelationScan] = []
    prefix: list[tuple[str, Series]] = []
    for label, series in ordered_basis_series:
        prefix.append((label, series))
        if len(prefix) < 2:
            continue

        basis_series_by_variable = {name: basis for name, basis in prefix}
        basis_names = tuple(basis_series_by_variable.keys())
        hits: list[TwoLayerFractionalLinearRelation] = []
        seen_signatures: set[str] = set()
        total_hits = 0
        factor_pairs = _two_layer_factor_index_pairs(len(basis_names))
        tuples_checked = len(factor_pairs)
        try:
            for factor_1, factor_2 in factor_pairs:
                numerator_variables = (
                    basis_names[factor_1[0]],
                    basis_names[factor_2[0]],
                )
                denominator_variables = (
                    basis_names[factor_1[1]],
                    basis_names[factor_2[1]],
                )
                try:
                    relation = search_two_layer_fractional_linear_relation(
                        target_series=target_series,
                        basis_series_by_variable=basis_series_by_variable,
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
                NamedTwoLayerFractionalLinearRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relations=tuple(hits),
                    total_hits=total_hits,
                    tuples_checked=tuples_checked,
                )
            )
        except ValueError as exc:
            scans.append(
                NamedTwoLayerFractionalLinearRelationScan(
                    basis_labels=tuple(name for name, _ in prefix),
                    relations=(),
                    total_hits=0,
                    tuples_checked=tuples_checked,
                    error=str(exc),
                )
            )
    return scans


def scan_named_polynomial_prefixes(
    *,
    target_series: Series,
    ordered_basis_series: tuple[tuple[str, Series], ...],
    order: int,
    degree_values: tuple[int, ...],
    required_variable: str | None = "F",
) -> list[NamedPolynomialRelationScan]:
    if not ordered_basis_series:
        return []

    scans: list[NamedPolynomialRelationScan] = []
    prefix: list[tuple[str, Series]] = []
    degrees = tuple(sorted({degree for degree in degree_values if degree >= 1}))
    for label, basis_series in ordered_basis_series:
        prefix.append((label, basis_series))
        variables = {"F": target_series}
        variables.update({name: series for name, series in prefix})
        for degree in degrees:
            try:
                relation = search_polynomial_relation(
                    series_by_variable=variables,
                    order=order,
                    max_total_degree=degree,
                    required_variable=required_variable,
                )
                scans.append(
                    NamedPolynomialRelationScan(
                        basis_labels=tuple(name for name, _ in prefix),
                        max_total_degree=degree,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                scans.append(
                    NamedPolynomialRelationScan(
                        basis_labels=tuple(name for name, _ in prefix),
                        max_total_degree=degree,
                        relation=None,
                        error=str(exc),
                    )
                )
    return scans


def scan_named_prefix_boxes(
    *,
    target_series: Series,
    ordered_basis_series: tuple[tuple[str, Series], ...],
    order: int,
    degree_values: tuple[int, ...] = (),
    required_variable: str | None = "F",
    max_abs_exponent: int = 8,
    solve_order: int | None = None,
    max_reported_two_layer_hits: int = 3,
    include_polynomial: bool = True,
    include_multiplicative: bool = True,
    include_fractional_linear: bool = True,
    include_two_layer: bool = True,
) -> NamedPrefixBoxScans:
    if not ordered_basis_series:
        return NamedPrefixBoxScans()

    if include_two_layer and max_reported_two_layer_hits < 1:
        raise ValueError("max_reported_two_layer_hits must be at least 1")

    polynomial_scans: list[NamedPolynomialRelationScan] = []
    multiplicative_scans: list[NamedMultiplicativeRelationScan] = []
    fractional_linear_scans: list[NamedFractionalLinearRelationScan] = []
    two_layer_fractional_linear_scans: list[NamedTwoLayerFractionalLinearRelationScan] = []
    prefix: list[tuple[str, Series]] = []
    degrees = tuple(sorted({degree for degree in degree_values if degree >= 1}))

    for label, basis_series in ordered_basis_series:
        prefix.append((label, basis_series))
        basis_labels = tuple(name for name, _ in prefix)
        basis_series_by_variable = {name: series for name, series in prefix}

        if include_polynomial:
            variables = {"F": target_series}
            variables.update(basis_series_by_variable)
            for degree in degrees:
                try:
                    relation = search_polynomial_relation(
                        series_by_variable=variables,
                        order=order,
                        max_total_degree=degree,
                        required_variable=required_variable,
                    )
                    polynomial_scans.append(
                        NamedPolynomialRelationScan(
                            basis_labels=basis_labels,
                            max_total_degree=degree,
                            relation=relation,
                        )
                    )
                except ValueError as exc:
                    polynomial_scans.append(
                        NamedPolynomialRelationScan(
                            basis_labels=basis_labels,
                            max_total_degree=degree,
                            relation=None,
                            error=str(exc),
                        )
                    )

        if include_multiplicative:
            try:
                relation = search_multiplicative_relation(
                    target_series=target_series,
                    basis_series_by_variable=basis_series_by_variable,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
                multiplicative_scans.append(
                    NamedMultiplicativeRelationScan(
                        basis_labels=basis_labels,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                multiplicative_scans.append(
                    NamedMultiplicativeRelationScan(
                        basis_labels=basis_labels,
                        relation=None,
                        error=str(exc),
                    )
                )

        if include_fractional_linear:
            try:
                relation = search_fractional_linear_relation(
                    target_series=target_series,
                    basis_series_by_variable=basis_series_by_variable,
                    order=order,
                )
                fractional_linear_scans.append(
                    NamedFractionalLinearRelationScan(
                        basis_labels=basis_labels,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                fractional_linear_scans.append(
                    NamedFractionalLinearRelationScan(
                        basis_labels=basis_labels,
                        relation=None,
                        error=str(exc),
                    )
                )

        if not include_two_layer or len(prefix) < 2:
            continue

        basis_names = tuple(basis_series_by_variable.keys())
        hits: list[TwoLayerFractionalLinearRelation] = []
        seen_signatures: set[str] = set()
        total_hits = 0
        factor_pairs = _two_layer_factor_index_pairs(len(basis_names))
        tuples_checked = len(factor_pairs)
        try:
            for factor_1, factor_2 in factor_pairs:
                numerator_variables = (
                    basis_names[factor_1[0]],
                    basis_names[factor_2[0]],
                )
                denominator_variables = (
                    basis_names[factor_1[1]],
                    basis_names[factor_2[1]],
                )
                try:
                    relation = search_two_layer_fractional_linear_relation(
                        target_series=target_series,
                        basis_series_by_variable=basis_series_by_variable,
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
                if len(hits) < max_reported_two_layer_hits:
                    hits.append(relation)

            two_layer_fractional_linear_scans.append(
                NamedTwoLayerFractionalLinearRelationScan(
                    basis_labels=basis_labels,
                    relations=tuple(hits),
                    total_hits=total_hits,
                    tuples_checked=tuples_checked,
                )
            )
        except ValueError as exc:
            two_layer_fractional_linear_scans.append(
                NamedTwoLayerFractionalLinearRelationScan(
                    basis_labels=basis_labels,
                    relations=(),
                    total_hits=0,
                    tuples_checked=tuples_checked,
                    error=str(exc),
                )
            )

    return NamedPrefixBoxScans(
        polynomial_scans=tuple(polynomial_scans),
        multiplicative_scans=tuple(multiplicative_scans),
        fractional_linear_scans=tuple(fractional_linear_scans),
        two_layer_fractional_linear_scans=tuple(two_layer_fractional_linear_scans),
    )


def scan_parameterized_source_family_power_boxes(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    order: int,
    degree_values: tuple[int, ...] = (1, 2),
    max_abs_exponent: int = 8,
    solve_order: int | None = None,
    max_reported_two_layer_hits: int = 3,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> list[ParameterizedSourceFamilyScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not ordered_base_families:
        return []

    scans: list[ParameterizedSourceFamilyScan] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        family_powers = set(unique_powers)
        if supplemental_powers_by_family is not None:
            family_powers.update(
                power
                for power in supplemental_powers_by_family.get(family_label, ())
                if power >= 2
            )
        ordered_basis_series: list[tuple[str, Series]] = [(family_label, base_series)]
        for power in tuple(sorted(family_powers)):
            ordered_basis_series.append(
                (
                    f"{family_label}{power}",
                    benchmark_power_substitution_series(base_series, power=power, order=order),
                )
            )
        ordered_basis_tuple = tuple(ordered_basis_series)

        quotient_basis_series: list[tuple[str, str, Series]] = []
        for label, basis_series in ordered_basis_series[1:]:
            power = int(label.removeprefix(family_label))
            quotient_basis_series.append(
                (
                    f"Q{power}",
                    f"{label} / {family_label}",
                    series_div(basis_series, base_series),
                )
            )
        quotient_basis_tuple = tuple(quotient_basis_series)
        quotient_ordered_basis_tuple = tuple((label, series) for label, _, series in quotient_basis_tuple)
        mixed_quotient_basis_series = [(family_label, family_label, base_series)]
        mixed_quotient_basis_series.extend(quotient_basis_series)
        mixed_quotient_basis_tuple = tuple(mixed_quotient_basis_series)
        mixed_quotient_ordered_basis_tuple = tuple((label, series) for label, _, series in mixed_quotient_basis_tuple)
        direct_scans = scan_named_prefix_boxes(
            target_series=target_series,
            ordered_basis_series=ordered_basis_tuple,
            order=order,
            degree_values=degree_values,
            max_abs_exponent=max_abs_exponent,
            solve_order=solve_order,
            max_reported_two_layer_hits=max_reported_two_layer_hits,
        )
        quotient_scans = scan_named_prefix_boxes(
            target_series=target_series,
            ordered_basis_series=quotient_ordered_basis_tuple,
            order=order,
            degree_values=degree_values,
            max_abs_exponent=max_abs_exponent,
            solve_order=solve_order,
            max_reported_two_layer_hits=max_reported_two_layer_hits,
        )
        mixed_quotient_scans = scan_named_prefix_boxes(
            target_series=target_series,
            ordered_basis_series=mixed_quotient_ordered_basis_tuple,
            order=order,
            degree_values=degree_values,
            max_abs_exponent=max_abs_exponent,
            solve_order=solve_order,
            max_reported_two_layer_hits=max_reported_two_layer_hits,
        )

        scans.append(
            ParameterizedSourceFamilyScan(
                family_label=family_label,
                benchmark_name=benchmark_name,
                ordered_basis_series=ordered_basis_tuple,
                polynomial_scans=direct_scans.polynomial_scans,
                multiplicative_scans=direct_scans.multiplicative_scans,
                fractional_linear_scans=direct_scans.fractional_linear_scans,
                two_layer_fractional_linear_scans=direct_scans.two_layer_fractional_linear_scans,
                quotient_basis_series=quotient_basis_tuple,
                quotient_polynomial_scans=quotient_scans.polynomial_scans,
                quotient_multiplicative_scans=quotient_scans.multiplicative_scans,
                quotient_fractional_linear_scans=quotient_scans.fractional_linear_scans,
                quotient_two_layer_fractional_linear_scans=quotient_scans.two_layer_fractional_linear_scans,
                mixed_quotient_basis_series=mixed_quotient_basis_tuple,
                mixed_quotient_polynomial_scans=mixed_quotient_scans.polynomial_scans,
                mixed_quotient_multiplicative_scans=mixed_quotient_scans.multiplicative_scans,
                mixed_quotient_fractional_linear_scans=mixed_quotient_scans.fractional_linear_scans,
                mixed_quotient_two_layer_fractional_linear_scans=mixed_quotient_scans.two_layer_fractional_linear_scans,
            )
        )
    return scans


def _explicit_source_family_ordered_basis_series(
    *,
    family_label: str,
    base_series: Series,
    powers: tuple[int, ...],
    order: int,
    supplemental_powers: tuple[int, ...] = (),
) -> tuple[tuple[str, Series], ...]:
    family_powers = set(power for power in powers if power >= 2)
    family_powers.update(power for power in supplemental_powers if power >= 2)
    ordered_basis_series: list[tuple[str, Series]] = [(family_label, base_series)]
    for power in tuple(sorted(family_powers)):
        ordered_basis_series.append(
            (
                f"{family_label}{power}",
                benchmark_power_substitution_series(base_series, power=power, order=order),
            )
        )
    return tuple(ordered_basis_series)


def _gg_modular_equation_ordered_basis_series(
    *,
    base_series: Series,
    order: int,
    supplemental_powers: tuple[int, ...] = (),
) -> tuple[tuple[str, str, Series], ...]:
    entries: list[tuple[str, str, Series]] = [
        ("GG", "GG(t)", base_series),
        ("GGneg", "GG(-t)", signed_argument_substitution_series(base_series, order=order)),
        ("GG2", "GG(t^2)", benchmark_power_substitution_series(base_series, power=2, order=order)),
        ("GG3", "GG(t^3)", benchmark_power_substitution_series(base_series, power=3, order=order)),
        ("GG4", "GG(t^4)", benchmark_power_substitution_series(base_series, power=4, order=order)),
    ]
    for power in tuple(sorted({value for value in supplemental_powers if value >= 5})):
        entries.append(
            (
                f"GG{power}",
                f"GG(t^{power})",
                benchmark_power_substitution_series(base_series, power=power, order=order),
            )
        )
    return tuple(entries)


@lru_cache(maxsize=None)
def _gg_exact_modular_relations(order: int) -> tuple[tuple[str, PolynomialRelation], ...]:
    u, v, t = sp.symbols("F G T")
    direct_level_3_expr = (
        u**4 * v**3 * t**5
        - 3 * u**3 * v**2 * t**3
        + u**3
        - 3 * u**2 * v**3 * t**4
        + 3 * u**2 * v * t
        - u * v**4 * t**5
        + 3 * u * v**2 * t**2
        - v
    )
    direct_level_3 = _polynomial_relation_from_sympy_expr(
        expr=direct_level_3_expr,
        variables=("F", "G", "T"),
        order=order,
    )
    w = sp.Symbol("G")
    direct_level_4_expr = (
        u**4 * w**3 * t**6
        + u**4 * w**2 * t**4
        + u**4 * w * t**2
        + u**4
        - 4 * u**2 * w**3 * t**5
        + 4 * u**2 * w * t
        + w**4 * t**6
        - w**3 * t**4
        + w**2 * t**2
        - w
    )
    direct_level_4 = _polynomial_relation_from_sympy_expr(
        expr=direct_level_4_expr,
        variables=("F", "G", "T"),
        order=order,
    )
    q = sp.Symbol("Q")
    quotient_level_3 = _polynomial_relation_from_sympy_expr(
        expr=sp.expand(direct_level_3_expr.subs(v, u * q)),
        variables=("F", "Q", "T"),
        order=order,
    )
    quotient_level_4 = _polynomial_relation_from_sympy_expr(
        expr=sp.expand(direct_level_4_expr.subs(w, u * q)),
        variables=("F", "Q", "T"),
        order=order,
    )
    return (
        ("Chan--Huang Cor. 3.2(i) on (F, G3)", direct_level_3),
        ("Chan--Huang Cor. 3.2(ii) on (F, G4)", direct_level_4),
        ("Chan--Huang Cor. 3.2(i) on (F, Q_3)", quotient_level_3),
        ("Chan--Huang Cor. 3.2(ii) on (F, Q_4)", quotient_level_4),
    )


def _gg_exact_modular_template_hits(
    *,
    target_series: Series,
    ordered_basis_entries: tuple[tuple[str, str, Series], ...],
    quotient_basis_entries: tuple[tuple[str, str, Series], ...],
    order: int,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    basis_by_label = {label: series for label, _, series in ordered_basis_entries}
    quotient_by_label = {label: series for label, _, series in quotient_basis_entries}
    relation_map = dict(_gg_exact_modular_relations(order))
    variable_series = [sp.Integer(0) for _ in range(order)]
    if order > 1:
        variable_series[1] = sp.Integer(1)

    exact_labels: list[str] = []
    exact_hits: list[str] = []
    if "GG3" in basis_by_label:
        label = "Chan--Huang Cor. 3.2(i) on (F, GG3)"
        exact_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(i) on (F, G3)"],
            series_by_variable={"F": target_series, "G": basis_by_label["GG3"], "T": variable_series},
            order=order,
        )
        if all(sp.simplify(value) == 0 for value in residual):
            exact_hits.append(label)
    if "GG4" in basis_by_label:
        label = "Chan--Huang Cor. 3.2(ii) on (F, GG4)"
        exact_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(ii) on (F, G4)"],
            series_by_variable={"F": target_series, "G": basis_by_label["GG4"], "T": variable_series},
            order=order,
        )
        if all(sp.simplify(value) == 0 for value in residual):
            exact_hits.append(label)

    quotient_labels: list[str] = []
    quotient_hits: list[str] = []
    if "Q_3" in quotient_by_label:
        label = "Chan--Huang Cor. 3.2(i) on (F, Q_3)"
        quotient_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(i) on (F, Q_3)"],
            series_by_variable={"F": target_series, "Q": quotient_by_label["Q_3"], "T": variable_series},
            order=order,
        )
        if all(sp.simplify(value) == 0 for value in residual):
            quotient_hits.append(label)
    if "Q_4" in quotient_by_label:
        label = "Chan--Huang Cor. 3.2(ii) on (F, Q_4)"
        quotient_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(ii) on (F, Q_4)"],
            series_by_variable={"F": target_series, "Q": quotient_by_label["Q_4"], "T": variable_series},
            order=order,
        )
        if all(sp.simplify(value) == 0 for value in residual):
            quotient_hits.append(label)

    return (
        tuple(exact_labels),
        tuple(exact_hits),
        tuple(quotient_labels),
        tuple(quotient_hits),
    )


def _explicit_source_family_template_series(
    ordered_basis_series: tuple[tuple[str, Series], ...],
) -> tuple[tuple[str, Series], ...]:
    templates: list[tuple[str, Series]] = []
    for label, basis_series in ordered_basis_series:
        templates.append((label, basis_series))
        templates.append((f"1 / {label}", series_invert(basis_series)))

    for numerator_label, numerator_series in ordered_basis_series:
        for denominator_label, denominator_series in ordered_basis_series:
            if numerator_label == denominator_label:
                continue
            templates.append(
                (
                    f"{numerator_label} / {denominator_label}",
                    series_div(numerator_series, denominator_series),
                )
            )
    return tuple(templates)


def scan_explicit_source_family_transform_templates(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    order: int,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> list[ExplicitSourceFamilyTransformScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    if not ordered_base_families:
        return []

    scans: list[ExplicitSourceFamilyTransformScan] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        ordered_basis_series = _explicit_source_family_ordered_basis_series(
            family_label=family_label,
            base_series=base_series,
            powers=unique_powers,
            order=order,
            supplemental_powers=()
            if supplemental_powers_by_family is None
            else supplemental_powers_by_family.get(family_label, ()),
        )

        checked_templates: list[str] = []
        hit_templates: list[str] = []
        for template_label, template_series in _explicit_source_family_template_series(ordered_basis_series):
            checked_templates.append(template_label)
            if _series_match(template_series, target_series, order=order):
                hit_templates.append(template_label)

        scans.append(
            ExplicitSourceFamilyTransformScan(
                family_label=family_label,
                benchmark_name=benchmark_name,
                ordered_basis_series=ordered_basis_series,
                checked_templates=tuple(checked_templates),
                hit_templates=tuple(hit_templates),
            )
        )
    return scans


def scan_explicit_source_family_eta_correction_templates(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> list[ExplicitSourceFamilyEtaCorrectionScan]:
    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not ordered_base_families or not normalized_eta_levels:
        return []

    scans: list[ExplicitSourceFamilyEtaCorrectionScan] = []
    for family_label, benchmark_name, base_series in ordered_base_families:
        ordered_basis_series = _explicit_source_family_ordered_basis_series(
            family_label=family_label,
            base_series=base_series,
            powers=unique_powers,
            order=order,
            supplemental_powers=()
            if supplemental_powers_by_family is None
            else supplemental_powers_by_family.get(family_label, ()),
        )
        checked_templates: list[str] = []
        hits: list[ExplicitTransformEtaCorrectionHit] = []
        for template_label, template_series in _explicit_source_family_template_series(ordered_basis_series):
            checked_templates.append(template_label)
            correction_series = series_div(target_series, template_series)
            for eta_scan in scan_ratio_eta_quotient_relations(
                ratio_series=correction_series,
                levels=normalized_eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            ):
                if eta_scan.relation is None:
                    continue
                hits.append(
                    ExplicitTransformEtaCorrectionHit(
                        template_label=template_label,
                        level=eta_scan.level,
                        relation=eta_scan.relation,
                    )
                )

        scans.append(
            ExplicitSourceFamilyEtaCorrectionScan(
                family_label=family_label,
                benchmark_name=benchmark_name,
                ordered_basis_series=ordered_basis_series,
                checked_templates=tuple(checked_templates),
                hits=tuple(hits),
            )
        )
    return scans


def scan_gg_modular_equation_box(
    *,
    target_series: Series,
    benchmark_name: str,
    gg_series: Series,
    order: int,
    degree_values: tuple[int, ...] = (1, 2),
    max_abs_exponent: int = 8,
    solve_order: int | None = None,
    supplemental_powers: tuple[int, ...] = (),
) -> GGModularEquationScan:
    ordered_basis_entries = _gg_modular_equation_ordered_basis_series(
        base_series=gg_series,
        order=order,
        supplemental_powers=supplemental_powers,
    )
    ordered_basis_series = tuple((label, series) for label, _, series in ordered_basis_entries)
    quotient_basis_series = _gg_modular_equation_quotient_basis_series(ordered_basis_entries)
    quotient_ordered_basis_series = tuple((label, series) for label, _, series in quotient_basis_series)
    mixed_quotient_basis_series = (
        ((ordered_basis_entries[0][0], ordered_basis_entries[0][1], ordered_basis_entries[0][2]),)
        + quotient_basis_series
    )
    mixed_quotient_ordered_basis_series = tuple(
        (label, series) for label, _, series in mixed_quotient_basis_series
    )
    direct_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=ordered_basis_series,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
    )
    quotient_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=quotient_ordered_basis_series,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
    )
    mixed_quotient_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=mixed_quotient_ordered_basis_series,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
    )
    (
        exact_polynomial_template_labels,
        exact_polynomial_template_hits,
        quotient_exact_polynomial_template_labels,
        quotient_exact_polynomial_template_hits,
    ) = _gg_exact_modular_template_hits(
        target_series=target_series,
        ordered_basis_entries=ordered_basis_entries,
        quotient_basis_entries=quotient_basis_series,
        order=order,
    )

    checked_templates: list[str] = []
    hit_templates: list[str] = []
    for template_label, template_series in _explicit_source_family_template_series(ordered_basis_series):
        checked_templates.append(template_label)
        if _series_match(template_series, target_series, order=order):
            hit_templates.append(template_label)

    return GGModularEquationScan(
        benchmark_name=benchmark_name,
        ordered_basis_series=ordered_basis_entries,
        checked_templates=tuple(checked_templates),
        hit_templates=tuple(hit_templates),
        exact_polynomial_template_labels=exact_polynomial_template_labels,
        exact_polynomial_template_hits=exact_polynomial_template_hits,
        polynomial_scans=direct_scans.polynomial_scans,
        multiplicative_scans=direct_scans.multiplicative_scans,
        fractional_linear_scans=direct_scans.fractional_linear_scans,
        two_layer_fractional_linear_scans=direct_scans.two_layer_fractional_linear_scans,
        quotient_basis_series=quotient_basis_series,
        quotient_exact_polynomial_template_labels=quotient_exact_polynomial_template_labels,
        quotient_exact_polynomial_template_hits=quotient_exact_polynomial_template_hits,
        quotient_polynomial_scans=quotient_scans.polynomial_scans,
        quotient_multiplicative_scans=quotient_scans.multiplicative_scans,
        quotient_fractional_linear_scans=quotient_scans.fractional_linear_scans,
        quotient_two_layer_fractional_linear_scans=quotient_scans.two_layer_fractional_linear_scans,
        mixed_quotient_basis_series=mixed_quotient_basis_series,
        mixed_quotient_polynomial_scans=mixed_quotient_scans.polynomial_scans,
        mixed_quotient_multiplicative_scans=mixed_quotient_scans.multiplicative_scans,
        mixed_quotient_fractional_linear_scans=mixed_quotient_scans.fractional_linear_scans,
        mixed_quotient_two_layer_fractional_linear_scans=mixed_quotient_scans.two_layer_fractional_linear_scans,
    )


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
        factor_pairs = _two_layer_factor_index_pairs(len(basis_names))
        tuples_checked = len(factor_pairs)
        try:
            for factor_1, factor_2 in factor_pairs:
                numerator_variables = (
                    basis_names[factor_1[0]],
                    basis_names[factor_2[0]],
                )
                denominator_variables = (
                    basis_names[factor_1[1]],
                    basis_names[factor_2[1]],
                )
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


def _closest_source_family_base_name(closest_benchmark: str) -> str | None:
    if closest_benchmark.startswith("rogers_ramanujan"):
        return "rogers_ramanujan_normalized"
    if closest_benchmark.startswith("ramanujan_cubic"):
        return "ramanujan_cubic_normalized"
    if closest_benchmark.startswith("gollnitz_gordon"):
        return "gollnitz_gordon_normalized"
    if closest_benchmark.startswith("hirschhorn_s"):
        return "hirschhorn_s_normalized"
    return None


def _source_family_basis_catalog(closest_benchmark: str) -> tuple[tuple[str, str], ...]:
    catalog: list[tuple[str, str]] = [
        ("RR", "rogers_ramanujan_normalized"),
        ("cubic", "ramanujan_cubic_normalized"),
        ("GG", "gollnitz_gordon_normalized"),
        ("S", "hirschhorn_s_normalized"),
    ]
    preferred = _closest_source_family_base_name(closest_benchmark)
    if preferred is not None:
        for index, (_, benchmark_name) in enumerate(catalog):
            if benchmark_name == preferred:
                catalog.insert(0, catalog.pop(index))
                break
    return tuple(catalog)


def _parameterized_source_family_powers(
    benchmark_powers: tuple[int, ...],
    *,
    smoke: bool,
) -> tuple[int, ...]:
    preferred = tuple(sorted({power for power in benchmark_powers if 2 <= power <= 4}))
    if preferred:
        return preferred[:2] if smoke else preferred
    return (2,) if smoke else (2, 3, 4)


def _supplemental_source_family_powers(*, smoke: bool) -> dict[str, tuple[int, ...]]:
    if smoke:
        return {}
    return {"GG": (5, 7, 11)}


def _eta_scan_levels(levels: tuple[int, ...]) -> tuple[int, ...]:
    normalized = {1}
    normalized.update(level for level in levels if level >= 1)
    return tuple(sorted(normalized))


def _rhs_uniqueness_moduli(levels: tuple[int, ...], *, smoke: bool) -> tuple[int, ...]:
    upper = 4 if smoke else 6
    preferred = tuple(sorted({level for level in levels if 2 <= level <= upper}))
    if preferred:
        return preferred
    return (2, 3, 4) if smoke else (2, 3, 4, 5, 6)


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
        profile_order = min(profile_order, 48)
        profile_degree = min(profile_degree, 3)
        profile_depth = min(profile_depth, 24)

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
    output_file = Path(output_path)
    source_family_scan_powers = _parameterized_source_family_powers(benchmark_powers, smoke=smoke)
    supplemental_source_family_powers = _supplemental_source_family_powers(smoke=smoke)
    eta_scan_levels = _eta_scan_levels(benchmark_powers)
    rhs_uniqueness_moduli = _rhs_uniqueness_moduli(benchmark_powers, smoke=smoke)
    progress_steps = [
        "series-and-benchmark-setup",
        "rhs-uniqueness-search",
        "source-family-scans",
        "cross-family-functional-scans",
        "explicit-gg-family-scans",
        "benchmark-tower-scans",
        "final-render",
    ]
    progress_status = {step_name: "pending" for step_name in progress_steps}
    build_started_at = perf_counter()
    stage_elapsed_seconds = {step_name: None for step_name in progress_steps}
    current_stage_start = build_started_at

    def advance_progress(*, completed_step: str | None = None, current_step: str) -> None:
        nonlocal current_stage_start
        now = perf_counter()
        if completed_step is not None:
            stage_elapsed_seconds[completed_step] = now - current_stage_start
            progress_status[completed_step] = "completed"
        progress_status[current_step] = "in_progress"
        current_stage_start = now
        write_progress(current_step=current_step)

    def write_progress(*, current_step: str) -> None:
        now = perf_counter()
        progress_lines = [
            f"# Identification Note Build In Progress: `{record.id}`",
            "",
            f"- Status: `in_progress`",
            f"- Current step: `{current_step}`",
            f"- Output target: `{output_path}`",
            f"- Elapsed build seconds: `{now - build_started_at:.2f}`",
            "",
            "## Progress",
            "",
        ]
        for step_name in progress_steps:
            elapsed = stage_elapsed_seconds[step_name]
            if progress_status[step_name] == "in_progress":
                elapsed_text = f"{now - current_stage_start:.2f}"
            elif elapsed is not None:
                elapsed_text = f"{elapsed:.2f}"
            else:
                elapsed_text = "-"
            progress_lines.append(
                f"- `{step_name}`: `{progress_status[step_name]}` (elapsed seconds: `{elapsed_text}`)"
            )
        progress_lines.append("")
        output_file.write_text("\n".join(progress_lines), encoding="utf-8")

    stage_elapsed_seconds["series-and-benchmark-setup"] = perf_counter() - build_started_at
    progress_status["series-and-benchmark-setup"] = "completed"
    progress_status["rhs-uniqueness-search"] = "in_progress"
    current_stage_start = perf_counter()
    write_progress(current_step="rhs-uniqueness-search")

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
    rhs_self_polynomial_scan = scan_self_polynomial_uniqueness_relations(
        target_series=ratio_series,
        moduli=rhs_uniqueness_moduli,
        order=profile_order,
        fg_degree_values=(1, 2) if smoke else (1, 2, 3),
        t_degree_values=(1, 2) if smoke else (1, 2, 3),
    )
    rhs_self_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
        target_series=ratio_series,
        moduli=rhs_uniqueness_moduli,
        order=profile_order,
        t_degree_values=(1, 2) if smoke else (1, 2, 3),
    )
    reduced_bridge_depth = min(profile_depth, 8 if smoke else 12)
    reduced_bridge_order = min(profile_order, 24 if smoke else 36)
    reduced_bridge_moduli = tuple(modulus for modulus in rhs_uniqueness_moduli if modulus <= (3 if smoke else 4))
    if not reduced_bridge_moduli:
        reduced_bridge_moduli = rhs_uniqueness_moduli[: min(len(rhs_uniqueness_moduli), 2)]
    reduced_bridge_error: str | None = None
    reduced_reciprocal_witness = None
    reduced_reciprocal_series: Series | None = None
    reduced_ratio_series: Series | None = None
    reduced_tail_transfer_equation: ReducedTailTransferEquation | None = None
    reduced_object_self_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_object_self_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_ratio_self_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_ratio_self_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_object_mahler_scan = SelfMahlerLinearScan((), (), (), ())
    reduced_ratio_mahler_scan = SelfMahlerLinearScan((), (), (), ())
    reduced_ratio_self_quotient_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_ratio_self_plus_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_ratio_self_signed_product_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_ratio_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_ratio_eta_quotient_scans: list[EtaQuotientRelationScan] = []
    try:
        reduced_reciprocal_witness, reduced_reciprocal_series = _reduced_reciprocal_bridge(
            template=reduced_candidate,
            symbol=sp.Symbol(series_symbol),
            depth=reduced_bridge_depth,
            order=reduced_bridge_order,
        )
        reduced_tail_transfer_equation = detect_reduced_tail_transfer_equation(
            reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
            symbol=sp.Symbol(series_symbol),
        )
        reduced_ratio_series = series_div(benchmark_recip[:reduced_bridge_order], reduced_reciprocal_series)
        reduced_object_self_polynomial_scan = scan_self_polynomial_uniqueness_relations(
            target_series=reduced_reciprocal_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            fg_degree_values=(1, 2),
            t_degree_values=(1, 2),
        )
        reduced_object_self_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
            target_series=reduced_reciprocal_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            t_degree_values=(1, 2),
        )
        reduced_object_mahler_scan = scan_self_mahler_linear_relations(
            target_series=reduced_reciprocal_series,
            moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 3),
            levels_checked=(2,),
            order=reduced_bridge_order,
            t_degree_values=(1, 2),
        )
        reduced_ratio_self_polynomial_scan = scan_self_polynomial_uniqueness_relations(
            target_series=reduced_ratio_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            fg_degree_values=(1, 2),
            t_degree_values=(1, 2),
        )
        reduced_ratio_self_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
            target_series=reduced_ratio_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            t_degree_values=(1, 2),
        )
        reduced_ratio_mahler_scan = scan_self_mahler_linear_relations(
            target_series=reduced_ratio_series,
            moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 3),
            levels_checked=(2,),
            order=reduced_bridge_order,
            t_degree_values=(1, 2),
        )
        reduced_ratio_self_quotient_product_scans = scan_ratio_self_quotient_product_relations(
            ratio_series=reduced_ratio_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
        reduced_ratio_self_plus_product_scans = scan_ratio_self_plus_product_relations(
            ratio_series=reduced_ratio_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
        reduced_ratio_self_signed_product_scans = scan_ratio_self_signed_product_relations(
            ratio_series=reduced_ratio_series,
            moduli=reduced_bridge_moduli,
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
        reduced_ratio_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
            ratio_series=reduced_ratio_series,
            moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 3),
            eta_levels=_eta_scan_levels(tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 3)),
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
        reduced_ratio_eta_quotient_scans = scan_ratio_eta_quotient_relations(
            ratio_series=reduced_ratio_series,
            levels=_eta_scan_levels(reduced_bridge_moduli),
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
    except Exception as exc:
        reduced_bridge_error = str(exc)
    advance_progress(
        completed_step="rhs-uniqueness-search",
        current_step="source-family-scans",
    )
    source_family_basis_catalog = _source_family_basis_catalog(record.closest_benchmark)
    source_family_base_series = tuple(
        (
            label,
            benchmark_name,
            continued_fraction_series_coeffs(
                get_benchmark(benchmark_name).canonical_template.normalized(),
                depth=profile_depth,
                order=profile_order,
            ),
        )
        for label, benchmark_name in source_family_basis_catalog
    )
    source_family_raw_basis_entries = _source_family_raw_basis_entries(
        ordered_base_families=source_family_base_series,
        powers=source_family_scan_powers,
        order=profile_order,
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    source_family_quotient_basis_entries = _source_family_quotient_basis_entries(
        ordered_base_families=source_family_base_series,
        powers=source_family_scan_powers,
        order=profile_order,
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    explicit_transform_family_base_series = tuple(
        item for item in source_family_base_series if item[0] in {"GG", "S"}
    )
    gg_base_family_entry = next(
        (item for item in source_family_base_series if item[0] == "GG"),
        None,
    )
    source_family_basis_series = tuple(
        (
            label,
            series,
        )
        for label, _, series in source_family_base_series
    )
    one_core_self_polynomial_scan = scan_source_correction_self_polynomial_uniqueness_relations(
        target_series=ratio_series,
        ordered_base_families=source_family_base_series,
        correction_size=1,
        moduli=rhs_uniqueness_moduli,
        order=profile_order,
        fg_degree_values=(1, 2) if smoke else (1, 2),
        t_degree_values=(1, 2) if smoke else (1, 2),
    )
    one_core_self_fractional_linear_scan = (
        scan_source_correction_self_fractional_linear_uniqueness_relations(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            correction_size=1,
            moduli=rhs_uniqueness_moduli,
            order=profile_order,
            t_degree_values=(1, 2) if smoke else (1, 2),
        )
    )
    two_core_self_polynomial_scan = scan_source_correction_self_polynomial_uniqueness_relations(
        target_series=ratio_series,
        ordered_base_families=source_family_base_series,
        correction_size=2,
        moduli=rhs_uniqueness_moduli,
        order=profile_order,
        fg_degree_values=(1, 2),
        t_degree_values=(1, 2),
    )
    two_core_self_fractional_linear_scan = (
        scan_source_correction_self_fractional_linear_uniqueness_relations(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            correction_size=2,
            moduli=rhs_uniqueness_moduli,
            order=profile_order,
            t_degree_values=(1, 2),
        )
    )
    source_family_prefix_scans = scan_named_prefix_boxes(
        target_series=ratio_series,
        ordered_basis_series=source_family_basis_series,
        order=profile_order,
        degree_values=(),
        include_polynomial=False,
        max_abs_exponent=4 if smoke else 6,
        solve_order=min(profile_order, 14 if smoke else 18),
    )
    source_family_multiplicative_scans = source_family_prefix_scans.multiplicative_scans
    source_family_fractional_linear_scans = source_family_prefix_scans.fractional_linear_scans
    source_family_two_layer_fractional_linear_scans = (
        source_family_prefix_scans.two_layer_fractional_linear_scans
    )
    parameterized_source_family_scans = scan_parameterized_source_family_power_boxes(
        target_series=ratio_series,
        ordered_base_families=source_family_base_series,
        powers=source_family_scan_powers,
        order=profile_order,
        degree_values=(1, 2),
        max_abs_exponent=4 if smoke else 6,
        solve_order=min(profile_order, 14 if smoke else 18),
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    advance_progress(
        completed_step="source-family-scans",
        current_step="cross-family-functional-scans",
    )
    two_core_source_family_eta_correction_scan = scan_two_core_source_family_eta_corrections(
        target_series=ratio_series,
        ordered_base_families=source_family_base_series,
        powers=source_family_scan_powers,
        eta_levels=eta_scan_levels,
        order=profile_order,
        max_abs_exponent=4 if smoke else 6,
        supplemental_powers_by_family=supplemental_source_family_powers,
        raw_basis_entries=source_family_raw_basis_entries,
    )
    quotient_core_source_family_eta_correction_scan = (
        scan_quotient_core_source_family_eta_corrections(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            eta_levels=eta_scan_levels,
            order=profile_order,
            max_abs_exponent=4 if smoke else 6,
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
        )
    )
    two_quotient_core_correction_entries = _two_quotient_core_correction_entries(
        target_series=ratio_series,
        quotient_basis_entries=source_family_quotient_basis_entries,
    )[1]
    two_quotient_core_source_family_eta_correction_scan = (
        scan_two_quotient_core_source_family_eta_corrections(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            eta_levels=eta_scan_levels,
            order=profile_order,
            max_abs_exponent=4 if smoke else 6,
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
        )
    )
    two_quotient_core_source_family_self_quotient_product_scan = (
        scan_two_quotient_core_source_family_self_quotient_products(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            moduli=source_family_scan_powers,
            order=profile_order,
            max_abs_exponent=4 if smoke else 6,
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
            correction_entries=two_quotient_core_correction_entries,
        )
    )
    two_quotient_core_source_family_self_polynomial_scan = (
        scan_two_quotient_core_source_family_self_polynomial_relations(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            moduli=source_family_scan_powers,
            order=profile_order,
            degree_values=(1, 2),
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
            correction_entries=two_quotient_core_correction_entries,
        )
    )
    two_quotient_core_source_family_self_eta_scan = (
        scan_two_quotient_core_source_family_self_eta_corrections(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            moduli=source_family_scan_powers,
            eta_levels=source_family_scan_powers,
            order=profile_order,
            max_abs_exponent=4 if smoke else 6,
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
            correction_entries=two_quotient_core_correction_entries,
        )
    )
    two_quotient_core_source_family_self_fractional_linear_scan = (
        scan_two_quotient_core_source_family_self_fractional_linear_relations(
            target_series=ratio_series,
            ordered_base_families=source_family_base_series,
            powers=source_family_scan_powers,
            moduli=source_family_scan_powers,
            eta_levels=source_family_scan_powers,
            order=profile_order,
            supplemental_powers_by_family=supplemental_source_family_powers,
            quotient_basis_entries=source_family_quotient_basis_entries,
            correction_entries=two_quotient_core_correction_entries,
        )
    )
    advance_progress(
        completed_step="cross-family-functional-scans",
        current_step="explicit-gg-family-scans",
    )
    source_family_eta_correction_scans = scan_source_family_eta_corrections(
        target_series=ratio_series,
        ordered_base_families=source_family_base_series,
        powers=source_family_scan_powers,
        eta_levels=eta_scan_levels,
        order=profile_order,
        max_abs_exponent=4 if smoke else 6,
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    explicit_source_family_transform_scans = scan_explicit_source_family_transform_templates(
        target_series=ratio_series,
        ordered_base_families=explicit_transform_family_base_series,
        powers=source_family_scan_powers,
        order=profile_order,
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    explicit_source_family_eta_correction_scans = scan_explicit_source_family_eta_correction_templates(
        target_series=ratio_series,
        ordered_base_families=explicit_transform_family_base_series,
        powers=source_family_scan_powers,
        eta_levels=eta_scan_levels,
        order=profile_order,
        max_abs_exponent=4 if smoke else 6,
        supplemental_powers_by_family=supplemental_source_family_powers,
    )
    gg_modular_equation_scan = (
        None
        if gg_base_family_entry is None
        else scan_gg_modular_equation_box(
            target_series=ratio_series,
            benchmark_name=gg_base_family_entry[1],
            gg_series=gg_base_family_entry[2],
            order=profile_order,
            degree_values=(1, 2),
            max_abs_exponent=4 if smoke else 6,
            solve_order=min(profile_order, 14 if smoke else 18),
            supplemental_powers=()
            if smoke
            else supplemental_source_family_powers.get("GG", ()),
        )
    )
    advance_progress(
        completed_step="explicit-gg-family-scans",
        current_step="benchmark-tower-scans",
    )
    power_tower_scans: list[BenchmarkPowerRelationScan] = []
    ratio_power_tower_scans: list[BenchmarkPowerRelationScan] = []
    ratio_self_quotient_product_scans: list[SelfQuotientProductRelationScan] = []
    ratio_eta_quotient_scans: list[EtaQuotientRelationScan] = []
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
        ratio_self_quotient_product_scans = scan_ratio_self_quotient_product_relations(
            ratio_series=ratio_series,
            moduli=tuple(benchmark_power_series),
            order=profile_order,
            max_abs_exponent=6 if smoke else 8,
        )
        ratio_eta_quotient_scans = scan_ratio_eta_quotient_relations(
            ratio_series=ratio_series,
            levels=_eta_scan_levels(tuple(benchmark_power_series)),
            order=profile_order,
            max_abs_exponent=6 if smoke else 8,
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
    advance_progress(
        completed_step="benchmark-tower-scans",
        current_step="final-render",
    )
    final_render_started_at = current_stage_start
    build_elapsed_before_render = final_render_started_at - build_started_at
    completed_stage_timing_lines = [
        f"- `{step_name}`: `{stage_elapsed_seconds[step_name]:.2f}`"
        for step_name in progress_steps[:-1]
        if stage_elapsed_seconds[step_name] is not None
    ]

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
        f"- Build elapsed seconds before final render: `{build_elapsed_before_render:.2f}`",
        "",
        "## Build Timing",
        "",
        *completed_stage_timing_lines,
        "- `final-render`: `in_progress`",
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

    lines.extend(
        [
            "",
            "## RHS Uniqueness Search",
            "",
            "We also ran a theorem-facing search directly on the ratio object, looking for a compact right-hand-side defining equation for:",
            "",
            f"- `F = candidate / {record.closest_benchmark}`",
            f"- Moduli checked: {', '.join(f'`m={modulus}`' for modulus in rhs_self_polynomial_scan.moduli_checked)}",
            f"- Self-polynomial boxes: {', '.join(f'`deg_(F,G) <= {degree}`' for degree in rhs_self_polynomial_scan.fg_degree_values)}",
            f"- `t`-degree boxes: {', '.join(f'`deg_t <= {degree}`' for degree in rhs_self_polynomial_scan.t_degree_values)}",
            f"- Self-fractional-linear `t`-degree boxes: {', '.join(f'`deg_t <= {degree}`' for degree in rhs_self_fractional_linear_scan.t_degree_values)}",
            "",
            "### Polynomial Functional Box",
            "",
            "```text",
            "P(t, F(t), F(t^m)) = 0",
            "```",
            "",
        ]
    )
    if not rhs_self_polynomial_scan.hits:
        lines.extend(
            [
                "No candidate-dependent self-polynomial uniqueness relation was found in the scanned box.",
                "",
            ]
        )
    for hit in rhs_self_polynomial_scan.hits:
        relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
        relation_expr = hit.relation.as_sympy(relation_symbols)
        relation_series_by_variable = {
            "T": _t_series(order=profile_order),
            "F": ratio_series,
            f"G{hit.modulus}": benchmark_power_substitution_series(
                ratio_series,
                power=hit.modulus,
                order=profile_order,
            ),
        }
        residual = _relation_residual_series(
            hit.relation,
            series_by_variable=relation_series_by_variable,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                f"- Modulus `m={hit.modulus}`, `deg_(F,G) <= {hit.max_fg_total_degree}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_expr(relation_expr),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    lines.extend(
        [
            "### Fractional-Linear Functional Box",
            "",
            "```text",
            "F(t) = (A(t) + B(t)*(F(t^m) - 1)) / (C(t) + D(t)*(F(t^m) - 1))",
            "```",
            "",
        ]
    )
    if not rhs_self_fractional_linear_scan.hits:
        lines.extend(
            [
                "No candidate-dependent self-fractional-linear uniqueness relation was found in the scanned box.",
                "",
            ]
        )
    for hit in rhs_self_fractional_linear_scan.hits:
        residual = _self_t_polynomial_fractional_linear_relation_residual_series(
            hit.relation,
            target_series=ratio_series,
            modulus=hit.modulus,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                f"- Modulus `m={hit.modulus}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_self_t_polynomial_fractional_linear_relation(
                    hit.relation,
                    modulus=hit.modulus,
                    target_variable="F",
                    series_symbol=series_symbol,
                ),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    lines.extend(
        [
            "### One-Source-Core Correction Objects",
            "",
            "We then stripped a single nearby source core and repeated the uniqueness search on the residual correction object:",
            "",
            "```text",
            "G = F / S",
            "P(t, G(t), G(t^m)) = 0",
            "G(t) = (A(t) + B(t)*(G(t^m) - 1)) / (C(t) + D(t)*(G(t^m) - 1))",
            "```",
            "",
            f"- Source cores checked: {', '.join(f'`{label}`' for label, _, _ in source_family_base_series)}",
            f"- Correction objects checked: `{one_core_self_polynomial_scan.total_corrections_checked}`",
            "",
        ]
    )
    if not one_core_self_polynomial_scan.hits:
        lines.append("No one-core self-polynomial correction hit was found in the scanned box.")
        lines.append("")
    for hit in one_core_self_polynomial_scan.hits:
        correction_series = ratio_series
        for label in hit.basis_labels:
            correction_series = series_div(
                correction_series,
                next(series for name, _, series in source_family_base_series if name == label),
            )
        relation_series_by_variable = {
            "T": _t_series(order=profile_order),
            "F": correction_series,
            f"G{hit.modulus}": benchmark_power_substitution_series(
                correction_series,
                power=hit.modulus,
                order=profile_order,
            ),
        }
        residual = _relation_residual_series(
            hit.relation,
            series_by_variable=relation_series_by_variable,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
        relation_expr = hit.relation.as_sympy(relation_symbols)
        lines.extend(
            [
                f"- One-core correction `{', '.join(hit.basis_labels)}` with `m={hit.modulus}`, `deg_(F,G) <= {hit.max_fg_total_degree}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_source_correction_expression(
                    basis_labels=hit.basis_labels,
                    basis_expressions=hit.basis_expressions,
                    target_variable="F",
                ),
                _format_expr(relation_expr),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    if not one_core_self_fractional_linear_scan.hits:
        lines.append("No one-core self-fractional-linear correction hit was found in the scanned box.")
        lines.append("")
    for hit in one_core_self_fractional_linear_scan.hits:
        correction_series = ratio_series
        for label in hit.basis_labels:
            correction_series = series_div(
                correction_series,
                next(series for name, _, series in source_family_base_series if name == label),
            )
        residual = _self_t_polynomial_fractional_linear_relation_residual_series(
            hit.relation,
            target_series=correction_series,
            modulus=hit.modulus,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                f"- One-core correction `{', '.join(hit.basis_labels)}` with `m={hit.modulus}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_source_correction_expression(
                    basis_labels=hit.basis_labels,
                    basis_expressions=hit.basis_expressions,
                    target_variable="F",
                ),
                _format_self_t_polynomial_fractional_linear_relation(
                    hit.relation,
                    modulus=hit.modulus,
                    target_variable="G",
                    series_symbol=series_symbol,
                ),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    lines.extend(
        [
            "### Two-Source-Core Correction Objects",
            "",
            "We also stripped products of two nearby source cores and repeated the same bounded uniqueness search on the residual object:",
            "",
            "```text",
            "H = F / (S1 * S2)",
            "P(t, H(t), H(t^m)) = 0",
            "H(t) = (A(t) + B(t)*(H(t^m) - 1)) / (C(t) + D(t)*(H(t^m) - 1))",
            "```",
            "",
            f"- Two-core correction objects checked: `{two_core_self_polynomial_scan.total_corrections_checked}`",
            "",
        ]
    )
    if not two_core_self_polynomial_scan.hits:
        lines.append("No two-core self-polynomial correction hit was found in the scanned box.")
        lines.append("")
    for hit in two_core_self_polynomial_scan.hits:
        correction_series = ratio_series
        for label in hit.basis_labels:
            correction_series = series_div(
                correction_series,
                next(series for name, _, series in source_family_base_series if name == label),
            )
        relation_series_by_variable = {
            "T": _t_series(order=profile_order),
            "F": correction_series,
            f"G{hit.modulus}": benchmark_power_substitution_series(
                correction_series,
                power=hit.modulus,
                order=profile_order,
            ),
        }
        residual = _relation_residual_series(
            hit.relation,
            series_by_variable=relation_series_by_variable,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
        relation_expr = hit.relation.as_sympy(relation_symbols)
        lines.extend(
            [
                f"- Two-core correction `{', '.join(hit.basis_labels)}` with `m={hit.modulus}`, `deg_(F,G) <= {hit.max_fg_total_degree}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_source_correction_expression(
                    basis_labels=hit.basis_labels,
                    basis_expressions=hit.basis_expressions,
                    target_variable="F",
                ),
                _format_expr(relation_expr),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    if not two_core_self_fractional_linear_scan.hits:
        lines.append("No two-core self-fractional-linear correction hit was found in the scanned box.")
        lines.append("")
    for hit in two_core_self_fractional_linear_scan.hits:
        correction_series = ratio_series
        for label in hit.basis_labels:
            correction_series = series_div(
                correction_series,
                next(series for name, _, series in source_family_base_series if name == label),
            )
        residual = _self_t_polynomial_fractional_linear_relation_residual_series(
            hit.relation,
            target_series=correction_series,
            modulus=hit.modulus,
            order=profile_order,
        )
        residual_ok = all(sp.simplify(value) == 0 for value in residual)
        lines.extend(
            [
                f"- Two-core correction `{', '.join(hit.basis_labels)}` with `m={hit.modulus}`, `deg_t <= {hit.max_t_degree}` produced:",
                "",
                "```text",
                _format_source_correction_expression(
                    basis_labels=hit.basis_labels,
                    basis_expressions=hit.basis_expressions,
                    target_variable="F",
                ),
                _format_self_t_polynomial_fractional_linear_relation(
                    hit.relation,
                    modulus=hit.modulus,
                    target_variable="G",
                    series_symbol=series_symbol,
                ),
                "```",
                "",
                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                "",
            ]
        )

    lines.extend(
        [
            "### Rational-Equivalence Reduced Object",
            "",
            "We also switched from the raw hero reciprocal to the reduced-by-factor object coming from the exact convergent-factor / rational-equivalence bridge:",
            "",
        ]
    )
    if reduced_bridge_error is not None:
        lines.extend(
            [
                "Reduced-object bridge construction failed:",
                "",
                "```text",
                reduced_bridge_error,
                "```",
                "",
            ]
        )
    elif reduced_reciprocal_witness is not None and reduced_reciprocal_series is not None and reduced_ratio_series is not None:
        reduced_coeffs = reduced_reciprocal_witness.reduction.reduced_coeffs
        lines.extend(
            [
                "```text",
                "R = reduced hero reciprocal object",
                f"F_red = B1 / R    where    B1 = 1 / {record.closest_benchmark}",
                f"a1_red = {_format_expr(reduced_coeffs.a_terms[1])}",
                f"b1_red = {_format_expr(reduced_coeffs.b_terms[1])}",
                f"a2_red = {_format_expr(reduced_coeffs.a_terms[2])}",
                f"b2_red = {_format_expr(reduced_coeffs.b_terms[2])}",
                f"a3_red = {_format_expr(reduced_coeffs.a_terms[3])}",
                f"b3_red = {_format_expr(reduced_coeffs.b_terms[3])}",
                f"r1 = {_format_expr(reduced_reciprocal_witness.scale_terms[1])}",
                f"r2 = {_format_expr(reduced_reciprocal_witness.scale_terms[2])}",
                f"r3 = {_format_expr(reduced_reciprocal_witness.scale_terms[3])}",
                "```",
                "",
                "- This is the first source-informed RHS lane: it uses the exact reduction/equivalence bridge rather than a generic low-degree guess on the raw ratio object.",
                f"- To keep this bridge tractable inside `identify`, the reduced-object lane is currently bounded at depth `{reduced_bridge_depth}` and order `{reduced_bridge_order}`.",
                "",
            ]
        )
        if reduced_tail_transfer_equation is not None:
            tail_header, tail_coeffs, tail_functional, tail_polynomial = _format_reduced_tail_transfer_equation(
                reduced_tail_transfer_equation
            )
            lines.extend(
                [
                    "- From stage `3` onward, the reduced coefficients collapse into one stationary tail family.",
                    "",
                    "```text",
                    tail_header,
                    tail_coeffs,
                    tail_functional,
                    tail_polynomial,
                    "```",
                    "",
                ]
            )

        if not reduced_object_self_polynomial_scan.hits:
            lines.append("No reduced-object self-polynomial hit was found in the scanned box.")
        else:
            lines.append("Reduced-object self-polynomial hits were found:")
            for hit in reduced_object_self_polynomial_scan.hits:
                relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
                relation_expr = hit.relation.as_sympy(relation_symbols)
                lines.append(
                    f"- `R`, `m={hit.modulus}`, `deg_(F,G) <= {hit.max_fg_total_degree}`, `deg_t <= {hit.max_t_degree}`: `{_format_expr(relation_expr)}`"
                )
        lines.append("")

        if not reduced_object_self_fractional_linear_scan.hits:
            lines.append("No reduced-object self-fractional-linear hit was found in the scanned box.")
        else:
            lines.append("Reduced-object self-fractional-linear hits were found:")
            for hit in reduced_object_self_fractional_linear_scan.hits:
                lines.append(
                    f"- `R`, `m={hit.modulus}`, `deg_t <= {hit.max_t_degree}`: `{_format_self_t_polynomial_fractional_linear_relation(hit.relation, modulus=hit.modulus, target_variable='R', series_symbol=series_symbol)}`"
                )
        lines.append("")

        if not reduced_object_mahler_scan.hits:
            lines.append("No reduced-object Mahler/transfer hit was found in the scanned box.")
        else:
            lines.append("Reduced-object Mahler/transfer hits were found:")
            for hit in reduced_object_mahler_scan.hits:
                relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
                relation_expr = hit.relation.as_sympy(relation_symbols)
                lines.append(
                    f"- `R`, `m={hit.modulus}`, `levels={hit.levels}`, `deg_t <= {hit.max_t_degree}`: `{_format_expr(relation_expr)}`"
                )
        lines.append("")

        if not reduced_ratio_self_polynomial_scan.hits:
            lines.append("No reduced-ratio self-polynomial hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self-polynomial hits were found:")
            for hit in reduced_ratio_self_polynomial_scan.hits:
                relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
                relation_expr = hit.relation.as_sympy(relation_symbols)
                lines.append(
                    f"- `F_red`, `m={hit.modulus}`, `deg_(F,G) <= {hit.max_fg_total_degree}`, `deg_t <= {hit.max_t_degree}`: `{_format_expr(relation_expr)}`"
                )
        lines.append("")

        if not reduced_ratio_self_fractional_linear_scan.hits:
            lines.append("No reduced-ratio self-fractional-linear hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self-fractional-linear hits were found:")
            for hit in reduced_ratio_self_fractional_linear_scan.hits:
                lines.append(
                    f"- `F_red`, `m={hit.modulus}`, `deg_t <= {hit.max_t_degree}`: `{_format_self_t_polynomial_fractional_linear_relation(hit.relation, modulus=hit.modulus, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        if not reduced_ratio_mahler_scan.hits:
            lines.append("No reduced-ratio Mahler/transfer hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio Mahler/transfer hits were found:")
            for hit in reduced_ratio_mahler_scan.hits:
                relation_symbols = tuple(sp.Symbol(name) for name in hit.relation.variables)
                relation_expr = hit.relation.as_sympy(relation_symbols)
                lines.append(
                    f"- `F_red`, `m={hit.modulus}`, `levels={hit.levels}`, `deg_t <= {hit.max_t_degree}`: `{_format_expr(relation_expr)}`"
                )
        lines.append("")

        reduced_ratio_self_product_hits = [scan for scan in reduced_ratio_self_quotient_product_scans if scan.relation is not None]
        if not reduced_ratio_self_product_hits:
            lines.append("No reduced-ratio self-quotient finite-product hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self-quotient finite-product hits were found:")
            for scan in reduced_ratio_self_product_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `{_format_self_quotient_product_relation(scan.relation, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        reduced_ratio_self_plus_hits = [scan for scan in reduced_ratio_self_plus_product_scans if scan.relation is not None]
        if not reduced_ratio_self_plus_hits:
            lines.append("No reduced-ratio self-quotient plus-product hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self-quotient plus-product hits were found:")
            for scan in reduced_ratio_self_plus_product_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `{_format_self_plus_product_relation(scan.relation, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        reduced_ratio_self_signed_hits = [scan for scan in reduced_ratio_self_signed_product_scans if scan.relation is not None]
        if not reduced_ratio_self_signed_hits:
            lines.append("No reduced-ratio self-quotient signed-product hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self-quotient signed-product hits were found:")
            for scan in reduced_ratio_self_signed_product_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `{_format_self_signed_product_relation(scan.relation, modulus=scan.modulus, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        reduced_ratio_self_signed_eta_hits = [scan for scan in reduced_ratio_self_signed_eta_scans if scan.relation is not None]
        if not reduced_ratio_self_signed_eta_hits:
            lines.append("No reduced-ratio self signed-eta transfer hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self signed-eta transfer hits were found:")
            for scan in reduced_ratio_self_signed_eta_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `m={scan.modulus}`, `N={scan.level}`: `{_format_self_signed_eta_relation(scan.relation, modulus=scan.modulus, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        reduced_ratio_eta_hits = [scan for scan in reduced_ratio_eta_quotient_scans if scan.relation is not None]
        if not reduced_ratio_eta_hits:
            lines.append("No reduced-ratio eta-quotient hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio eta-quotient hits were found:")
            for scan in reduced_ratio_eta_quotient_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- Level `N={scan.level}`: `{_format_eta_quotient_relation(scan.relation, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

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

    if source_family_multiplicative_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Source-Family Multiplicative Scan",
                "",
                "We also searched for exact multiplicative corrections built from nearby named source families:",
                "",
                "```text",
                "F = prod_i S_i^e_i",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- The source-family bases are evaluated in the same variable view used above.",
            ]
        )
        for label, benchmark_name in source_family_basis_catalog:
            lines.append(f"- `{label} = {benchmark_name}`")
        lines.extend(
            [
                "",
                "- Prefixes are scanned in that order, solving exact integer exponents from the log-series constraints and then verifying by exact series re-expansion.",
                "",
            ]
        )

        any_source_hit = any(scan.relation is not None for scan in source_family_multiplicative_scans)
        if not any_source_hit:
            lines.append("No source-family multiplicative relation was found in any scanned prefix box.")
            lines.append("")

        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in source_family_multiplicative_scans
            if scan.error is None
        ]
        if not any_source_hit and no_hit_labels:
            lines.append(
                f"- No hit for source-family prefixes ending at {', '.join(no_hit_labels)}."
            )

        for scan in source_family_multiplicative_scans:
            if scan.error is not None:
                lines.append(
                    f"- Source-family prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )
            elif scan.relation is not None:
                residual = _multiplicative_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        label: series
                        for label, series in source_family_basis_series
                        if label in scan.basis_labels
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Source-family prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_multiplicative_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

    if source_family_fractional_linear_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Source-Family Fractional-Linear Scan",
                "",
                "We also searched for low-complexity fractional-linear corrections built from nearby named source families:",
                "",
                "```text",
                "F = (1 + sum a_i*(S_i - 1)) / (1 + sum b_i*(S_i - 1))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- The source-family bases are evaluated in the same variable view used above.",
            ]
        )
        for label, benchmark_name in source_family_basis_catalog:
            lines.append(f"- `{label} = {benchmark_name}`")
        lines.extend(
            [
                "",
                "- Prefixes are scanned in that order, solving an exact linear system for the numerator and denominator correction coefficients in each source-family box.",
                "",
            ]
        )

        any_source_fractional_hit = any(
            scan.relation is not None for scan in source_family_fractional_linear_scans
        )
        if not any_source_fractional_hit:
            lines.append("No source-family fractional-linear relation was found in any scanned prefix box.")
            lines.append("")

        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in source_family_fractional_linear_scans
            if scan.error is None
        ]
        if not any_source_fractional_hit and no_hit_labels:
            lines.append(
                f"- No hit for source-family fractional-linear prefixes ending at {', '.join(no_hit_labels)}."
            )

        for scan in source_family_fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Source-family fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )
            elif scan.relation is not None:
                residual = _fractional_linear_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        label: series
                        for label, series in source_family_basis_series
                        if label in scan.basis_labels
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Source-family fractional-linear prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_fractional_linear_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

    if source_family_two_layer_fractional_linear_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Source-Family Two-Layer Fractional-Linear Scan",
                "",
                "We then expanded to a second-ring nonlinear box built from two single-basis factors drawn from the named source-family prefixes:",
                "",
                "```text",
                "F = ((1 + a0*(X1 - 1)) / (1 + b0*(Y1 - 1))) * ((1 + a1*(X2 - 1)) / (1 + b1*(Y2 - 1)))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- The source-family bases are evaluated in the same variable view used above.",
            ]
        )
        for label, benchmark_name in source_family_basis_catalog:
            lines.append(f"- `{label} = {benchmark_name}`")
        lines.extend(
            [
                "",
                "- Prefixes checked: `(RR, cubic)`, then `(RR, cubic, GG)`, and so on through the final listed source-family basis.",
                "- Each prefix scans low-complexity two-factor templates and verifies any candidate hit by exact series re-expansion.",
                "",
            ]
        )

        any_source_two_layer_hit = any(
            scan.total_hits > 0 for scan in source_family_two_layer_fractional_linear_scans
        )
        if not any_source_two_layer_hit:
            lines.append("No source-family two-layer fractional-linear relation was found in any scanned prefix box.")
            lines.append("")

        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in source_family_two_layer_fractional_linear_scans
            if scan.error is None and scan.total_hits == 0
        ]
        if not any_source_two_layer_hit and no_hit_labels:
            lines.append(
                f"- No hit for source-family two-layer prefixes ending at {', '.join(no_hit_labels)}."
            )

        for scan in source_family_two_layer_fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Source-family two-layer prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )
                continue
            if scan.total_hits == 0:
                continue
            lines.append(
                f"- Source-family two-layer prefix ending at `{scan.basis_labels[-1]}` found `{scan.total_hits}` exact hit(s) after checking `{scan.tuples_checked}` template(s):"
            )
            lines.append("")
            basis_series_by_variable = {
                label: series
                for label, series in source_family_basis_series
                if label in scan.basis_labels
            }
            for relation in scan.relations:
                residual = _two_layer_fractional_linear_relation_residual_series(
                    relation,
                    target_series=ratio_series,
                    basis_series_by_variable=basis_series_by_variable,
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

    if parameterized_source_family_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Parameterized Source-Family Power Scan",
                "",
                "We also scanned short power ladders inside each named source family so the family meaning stays explicit:",
                "",
                "```text",
                "P(F, T_i) = 0",
                "F = prod_i T_i^e_i",
                "F = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))",
                "F = prod_j (1 + a_j*(T_r(j) - 1)) / (1 + b_j*(T_s(j) - 1))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Each family is scanned separately, using the base object together with powered substitutions in the same variable view.",
                "- This keeps the Gordon/Hirschhorn family labels explicit instead of collapsing them into one anonymous mixed basis.",
                "- The low-degree polynomial box is motivated by literature where `GG` / Hirschhorn-type objects can satisfy nontrivial power-substitution identities without reducing to a pure product or a simple quotient.",
                "- We now also include a family-preserving two-layer fractional-linear box, so simple nonlinear corrections stay inside one literature family instead of mixing labels.",
                "- We also scan a within-family quotient ladder `Qk = Tk / T1`, which is often a more natural coordinate for power-substitution identities than the raw powered objects themselves.",
                "- That quotient ladder now also gets its own two-layer fractional-linear pass, so the first nonlinear corrections can stay in quotient coordinates without crossing families.",
                "- We also scan a mixed quotient basis built from the family base object together with the quotient ladder, so relations of the form `T1 * correction(Q2, Q3, ...)` can surface without mixing literature families.",
                "- That mixed quotient basis now also gets its own two-layer fractional-linear pass, so the first nonlinear corrections can use both the base object and quotient coordinates while still staying in one family.",
            ]
        )
        lines.append("- The exact powered labels are listed separately inside each family subsection, because the literature-motivated ladders are now family-specific.")
        lines.append("")

        for family_scan in parameterized_source_family_scans:
            basis_series_by_variable = dict(family_scan.ordered_basis_series)
            lines.extend(
                [
                    f"### `{family_scan.family_label}` Family",
                    "",
                    f"- Base benchmark: `{family_scan.benchmark_name}`",
                ]
            )
            basis_descriptions = [f"`{family_scan.family_label} = {family_scan.family_label}({series_symbol})`"]
            for label, _ in family_scan.ordered_basis_series[1:]:
                power = int(label.removeprefix(family_scan.family_label))
                basis_descriptions.append(
                    f"`{label} = {family_scan.family_label}({series_symbol}^{power})`"
                )
            lines.append(f"- Basis ladder: {', '.join(basis_descriptions)}")

            grouped_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
            for scan in family_scan.polynomial_scans:
                grouped_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

            any_polynomial_hit = any(scan.relation is not None for scan in family_scan.polynomial_scans)
            if not any_polynomial_hit:
                lines.append("- Polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes.")
            for degree in sorted(grouped_polynomial_scans):
                scans = grouped_polynomial_scans[degree]
                no_hit_labels = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in scans
                    if scan.error is None and scan.relation is None
                ]
                if not any_polynomial_hit and no_hit_labels:
                    lines.append(
                        f"- Polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                    )
                for scan in scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.relation is None:
                        continue
                    relation_series = {"F": ratio_series}
                    relation_series.update(
                        {
                            label: basis_series_by_variable[label]
                            for label in scan.basis_labels
                        }
                    )
                    residual = _relation_residual_series(
                        scan.relation,
                        series_by_variable=relation_series,
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    sym_map = {name: sp.Symbol(name) for name in scan.relation.variables}
                    poly = scan.relation.as_sympy(tuple(sym_map[name] for name in scan.relation.variables))
                    lines.extend(
                        [
                            f"- Polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_expr(poly),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

            multiplicative_hits = [scan for scan in family_scan.multiplicative_scans if scan.relation is not None]
            multiplicative_no_hits = [
                f"`{scan.basis_labels[-1]}`"
                for scan in family_scan.multiplicative_scans
                if scan.error is None and scan.relation is None
            ]
            if not multiplicative_hits and multiplicative_no_hits:
                lines.append(
                    f"- Multiplicative scan: no hit for prefixes ending at {', '.join(multiplicative_no_hits)}."
                )
            for scan in family_scan.multiplicative_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )
                    continue
                if scan.relation is None:
                    continue
                residual = _multiplicative_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        label: basis_series_by_variable[label]
                        for label in scan.basis_labels
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Multiplicative prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_multiplicative_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

            fractional_hits = [scan for scan in family_scan.fractional_linear_scans if scan.relation is not None]
            fractional_no_hits = [
                f"`{scan.basis_labels[-1]}`"
                for scan in family_scan.fractional_linear_scans
                if scan.error is None and scan.relation is None
            ]
            if not fractional_hits and fractional_no_hits:
                lines.append(
                    f"- Fractional-linear scan: no hit for prefixes ending at {', '.join(fractional_no_hits)}."
                )
            for scan in family_scan.fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )
                    continue
                if scan.relation is None:
                    continue
                residual = _fractional_linear_relation_residual_series(
                    scan.relation,
                    target_series=ratio_series,
                    basis_series_by_variable={
                        label: basis_series_by_variable[label]
                        for label in scan.basis_labels
                    },
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Fractional-linear prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                        "",
                        "```text",
                        _format_fractional_linear_relation(scan.relation, target_variable="F"),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                        "",
                    ]
                )

            two_layer_hits = [
                scan for scan in family_scan.two_layer_fractional_linear_scans if scan.total_hits > 0
            ]
            two_layer_no_hits = [
                f"`{scan.basis_labels[-1]}`"
                for scan in family_scan.two_layer_fractional_linear_scans
                if scan.error is None and scan.total_hits == 0
            ]
            if not two_layer_hits and two_layer_no_hits:
                lines.append(
                    f"- Two-layer fractional-linear scan: no hit for prefixes ending at {', '.join(two_layer_no_hits)}."
                )
            for scan in family_scan.two_layer_fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )
                    continue
                if scan.total_hits == 0:
                    continue
                lines.append(
                    f"- Two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` found `{scan.total_hits}` exact hit(s) after checking `{scan.tuples_checked}` template(s):"
                )
                lines.append("")
                relation_basis = {
                    label: basis_series_by_variable[label]
                    for label in scan.basis_labels
                }
                for relation in scan.relations:
                    residual = _two_layer_fractional_linear_relation_residual_series(
                        relation,
                        target_series=ratio_series,
                        basis_series_by_variable=relation_basis,
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            "```text",
                            _format_two_layer_fractional_linear_relation(
                                relation,
                                target_variable="F",
                            ),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

            if family_scan.quotient_basis_series:
                quotient_basis_series_by_variable = {
                    label: series for label, _, series in family_scan.quotient_basis_series
                }
                quotient_descriptions = [
                    f"`{label} = {expr}`"
                    for label, expr, _ in family_scan.quotient_basis_series
                ]
                lines.append(f"- Quotient ladder: {', '.join(quotient_descriptions)}")

                grouped_quotient_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
                for scan in family_scan.quotient_polynomial_scans:
                    grouped_quotient_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

                any_quotient_polynomial_hit = any(
                    scan.relation is not None for scan in family_scan.quotient_polynomial_scans
                )
                if not any_quotient_polynomial_hit:
                    lines.append(
                        "- Quotient-ladder polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes."
                    )
                for degree in sorted(grouped_quotient_polynomial_scans):
                    scans = grouped_quotient_polynomial_scans[degree]
                    no_hit_labels = [
                        f"`{scan.basis_labels[-1]}`"
                        for scan in scans
                        if scan.error is None and scan.relation is None
                    ]
                    if not any_quotient_polynomial_hit and no_hit_labels:
                        lines.append(
                            f"- Quotient-ladder polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                        )
                    for scan in scans:
                        if scan.error is not None:
                            lines.append(
                                f"- Quotient-ladder polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                            )
                            continue
                        if scan.relation is None:
                            continue
                        relation_series = {"F": ratio_series}
                        relation_series.update(
                            {
                                label: quotient_basis_series_by_variable[label]
                                for label in scan.basis_labels
                            }
                        )
                        residual = _relation_residual_series(
                            scan.relation,
                            series_by_variable=relation_series,
                            order=profile_order,
                        )
                        residual_ok = all(sp.simplify(value) == 0 for value in residual)
                        sym_map = {name: sp.Symbol(name) for name in scan.relation.variables}
                        poly = scan.relation.as_sympy(tuple(sym_map[name] for name in scan.relation.variables))
                        lines.extend(
                            [
                                f"- Quotient-ladder polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                                "",
                                "```text",
                                _format_expr(poly),
                                "```",
                                "",
                                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                                "",
                            ]
                        )

                quotient_multiplicative_hits = [
                    scan for scan in family_scan.quotient_multiplicative_scans if scan.relation is not None
                ]
                quotient_multiplicative_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.quotient_multiplicative_scans
                    if scan.error is None and scan.relation is None
                ]
                if not quotient_multiplicative_hits and quotient_multiplicative_no_hits:
                    lines.append(
                        f"- Quotient-ladder multiplicative scan: no hit for prefixes ending at {', '.join(quotient_multiplicative_no_hits)}."
                    )
                for scan in family_scan.quotient_multiplicative_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Quotient-ladder multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.relation is None:
                        continue
                    residual = _multiplicative_relation_residual_series(
                        scan.relation,
                        target_series=ratio_series,
                        basis_series_by_variable={
                            label: quotient_basis_series_by_variable[label]
                            for label in scan.basis_labels
                        },
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- Quotient-ladder multiplicative prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_multiplicative_relation(scan.relation, target_variable="F"),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

                quotient_fractional_hits = [
                    scan for scan in family_scan.quotient_fractional_linear_scans if scan.relation is not None
                ]
                quotient_fractional_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.quotient_fractional_linear_scans
                    if scan.error is None and scan.relation is None
                ]
                if not quotient_fractional_hits and quotient_fractional_no_hits:
                    lines.append(
                        f"- Quotient-ladder fractional-linear scan: no hit for prefixes ending at {', '.join(quotient_fractional_no_hits)}."
                    )
                for scan in family_scan.quotient_fractional_linear_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Quotient-ladder fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.relation is None:
                        continue
                    residual = _fractional_linear_relation_residual_series(
                        scan.relation,
                        target_series=ratio_series,
                        basis_series_by_variable={
                            label: quotient_basis_series_by_variable[label]
                            for label in scan.basis_labels
                        },
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- Quotient-ladder fractional-linear prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_fractional_linear_relation(scan.relation, target_variable="F"),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

                quotient_two_layer_hits = [
                    scan
                    for scan in family_scan.quotient_two_layer_fractional_linear_scans
                    if scan.total_hits > 0
                ]
                quotient_two_layer_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.quotient_two_layer_fractional_linear_scans
                    if scan.error is None and scan.total_hits == 0
                ]
                if not quotient_two_layer_hits and quotient_two_layer_no_hits:
                    lines.append(
                        f"- Quotient-ladder two-layer fractional-linear scan: no hit for prefixes ending at {', '.join(quotient_two_layer_no_hits)}."
                    )
                for scan in family_scan.quotient_two_layer_fractional_linear_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Quotient-ladder two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.total_hits == 0:
                        continue
                    lines.append(
                        f"- Quotient-ladder two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` found `{scan.total_hits}` exact hit(s) after checking `{scan.tuples_checked}` template(s):"
                    )
                    lines.append("")
                    relation_basis = {
                        label: quotient_basis_series_by_variable[label]
                        for label in scan.basis_labels
                    }
                    for relation in scan.relations:
                        residual = _two_layer_fractional_linear_relation_residual_series(
                            relation,
                            target_series=ratio_series,
                            basis_series_by_variable=relation_basis,
                            order=profile_order,
                        )
                        residual_ok = all(sp.simplify(value) == 0 for value in residual)
                        lines.extend(
                            [
                                "```text",
                                _format_two_layer_fractional_linear_relation(
                                    relation,
                                    target_variable="F",
                                ),
                                "```",
                                "",
                                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                                "",
                            ]
                        )

                mixed_quotient_basis_series_by_variable = {
                    label: series for label, _, series in family_scan.mixed_quotient_basis_series
                }
                mixed_quotient_descriptions = []
                for label, expr, _ in family_scan.mixed_quotient_basis_series:
                    if label == family_scan.family_label and expr == family_scan.family_label:
                        mixed_quotient_descriptions.append(
                            f"`{label} = {family_scan.family_label}({series_symbol})`"
                        )
                    else:
                        mixed_quotient_descriptions.append(f"`{label} = {expr}`")
                lines.append(f"- Mixed quotient basis: {', '.join(mixed_quotient_descriptions)}")

                grouped_mixed_quotient_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
                for scan in family_scan.mixed_quotient_polynomial_scans:
                    grouped_mixed_quotient_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

                any_mixed_quotient_polynomial_hit = any(
                    scan.relation is not None for scan in family_scan.mixed_quotient_polynomial_scans
                )
                if not any_mixed_quotient_polynomial_hit:
                    lines.append(
                        "- Mixed-quotient polynomial scan: no candidate-dependent hit was found in the checked low-degree boxes."
                    )
                for degree in sorted(grouped_mixed_quotient_polynomial_scans):
                    scans = grouped_mixed_quotient_polynomial_scans[degree]
                    no_hit_labels = [
                        f"`{scan.basis_labels[-1]}`"
                        for scan in scans
                        if scan.error is None and scan.relation is None
                    ]
                    if not any_mixed_quotient_polynomial_hit and no_hit_labels:
                        lines.append(
                            f"- Mixed-quotient polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                        )
                    for scan in scans:
                        if scan.error is not None:
                            lines.append(
                                f"- Mixed-quotient polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                            )
                            continue
                        if scan.relation is None:
                            continue
                        relation_series = {"F": ratio_series}
                        relation_series.update(
                            {
                                label: mixed_quotient_basis_series_by_variable[label]
                                for label in scan.basis_labels
                            }
                        )
                        residual = _relation_residual_series(
                            scan.relation,
                            series_by_variable=relation_series,
                            order=profile_order,
                        )
                        residual_ok = all(sp.simplify(value) == 0 for value in residual)
                        sym_map = {name: sp.Symbol(name) for name in scan.relation.variables}
                        poly = scan.relation.as_sympy(tuple(sym_map[name] for name in scan.relation.variables))
                        lines.extend(
                            [
                                f"- Mixed-quotient polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                                "",
                                "```text",
                                _format_expr(poly),
                                "```",
                                "",
                                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                                "",
                            ]
                        )

                mixed_quotient_multiplicative_hits = [
                    scan for scan in family_scan.mixed_quotient_multiplicative_scans if scan.relation is not None
                ]
                mixed_quotient_multiplicative_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.mixed_quotient_multiplicative_scans
                    if scan.error is None and scan.relation is None
                ]
                if not mixed_quotient_multiplicative_hits and mixed_quotient_multiplicative_no_hits:
                    lines.append(
                        f"- Mixed-quotient multiplicative scan: no hit for prefixes ending at {', '.join(mixed_quotient_multiplicative_no_hits)}."
                    )
                for scan in family_scan.mixed_quotient_multiplicative_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Mixed-quotient multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.relation is None:
                        continue
                    residual = _multiplicative_relation_residual_series(
                        scan.relation,
                        target_series=ratio_series,
                        basis_series_by_variable={
                            label: mixed_quotient_basis_series_by_variable[label]
                            for label in scan.basis_labels
                        },
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- Mixed-quotient multiplicative prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_multiplicative_relation(scan.relation, target_variable="F"),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

                mixed_quotient_fractional_hits = [
                    scan for scan in family_scan.mixed_quotient_fractional_linear_scans if scan.relation is not None
                ]
                mixed_quotient_fractional_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.mixed_quotient_fractional_linear_scans
                    if scan.error is None and scan.relation is None
                ]
                if not mixed_quotient_fractional_hits and mixed_quotient_fractional_no_hits:
                    lines.append(
                        f"- Mixed-quotient fractional-linear scan: no hit for prefixes ending at {', '.join(mixed_quotient_fractional_no_hits)}."
                    )
                for scan in family_scan.mixed_quotient_fractional_linear_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Mixed-quotient fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.relation is None:
                        continue
                    residual = _fractional_linear_relation_residual_series(
                        scan.relation,
                        target_series=ratio_series,
                        basis_series_by_variable={
                            label: mixed_quotient_basis_series_by_variable[label]
                            for label in scan.basis_labels
                        },
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- Mixed-quotient fractional-linear prefix ending at `{scan.basis_labels[-1]}` produced a candidate relation:",
                            "",
                            "```text",
                            _format_fractional_linear_relation(scan.relation, target_variable="F"),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

                mixed_quotient_two_layer_hits = [
                    scan
                    for scan in family_scan.mixed_quotient_two_layer_fractional_linear_scans
                    if scan.total_hits > 0
                ]
                mixed_quotient_two_layer_no_hits = [
                    f"`{scan.basis_labels[-1]}`"
                    for scan in family_scan.mixed_quotient_two_layer_fractional_linear_scans
                    if scan.error is None and scan.total_hits == 0
                ]
                if not mixed_quotient_two_layer_hits and mixed_quotient_two_layer_no_hits:
                    lines.append(
                        f"- Mixed-quotient two-layer fractional-linear scan: no hit for prefixes ending at {', '.join(mixed_quotient_two_layer_no_hits)}."
                    )
                for scan in family_scan.mixed_quotient_two_layer_fractional_linear_scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Mixed-quotient two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )
                        continue
                    if scan.total_hits == 0:
                        continue
                    lines.append(
                        f"- Mixed-quotient two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` found `{scan.total_hits}` exact hit(s) after checking `{scan.tuples_checked}` template(s):"
                    )
                    lines.append("")
                    relation_basis = {
                        label: mixed_quotient_basis_series_by_variable[label]
                        for label in scan.basis_labels
                    }
                    for relation in scan.relations:
                        residual = _two_layer_fractional_linear_relation_residual_series(
                            relation,
                            target_series=ratio_series,
                            basis_series_by_variable=relation_basis,
                            order=profile_order,
                        )
                        residual_ok = all(sp.simplify(value) == 0 for value in residual)
                        lines.extend(
                            [
                                "```text",
                                _format_two_layer_fractional_linear_relation(
                                    relation,
                                    target_variable="F",
                                ),
                                "```",
                                "",
                                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                                "",
                            ]
                        )

    if source_family_eta_correction_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Source-Family Eta-Correction Scan",
                "",
                "We also checked whether the ratio object can be written as one nearby source-family basis object times a small eta-quotient correction:",
                "",
                "```text",
                "F = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a more direct closed-form recognition lane than the polynomial / fractional-linear boxes above.",
                f"- Eta levels checked: {', '.join(f'`N={level}`' for level in _eta_scan_levels(benchmark_powers))}",
                "",
            ]
        )

        for family_scan in source_family_eta_correction_scans:
            lines.extend(
                [
                    f"### `{family_scan.family_label}` Eta-Correction Box",
                    "",
                    f"- Base benchmark: `{family_scan.benchmark_name}`",
                ]
            )

            direct_labels = ", ".join(
                f"`{scan.basis_label}`" for scan in family_scan.direct_basis_scans
            )
            lines.append(f"- Raw basis choices: {direct_labels}")

            direct_hit_scans = [
                scan
                for scan in family_scan.direct_basis_scans
                if any(eta_scan.relation is not None for eta_scan in scan.eta_scans)
            ]
            direct_no_hits = [
                f"`{scan.basis_label}`"
                for scan in family_scan.direct_basis_scans
                if all(
                    eta_scan.error is None and eta_scan.relation is None
                    for eta_scan in scan.eta_scans
                )
            ]
            if not direct_hit_scans and direct_no_hits:
                lines.append(
                    f"- Raw-basis eta-correction scan: no hit for basis choices {', '.join(direct_no_hits)}."
                )

            for basis_scan in direct_hit_scans:
                correction_series = series_div(ratio_series, basis_scan.basis_series)
                for eta_scan in basis_scan.eta_scans:
                    if eta_scan.relation is None:
                        continue
                    residual = _multiplicative_relation_residual_series(
                        eta_scan.relation,
                        target_series=correction_series,
                        basis_series_by_variable=_eta_quotient_basis_series(
                            level=eta_scan.level,
                            order=profile_order,
                        ),
                        order=profile_order,
                    )
                    residual_ok = all(sp.simplify(value) == 0 for value in residual)
                    lines.extend(
                        [
                            f"- Raw basis `{basis_scan.basis_label}` with eta level `N={eta_scan.level}` produced a candidate correction:",
                            "",
                            "```text",
                            _format_source_family_eta_correction(
                                basis_expression=basis_scan.basis_expression,
                                relation=eta_scan.relation,
                                target_variable="F",
                                series_symbol=series_symbol,
                            ),
                            "```",
                            "",
                            f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                            "",
                        ]
                    )

            if family_scan.quotient_basis_scans:
                quotient_descriptions = ", ".join(
                    f"`{scan.basis_label} = {scan.basis_expression}`"
                    for scan in family_scan.quotient_basis_scans
                )
                lines.append(f"- Quotient basis choices: {quotient_descriptions}")

                quotient_hit_scans = [
                    scan
                    for scan in family_scan.quotient_basis_scans
                    if any(eta_scan.relation is not None for eta_scan in scan.eta_scans)
                ]
                quotient_no_hits = [
                    f"`{scan.basis_label}`"
                    for scan in family_scan.quotient_basis_scans
                    if all(
                        eta_scan.error is None and eta_scan.relation is None
                        for eta_scan in scan.eta_scans
                    )
                ]
                if not quotient_hit_scans and quotient_no_hits:
                    lines.append(
                        f"- Quotient-basis eta-correction scan: no hit for basis choices {', '.join(quotient_no_hits)}."
                    )

                for basis_scan in quotient_hit_scans:
                    correction_series = series_div(ratio_series, basis_scan.basis_series)
                    for eta_scan in basis_scan.eta_scans:
                        if eta_scan.relation is None:
                            continue
                        residual = _multiplicative_relation_residual_series(
                            eta_scan.relation,
                            target_series=correction_series,
                            basis_series_by_variable=_eta_quotient_basis_series(
                                level=eta_scan.level,
                                order=profile_order,
                            ),
                            order=profile_order,
                        )
                        residual_ok = all(sp.simplify(value) == 0 for value in residual)
                        lines.extend(
                            [
                                f"- Quotient basis `{basis_scan.basis_label}` with eta level `N={eta_scan.level}` produced a candidate correction:",
                                "",
                                "```text",
                                _format_source_family_eta_correction(
                                    basis_expression=basis_scan.basis_expression,
                                    relation=eta_scan.relation,
                                    target_variable="F",
                                    series_symbol=series_symbol,
                                ),
                                "```",
                                "",
                                f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                                "",
                            ]
                        )

    if two_core_source_family_eta_correction_scan.total_boxes_checked > 0:
        two_core_basis_series_by_label = {
            label: series
            for _, _, label, series in _source_family_raw_basis_entries(
                ordered_base_families=source_family_base_series,
                powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
                order=profile_order,
                supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
            )
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Core Source-Family Eta Scan",
                "",
                "We also checked a low-complexity hybrid source box built from two raw basis objects from different nearby families together with a small eta tail:",
                "",
                "```text",
                "F = T1 * T2 * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Here `T1` and `T2` come from distinct named-family raw ladders, and each source-core exponent is restricted to `±1`.",
                f"- Basis pairs checked: `{two_core_source_family_eta_correction_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_core_source_family_eta_correction_scan.total_boxes_checked}`",
                (
                    "- Family-pair split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_core_source_family_eta_correction_scan.family_pair_basis_counts
                    )
                    if two_core_source_family_eta_correction_scan.family_pair_basis_counts
                    else "- Family-pair split: none"
                ),
                "",
            ]
        )
        if not two_core_source_family_eta_correction_scan.hits:
            lines.append("No cross-family two-core eta-correction hit was found in the scanned box.")
            lines.append("")
        for hit in two_core_source_family_eta_correction_scan.hits:
            residual = _multiplicative_relation_residual_series(
                hit.relation,
                target_series=ratio_series,
                basis_series_by_variable={
                    hit.basis_labels[0]: two_core_basis_series_by_label[hit.basis_labels[0]],
                    hit.basis_labels[1]: two_core_basis_series_by_label[hit.basis_labels[1]],
                    **_eta_quotient_basis_series(level=hit.level, order=profile_order),
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Basis pair `{hit.basis_labels[0]}`, `{hit.basis_labels[1]}` with eta level `N={hit.level}` produced a candidate correction:",
                    "",
                    "```text",
                    _format_two_core_source_family_eta_correction(
                        relation=hit.relation,
                        target_variable="F",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if quotient_core_source_family_eta_correction_scan.total_boxes_checked > 0:
        quotient_core_basis_series_by_label = {
            label: series
            for _, _, label, _, series in _source_family_quotient_basis_entries(
                ordered_base_families=source_family_base_series,
                powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
                order=profile_order,
                supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
            )
        }
        raw_basis_series_by_label = {
            label: series
            for _, _, label, series in _source_family_raw_basis_entries(
                ordered_base_families=source_family_base_series,
                powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
                order=profile_order,
                supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
            )
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Quotient-Core Source-Family Eta Scan",
                "",
                "We also checked a hybrid source box where one nearby family contributes a quotient core and a second family contributes one raw basis object, again with a small eta tail:",
                "",
                "```text",
                "F = Q * T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Here `Q = T_k / T_1` comes from one named family, `T` comes from a different family's raw ladder, and both source-core exponents are restricted to `±1`.",
                f"- Quotient/raw basis pairs checked: `{quotient_core_source_family_eta_correction_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{quotient_core_source_family_eta_correction_scan.total_boxes_checked}`",
                (
                    "- Quotient/raw family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in quotient_core_source_family_eta_correction_scan.family_pair_basis_counts
                    )
                    if quotient_core_source_family_eta_correction_scan.family_pair_basis_counts
                    else "- Quotient/raw family split: none"
                ),
                "",
            ]
        )
        if not quotient_core_source_family_eta_correction_scan.hits:
            lines.append("No cross-family quotient-core eta-correction hit was found in the scanned box.")
            lines.append("")
        for hit in quotient_core_source_family_eta_correction_scan.hits:
            residual = _multiplicative_relation_residual_series(
                hit.relation,
                target_series=ratio_series,
                basis_series_by_variable={
                    hit.quotient_label: quotient_core_basis_series_by_label[hit.quotient_label],
                    hit.raw_label: raw_basis_series_by_label[hit.raw_label],
                    **_eta_quotient_basis_series(level=hit.level, order=profile_order),
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Quotient core `{hit.quotient_label} = {hit.quotient_expression}` with raw basis `{hit.raw_label}` and eta level `N={hit.level}` produced a candidate correction:",
                    "",
                    "```text",
                    _format_quotient_core_source_family_eta_correction(
                        quotient_label=hit.quotient_label,
                        quotient_expression=hit.quotient_expression,
                        raw_label=hit.raw_label,
                        raw_expression=hit.raw_expression,
                        relation=hit.relation,
                        target_variable="F",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if two_quotient_core_source_family_eta_correction_scan.total_boxes_checked > 0:
        quotient_core_basis_series_by_label = {
            label: series
            for _, _, label, _, series in _source_family_quotient_basis_entries(
                ordered_base_families=source_family_base_series,
                powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
                order=profile_order,
                supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
            )
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Quotient-Core Source-Family Eta Scan",
                "",
                "We also checked a quotient-only hybrid source box where two distinct nearby families each contribute one quotient core, again with a small eta tail:",
                "",
                "```text",
                "F = Q1 * Q2 * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Here `Q1 = T_i / T_1` and `Q2 = U_j / U_1` come from distinct named families, and both quotient-core exponents are restricted to `±1`.",
                f"- Quotient-pair basis pairs checked: `{two_quotient_core_source_family_eta_correction_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_quotient_core_source_family_eta_correction_scan.total_boxes_checked}`",
                (
                    "- Quotient-pair family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_quotient_core_source_family_eta_correction_scan.family_pair_basis_counts
                    )
                    if two_quotient_core_source_family_eta_correction_scan.family_pair_basis_counts
                    else "- Quotient-pair family split: none"
                ),
                "",
            ]
        )
        if not two_quotient_core_source_family_eta_correction_scan.hits:
            lines.append("No cross-family two-quotient-core eta-correction hit was found in the scanned box.")
            lines.append("")
        for hit in two_quotient_core_source_family_eta_correction_scan.hits:
            residual = _multiplicative_relation_residual_series(
                hit.relation,
                target_series=ratio_series,
                basis_series_by_variable={
                    hit.quotient_labels[0]: quotient_core_basis_series_by_label[hit.quotient_labels[0]],
                    hit.quotient_labels[1]: quotient_core_basis_series_by_label[hit.quotient_labels[1]],
                    **_eta_quotient_basis_series(level=hit.level, order=profile_order),
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Quotient pair `{hit.quotient_labels[0]}`, `{hit.quotient_labels[1]}` with eta level `N={hit.level}` produced a candidate correction:",
                    "",
                    "```text",
                    _format_two_quotient_core_source_family_eta_correction(
                        quotient_labels=hit.quotient_labels,
                        quotient_expressions=hit.quotient_expressions,
                        relation=hit.relation,
                        target_variable="F",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if two_quotient_core_source_family_self_quotient_product_scan.total_boxes_checked > 0:
        quotient_core_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=source_family_base_series,
            powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
            order=profile_order,
            supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
        )
        quotient_core_basis_series_by_label = {
            label: series for _, _, label, _, series in quotient_core_basis_entries
        }
        quotient_core_expression_by_label = {
            label: expression for _, _, label, expression, _ in quotient_core_basis_entries
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Quotient-Core Self-Quotient Finite-Product Scan",
                "",
                "We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object with a compact finite-product self-quotient equation:",
                "",
                "```text",
                "G = F / (Q1 * Q2)",
                "G(t) / G(t^m) = prod_{r=1}^{m-1} (1 - t^r)^{e_r}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a source-aware Mahler-style lane: a hit would point to a recursive product correction after factoring out two nearby quotient cores.",
                f"- Moduli checked: {', '.join(f'`m={modulus}`' for modulus in two_quotient_core_source_family_self_quotient_product_scan.moduli_checked)}",
                f"- Quotient-pair basis pairs checked: `{two_quotient_core_source_family_self_quotient_product_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_quotient_core_source_family_self_quotient_product_scan.total_boxes_checked}`",
                (
                    "- Quotient-pair family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_quotient_core_source_family_self_quotient_product_scan.family_pair_basis_counts
                    )
                    if two_quotient_core_source_family_self_quotient_product_scan.family_pair_basis_counts
                    else "- Quotient-pair family split: none"
                ),
                "",
            ]
        )
        if not two_quotient_core_source_family_self_quotient_product_scan.hits:
            lines.append(
                "No cross-family two-quotient-core finite-product self-quotient hit was found in the scanned box."
            )
            lines.append("")
        for hit in two_quotient_core_source_family_self_quotient_product_scan.hits:
            correction_series = series_div(
                ratio_series,
                series_mul(
                    quotient_core_basis_series_by_label[hit.quotient_labels[0]],
                    quotient_core_basis_series_by_label[hit.quotient_labels[1]],
                ),
            )
            residual = _self_quotient_product_relation_residual_series(
                hit.relation,
                target_series=correction_series,
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Quotient pair `{hit.quotient_labels[0]}`, `{hit.quotient_labels[1]}` with self-quotient modulus `m={hit.modulus}` produced a candidate finite-product correction:",
                    "",
                    "```text",
                    "G = F / "
                    + f"(({quotient_core_expression_by_label[hit.quotient_labels[0]]}) * ({quotient_core_expression_by_label[hit.quotient_labels[1]]}))",
                    _format_self_quotient_product_relation(
                        hit.relation,
                        target_variable="G",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if two_quotient_core_source_family_self_eta_scan.total_boxes_checked > 0:
        quotient_core_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=source_family_base_series,
            powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
            order=profile_order,
            supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
        )
        quotient_core_basis_series_by_label = {
            label: series for _, _, label, _, series in quotient_core_basis_entries
        }
        quotient_core_expression_by_label = {
            label: expression for _, _, label, expression, _ in quotient_core_basis_entries
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Quotient-Core Self-Eta Functional Scan",
                "",
                "We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-complexity self-eta functional equation:",
                "",
                "```text",
                "G = F / (Q1 * Q2)",
                "G(t) = G(t^m)^a * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a product-theoretic theorem-facing lane: a hit would identify the residual correction through a recursive self-eta equation after factoring out two nearby quotient cores.",
                f"- Moduli checked: {', '.join(f'`m={modulus}`' for modulus in two_quotient_core_source_family_self_eta_scan.moduli_checked)}",
                f"- Eta levels checked: {', '.join(f'`N={level}`' for level in two_quotient_core_source_family_self_eta_scan.levels_checked)}",
                f"- Quotient-pair basis pairs checked: `{two_quotient_core_source_family_self_eta_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_quotient_core_source_family_self_eta_scan.total_boxes_checked}`",
                (
                    "- Quotient-pair family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_quotient_core_source_family_self_eta_scan.family_pair_basis_counts
                    )
                    if two_quotient_core_source_family_self_eta_scan.family_pair_basis_counts
                    else "- Quotient-pair family split: none"
                ),
                "",
            ]
        )
        if not two_quotient_core_source_family_self_eta_scan.hits:
            lines.append(
                "No cross-family two-quotient-core self-eta functional-equation hit was found in the scanned box."
            )
            lines.append("")
        for hit in two_quotient_core_source_family_self_eta_scan.hits:
            correction_series = series_div(
                ratio_series,
                series_mul(
                    quotient_core_basis_series_by_label[hit.quotient_labels[0]],
                    quotient_core_basis_series_by_label[hit.quotient_labels[1]],
                ),
            )
            powered_correction = benchmark_power_substitution_series(
                correction_series,
                power=hit.modulus,
                order=profile_order,
            )
            residual = _multiplicative_relation_residual_series(
                hit.relation,
                target_series=correction_series,
                basis_series_by_variable={
                    f"G{hit.modulus}": powered_correction,
                    **_eta_quotient_basis_series(level=hit.level, order=profile_order),
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Quotient pair `{hit.quotient_labels[0]}`, `{hit.quotient_labels[1]}` with self-functional modulus `m={hit.modulus}` and eta level `N={hit.level}` produced a candidate relation:",
                    "",
                    "```text",
                    "G = F / "
                    + f"(({quotient_core_expression_by_label[hit.quotient_labels[0]]}) * ({quotient_core_expression_by_label[hit.quotient_labels[1]]}))",
                    _format_self_eta_correction(
                        relation=hit.relation,
                        modulus=hit.modulus,
                        target_variable="G",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if two_quotient_core_source_family_self_fractional_linear_scan.total_boxes_checked > 0:
        quotient_core_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=source_family_base_series,
            powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
            order=profile_order,
            supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
        )
        quotient_core_basis_series_by_label = {
            label: series for _, _, label, _, series in quotient_core_basis_entries
        }
        quotient_core_expression_by_label = {
            label: expression for _, _, label, expression, _ in quotient_core_basis_entries
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Quotient-Core Self-Fractional-Linear Scan",
                "",
                "We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-complexity self-fractional-linear equation with a small eta tail:",
                "",
                "```text",
                "G = F / (Q1 * Q2)",
                "G(t) = (1 + a*(G(t^m) - 1) + ... ) / (1 + b*(G(t^m) - 1) + ... )",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a theorem-facing nonlinear lane: a hit would give a compact recursive rational equation for the residual correction after factoring out two nearby quotient cores.",
                f"- Moduli checked: {', '.join(f'`m={modulus}`' for modulus in two_quotient_core_source_family_self_fractional_linear_scan.moduli_checked)}",
                f"- Eta levels checked: {', '.join(f'`N={level}`' for level in two_quotient_core_source_family_self_fractional_linear_scan.levels_checked)}",
                f"- Quotient-pair basis pairs checked: `{two_quotient_core_source_family_self_fractional_linear_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_quotient_core_source_family_self_fractional_linear_scan.total_boxes_checked}`",
                (
                    "- Quotient-pair family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_quotient_core_source_family_self_fractional_linear_scan.family_pair_basis_counts
                    )
                    if two_quotient_core_source_family_self_fractional_linear_scan.family_pair_basis_counts
                    else "- Quotient-pair family split: none"
                ),
                "",
            ]
        )
        if not two_quotient_core_source_family_self_fractional_linear_scan.hits:
            lines.append(
                "No cross-family two-quotient-core self-fractional-linear hit was found in the scanned box."
            )
            lines.append("")
        for hit in two_quotient_core_source_family_self_fractional_linear_scan.hits:
            correction_series = series_div(
                ratio_series,
                series_mul(
                    quotient_core_basis_series_by_label[hit.quotient_labels[0]],
                    quotient_core_basis_series_by_label[hit.quotient_labels[1]],
                ),
            )
            residual = _fractional_linear_relation_residual_series(
                hit.relation,
                target_series=correction_series,
                basis_series_by_variable={
                    f"G{hit.modulus}": benchmark_power_substitution_series(
                        correction_series,
                        power=hit.modulus,
                        order=profile_order,
                    ),
                    **_eta_quotient_basis_series(level=hit.level, order=profile_order),
                },
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Quotient pair `{hit.quotient_labels[0]}`, `{hit.quotient_labels[1]}` with self-functional modulus `m={hit.modulus}` and eta level `N={hit.level}` produced a candidate relation:",
                    "",
                    "```text",
                    "G = F / "
                    + f"(({quotient_core_expression_by_label[hit.quotient_labels[0]]}) * ({quotient_core_expression_by_label[hit.quotient_labels[1]]}))",
                    _format_fractional_linear_relation(hit.relation, target_variable="G"),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if two_quotient_core_source_family_self_polynomial_scan.total_boxes_checked > 0:
        quotient_core_basis_entries = _source_family_quotient_basis_entries(
            ordered_base_families=source_family_base_series,
            powers=_parameterized_source_family_powers(benchmark_powers, smoke=smoke),
            order=profile_order,
            supplemental_powers_by_family=_supplemental_source_family_powers(smoke=smoke),
        )
        quotient_core_basis_series_by_label = {
            label: series for _, _, label, _, series in quotient_core_basis_entries
        }
        quotient_core_expression_by_label = {
            label: expression for _, _, label, expression, _ in quotient_core_basis_entries
        }
        lines.extend(
            [
                "",
                "## Ratio-Object Two-Quotient-Core Self-Polynomial Functional Scan",
                "",
                "We also checked whether dividing by one quotient core from each of two distinct nearby families leaves a correction object satisfying a low-degree algebraic self-functional equation:",
                "",
                "```text",
                "G = F / (Q1 * Q2)",
                "P(G(t), G(t^m)) = 0",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a theorem-facing lane: a hit would suggest a compact defining functional equation for the residual correction after factoring out two nearby quotient cores.",
                f"- Moduli checked: {', '.join(f'`m={modulus}`' for modulus in two_quotient_core_source_family_self_polynomial_scan.moduli_checked)}",
                f"- Degrees checked: {', '.join(f'`total degree <= {degree}`' for degree in two_quotient_core_source_family_self_polynomial_scan.degree_values)}",
                f"- Quotient-pair basis pairs checked: `{two_quotient_core_source_family_self_polynomial_scan.total_basis_pairs_checked}`",
                f"- Total pair-level boxes checked: `{two_quotient_core_source_family_self_polynomial_scan.total_boxes_checked}`",
                (
                    "- Quotient-pair family split: "
                    + ", ".join(
                        f"`{label}` -> `{count}` pair(s)"
                        for label, count in two_quotient_core_source_family_self_polynomial_scan.family_pair_basis_counts
                    )
                    if two_quotient_core_source_family_self_polynomial_scan.family_pair_basis_counts
                    else "- Quotient-pair family split: none"
                ),
                "",
            ]
        )
        if not two_quotient_core_source_family_self_polynomial_scan.hits:
            lines.append(
                "No cross-family two-quotient-core self-polynomial functional-equation hit was found in the scanned box."
            )
            lines.append("")
        for hit in two_quotient_core_source_family_self_polynomial_scan.hits:
            correction_series = series_div(
                ratio_series,
                series_mul(
                    quotient_core_basis_series_by_label[hit.quotient_labels[0]],
                    quotient_core_basis_series_by_label[hit.quotient_labels[1]],
                ),
            )
            series_by_variable = {
                "G": correction_series,
                f"G{hit.modulus}": benchmark_power_substitution_series(
                    correction_series,
                    power=hit.modulus,
                    order=profile_order,
                ),
            }
            residual = _relation_residual_series(
                hit.relation,
                series_by_variable=series_by_variable,
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            sym_map = {name: sp.Symbol(name) for name in hit.relation.variables}
            poly = hit.relation.as_sympy(tuple(sym_map[name] for name in hit.relation.variables))
            lines.extend(
                [
                    f"- Quotient pair `{hit.quotient_labels[0]}`, `{hit.quotient_labels[1]}` with self-functional modulus `m={hit.modulus}` and polynomial `total degree <= {hit.max_total_degree}` produced a candidate relation:",
                    "",
                    "```text",
                    "G = F / "
                    + f"(({quotient_core_expression_by_label[hit.quotient_labels[0]]}) * ({quotient_core_expression_by_label[hit.quotient_labels[1]]}))",
                    f"{_format_expr(poly)} = 0",
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if explicit_source_family_transform_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Explicit GG/S Transform Template Scan",
                "",
                "We also checked a smaller family-meaning-preserving box tailored to the Gordon/Hirschhorn orbit:",
                "",
                "```text",
                "F = T",
                "F = 1 / T",
                "F = T_i / T_j",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Families checked here are the literature-family ladders `GG` and `S`.",
                "- This does not enlarge the algebraic search box much; it makes the reciprocal / quotient interpretations explicit in the note.",
                "",
            ]
        )

        for family_scan in explicit_source_family_transform_scans:
            basis_descriptions = [f"`{family_scan.family_label} = {family_scan.family_label}({series_symbol})`"]
            for label, _ in family_scan.ordered_basis_series[1:]:
                power = int(label.removeprefix(family_scan.family_label))
                basis_descriptions.append(
                    f"`{label} = {family_scan.family_label}({series_symbol}^{power})`"
                )
            lines.extend(
                [
                    f"### `{family_scan.family_label}` Explicit Transform Box",
                    "",
                    f"- Base benchmark: `{family_scan.benchmark_name}`",
                    f"- Basis ladder: {', '.join(basis_descriptions)}",
                    f"- Templates checked: `{len(family_scan.checked_templates)}` exact direct / reciprocal / quotient templates.",
                ]
            )
            if family_scan.hit_templates:
                lines.append(
                    f"- Exact hit(s): {', '.join(f'`{label}`' for label in family_scan.hit_templates)}."
                )
            else:
                lines.append("- No exact direct / reciprocal / quotient template hit was found in this family.")
            lines.append("")

    if gg_modular_equation_scan is not None:
        lines.extend(
            [
                "",
                "## Ratio-Object GG Modular-Equation Template Scan",
                "",
                "We also checked a narrower literature-driven `GG` box motivated by the modular-equation papers of Chan--Huang and Cho--Koo--Park:",
                "",
                "```text",
                "F = T",
                "F = 1 / T",
                "F = T_i / T_j",
                "P(F, T_i) = 0",
                "F = prod_i T_i^e_i",
                "F = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                f"- Base benchmark: `{gg_modular_equation_scan.benchmark_name}`",
                "- This lane keeps the sign and substitution objects explicit instead of flattening them into a larger anonymous basis box.",
                "- The literature-motivated basis here starts with `GG(t)`, `GG(-t)`, `GG(t^2)`, `GG(t^3)`, and `GG(t^4)`, and in the full profile it also includes the odd-prime descendants suggested by the GG modular-equation papers.",
                "- We also run a mixed quotient-coordinate pass that keeps `GG(t)` explicit while letting the correction move in quotient coordinates such as `GG(-t)/GG(t)` and `GG(t^p)/GG(t)`.",
            ]
        )
        basis_descriptions = [
            f"`{label} = {expression}`"
            for label, expression, _ in gg_modular_equation_scan.ordered_basis_series
        ]
        lines.append(f"- Basis ladder: {', '.join(basis_descriptions)}")
        lines.append(
            f"- Exact direct / reciprocal / quotient templates checked: `{len(gg_modular_equation_scan.checked_templates)}`."
        )
        if gg_modular_equation_scan.hit_templates:
            lines.append(
                f"- Exact template hit(s): {', '.join(f'`{label}`' for label in gg_modular_equation_scan.hit_templates)}."
            )
        else:
            lines.append("- No exact direct / reciprocal / quotient template hit was found in this modular-equation box.")
        if gg_modular_equation_scan.exact_polynomial_template_labels:
            lines.append(
                "- Exact literature polynomial templates checked: "
                + ", ".join(f"`{label}`" for label in gg_modular_equation_scan.exact_polynomial_template_labels)
                + "."
            )
            if gg_modular_equation_scan.exact_polynomial_template_hits:
                lines.append(
                    "- Exact literature polynomial hit(s): "
                    + ", ".join(f"`{label}`" for label in gg_modular_equation_scan.exact_polynomial_template_hits)
                    + "."
                )
            else:
                lines.append("- No exact Chan--Huang direct modular-equation polynomial hit was found.")
        lines.append("")

        grouped_gg_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
        for scan in gg_modular_equation_scan.polynomial_scans:
            grouped_gg_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

        any_gg_polynomial_hit = any(
            scan.relation is not None for scan in gg_modular_equation_scan.polynomial_scans
        )
        if not any_gg_polynomial_hit:
            lines.append("- Polynomial scan: no candidate-dependent hit was found in the checked modular-equation prefixes.")
        for degree in sorted(grouped_gg_polynomial_scans):
            scans = grouped_gg_polynomial_scans[degree]
            prefix_labels = [f"`{scan.basis_labels[-1]}`" for scan in scans if scan.error is None]
            if not any_gg_polynomial_hit and prefix_labels:
                lines.append(
                    f"- Polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(prefix_labels)}."
                )
            for scan in scans:
                if scan.error is not None:
                    lines.append(
                        f"- Polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )

        any_gg_multiplicative_hit = any(
            scan.relation is not None for scan in gg_modular_equation_scan.multiplicative_scans
        )
        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in gg_modular_equation_scan.multiplicative_scans
            if scan.error is None and scan.relation is None
        ]
        if not any_gg_multiplicative_hit and no_hit_labels:
            lines.append(
                f"- Multiplicative scan: no hit for modular-equation prefixes ending at {', '.join(no_hit_labels)}."
            )
        for scan in gg_modular_equation_scan.multiplicative_scans:
            if scan.error is not None:
                lines.append(
                    f"- Multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )

        any_gg_fractional_hit = any(
            scan.relation is not None for scan in gg_modular_equation_scan.fractional_linear_scans
        )
        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in gg_modular_equation_scan.fractional_linear_scans
            if scan.error is None and scan.relation is None
        ]
        if not any_gg_fractional_hit and no_hit_labels:
            lines.append(
                f"- Fractional-linear scan: no hit for modular-equation prefixes ending at {', '.join(no_hit_labels)}."
            )
        for scan in gg_modular_equation_scan.fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )

        any_gg_two_layer_hit = any(
            scan.total_hits > 0 for scan in gg_modular_equation_scan.two_layer_fractional_linear_scans
        )
        no_hit_labels = [
            f"`{scan.basis_labels[-1]}`"
            for scan in gg_modular_equation_scan.two_layer_fractional_linear_scans
            if scan.error is None and scan.total_hits == 0
        ]
        if not any_gg_two_layer_hit and no_hit_labels:
            lines.append(
                f"- Two-layer fractional-linear scan: no hit for modular-equation prefixes ending at {', '.join(no_hit_labels)}."
            )
        for scan in gg_modular_equation_scan.two_layer_fractional_linear_scans:
            if scan.error is not None:
                lines.append(
                    f"- Two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                )
        if gg_modular_equation_scan.quotient_basis_series:
            quotient_descriptions = [
                f"`{label} = {expression}`"
                for label, expression, _ in gg_modular_equation_scan.quotient_basis_series
            ]
            lines.append(f"- Quotient basis: {', '.join(quotient_descriptions)}")
            if gg_modular_equation_scan.quotient_exact_polynomial_template_labels:
                lines.append(
                    "- Exact quotient-coordinate literature templates checked: "
                    + ", ".join(
                        f"`{label}`"
                        for label in gg_modular_equation_scan.quotient_exact_polynomial_template_labels
                    )
                    + "."
                )
                if gg_modular_equation_scan.quotient_exact_polynomial_template_hits:
                    lines.append(
                        "- Exact quotient-coordinate literature hit(s): "
                        + ", ".join(
                            f"`{label}`"
                            for label in gg_modular_equation_scan.quotient_exact_polynomial_template_hits
                        )
                        + "."
                    )
                else:
                    lines.append(
                        "- No exact Chan--Huang quotient-coordinate modular-equation polynomial hit was found."
                    )

            grouped_gg_quotient_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
            for scan in gg_modular_equation_scan.quotient_polynomial_scans:
                grouped_gg_quotient_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

            any_gg_quotient_polynomial_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.quotient_polynomial_scans
            )
            if not any_gg_quotient_polynomial_hit:
                lines.append("- Quotient-coordinate polynomial scan: no candidate-dependent hit was found in the checked quotient prefixes.")
            for degree in sorted(grouped_gg_quotient_polynomial_scans):
                scans = grouped_gg_quotient_polynomial_scans[degree]
                prefix_labels = [f"`{scan.basis_labels[-1]}`" for scan in scans if scan.error is None]
                if not any_gg_quotient_polynomial_hit and prefix_labels:
                    lines.append(
                        f"- Quotient-coordinate polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(prefix_labels)}."
                    )
                for scan in scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Quotient-coordinate polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )

            any_gg_quotient_multiplicative_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.quotient_multiplicative_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.quotient_multiplicative_scans
                if scan.error is None and scan.relation is None
            ]
            if not any_gg_quotient_multiplicative_hit and no_hit_labels:
                lines.append(
                    f"- Quotient-coordinate multiplicative scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.quotient_multiplicative_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Quotient-coordinate multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )

            any_gg_quotient_fractional_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.quotient_fractional_linear_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.quotient_fractional_linear_scans
                if scan.error is None and scan.relation is None
            ]
            if not any_gg_quotient_fractional_hit and no_hit_labels:
                lines.append(
                    f"- Quotient-coordinate fractional-linear scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.quotient_fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Quotient-coordinate fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )

            any_gg_quotient_two_layer_hit = any(
                scan.total_hits > 0 for scan in gg_modular_equation_scan.quotient_two_layer_fractional_linear_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.quotient_two_layer_fractional_linear_scans
                if scan.error is None and scan.total_hits == 0
            ]
            if not any_gg_quotient_two_layer_hit and no_hit_labels:
                lines.append(
                    f"- Quotient-coordinate two-layer fractional-linear scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.quotient_two_layer_fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Quotient-coordinate two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )
            mixed_descriptions = [
                f"`{label} = {expression}`"
                for label, expression, _ in gg_modular_equation_scan.mixed_quotient_basis_series
            ]
            lines.append(f"- Mixed quotient basis: {', '.join(mixed_descriptions)}")

            grouped_gg_mixed_polynomial_scans: dict[int, list[NamedPolynomialRelationScan]] = {}
            for scan in gg_modular_equation_scan.mixed_quotient_polynomial_scans:
                grouped_gg_mixed_polynomial_scans.setdefault(scan.max_total_degree, []).append(scan)

            any_gg_mixed_polynomial_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.mixed_quotient_polynomial_scans
            )
            if not any_gg_mixed_polynomial_hit:
                lines.append("- Mixed quotient-coordinate polynomial scan: no candidate-dependent hit was found in the checked prefixes.")
            for degree in sorted(grouped_gg_mixed_polynomial_scans):
                scans = grouped_gg_mixed_polynomial_scans[degree]
                prefix_labels = [f"`{scan.basis_labels[-1]}`" for scan in scans if scan.error is None]
                if not any_gg_mixed_polynomial_hit and prefix_labels:
                    lines.append(
                        f"- Mixed quotient-coordinate polynomial `total degree <= {degree}`: no hit for prefixes ending at {', '.join(prefix_labels)}."
                    )
                for scan in scans:
                    if scan.error is not None:
                        lines.append(
                            f"- Mixed quotient-coordinate polynomial `total degree <= {degree}` prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                        )

            any_gg_mixed_multiplicative_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.mixed_quotient_multiplicative_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.mixed_quotient_multiplicative_scans
                if scan.error is None and scan.relation is None
            ]
            if not any_gg_mixed_multiplicative_hit and no_hit_labels:
                lines.append(
                    f"- Mixed quotient-coordinate multiplicative scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.mixed_quotient_multiplicative_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Mixed quotient-coordinate multiplicative prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )

            any_gg_mixed_fractional_hit = any(
                scan.relation is not None for scan in gg_modular_equation_scan.mixed_quotient_fractional_linear_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.mixed_quotient_fractional_linear_scans
                if scan.error is None and scan.relation is None
            ]
            if not any_gg_mixed_fractional_hit and no_hit_labels:
                lines.append(
                    f"- Mixed quotient-coordinate fractional-linear scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.mixed_quotient_fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Mixed quotient-coordinate fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )

            any_gg_mixed_two_layer_hit = any(
                scan.total_hits > 0 for scan in gg_modular_equation_scan.mixed_quotient_two_layer_fractional_linear_scans
            )
            no_hit_labels = [
                f"`{scan.basis_labels[-1]}`"
                for scan in gg_modular_equation_scan.mixed_quotient_two_layer_fractional_linear_scans
                if scan.error is None and scan.total_hits == 0
            ]
            if not any_gg_mixed_two_layer_hit and no_hit_labels:
                lines.append(
                    f"- Mixed quotient-coordinate two-layer fractional-linear scan: no hit for prefixes ending at {', '.join(no_hit_labels)}."
                )
            for scan in gg_modular_equation_scan.mixed_quotient_two_layer_fractional_linear_scans:
                if scan.error is not None:
                    lines.append(
                        f"- Mixed quotient-coordinate two-layer fractional-linear prefix ending at `{scan.basis_labels[-1]}` skipped: {scan.error}"
                    )
        lines.append("")

    if explicit_source_family_eta_correction_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Explicit GG/S Template Eta-Correction Scan",
                "",
                "We also checked whether one explicit Gordon/Hirschhorn-orbit template times a small eta tail explains the ratio object:",
                "",
                "```text",
                "F = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- Here `T` ranges over the exact direct / reciprocal / quotient templates from the preceding GG/S transform box.",
                f"- Eta levels checked: {', '.join(f'`N={level}`' for level in _eta_scan_levels(benchmark_powers))}",
                "",
            ]
        )

        for family_scan in explicit_source_family_eta_correction_scans:
            lines.extend(
                [
                    f"### `{family_scan.family_label}` Explicit Eta-Correction Box",
                    "",
                    f"- Base benchmark: `{family_scan.benchmark_name}`",
                    f"- Templates checked: `{len(family_scan.checked_templates)}` explicit direct / reciprocal / quotient templates.",
                ]
            )
            if not family_scan.hits:
                lines.append("- No explicit-template eta-correction hit was found in this family.")
                lines.append("")
                continue

            for hit in family_scan.hits:
                template_series = next(
                    series
                    for label, series in _explicit_source_family_template_series(family_scan.ordered_basis_series)
                    if label == hit.template_label
                )
                correction_series = series_div(ratio_series, template_series)
                residual = _multiplicative_relation_residual_series(
                    hit.relation,
                    target_series=correction_series,
                    basis_series_by_variable=_eta_quotient_basis_series(
                        level=hit.level,
                        order=profile_order,
                    ),
                    order=profile_order,
                )
                residual_ok = all(sp.simplify(value) == 0 for value in residual)
                lines.extend(
                    [
                        f"- Template `{hit.template_label}` with eta level `N={hit.level}` produced a candidate correction:",
                        "",
                        "```text",
                        _format_source_family_eta_correction(
                            basis_expression=hit.template_label,
                            relation=hit.relation,
                            target_variable="F",
                            series_symbol=series_symbol,
                        ),
                        "```",
                        "",
                        f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
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

    if ratio_self_quotient_product_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Self-Quotient Finite-Product Scan",
                "",
                "We also checked a simple finite-product self-quotient box for the ratio object:",
                "",
                "```text",
                "F(t) / F(t^m) = prod_{r=1}^{m-1} (1 - t^r)^{e_r}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a Mahler-style finite-product functional equation: a hit would give a compact recursive product description, but a miss does not rule out general q-Pochhammer products.",
                "",
            ]
        )

        any_self_quotient_hit = any(scan.relation is not None for scan in ratio_self_quotient_product_scans)
        if not any_self_quotient_hit:
            lines.append("No finite-product self-quotient relation was found in any scanned modulus.")
            lines.append("")

        no_hit_labels = [f"`m={scan.modulus}`" for scan in ratio_self_quotient_product_scans if scan.error is None and scan.relation is None]
        if not any_self_quotient_hit and no_hit_labels:
            lines.append(
                f"- No hit for moduli {', '.join(no_hit_labels)}."
            )

        for scan in ratio_self_quotient_product_scans:
            if scan.error is not None:
                lines.append(
                    f"- Self-quotient modulus `m={scan.modulus}` skipped: {scan.error}"
                )
                continue
            if scan.relation is None:
                continue
            residual = _self_quotient_product_relation_residual_series(
                scan.relation,
                target_series=ratio_series,
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Self-quotient modulus `m={scan.modulus}` produced a candidate finite-product relation:",
                    "",
                    "```text",
                    _format_self_quotient_product_relation(
                        scan.relation,
                        target_variable="F",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if ratio_eta_quotient_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Eta-Quotient Scan",
                "",
                "We also checked whether the ratio object itself is already a small-level eta-quotient:",
                "",
                "```text",
                "F = prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a direct closed-form recognition lane rather than another transform-elimination box.",
                "",
            ]
        )

        any_eta_hit = any(scan.relation is not None for scan in ratio_eta_quotient_scans)
        if not any_eta_hit:
            lines.append("No eta-quotient relation was found in any scanned level.")
            lines.append("")

        no_hit_labels = [
            f"`N={scan.level}`"
            for scan in ratio_eta_quotient_scans
            if scan.error is None and scan.relation is None
        ]
        if not any_eta_hit and no_hit_labels:
            lines.append(f"- No hit for eta levels {', '.join(no_hit_labels)}.")

        for scan in ratio_eta_quotient_scans:
            if scan.error is not None:
                lines.append(
                    f"- Eta-quotient level `N={scan.level}` skipped: {scan.error}"
                )
                continue
            if scan.relation is None:
                continue
            residual = _multiplicative_relation_residual_series(
                scan.relation,
                target_series=ratio_series,
                basis_series_by_variable=_eta_quotient_basis_series(
                    level=scan.level,
                    order=profile_order,
                ),
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Eta-quotient level `N={scan.level}` produced a candidate relation:",
                    "",
                    "```text",
                    _format_eta_quotient_relation(
                        scan.relation,
                        target_variable="F",
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{profile_order}`: `{residual_ok}`",
                    "",
                ]
            )

    if ratio_multiplicative_scans:
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
    final_render_elapsed = perf_counter() - final_render_started_at
    stage_elapsed_seconds["final-render"] = final_render_elapsed
    for index, line in enumerate(lines):
        if line == "- `final-render`: `in_progress`":
            lines[index] = f"- `final-render`: `{final_render_elapsed:.2f}`"
            break

    progress_status["final-render"] = "completed"
    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
