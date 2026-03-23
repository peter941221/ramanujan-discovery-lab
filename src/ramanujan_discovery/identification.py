from __future__ import annotations

from dataclasses import dataclass, replace
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


@lru_cache(maxsize=None)
def _canonical_benchmark_series(
    benchmark_name: str,
    *,
    depth: int,
    order: int,
) -> tuple[sp.Expr, ...]:
    return tuple(
        continued_fraction_series_coeffs(
            get_benchmark(benchmark_name).canonical_template.normalized(),
            depth=depth,
            order=order,
        )
    )


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
class SelfPlusPochhammerRelationScan:
    """Outcome for one modulus scan of a periodic plus-Pochhammer self-transfer box."""

    modulus: int
    relation: MultiplicativeRelation | None
    error: str | None = None


@dataclass(frozen=True)
class SelfPlusPochhammerEtaRelationScan:
    """Outcome for one modulus/level scan of a plus-Pochhammer + eta self-transfer box."""

    modulus: int
    level: int
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
class ReducedTailAnchor:
    """A concrete tail object obtained by anchoring the stationary tail law at one stage."""

    start_stage: int
    state_expr: sp.Expr
    tail_series: tuple[sp.Expr, ...]
    normalized_series: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class GapNormalizedSeries:
    """A constant-1 series normalized by its first nonzero gap after the constant term."""

    shift: int
    leading_coefficient: sp.Expr
    normalized_series: tuple[sp.Expr, ...]


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
class ModularUnitEtaRelationScan:
    """Outcome for one modular-unit / eta scan."""

    modulus: int
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
class SelfPlusPochhammerEtaCorrectionBasisScan:
    """Self plus-Pochhammer + eta scans for one fixed source-family basis choice."""

    basis_label: str
    basis_expression: str
    basis_series: Series
    self_scans: tuple[SelfPlusPochhammerEtaRelationScan, ...]


@dataclass(frozen=True)
class SourceFamilySelfPlusPochhammerEtaCorrectionScan:
    """Per-family source-core self plus-Pochhammer + eta scans over raw and quotient basis choices."""

    family_label: str
    benchmark_name: str
    direct_basis_scans: tuple[SelfPlusPochhammerEtaCorrectionBasisScan, ...]
    quotient_basis_scans: tuple[SelfPlusPochhammerEtaCorrectionBasisScan, ...]


@dataclass(frozen=True)
class TailFamilySourceEtaSample:
    """One tail-family specialization/gap-depth sample with nearby source-core eta scans."""

    label: str
    expression: str
    start_stage: int
    gap_depth: int
    state_expr: sp.Expr
    series: tuple[sp.Expr, ...]
    source_family_eta_scans: tuple[SourceFamilyEtaCorrectionScan, ...]
    direct_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    direct_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    gg_modular_equation_scan: GGModularEquationScan | None = None
    morton_periodic_point_scan: MortonPeriodicPointScan | None = None
    weber_g_class_invariant_scan: WeberClassInvariantScan | None = None
    weber_p_class_invariant_scan: WeberClassInvariantScan | None = None
    weber_residual_bridge_scan: WeberResidualBridgeScan | None = None
    weber_j_pb_bridge_scan: "ConstantOnePairBridgeScan | None" = None
    weber_j_lift_pivot_bridge_scans: tuple[ConstantOnePairBridgeLiteScan, ...] = ()


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
class NamedCoordinateOrbitScan:
    """A focused one-coordinate source ladder with direct, quotient, and mixed scans."""

    family_label: str
    base_label: str
    base_expression: str
    ordered_basis_series: tuple[tuple[str, str, Series], ...]
    direct_scans: NamedPrefixBoxScans
    quotient_basis_series: tuple[tuple[str, str, Series], ...]
    quotient_scans: NamedPrefixBoxScans
    mixed_quotient_basis_series: tuple[tuple[str, str, Series], ...]
    mixed_quotient_scans: NamedPrefixBoxScans


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
class GGWeightedCoordinateDiagnostic:
    """A weighted quotient-coordinate diagnostic inside the GG modular-equation lane."""

    label: str
    expression: str
    log_expression: str
    correction_expression: str
    first_difference_power: int | None
    first_difference_coeff: sp.Expr | None
    first_log_difference_power: int | None
    first_log_difference_coeff: sp.Expr | None
    correction_first_gap_power: int | None
    correction_first_gap_coeff: sp.Expr | None
    polynomial_degree1_relation: PolynomialRelation | None
    polynomial_degree2_relation: PolynomialRelation | None
    fractional_linear_relation: FractionalLinearRelation | None
    correction_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    correction_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    normalized_correction_label: str | None = None
    normalized_correction_gap: GapNormalizedSeries | None = None
    normalized_correction_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    normalized_correction_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    normalized_correction_source_family_eta_scans: tuple[SourceFamilyEtaCorrectionScan, ...] = ()
    second_normalized_correction_label: str | None = None
    second_normalized_correction_gap: GapNormalizedSeries | None = None
    second_normalized_correction_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    second_normalized_correction_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    second_normalized_correction_source_family_eta_scans: tuple[SourceFamilyEtaCorrectionScan, ...] = ()
    second_normalized_correction_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...] = ()
    second_normalized_correction_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...] = ()
    second_normalized_correction_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_explicit_transform_eta_scans: tuple[ExplicitSourceFamilyEtaCorrectionScan, ...] = ()


@dataclass(frozen=True)
class GGDescendantPreview:
    """A lightweight preview of extra named GG descendants worth checking next."""

    direct_labels: tuple[str, ...]
    quotient_labels: tuple[str, ...]


@dataclass(frozen=True)
class GGDescendantFocusedScan:
    """A very small direct/quotient odd-prime descendant scan for one target object."""

    direct_labels: tuple[str, ...]
    quotient_labels: tuple[str, ...]
    order_checked: int
    degree_values: tuple[int, ...]
    max_abs_exponent: int
    direct_scans: NamedPrefixBoxScans
    quotient_scans: NamedPrefixBoxScans


@dataclass(frozen=True)
class GGModularEquationScan:
    """A literature-driven GG box built from sign and low-power substitutions."""

    benchmark_name: str
    ordered_basis_series: tuple[tuple[str, str, Series], ...]
    checked_templates: tuple[str, ...]
    hit_templates: tuple[str, ...]
    exact_polynomial_template_labels: tuple[str, ...]
    exact_polynomial_template_hits: tuple[str, ...]
    exact_polynomial_template_obstructions: tuple[tuple[str, int | None, sp.Expr | None], ...]
    polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    quotient_basis_series: tuple[tuple[str, str, Series], ...]
    quotient_exact_polynomial_template_labels: tuple[str, ...]
    quotient_exact_polynomial_template_hits: tuple[str, ...]
    quotient_exact_polynomial_template_obstructions: tuple[tuple[str, int | None, sp.Expr | None], ...]
    quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]
    weighted_coordinate_diagnostics: tuple[GGWeightedCoordinateDiagnostic, ...]
    mixed_quotient_basis_series: tuple[tuple[str, str, Series], ...]
    mixed_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...]
    mixed_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...]
    mixed_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...]
    mixed_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...]


@dataclass(frozen=True)
class MortonPeriodicPointTemplateResult:
    """One Morton-inspired exact periodic-point / algebraic-function template outcome."""

    label: str
    first_failure_power: int | None
    first_failure_coeff: sp.Expr | None
    hit: bool


@dataclass(frozen=True)
class MortonNamedCoordinateScan:
    """A named Morton-source coordinate scan inside the periodic-point lane."""

    family_label: str
    label: str
    expression: str
    template_results: tuple[MortonPeriodicPointTemplateResult, ...]
    leading_normalized_scan: "LeadingNormalizedCoordinateScan | None" = None


@dataclass(frozen=True)
class LeadingNormalizedCoordinateScan:
    """A constant-1 micro-scan obtained by dividing a coordinate by its leading term."""

    label: str
    expression: str
    shift: int
    leading_coefficient: sp.Expr
    normalized_series: tuple[sp.Expr, ...]
    scan: "ConstantOneSeriesScan"


@dataclass(frozen=True)
class WeberLeadingBridgeReferenceLadder:
    """Focused true-GG reference objects for the Weber leading-normalized bridge lane."""

    q_pb_series: tuple[sp.Expr, ...]
    k_pb_series: tuple[sp.Expr, ...]
    q_pk_series: tuple[sp.Expr, ...]
    l_pk_series: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class MortonPeriodicPointScan:
    """A small exact template box inspired by Morton's periodic-point algebraic functions."""

    template_results: tuple[MortonPeriodicPointTemplateResult, ...]
    named_coordinate_scans: tuple[MortonNamedCoordinateScan, ...] = ()
    leading_normalized_bridge_scans: tuple["ConstantOnePairBridgeScan", ...] = ()


@dataclass(frozen=True)
class WeberClassInvariantScan:
    """A literature-backed Weber class-invariant coordinate scan."""

    label: str
    expression: str
    template_label: str
    template_expression: str
    template_hit: bool
    template_first_failure_power: int | None
    template_first_failure_coeff: sp.Expr | None
    correction_label: str
    correction_expression: str
    direct_eta_scans: tuple[EtaQuotientRelationScan, ...]
    direct_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...]
    correction_self_plus_pochhammer_scans: tuple[SelfPlusPochhammerRelationScan, ...]
    correction_self_plus_pochhammer_eta_scans: tuple[SelfPlusPochhammerEtaRelationScan, ...]


@dataclass(frozen=True)
class WeberResidualBridgeScan:
    """A focused bridge between the two Weber class-invariant residuals."""

    primary_label: str
    primary_expression: str
    primary_reason: str
    companion_label: str
    companion_expression: str
    exact_bridge_expression: str
    exact_bridge_holds: bool
    exact_bridge_first_failure_power: int | None
    exact_bridge_first_failure_coeff: sp.Expr | None
    residual_bridge_expression: str
    classical_product_coordinate_label: str
    classical_product_coordinate_expression: str
    classical_product_coordinate_bridge_expression: str
    classical_product_coordinate_scan: "ConstantOneSeriesScan"
    canonical_j_coordinate_label: str
    canonical_j_coordinate_expression: str
    canonical_j_coordinate_reason: str
    canonical_j_coordinate_bridge_expression: str
    canonical_j_coordinate_scan: "ConstantOneSeriesScan"
    quotient_coordinate_label: str
    quotient_coordinate_expression: str
    quotient_coordinate_bridge_expression: str
    quotient_coordinate_bridge_holds: bool
    quotient_coordinate_bridge_first_failure_power: int | None
    quotient_coordinate_bridge_first_failure_coeff: sp.Expr | None
    quotient_coordinate_template_bridge_expression: str
    quotient_coordinate_template_scan: "ConstantOneSeriesScan"
    anchor_canonical_j_coordinate_label: str
    anchor_canonical_j_coordinate_expression: str
    anchor_canonical_j_coordinate_reason: str
    anchor_canonical_j_coordinate_bridge_expression: str
    anchor_canonical_j_coordinate_scan: "ConstantOneSeriesScan"
    primary_named_gg_modular_equation_scan: GGModularEquationScan | None
    quotient_coordinate_template_named_gg_modular_equation_scan: GGModularEquationScan | None
    quotient_label: str
    quotient_expression: str
    quotient_first_failure_power: int | None
    quotient_first_failure_coeff: sp.Expr | None
    quotient_self_polynomial_scan: SelfPolynomialUniquenessScan
    quotient_self_fractional_linear_scan: SelfFractionalLinearUniquenessScan
    quotient_self_quotient_product_scans: tuple[SelfQuotientProductRelationScan, ...]
    quotient_eta_scans: tuple[EtaQuotientRelationScan, ...]
    quotient_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...]
    quotient_self_plus_pochhammer_scans: tuple[SelfPlusPochhammerRelationScan, ...]
    quotient_self_plus_pochhammer_eta_scans: tuple[SelfPlusPochhammerEtaRelationScan, ...]
    canonical_j_anchor_bridge_scan: "ConstantOnePairBridgeScan | None" = None
    canonical_j_lift_bridge_scan: "ConstantOnePairBridgeLiteScan | None" = None
    alternate_anchor_canonical_j_coordinate_label: str | None = None
    alternate_anchor_canonical_j_coordinate_expression: str | None = None
    alternate_anchor_canonical_j_coordinate_reason: str | None = None
    alternate_anchor_canonical_j_coordinate_bridge_expression: str | None = None
    alternate_anchor_canonical_j_coordinate_first_failure_power: int | None = None
    alternate_anchor_canonical_j_coordinate_first_failure_coeff: sp.Expr | None = None
    canonical_j_alt_lift_bridge_scan: "ConstantOnePairBridgeLiteScan | None" = None
    normalized_followup: "NormalizedResidualFollowupScan | None" = None
    followup_bridge_scan: "ConstantOnePairBridgeScan | None" = None


@dataclass(frozen=True)
class NormalizedResidualFollowupScan:
    """A gap-normalized follow-up object after the first residual failure term."""

    label: str
    expression: str
    first_failure_power: int | None
    first_failure_coeff: sp.Expr | None
    self_polynomial_scan: SelfPolynomialUniquenessScan
    self_fractional_linear_scan: SelfFractionalLinearUniquenessScan
    self_quotient_product_scans: tuple[SelfQuotientProductRelationScan, ...]
    eta_scans: tuple[EtaQuotientRelationScan, ...]
    modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...]
    self_plus_pochhammer_scans: tuple[SelfPlusPochhammerRelationScan, ...]
    self_plus_pochhammer_eta_scans: tuple[SelfPlusPochhammerEtaRelationScan, ...]
    named_gg_modular_equation_scan: GGModularEquationScan | None = None
    named_gg_descendant_preview: GGDescendantPreview | None = None
    named_gg_descendant_focused_scan: GGDescendantFocusedScan | None = None


@dataclass(frozen=True)
class ConstantOneSeriesScan:
    """A constant-1 diagnostic object together with its first gap-normalized follow-up."""

    label: str
    expression: str
    first_failure_power: int | None
    first_failure_coeff: sp.Expr | None
    self_polynomial_scan: SelfPolynomialUniquenessScan
    self_fractional_linear_scan: SelfFractionalLinearUniquenessScan
    self_quotient_product_scans: tuple[SelfQuotientProductRelationScan, ...]
    eta_scans: tuple[EtaQuotientRelationScan, ...]
    modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...]
    self_plus_pochhammer_scans: tuple[SelfPlusPochhammerRelationScan, ...]
    self_plus_pochhammer_eta_scans: tuple[SelfPlusPochhammerEtaRelationScan, ...]
    named_gg_modular_equation_scan: GGModularEquationScan | None = None
    normalized_followup: "NormalizedResidualFollowupScan | None" = None


@dataclass(frozen=True)
class PolynomialBridgeRelationScan:
    """A bounded polynomial bridge search between two constant-1 series."""

    degree: int
    relation: PolynomialRelation | None
    error: str | None = None


@dataclass(frozen=True)
class ConstantOnePairBridgeScan:
    """A focused comparison between two normalized constant-1 follow-up objects."""

    left_label: str
    right_label: str
    difference_label: str
    difference_expression: str
    difference_first_failure_power: int | None
    difference_first_failure_coeff: sp.Expr | None
    quotient_label: str
    quotient_expression: str
    quotient_first_failure_power: int | None
    quotient_first_failure_coeff: sp.Expr | None
    quotient_scan: ConstantOneSeriesScan
    polynomial_scans: tuple[PolynomialBridgeRelationScan, ...]
    fractional_linear_relation: FractionalLinearRelation | None
    fractional_linear_error: str | None = None
    quotient_named_coordinate_orbit_scan: NamedCoordinateOrbitScan | None = None
    quotient_followup_named_coordinate_orbit_scan: NamedCoordinateOrbitScan | None = None
    quotient_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_followup_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_followup_bridge_scan: "ConstantOnePairBridgeScan | None" = None


@dataclass(frozen=True)
class ConstantOnePairBridgeLiteScan:
    """A lightweight comparison between two constant-1 series.

    This keeps the first-failure fingerprint together with the smallest bridge boxes
    (bounded polynomial and one-coordinate fractional-linear comparisons), but avoids
    the heavier closure scans on the bridge quotient.
    """

    left_label: str
    right_label: str
    difference_label: str
    difference_expression: str
    difference_first_failure_power: int | None
    difference_first_failure_coeff: sp.Expr | None
    quotient_label: str
    quotient_expression: str
    quotient_first_failure_power: int | None
    quotient_first_failure_coeff: sp.Expr | None
    polynomial_scans: tuple[PolynomialBridgeRelationScan, ...]
    fractional_linear_relation: FractionalLinearRelation | None
    fractional_linear_error: str | None = None
    quotient_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_eta_scans: tuple[EtaQuotientRelationScan, ...] | None = None
    quotient_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] | None = None
    quotient_followup_label: str | None = None
    quotient_followup_expression: str | None = None
    quotient_followup_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_followup_eta_scans: tuple[EtaQuotientRelationScan, ...] | None = None
    quotient_followup_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] | None = None


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


def search_polynomial_bridge_relation(
    *,
    series_by_variable: dict[str, Series],
    order: int,
    max_total_degree: int,
    required_variables: tuple[str, ...],
) -> PolynomialRelation | None:
    """Find a bounded polynomial relation that genuinely uses each required variable."""
    if order < 2:
        raise ValueError("order must be at least 2")
    if max_total_degree < 1:
        raise ValueError("max_total_degree must be at least 1")
    if len(series_by_variable) < 2:
        raise ValueError("need at least two variables for a relation search")
    if any(len(series) < order for series in series_by_variable.values()):
        raise ValueError("series are shorter than requested order")
    if not required_variables:
        raise ValueError("required_variables must be non-empty")
    if any(required not in series_by_variable for required in required_variables):
        raise ValueError("required_variables must be series variable names")

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
        required_variables=required_variables,
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


def build_reduced_tail_anchor(
    *,
    reduced_coeffs: ContinuedFractionCoeffs,
    symbol: sp.Symbol,
    start_stage: int = 3,
    order: int,
) -> ReducedTailAnchor | None:
    """Build the anchored tail object T(state_expr) defined by T_n = b_n + a_n / T_{n+1}."""
    if start_stage < 1:
        raise ValueError("start_stage must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(reduced_coeffs.a_terms) != len(reduced_coeffs.b_terms):
        raise ValueError("reduced coefficient lists must have the same length")
    if len(reduced_coeffs.a_terms) <= start_stage:
        return None

    tail_series = _expr_to_series(reduced_coeffs.b_terms[-1], symbol=symbol, order=order)
    for stage in range(len(reduced_coeffs.a_terms) - 2, start_stage - 1, -1):
        denominator_series = tail_series
        numerator_series = _expr_to_series(reduced_coeffs.a_terms[stage], symbol=symbol, order=order)
        tail_series = series_add(
            _expr_to_series(reduced_coeffs.b_terms[stage], symbol=symbol, order=order),
            series_mul(numerator_series, series_invert(denominator_series)),
        )

    state_expr = sp.expand(reduced_coeffs.b_terms[start_stage] - 1)
    normalization_series = _expr_to_series(1 + state_expr, symbol=symbol, order=order)
    normalized_series = series_div(tail_series, normalization_series)
    return ReducedTailAnchor(
        start_stage=start_stage,
        state_expr=state_expr,
        tail_series=tuple(tail_series),
        normalized_series=tuple(normalized_series),
    )


def build_gap_normalized_series(
    *,
    target_series: Series,
) -> GapNormalizedSeries | None:
    """Normalize a constant-1 series by its first nonzero term after the constant coefficient."""
    if len(target_series) < 2:
        raise ValueError("target_series must have length at least 2")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for gap normalization")

    residual = [sp.simplify(value - (1 if index == 0 else 0)) for index, value in enumerate(target_series)]
    shift: int | None = None
    leading_coefficient: sp.Expr | None = None
    for index in range(1, len(residual)):
        if sp.simplify(residual[index]) != 0:
            shift = index
            leading_coefficient = sp.simplify(residual[index])
            break
    if shift is None or leading_coefficient is None:
        return None

    normalized: Series = [sp.Integer(0) for _ in range(len(target_series))]
    for index in range(len(target_series) - shift):
        normalized[index] = sp.simplify(residual[index + shift] / leading_coefficient)
    return GapNormalizedSeries(
        shift=shift,
        leading_coefficient=leading_coefficient,
        normalized_series=tuple(normalized),
    )


def _leading_term_scale_expr(*, shift: int, leading_coefficient: sp.Expr) -> str:
    coeff = sp.simplify(leading_coefficient)
    if shift == 0:
        return _format_expr(coeff)
    power_term = "t" if shift == 1 else f"t^{shift}"
    if coeff == 1:
        return power_term
    if coeff == -1:
        return f"-{power_term}"
    return f"{_format_expr(coeff)}*{power_term}"


def _build_leading_normalized_coordinate_scan(
    *,
    source_label: str,
    target_series: Series,
    order: int,
    max_abs_exponent: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
) -> LeadingNormalizedCoordinateScan | None:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    shift, leading_coefficient = _first_nonzero_residual_term(target_series[:order])
    if shift is None or leading_coefficient is None:
        return None

    normalized: Series = [sp.Integer(0) for _ in range(order)]
    for index in range(order - shift):
        normalized[index] = sp.simplify(target_series[index + shift] / leading_coefficient)

    label = f"N_{source_label}"
    scale_expr = _leading_term_scale_expr(
        shift=shift,
        leading_coefficient=leading_coefficient,
    )
    expression = f"{label} = {source_label} / ({scale_expr})"
    scan = _scan_constant_one_series(
        label=label,
        expression=expression,
        target_series=normalized,
        order=order,
        eta_levels=eta_levels,
        moduli=moduli,
        max_abs_exponent=max_abs_exponent,
    )
    return LeadingNormalizedCoordinateScan(
        label=label,
        expression=expression,
        shift=shift,
        leading_coefficient=leading_coefficient,
        normalized_series=tuple(normalized),
        scan=scan,
    )


def _constant_one_residual_series(
    *,
    target_series: Series,
    order: int,
) -> tuple[Series, int | None, sp.Expr | None]:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1")

    residual = [
        sp.simplify(value - (sp.Integer(1) if index == 0 else sp.Integer(0)))
        for index, value in enumerate(target_series[:order])
    ]
    first_failure_power, first_failure_coeff = _first_nonzero_residual_term(residual)
    return residual, first_failure_power, first_failure_coeff


def _normalized_constant_one_followup_series(
    *,
    target_series: Series,
    order: int,
) -> Series | None:
    residual, first_failure_power, first_failure_coeff = _constant_one_residual_series(
        target_series=target_series,
        order=order,
    )
    if first_failure_power is None or first_failure_coeff is None:
        return None

    followup: Series = [sp.Integer(0) for _ in range(order)]
    for index in range(order - first_failure_power):
        followup[index] = sp.simplify(
            residual[index + first_failure_power] / first_failure_coeff
        )
    return followup


def _named_coordinate_power_entries(
    *,
    base_label: str,
    base_expression: str,
    base_series: Series,
    order: int,
    powers: tuple[int, ...],
) -> tuple[tuple[str, str, Series], ...]:
    if len(base_series) < order:
        raise ValueError("base_series is shorter than requested order")

    entries: list[tuple[str, str, Series]] = [
        (base_label, base_expression, list(base_series[:order]))
    ]
    for power in tuple(sorted({value for value in powers if value >= 2})):
        entries.append(
            (
                f"{base_label}_{power}",
                f"{base_label}(t^{power})",
                benchmark_power_substitution_series(base_series, power=power, order=order),
            )
        )
    return tuple(entries)


def scan_named_coordinate_orbit_box(
    *,
    target_series: Series,
    family_label: str,
    base_label: str,
    base_expression: str,
    base_series: Series,
    order: int,
    powers: tuple[int, ...] = (2, 3, 4),
    degree_values: tuple[int, ...] = (1, 2),
    max_abs_exponent: int = 2,
    solve_order: int | None = None,
) -> NamedCoordinateOrbitScan:
    direct_entries = _named_coordinate_power_entries(
        base_label=base_label,
        base_expression=base_expression,
        base_series=base_series,
        order=order,
        powers=powers,
    )
    direct_ordered_basis = tuple((label, series) for label, _, series in direct_entries)
    direct_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=direct_ordered_basis,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
        include_two_layer=False,
    )

    base_series_trunc = direct_entries[0][2]
    quotient_entries: list[tuple[str, str, Series]] = []
    for label, _, series in direct_entries[1:]:
        suffix = label.removeprefix(f"{base_label}_")
        quotient_entries.append(
            (
                f"{base_label}_q{suffix}",
                f"{label} / {base_label}",
                series_div(series, base_series_trunc),
            )
        )
    quotient_ordered_basis = tuple((label, series) for label, _, series in quotient_entries)
    quotient_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=quotient_ordered_basis,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
        include_two_layer=False,
    )

    mixed_entries = ((direct_entries[0][0], direct_entries[0][1], direct_entries[0][2]),) + tuple(
        quotient_entries
    )
    mixed_ordered_basis = tuple((label, series) for label, _, series in mixed_entries)
    mixed_scans = scan_named_prefix_boxes(
        target_series=target_series,
        ordered_basis_series=mixed_ordered_basis,
        order=order,
        degree_values=degree_values,
        required_variable="F",
        max_abs_exponent=max_abs_exponent,
        solve_order=solve_order,
        include_two_layer=False,
    )

    return NamedCoordinateOrbitScan(
        family_label=family_label,
        base_label=base_label,
        base_expression=base_expression,
        ordered_basis_series=direct_entries,
        direct_scans=direct_scans,
        quotient_basis_series=tuple(quotient_entries),
        quotient_scans=quotient_scans,
        mixed_quotient_basis_series=mixed_entries,
        mixed_quotient_scans=mixed_scans,
    )


def _scan_constant_one_series(
    *,
    label: str,
    expression: str,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...],
    moduli: tuple[int, ...],
    max_abs_exponent: int,
    followup_label: str | None = None,
    named_gg_benchmark_name: str | None = None,
    named_gg_series: Series | None = None,
    named_gg_degree_values: tuple[int, ...] = (1, 2),
    named_gg_solve_order: int | None = None,
    named_gg_supplemental_powers: tuple[int, ...] = (),
    named_gg_include_weighted_coordinate_diagnostics: bool = False,
    followup_named_gg_benchmark_name: str | None = None,
    followup_named_gg_series: Series | None = None,
    followup_named_gg_degree_values: tuple[int, ...] = (1, 2),
    followup_named_gg_solve_order: int | None = None,
    followup_named_gg_supplemental_powers: tuple[int, ...] = (),
    followup_named_gg_descendant_preview_powers: tuple[int, ...] = (),
    followup_named_gg_descendant_focus_powers: tuple[int, ...] = (),
    followup_named_gg_descendant_focus_degree_values: tuple[int, ...] = (1,),
    followup_named_gg_descendant_focus_max_abs_exponent: int = 2,
    followup_named_gg_include_weighted_coordinate_diagnostics: bool = False,
) -> ConstantOneSeriesScan:
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1")

    residual = [sp.simplify(value - (1 if index == 0 else 0)) for index, value in enumerate(target_series[:order])]
    first_failure_power, first_failure_coeff = _first_nonzero_residual_term(residual)
    focused_named_gg_order = min(order, 12 if max_abs_exponent <= 4 else 16)
    focused_named_gg_solve_order = min(focused_named_gg_order, 12)

    named_gg_modular_equation_scan: GGModularEquationScan | None = None
    if (
        named_gg_benchmark_name is not None
        and named_gg_series is not None
        and max_abs_exponent > 4
    ):
        named_gg_modular_equation_scan = scan_gg_modular_equation_box(
            target_series=target_series[:focused_named_gg_order],
            benchmark_name=named_gg_benchmark_name,
            gg_series=named_gg_series[:focused_named_gg_order],
            order=focused_named_gg_order,
            degree_values=named_gg_degree_values,
            max_abs_exponent=max_abs_exponent,
            solve_order=named_gg_solve_order or focused_named_gg_solve_order,
            supplemental_powers=named_gg_supplemental_powers,
            include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
        )

    normalized_followup: NormalizedResidualFollowupScan | None = None
    if (
        followup_label is not None
        and first_failure_power is not None
        and first_failure_coeff is not None
    ):
        followup_series: Series = [sp.Integer(0) for _ in range(order)]
        for index in range(order):
            source_index = index + first_failure_power
            if source_index >= order:
                break
            followup_series[index] = sp.simplify(
                residual[source_index] / first_failure_coeff
            )
        followup_residual = [sp.simplify(value) for value in followup_series]
        followup_residual[0] = sp.simplify(followup_residual[0] - 1)
        followup_first_failure_power, followup_first_failure_coeff = _first_nonzero_residual_term(
            followup_residual
        )
        followup_named_gg_modular_equation_scan: GGModularEquationScan | None = None
        followup_named_gg_descendant_preview: GGDescendantPreview | None = None
        followup_named_gg_descendant_focused_scan: GGDescendantFocusedScan | None = None
        if (
            followup_named_gg_benchmark_name is not None
            and followup_named_gg_series is not None
            and max_abs_exponent > 4
        ):
            followup_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
                target_series=followup_series[:focused_named_gg_order],
                benchmark_name=followup_named_gg_benchmark_name,
                gg_series=followup_named_gg_series[:focused_named_gg_order],
                order=focused_named_gg_order,
                degree_values=followup_named_gg_degree_values,
                max_abs_exponent=max_abs_exponent,
                solve_order=followup_named_gg_solve_order or focused_named_gg_solve_order,
                supplemental_powers=followup_named_gg_supplemental_powers,
                include_weighted_coordinate_diagnostics=followup_named_gg_include_weighted_coordinate_diagnostics,
            )
        if followup_named_gg_series is not None and max_abs_exponent > 4:
            followup_named_gg_descendant_preview = _gg_descendant_preview(
                base_series=followup_named_gg_series[:focused_named_gg_order],
                order=focused_named_gg_order,
                supplemental_powers=followup_named_gg_descendant_preview_powers,
            )
            followup_named_gg_descendant_focused_scan = _scan_gg_descendant_focus_box(
                target_series=followup_series,
                gg_series=followup_named_gg_series,
                order=order,
                supplemental_powers=followup_named_gg_descendant_focus_powers,
                degree_values=followup_named_gg_descendant_focus_degree_values,
                max_abs_exponent=followup_named_gg_descendant_focus_max_abs_exponent,
            )
        normalized_followup = NormalizedResidualFollowupScan(
            label=followup_label,
            expression=(
                f"{followup_label} = ({label} - 1) / "
                f"({_format_expr(first_failure_coeff)}*t^{first_failure_power})"
            ),
            first_failure_power=followup_first_failure_power,
            first_failure_coeff=followup_first_failure_coeff,
            self_polynomial_scan=scan_self_polynomial_uniqueness_relations(
                target_series=followup_series,
                moduli=moduli,
                order=order,
                fg_degree_values=(1, 2),
                t_degree_values=(1, 2, 3),
            ),
            self_fractional_linear_scan=scan_self_fractional_linear_uniqueness_relations(
                target_series=followup_series,
                moduli=moduli,
                order=order,
                t_degree_values=(1, 2, 3),
            ),
            self_quotient_product_scans=tuple(
                scan_ratio_self_quotient_product_relations(
                    ratio_series=followup_series,
                    moduli=moduli,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            eta_scans=tuple(
                scan_ratio_eta_quotient_relations(
                    ratio_series=followup_series,
                    levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            modular_unit_eta_scans=tuple(
                scan_ratio_modular_unit_eta_relations(
                    ratio_series=followup_series,
                    moduli=moduli,
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            self_plus_pochhammer_scans=tuple(
                scan_ratio_self_plus_pochhammer_relations(
                    ratio_series=followup_series,
                    moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            self_plus_pochhammer_eta_scans=tuple(
                scan_ratio_self_plus_pochhammer_eta_relations(
                    ratio_series=followup_series,
                    moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            named_gg_modular_equation_scan=followup_named_gg_modular_equation_scan,
            named_gg_descendant_preview=followup_named_gg_descendant_preview,
            named_gg_descendant_focused_scan=followup_named_gg_descendant_focused_scan,
        )

    return ConstantOneSeriesScan(
        label=label,
        expression=expression,
        first_failure_power=first_failure_power,
        first_failure_coeff=first_failure_coeff,
        self_polynomial_scan=scan_self_polynomial_uniqueness_relations(
            target_series=target_series[:order],
            moduli=moduli,
            order=order,
            fg_degree_values=(1, 2),
            t_degree_values=(1, 2, 3),
        ),
        self_fractional_linear_scan=scan_self_fractional_linear_uniqueness_relations(
            target_series=target_series[:order],
            moduli=moduli,
            order=order,
            t_degree_values=(1, 2, 3),
        ),
        self_quotient_product_scans=tuple(
            scan_ratio_self_quotient_product_relations(
                ratio_series=target_series[:order],
                moduli=moduli,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        eta_scans=tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=target_series[:order],
                levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        modular_unit_eta_scans=tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=target_series[:order],
                moduli=moduli,
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        self_plus_pochhammer_scans=tuple(
            scan_ratio_self_plus_pochhammer_relations(
                ratio_series=target_series[:order],
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        self_plus_pochhammer_eta_scans=tuple(
            scan_ratio_self_plus_pochhammer_eta_relations(
                ratio_series=target_series[:order],
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        named_gg_modular_equation_scan=named_gg_modular_equation_scan,
        normalized_followup=normalized_followup,
    )


def _scan_constant_one_series_pair_bridge(
    *,
    left_label: str,
    left_series: Series,
    right_label: str,
    right_series: Series,
    order: int,
    difference_label: str = "D_XR_ws",
    quotient_label: str = "Q_XR_ws",
    quotient_followup_label: str = "K_XR_ws",
    polynomial_degree_values: tuple[int, ...] = (1, 2, 3),
    solve_order: int | None = 24,
    quotient_followup_bridge_left_label: str | None = None,
    quotient_followup_bridge_difference_label: str | None = None,
    quotient_followup_bridge_quotient_label: str | None = None,
    quotient_followup_bridge_quotient_followup_label: str | None = None,
    named_gg_benchmark_name: str | None = None,
    named_gg_series: Series | None = None,
    named_gg_degree_values: tuple[int, ...] = (1, 2),
    named_gg_max_abs_exponent: int = 8,
    named_gg_solve_order: int | None = None,
    named_gg_supplemental_powers: tuple[int, ...] = (),
    named_gg_target_pairs: tuple[tuple[str, str], ...] = (),
    quotient_followup_named_gg_series: Series | None = None,
    quotient_followup_named_gg_descendant_preview_powers: tuple[int, ...] = (),
    quotient_followup_named_gg_descendant_focus_powers: tuple[int, ...] = (),
    quotient_followup_named_gg_descendant_focus_degree_values: tuple[int, ...] = (1,),
    quotient_followup_named_gg_descendant_focus_max_abs_exponent: int = 2,
    nested_quotient_followup_named_gg_series: Series | None = None,
    nested_quotient_followup_named_gg_descendant_preview_powers: tuple[int, ...] = (),
    nested_quotient_followup_named_gg_descendant_focus_powers: tuple[int, ...] = (),
    nested_quotient_followup_named_gg_descendant_focus_degree_values: tuple[int, ...] = (1,),
    nested_quotient_followup_named_gg_descendant_focus_max_abs_exponent: int = 2,
    named_gg_include_weighted_coordinate_diagnostics: bool = False,
) -> ConstantOnePairBridgeScan:
    if len(left_series) < order or len(right_series) < order:
        raise ValueError("series are shorter than requested order")
    if sp.simplify(left_series[0] - 1) != 0 or sp.simplify(right_series[0] - 1) != 0:
        raise ValueError("pair-bridge scan requires constant-1 series")

    difference_series = [
        sp.simplify(right_series[index] - left_series[index])
        for index in range(order)
    ]
    (
        difference_first_failure_power,
        difference_first_failure_coeff,
    ) = _first_nonzero_residual_term(difference_series)
    quotient_series = series_div(right_series, left_series)
    quotient_residual = [sp.simplify(value) for value in quotient_series]
    quotient_residual[0] = sp.simplify(quotient_residual[0] - 1)
    quotient_first_failure_power, quotient_first_failure_coeff = _first_nonzero_residual_term(
        quotient_residual
    )
    quotient_scan = _scan_constant_one_series(
        label=quotient_label,
        expression=f"{quotient_label} = {right_label} / {left_label}",
        target_series=quotient_series,
        order=order,
        eta_levels=(1, 2, 4),
        moduli=(2, 3, 4),
        max_abs_exponent=8,
        followup_label=quotient_followup_label,
        followup_named_gg_series=quotient_followup_named_gg_series,
        followup_named_gg_descendant_preview_powers=quotient_followup_named_gg_descendant_preview_powers,
        followup_named_gg_descendant_focus_powers=quotient_followup_named_gg_descendant_focus_powers,
        followup_named_gg_descendant_focus_degree_values=quotient_followup_named_gg_descendant_focus_degree_values,
        followup_named_gg_descendant_focus_max_abs_exponent=quotient_followup_named_gg_descendant_focus_max_abs_exponent,
    )

    checked_order = min(order, solve_order or order)
    variable_series: Series = [sp.Integer(0) for _ in range(checked_order)]
    if checked_order > 1:
        variable_series[1] = sp.Integer(1)

    polynomial_scans: list[PolynomialBridgeRelationScan] = []
    for degree in polynomial_degree_values:
        try:
            relation = search_polynomial_bridge_relation(
                series_by_variable={
                    "HX": left_series[:checked_order],
                    "HR": right_series[:checked_order],
                    "T": variable_series,
                },
                order=checked_order,
                max_total_degree=degree,
                required_variables=("HX", "HR"),
            )
            polynomial_scans.append(
                PolynomialBridgeRelationScan(
                    degree=degree,
                    relation=relation,
                )
            )
        except ValueError as exc:
            polynomial_scans.append(
                PolynomialBridgeRelationScan(
                    degree=degree,
                    relation=None,
                    error=str(exc),
                )
            )

    fractional_linear_relation = None
    fractional_linear_error = None
    try:
        fractional_linear_relation = search_fractional_linear_relation(
            target_series=right_series[:checked_order],
            basis_series_by_variable={"HX": left_series[:checked_order]},
            order=checked_order,
        )
    except ValueError as exc:
        fractional_linear_error = str(exc)

    quotient_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_followup_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    named_gg_matches_target = (
        named_gg_benchmark_name is not None
        and named_gg_series is not None
        and named_gg_max_abs_exponent > 4
        and any(
            quotient_label == target_quotient_label
            and quotient_followup_label == target_followup_label
            for target_quotient_label, target_followup_label in named_gg_target_pairs
        )
    )
    if named_gg_matches_target:
        quotient_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
            target_series=quotient_series,
            benchmark_name=named_gg_benchmark_name,
            gg_series=named_gg_series,
            order=order,
            degree_values=named_gg_degree_values,
            max_abs_exponent=named_gg_max_abs_exponent,
            solve_order=named_gg_solve_order,
            supplemental_powers=named_gg_supplemental_powers,
            include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
        )
        if quotient_scan.normalized_followup is not None:
            quotient_followup_series: Series = [sp.Integer(0) for _ in range(order)]
            quotient_followup_series[0] = sp.Integer(1)
            for index in range(order):
                source_index = index + quotient_first_failure_power
                if source_index >= order:
                    break
                quotient_followup_series[index] = sp.simplify(
                    quotient_residual[source_index] / quotient_first_failure_coeff
                )
            quotient_followup_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
                target_series=quotient_followup_series,
                benchmark_name=named_gg_benchmark_name,
                gg_series=named_gg_series,
                order=order,
                degree_values=named_gg_degree_values,
                max_abs_exponent=named_gg_max_abs_exponent,
                solve_order=named_gg_solve_order,
                supplemental_powers=named_gg_supplemental_powers,
                include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
            )

    quotient_followup_bridge_scan: ConstantOnePairBridgeScan | None = None
    if (
        quotient_followup_bridge_left_label is not None
        and quotient_followup_bridge_difference_label is not None
        and quotient_followup_bridge_quotient_label is not None
        and quotient_followup_bridge_quotient_followup_label is not None
        and quotient_scan.normalized_followup is not None
        and quotient_first_failure_power is not None
        and quotient_first_failure_coeff is not None
    ):
        quotient_followup_series: Series = [sp.Integer(0) for _ in range(order)]
        quotient_followup_series[0] = sp.Integer(1)
        for index in range(order):
            source_index = index + quotient_first_failure_power
            if source_index >= order:
                break
            quotient_followup_series[index] = sp.simplify(
                quotient_residual[source_index] / quotient_first_failure_coeff
            )
        quotient_followup_bridge_scan = _scan_constant_one_series_pair_bridge(
            left_label=quotient_followup_bridge_left_label,
            left_series=left_series,
            right_label=quotient_scan.normalized_followup.label,
            right_series=quotient_followup_series,
            order=order,
            difference_label=quotient_followup_bridge_difference_label,
            quotient_label=quotient_followup_bridge_quotient_label,
            quotient_followup_label=quotient_followup_bridge_quotient_followup_label,
            named_gg_benchmark_name=named_gg_benchmark_name,
            named_gg_series=named_gg_series,
            named_gg_degree_values=named_gg_degree_values,
            named_gg_max_abs_exponent=named_gg_max_abs_exponent,
            named_gg_solve_order=named_gg_solve_order,
            named_gg_supplemental_powers=named_gg_supplemental_powers,
            named_gg_target_pairs=named_gg_target_pairs,
            named_gg_include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
            quotient_followup_named_gg_series=nested_quotient_followup_named_gg_series,
            quotient_followup_named_gg_descendant_preview_powers=nested_quotient_followup_named_gg_descendant_preview_powers,
            quotient_followup_named_gg_descendant_focus_powers=nested_quotient_followup_named_gg_descendant_focus_powers,
            quotient_followup_named_gg_descendant_focus_degree_values=nested_quotient_followup_named_gg_descendant_focus_degree_values,
            quotient_followup_named_gg_descendant_focus_max_abs_exponent=nested_quotient_followup_named_gg_descendant_focus_max_abs_exponent,
        )

    return ConstantOnePairBridgeScan(
        left_label=left_label,
        right_label=right_label,
        difference_label=difference_label,
        difference_expression=f"{difference_label} = {right_label} - {left_label}",
        difference_first_failure_power=difference_first_failure_power,
        difference_first_failure_coeff=difference_first_failure_coeff,
        quotient_label=quotient_label,
        quotient_expression=f"{quotient_label} = {right_label} / {left_label}",
        quotient_first_failure_power=quotient_first_failure_power,
        quotient_first_failure_coeff=quotient_first_failure_coeff,
        quotient_scan=quotient_scan,
        polynomial_scans=tuple(polynomial_scans),
        fractional_linear_relation=fractional_linear_relation,
        fractional_linear_error=fractional_linear_error,
        quotient_named_gg_modular_equation_scan=quotient_named_gg_modular_equation_scan,
        quotient_followup_named_gg_modular_equation_scan=quotient_followup_named_gg_modular_equation_scan,
        quotient_followup_bridge_scan=quotient_followup_bridge_scan,
    )


def _scan_constant_one_series_pair_bridge_lite(
    *,
    left_label: str,
    left_series: Series,
    right_label: str,
    right_series: Series,
    order: int,
    difference_label: str = "D_XR_ws",
    quotient_label: str = "Q_XR_ws",
    polynomial_degree_values: tuple[int, ...] = (1, 2, 3),
    solve_order: int | None = 24,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    named_gg_benchmark_name: str | None = None,
    named_gg_series: Series | None = None,
    named_gg_order: int | None = None,
    named_gg_degree_values: tuple[int, ...] = (1, 2),
    named_gg_max_abs_exponent: int = 8,
    named_gg_solve_order: int | None = None,
    named_gg_supplemental_powers: tuple[int, ...] = (),
    named_gg_include_weighted_coordinate_diagnostics: bool = False,
    quotient_followup_label: str | None = None,
) -> ConstantOnePairBridgeLiteScan:
    if len(left_series) < order or len(right_series) < order:
        raise ValueError("series are shorter than requested order")
    if sp.simplify(left_series[0] - 1) != 0 or sp.simplify(right_series[0] - 1) != 0:
        raise ValueError("pair-bridge scan requires constant-1 series")

    difference_series = [
        sp.simplify(right_series[index] - left_series[index])
        for index in range(order)
    ]
    difference_first_failure_power, difference_first_failure_coeff = _first_nonzero_residual_term(
        difference_series
    )
    quotient_series = series_div(right_series, left_series)
    quotient_residual = [sp.simplify(value) for value in quotient_series]
    quotient_residual[0] = sp.simplify(quotient_residual[0] - 1)
    quotient_first_failure_power, quotient_first_failure_coeff = _first_nonzero_residual_term(
        quotient_residual
    )

    quotient_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_eta_scans: tuple[EtaQuotientRelationScan, ...] | None = None
    quotient_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] | None = None
    quotient_followup_expression: str | None = None
    quotient_followup_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    quotient_followup_eta_scans: tuple[EtaQuotientRelationScan, ...] | None = None
    quotient_followup_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] | None = None
    if (
        named_gg_benchmark_name is not None
        and named_gg_series is not None
        and named_gg_max_abs_exponent > 4
        and quotient_first_failure_power is not None
        and quotient_first_failure_coeff is not None
    ):
        gg_scan_order = order if named_gg_order is None else max(2, min(named_gg_order, order))
        quotient_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
            target_series=quotient_series[:gg_scan_order],
            benchmark_name=named_gg_benchmark_name,
            gg_series=named_gg_series[:gg_scan_order],
            order=gg_scan_order,
            degree_values=named_gg_degree_values,
            max_abs_exponent=named_gg_max_abs_exponent,
            solve_order=named_gg_solve_order,
            supplemental_powers=named_gg_supplemental_powers,
            include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
        )
        quotient_eta_scans = tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=quotient_series[:gg_scan_order],
                levels=eta_levels,
                order=gg_scan_order,
                max_abs_exponent=min(named_gg_max_abs_exponent, 4),
            )
        )
        quotient_modular_unit_eta_scans = tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=quotient_series[:gg_scan_order],
                moduli=moduli,
                eta_levels=eta_levels,
                order=gg_scan_order,
                max_abs_exponent=min(named_gg_max_abs_exponent, 4),
            )
        )

        if quotient_followup_label is not None:
            quotient_followup_series: Series = [sp.Integer(0) for _ in range(order)]
            for index in range(order):
                source_index = index + quotient_first_failure_power
                if source_index >= order:
                    break
                quotient_followup_series[index] = sp.simplify(
                    quotient_residual[source_index] / quotient_first_failure_coeff
                )
            quotient_followup_expression = (
                f"{quotient_followup_label} = ({quotient_label} - 1) / "
                f"({_format_expr(quotient_first_failure_coeff)}*t^{quotient_first_failure_power})"
            )
            quotient_followup_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
                target_series=quotient_followup_series[:gg_scan_order],
                benchmark_name=named_gg_benchmark_name,
                gg_series=named_gg_series[:gg_scan_order],
                order=gg_scan_order,
                degree_values=named_gg_degree_values,
                max_abs_exponent=named_gg_max_abs_exponent,
                solve_order=named_gg_solve_order,
                supplemental_powers=named_gg_supplemental_powers,
                include_weighted_coordinate_diagnostics=named_gg_include_weighted_coordinate_diagnostics,
            )
            quotient_followup_eta_scans = tuple(
                scan_ratio_eta_quotient_relations(
                    ratio_series=quotient_followup_series[:gg_scan_order],
                    levels=eta_levels,
                    order=gg_scan_order,
                    max_abs_exponent=min(named_gg_max_abs_exponent, 4),
                )
            )
            quotient_followup_modular_unit_eta_scans = tuple(
                scan_ratio_modular_unit_eta_relations(
                    ratio_series=quotient_followup_series[:gg_scan_order],
                    moduli=moduli,
                    eta_levels=eta_levels,
                    order=gg_scan_order,
                    max_abs_exponent=min(named_gg_max_abs_exponent, 4),
                )
            )

    checked_order = min(order, solve_order or order)
    variable_series: Series = [sp.Integer(0) for _ in range(checked_order)]
    if checked_order > 1:
        variable_series[1] = sp.Integer(1)

    polynomial_scans: list[PolynomialBridgeRelationScan] = []
    for degree in polynomial_degree_values:
        try:
            relation = search_polynomial_bridge_relation(
                series_by_variable={
                    "HX": left_series[:checked_order],
                    "HR": right_series[:checked_order],
                    "T": variable_series,
                },
                order=checked_order,
                max_total_degree=degree,
                required_variables=("HX", "HR"),
            )
            polynomial_scans.append(
                PolynomialBridgeRelationScan(
                    degree=degree,
                    relation=relation,
                )
            )
        except ValueError as exc:
            polynomial_scans.append(
                PolynomialBridgeRelationScan(
                    degree=degree,
                    relation=None,
                    error=str(exc),
                )
            )

    fractional_linear_relation = None
    fractional_linear_error = None
    try:
        fractional_linear_relation = search_fractional_linear_relation(
            target_series=right_series[:checked_order],
            basis_series_by_variable={"HX": left_series[:checked_order]},
            order=checked_order,
        )
    except ValueError as exc:
        fractional_linear_error = str(exc)

    return ConstantOnePairBridgeLiteScan(
        left_label=left_label,
        right_label=right_label,
        difference_label=difference_label,
        difference_expression=f"{difference_label} = {right_label} - {left_label}",
        difference_first_failure_power=difference_first_failure_power,
        difference_first_failure_coeff=difference_first_failure_coeff,
        quotient_label=quotient_label,
        quotient_expression=f"{quotient_label} = {right_label} / {left_label}",
        quotient_first_failure_power=quotient_first_failure_power,
        quotient_first_failure_coeff=quotient_first_failure_coeff,
        polynomial_scans=tuple(polynomial_scans),
        fractional_linear_relation=fractional_linear_relation,
        fractional_linear_error=fractional_linear_error,
        quotient_named_gg_modular_equation_scan=quotient_named_gg_modular_equation_scan,
        quotient_eta_scans=quotient_eta_scans,
        quotient_modular_unit_eta_scans=quotient_modular_unit_eta_scans,
        quotient_followup_label=quotient_followup_label if quotient_followup_expression is not None else None,
        quotient_followup_expression=quotient_followup_expression,
        quotient_followup_named_gg_modular_equation_scan=quotient_followup_named_gg_modular_equation_scan,
        quotient_followup_eta_scans=quotient_followup_eta_scans,
        quotient_followup_modular_unit_eta_scans=quotient_followup_modular_unit_eta_scans,
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
def _plus_residue_pochhammer_series_tuple(*, residue: int, modulus: int, order: int) -> tuple[sp.Expr, ...]:
    if residue < 1:
        raise ValueError("residue must be positive")
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if residue >= modulus:
        raise ValueError("residue must be strictly smaller than modulus")
    if order < 2:
        raise ValueError("order must be at least 2")

    series: Series = [sp.Integer(0) for _ in range(order)]
    series[0] = sp.Integer(1)
    for power in range(residue, order, modulus):
        series = series_mul(series, list(_one_plus_power_series_tuple(power=power, order=order)))
    return tuple(series)


def _plus_residue_pochhammer_series(*, residue: int, modulus: int, order: int) -> Series:
    return list(_plus_residue_pochhammer_series_tuple(residue=residue, modulus=modulus, order=order))


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


@lru_cache(maxsize=None)
def _modular_unit_eta_basis_series_tuple(
    *,
    modulus: int,
    level: int,
    order: int,
) -> tuple[tuple[str, tuple[sp.Expr, ...]], ...]:
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    basis_entries: list[tuple[str, tuple[sp.Expr, ...]]] = []
    for residue in range(1, modulus):
        basis_entries.append((f"M{residue}", _one_minus_power_series_tuple(power=residue, order=order)))
    for residue in range(1, modulus):
        basis_entries.append((f"P{residue}", _one_plus_power_series_tuple(power=residue, order=order)))
    basis_entries.extend(_eta_quotient_basis_series_tuple(level=level, order=order))
    return tuple(basis_entries)


def _modular_unit_eta_basis_series(*, modulus: int, level: int, order: int) -> dict[str, Series]:
    return {
        name: list(series)
        for name, series in _modular_unit_eta_basis_series_tuple(
            modulus=modulus,
            level=level,
            order=order,
        )
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


def search_self_plus_pochhammer_relation(
    *,
    target_series: Series,
    modulus: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = F(t^m)^a * prod_r (-t^r; t^m)_inf^e_r with bounded exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a self plus-Pochhammer search")

    g_power_label = f"G{modulus}"
    basis_series_by_variable = {
        g_power_label: benchmark_power_substitution_series(target_series, power=modulus, order=order),
        **{
            f"PP{residue}": _plus_residue_pochhammer_series(residue=residue, modulus=modulus, order=order)
            for residue in range(1, modulus)
        },
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


def scan_ratio_self_plus_pochhammer_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SelfPlusPochhammerRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    if not unique_moduli:
        return []

    scans: list[SelfPlusPochhammerRelationScan] = []
    for modulus in unique_moduli:
        try:
            relation = search_self_plus_pochhammer_relation(
                target_series=ratio_series,
                modulus=modulus,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
            scans.append(SelfPlusPochhammerRelationScan(modulus=modulus, relation=relation))
        except ValueError as exc:
            scans.append(
                SelfPlusPochhammerRelationScan(
                    modulus=modulus,
                    relation=None,
                    error=str(exc),
                )
            )
    return scans


def search_self_plus_pochhammer_eta_relation(
    *,
    target_series: Series,
    modulus: int,
    level: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = F(t^m)^a * prod_r (-t^r; t^m)_inf^u_r * eta_tail with bounded exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if level < 1:
        raise ValueError("level must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a self plus-Pochhammer eta search")

    g_power_label = f"G{modulus}"
    basis_series_by_variable = {
        g_power_label: benchmark_power_substitution_series(target_series, power=modulus, order=order),
        **{
            f"PP{residue}": _plus_residue_pochhammer_series(residue=residue, modulus=modulus, order=order)
            for residue in range(1, modulus)
        },
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


def scan_ratio_self_plus_pochhammer_eta_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[SelfPlusPochhammerEtaRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    unique_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not unique_moduli or not unique_levels:
        return []

    scans: list[SelfPlusPochhammerEtaRelationScan] = []
    for modulus in unique_moduli:
        for level in unique_levels:
            try:
                relation = search_self_plus_pochhammer_eta_relation(
                    target_series=ratio_series,
                    modulus=modulus,
                    level=level,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
                scans.append(SelfPlusPochhammerEtaRelationScan(modulus=modulus, level=level, relation=relation))
            except ValueError as exc:
                scans.append(
                    SelfPlusPochhammerEtaRelationScan(
                        modulus=modulus,
                        level=level,
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


def search_modular_unit_eta_relation(
    *,
    target_series: Series,
    modulus: int,
    level: int,
    order: int,
    max_abs_exponent: int = 8,
) -> MultiplicativeRelation | None:
    """Search F = prod (1-t^r)^a_r (1+t^r)^b_r * eta_tail with bounded integer exponents."""
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    if level < 1:
        raise ValueError("level must be positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        raise ValueError("target series must have constant term 1 for a modular-unit / eta search")

    basis_series_by_variable = _modular_unit_eta_basis_series(
        modulus=modulus,
        level=level,
        order=order,
    )
    return search_multiplicative_relation(
        target_series=target_series,
        basis_series_by_variable=basis_series_by_variable,
        order=order,
        max_abs_exponent=max_abs_exponent,
    )


def scan_ratio_modular_unit_eta_relations(
    *,
    ratio_series: Series,
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
) -> list[ModularUnitEtaRelationScan]:
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    unique_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not unique_moduli or not unique_levels:
        return []

    scans: list[ModularUnitEtaRelationScan] = []
    for modulus in unique_moduli:
        for level in unique_levels:
            try:
                relation = search_modular_unit_eta_relation(
                    target_series=ratio_series,
                    modulus=modulus,
                    level=level,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
                scans.append(
                    ModularUnitEtaRelationScan(
                        modulus=modulus,
                        level=level,
                        relation=relation,
                    )
                )
            except ValueError as exc:
                scans.append(
                    ModularUnitEtaRelationScan(
                        modulus=modulus,
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


def scan_source_family_self_plus_pochhammer_eta_corrections(
    *,
    target_series: Series,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    powers: tuple[int, ...],
    moduli: tuple[int, ...],
    eta_levels: tuple[int, ...],
    order: int,
    max_abs_exponent: int = 8,
    supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> list[SourceFamilySelfPlusPochhammerEtaCorrectionScan]:
    if not ordered_base_families:
        return []

    unique_powers = tuple(sorted({power for power in powers if power >= 2}))
    unique_moduli = tuple(sorted({modulus for modulus in moduli if modulus >= 2}))
    normalized_eta_levels = tuple(sorted({level for level in eta_levels if level >= 1}))
    if not unique_moduli or not normalized_eta_levels:
        return []

    scans: list[SourceFamilySelfPlusPochhammerEtaCorrectionScan] = []
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
                    f"{family_label}_Q{int(label.removeprefix(family_label))}",
                    f"{label} / {family_label}",
                    series_div(basis_series, base_series),
                )
            )

        direct_basis_scans = tuple(
            SelfPlusPochhammerEtaCorrectionBasisScan(
                basis_label=label,
                basis_expression=expr,
                basis_series=basis_series,
                self_scans=tuple(
                    scan_ratio_self_plus_pochhammer_eta_relations(
                        ratio_series=series_div(target_series, basis_series),
                        moduli=unique_moduli,
                        eta_levels=normalized_eta_levels,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                ),
            )
            for label, expr, basis_series in direct_basis_entries
        )
        quotient_basis_scans = tuple(
            SelfPlusPochhammerEtaCorrectionBasisScan(
                basis_label=label,
                basis_expression=expr,
                basis_series=basis_series,
                self_scans=tuple(
                    scan_ratio_self_plus_pochhammer_eta_relations(
                        ratio_series=series_div(target_series, basis_series),
                        moduli=unique_moduli,
                        eta_levels=normalized_eta_levels,
                        order=order,
                        max_abs_exponent=max_abs_exponent,
                    )
                ),
            )
            for label, expr, basis_series in quotient_basis_entries
        )

        scans.append(
            SourceFamilySelfPlusPochhammerEtaCorrectionScan(
                family_label=family_label,
                benchmark_name=benchmark_name,
                direct_basis_scans=direct_basis_scans,
                quotient_basis_scans=quotient_basis_scans,
            )
        )
    return scans


def scan_tail_family_source_eta_ladder(
    *,
    reduced_coeffs: ContinuedFractionCoeffs,
    symbol: sp.Symbol,
    ordered_base_families: tuple[tuple[str, str, Series], ...],
    start_stages: tuple[int, ...] = (3, 4, 5),
    max_gap_depth: int = 3,
    order: int,
    powers: tuple[int, ...] = (2, 3, 4),
    eta_levels: tuple[int, ...] = (1, 2, 4),
    max_abs_exponent: int = 8,
    gg_order: int | None = None,
    gg_benchmark_name: str | None = None,
    gg_base_series: Series | None = None,
    gg_degree_values: tuple[int, ...] = (1, 2),
    gg_max_abs_exponent: int = 8,
    gg_solve_order: int | None = None,
    gg_supplemental_powers: tuple[int, ...] = (),
    morton_order: int | None = None,
) -> tuple[TailFamilySourceEtaSample, ...]:
    if order < 2:
        raise ValueError("order must be at least 2")
    if max_gap_depth < 0:
        raise ValueError("max_gap_depth must be non-negative")

    normalized_stages = tuple(sorted({stage for stage in start_stages if stage >= 1}))
    if not normalized_stages:
        return ()
    gg_scan_order = order if gg_order is None else max(2, min(gg_order, order))
    morton_scan_order = order if morton_order is None else max(2, min(morton_order, order))

    series_symbol = str(symbol)
    samples: list[TailFamilySourceEtaSample] = []
    for start_stage in normalized_stages:
        anchor = build_reduced_tail_anchor(
            reduced_coeffs=reduced_coeffs,
            symbol=symbol,
            start_stage=start_stage,
            order=order,
        )
        if anchor is None:
            continue

        state_expr = sp.expand(anchor.state_expr)
        state_text = _format_expr(state_expr)
        state_suffix = "".join(ch for ch in state_text if ch.isalnum()) or f"s{start_stage}"
        base_label = f"U_{state_suffix}"
        current_label = base_label
        current_expression = f"{current_label} = T({state_text}) / (1 + {state_text})"
        current_series = list(anchor.normalized_series)

        for gap_depth in range(max_gap_depth + 1):
            truncated_base_families = tuple(
                (label, benchmark_name, basis_series[:order])
                for label, benchmark_name, basis_series in ordered_base_families
            )
            direct_eta_scans = tuple(
                scan_ratio_eta_quotient_relations(
                    ratio_series=current_series,
                    levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            )
            direct_modular_unit_eta_scans = tuple(
                scan_ratio_modular_unit_eta_relations(
                    ratio_series=current_series,
                    moduli=powers,
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            )
            source_scans = tuple(
                scan_source_family_eta_corrections(
                    target_series=current_series,
                    ordered_base_families=truncated_base_families,
                    powers=powers,
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            )
            weighted_correction_source_families = tuple(
                entry for entry in truncated_base_families if entry[0] in {"RR", "GG"}
            )
            gg_scan = (
                None
                if gg_benchmark_name is None or gg_base_series is None
                else scan_gg_modular_equation_box(
                    target_series=current_series[:gg_scan_order],
                    benchmark_name=gg_benchmark_name,
                    gg_series=gg_base_series[:gg_scan_order],
                    order=gg_scan_order,
                    degree_values=gg_degree_values,
                    max_abs_exponent=gg_max_abs_exponent,
                    solve_order=gg_solve_order,
                    supplemental_powers=gg_supplemental_powers,
                    weighted_correction_eta_levels=eta_levels,
                    weighted_correction_moduli=tuple(modulus for modulus in powers if modulus <= 4) or (2, 3, 4),
                    weighted_correction_max_abs_exponent=max_abs_exponent,
                    weighted_correction_source_families=weighted_correction_source_families,
                    weighted_correction_source_powers=powers,
                )
            )
            morton_scan = scan_morton_periodic_point_box(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                max_abs_exponent=max_abs_exponent,
            )
            weber_g_class_invariant_scan = scan_weber_class_invariant_box(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                eta_levels=eta_levels,
                moduli=powers,
                max_abs_exponent=max_abs_exponent,
            )
            weber_p_class_invariant_scan = scan_weber_p_class_invariant_box(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                eta_levels=eta_levels,
                moduli=powers,
                max_abs_exponent=max_abs_exponent,
            )
            weber_residual_bridge_scan = scan_weber_class_invariant_bridge_box(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                eta_levels=eta_levels,
                moduli=powers,
                max_abs_exponent=max_abs_exponent,
            )
            weber_j_pb_bridge_scan = _build_weber_j_pb_bridge_scan(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                eta_levels=eta_levels,
                moduli=powers,
                max_abs_exponent=max_abs_exponent,
            )
            weber_j_lift_pivot_bridge_scans = _build_weber_j_lift_pivot_bridge_scans(
                target_series=current_series[:morton_scan_order],
                order=morton_scan_order,
                eta_levels=eta_levels,
                moduli=powers,
                max_abs_exponent=max_abs_exponent,
            )
            samples.append(
                TailFamilySourceEtaSample(
                    label=current_label,
                    expression=current_expression,
                    start_stage=start_stage,
                    gap_depth=gap_depth,
                    state_expr=state_expr,
                    series=tuple(current_series),
                    direct_eta_scans=direct_eta_scans,
                    direct_modular_unit_eta_scans=direct_modular_unit_eta_scans,
                    source_family_eta_scans=source_scans,
                    gg_modular_equation_scan=gg_scan,
                    morton_periodic_point_scan=morton_scan,
                    weber_g_class_invariant_scan=weber_g_class_invariant_scan,
                    weber_p_class_invariant_scan=weber_p_class_invariant_scan,
                    weber_residual_bridge_scan=weber_residual_bridge_scan,
                    weber_j_pb_bridge_scan=weber_j_pb_bridge_scan,
                    weber_j_lift_pivot_bridge_scans=weber_j_lift_pivot_bridge_scans,
                )
            )
            if gap_depth == max_gap_depth:
                break
            gap = build_gap_normalized_series(target_series=current_series)
            if gap is None:
                break
            next_label = f"{base_label}_g{gap_depth + 1}"
            current_expression = _format_gap_normalization_formula(
                source_variable=current_label,
                target_variable=next_label,
                gap=gap,
                series_symbol=series_symbol,
            )
            current_label = next_label
            current_series = list(gap.normalized_series)

    return tuple(samples)


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


def _format_self_plus_pochhammer_relation(
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
        if name.startswith("PP"):
            residue = int(name.removeprefix("PP"))
            base = (
                f"(-{series_symbol}; {series_symbol}^{modulus})_inf"
                if residue == 1
                else f"(-{series_symbol}^{residue}; {series_symbol}^{modulus})_inf"
            )
            terms.append(base if exponent == 1 else f"{base}^{exponent}")
            continue
        terms.append(name if exponent == 1 else f"{name}^{exponent}")
    rhs = " * ".join(terms) if terms else "1"
    return f"{target_variable} = {rhs}"


def _format_self_plus_pochhammer_eta_relation(
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
        if name.startswith("PP"):
            residue = int(name.removeprefix("PP"))
            base = (
                f"(-{series_symbol}; {series_symbol}^{modulus})_inf"
                if residue == 1
                else f"(-{series_symbol}^{residue}; {series_symbol}^{modulus})_inf"
            )
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


def _format_self_mahler_linear_relation(
    relation: PolynomialRelation,
    *,
    target_variable: str,
    series_symbol: str,
) -> str:
    symbol_map: dict[str, sp.Symbol] = {
        "T": sp.Symbol(series_symbol),
        "F": sp.Symbol(target_variable),
    }
    for variable in relation.variables:
        if variable.startswith("G"):
            symbol_map[variable] = sp.Symbol(variable)
    expr = relation.as_sympy(tuple(symbol_map[name] for name in relation.variables))
    return f"{_format_expr(sp.expand(expr))} = 0"


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


def _format_gap_normalization_formula(
    *,
    source_variable: str,
    target_variable: str,
    gap: GapNormalizedSeries,
    series_symbol: str,
) -> str:
    leading_coeff = _format_expr(gap.leading_coefficient)
    shift = gap.shift
    if leading_coeff == "1":
        return f"{target_variable} = ({source_variable} - 1) / {series_symbol}^{shift}"
    if leading_coeff == "-1":
        return f"{target_variable} = (1 - {source_variable}) / {series_symbol}^{shift}"
    return f"{target_variable} = ({source_variable} - 1) / ({leading_coeff}*{series_symbol}^{shift})"


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


def _format_modular_unit_eta_relation(
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
        if name.startswith("M"):
            residue = int(name.removeprefix("M"))
            base = f"(1 - {series_symbol})" if residue == 1 else f"(1 - {series_symbol}^{residue})"
        elif name.startswith("P"):
            residue = int(name.removeprefix("P"))
            base = f"(1 + {series_symbol})" if residue == 1 else f"(1 + {series_symbol}^{residue})"
        elif name.startswith("E"):
            divisor = int(name.removeprefix("E"))
            base = (
                f"({series_symbol}; {series_symbol})_inf"
                if divisor == 1
                else f"({series_symbol}^{divisor}; {series_symbol}^{divisor})_inf"
            )
        else:
            base = name
        terms.append(base if exponent == 1 else f"{base}^{exponent}")
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


def _flatten_source_family_eta_hits(
    scans: list[SourceFamilyEtaCorrectionScan] | tuple[SourceFamilyEtaCorrectionScan, ...],
) -> list[tuple[str, str, str, int, MultiplicativeRelation]]:
    hits: list[tuple[str, str, str, int, MultiplicativeRelation]] = []
    for family_scan in scans:
        for basis_scan in family_scan.direct_basis_scans:
            for eta_scan in basis_scan.eta_scans:
                if eta_scan.relation is not None:
                    hits.append(
                        (
                            basis_scan.basis_label,
                            basis_scan.basis_expression,
                            "raw",
                            eta_scan.level,
                            eta_scan.relation,
                        )
                    )
        for basis_scan in family_scan.quotient_basis_scans:
            for eta_scan in basis_scan.eta_scans:
                if eta_scan.relation is not None:
                    hits.append(
                        (
                            basis_scan.basis_label,
                            basis_scan.basis_expression,
                            "quotient",
                            eta_scan.level,
                            eta_scan.relation,
                        )
                    )
    return hits


def _named_prefix_scan_hit_labels(scans, *, hit_predicate) -> tuple[str, ...]:
    return tuple(
        scan.basis_labels[-1]
        for scan in scans
        if scan.error is None and hit_predicate(scan)
    )


def _named_prefix_scan_checked_count(scans) -> int:
    return sum(1 for scan in scans if scan.error is None)


def _named_prefix_scan_skipped_count(scans) -> int:
    return sum(1 for scan in scans if scan.error is not None)


def _format_tail_prefix_summary(scans, *, label: str, hit_predicate) -> str:
    hit_labels = _named_prefix_scan_hit_labels(scans, hit_predicate=hit_predicate)
    checked_count = _named_prefix_scan_checked_count(scans)
    skipped_count = _named_prefix_scan_skipped_count(scans)
    summary = f"{label} `{len(hit_labels)}` / `{checked_count}` hit prefixes"
    if hit_labels:
        summary += f" ({', '.join(f'`{value}`' for value in hit_labels)})"
    if skipped_count:
        summary += f"; skipped `{skipped_count}`"
    return summary


def _gg_modular_equation_scan_has_hit(scan: GGModularEquationScan) -> bool:
    if scan.hit_templates:
        return True
    if scan.exact_polynomial_template_hits or scan.quotient_exact_polynomial_template_hits:
        return True
    if any(_gg_weighted_coordinate_diagnostic_has_hit(item) for item in scan.weighted_coordinate_diagnostics):
        return True
    if any(item.relation is not None for item in scan.polynomial_scans):
        return True
    if any(item.relation is not None for item in scan.multiplicative_scans):
        return True
    if any(item.relation is not None for item in scan.fractional_linear_scans):
        return True
    if any(item.total_hits > 0 for item in scan.two_layer_fractional_linear_scans):
        return True
    if any(item.relation is not None for item in scan.quotient_polynomial_scans):
        return True
    if any(item.relation is not None for item in scan.quotient_multiplicative_scans):
        return True
    if any(item.relation is not None for item in scan.quotient_fractional_linear_scans):
        return True
    if any(item.total_hits > 0 for item in scan.quotient_two_layer_fractional_linear_scans):
        return True
    if any(item.relation is not None for item in scan.mixed_quotient_polynomial_scans):
        return True
    if any(item.relation is not None for item in scan.mixed_quotient_multiplicative_scans):
        return True
    if any(item.relation is not None for item in scan.mixed_quotient_fractional_linear_scans):
        return True
    if any(item.total_hits > 0 for item in scan.mixed_quotient_two_layer_fractional_linear_scans):
        return True
    return False


def _gg_modular_equation_scan_has_exact_quotient_hit(scan: GGModularEquationScan) -> bool:
    return bool(scan.quotient_exact_polynomial_template_hits)


def _morton_periodic_point_scan_has_hit(scan: MortonPeriodicPointScan) -> bool:
    return any(item.hit for item in scan.template_results) or any(
        item.hit
        for coordinate_scan in scan.named_coordinate_scans
        for item in coordinate_scan.template_results
    )


def _weber_class_invariant_scan_has_hit(scan: WeberClassInvariantScan) -> bool:
    return (
        scan.template_hit
        or any(item.relation is not None for item in scan.direct_eta_scans)
        or any(item.relation is not None for item in scan.direct_modular_unit_eta_scans)
        or any(item.relation is not None for item in scan.correction_self_plus_pochhammer_scans)
        or any(item.relation is not None for item in scan.correction_self_plus_pochhammer_eta_scans)
    )


def _constant_one_series_scan_has_hit(scan: ConstantOneSeriesScan) -> bool:
    return (
        bool(scan.self_polynomial_scan.hits)
        or bool(scan.self_fractional_linear_scan.hits)
        or any(item.relation is not None for item in scan.self_quotient_product_scans)
        or any(item.relation is not None for item in scan.eta_scans)
        or any(item.relation is not None for item in scan.modular_unit_eta_scans)
        or any(item.relation is not None for item in scan.self_plus_pochhammer_scans)
        or any(item.relation is not None for item in scan.self_plus_pochhammer_eta_scans)
    )


def _gg_weighted_coordinate_diagnostic_has_hit(diagnostic: GGWeightedCoordinateDiagnostic) -> bool:
    return (
        (diagnostic.first_difference_power is None and diagnostic.first_difference_coeff is None)
        or (diagnostic.first_log_difference_power is None and diagnostic.first_log_difference_coeff is None)
        or any(scan.relation is not None for scan in diagnostic.correction_eta_scans)
        or any(scan.relation is not None for scan in diagnostic.correction_modular_unit_eta_scans)
        or any(scan.relation is not None for scan in diagnostic.normalized_correction_eta_scans)
        or any(scan.relation is not None for scan in diagnostic.normalized_correction_modular_unit_eta_scans)
        or bool(_flatten_source_family_eta_hits(diagnostic.normalized_correction_source_family_eta_scans))
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_eta_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_modular_unit_eta_scans)
        or bool(_flatten_source_family_eta_hits(diagnostic.second_normalized_correction_source_family_eta_scans))
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_quotient_polynomial_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_quotient_multiplicative_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_quotient_fractional_linear_scans)
        or any(scan.total_hits > 0 for scan in diagnostic.second_normalized_correction_quotient_two_layer_fractional_linear_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_mixed_quotient_polynomial_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_mixed_quotient_multiplicative_scans)
        or any(scan.relation is not None for scan in diagnostic.second_normalized_correction_mixed_quotient_fractional_linear_scans)
        or any(scan.total_hits > 0 for scan in diagnostic.second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans)
        or any(scan.hits for scan in diagnostic.second_normalized_correction_explicit_transform_eta_scans)
        or diagnostic.polynomial_degree1_relation is not None
        or diagnostic.polynomial_degree2_relation is not None
        or diagnostic.fractional_linear_relation is not None
    )


def _format_exact_polynomial_obstruction(
    obstruction: tuple[str, int | None, sp.Expr | None],
    *,
    series_symbol: str,
) -> str:
    label, power, coeff = obstruction
    if power is None or coeff is None:
        return f"`{label}` matches through the checked truncation"
    return f"`{label}` first fails at `{series_symbol}^{power}` with coefficient `{_format_expr(coeff)}`"


def _format_weber_class_invariant_obstruction(
    scan: WeberClassInvariantScan,
    *,
    series_symbol: str,
) -> str:
    power = scan.template_first_failure_power
    coeff = scan.template_first_failure_coeff
    if power is None or coeff is None:
        return f"`{scan.template_label}` matches through the checked truncation"
    return (
        f"`{scan.template_label}` first differs from `{scan.template_expression}` "
        f"at `{series_symbol}^{power}` with coefficient `{_format_expr(coeff)}`"
    )


def _format_weighted_coordinate_obstruction(
    *,
    power: int | None,
    coeff: sp.Expr | None,
    lhs: str,
    series_symbol: str,
) -> str:
    if power is None or coeff is None:
        return f"`{lhs}` matches through the checked truncation"
    return f"`{lhs}` first fails at `{series_symbol}^{power}` with coefficient `{_format_expr(coeff)}`"


def _append_named_gg_bridge_lines(
    lines: list[str],
    *,
    prefix: str,
    gg_scan: GGModularEquationScan,
    series_symbol: str,
) -> None:
    direct_prefix_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                gg_scan.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.two_layer_fractional_linear_scans,
                label="two-layer fractional-linear",
                hit_predicate=lambda item: item.total_hits > 0,
            ),
        )
    )
    quotient_prefix_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                gg_scan.quotient_polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.quotient_multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.quotient_fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.quotient_two_layer_fractional_linear_scans,
                label="two-layer fractional-linear",
                hit_predicate=lambda item: item.total_hits > 0,
            ),
        )
    )
    mixed_prefix_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                gg_scan.mixed_quotient_polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.mixed_quotient_multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.mixed_quotient_fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                gg_scan.mixed_quotient_two_layer_fractional_linear_scans,
                label="two-layer fractional-linear",
                hit_predicate=lambda item: item.total_hits > 0,
            ),
        )
    )
    exact_template_summary = (
        f"`{len(gg_scan.hit_templates)}` / `{len(gg_scan.checked_templates)}` exact template hits"
    )
    if gg_scan.hit_templates:
        exact_template_summary += (
            f" ({', '.join(f'`{label}`' for label in gg_scan.hit_templates)})"
        )
    direct_exact_summary = (
        f"`{len(gg_scan.exact_polynomial_template_hits)}` / `{len(gg_scan.exact_polynomial_template_labels)}` exact direct Chan--Huang hits"
    )
    if gg_scan.exact_polynomial_template_hits:
        direct_exact_summary += (
            f" ({', '.join(f'`{label}`' for label in gg_scan.exact_polynomial_template_hits)})"
        )
    quotient_exact_summary = (
        f"`{len(gg_scan.quotient_exact_polynomial_template_hits)}` / `{len(gg_scan.quotient_exact_polynomial_template_labels)}` exact quotient-coordinate Chan--Huang hits"
    )
    if gg_scan.quotient_exact_polynomial_template_hits:
        quotient_exact_summary += (
            f" ({', '.join(f'`{label}`' for label in gg_scan.quotient_exact_polynomial_template_hits)})"
        )
    lines.extend(
        [
            f"- {prefix} named `GG` exact templates: {exact_template_summary}.",
            f"- {prefix} named `GG` direct prefixes: {direct_prefix_summary}.",
            f"- {prefix} named `GG` quotient prefixes: {quotient_prefix_summary}.",
            f"- {prefix} named `GG` mixed quotient prefixes: {mixed_prefix_summary}.",
            f"- {prefix} named `GG` direct exact modular-equation templates: {direct_exact_summary}.",
            f"- {prefix} named `GG` quotient exact modular-equation templates: {quotient_exact_summary}.",
            f"- {prefix} named `GG` direct exact obstruction witnesses: "
            + "; ".join(
                _format_exact_polynomial_obstruction(
                    obstruction,
                    series_symbol=series_symbol,
                )
                for obstruction in gg_scan.exact_polynomial_template_obstructions
            )
            + ".",
            f"- {prefix} named `GG` quotient exact obstruction witnesses: "
            + "; ".join(
                _format_exact_polynomial_obstruction(
                    obstruction,
                    series_symbol=series_symbol,
                )
                for obstruction in gg_scan.quotient_exact_polynomial_template_obstructions
            )
            + ".",
        ]
    )


def _append_constant_one_scan_summary_lines(
    lines: list[str],
    *,
    prefix: str,
    scan: ConstantOneSeriesScan,
    series_symbol: str,
) -> None:
    self_product_hits = [
        item for item in scan.self_quotient_product_scans if item.relation is not None
    ]
    eta_hits = [item for item in scan.eta_scans if item.relation is not None]
    modular_hits = [item for item in scan.modular_unit_eta_scans if item.relation is not None]
    plus_hits = [item for item in scan.self_plus_pochhammer_scans if item.relation is not None]
    plus_eta_hits = [
        item for item in scan.self_plus_pochhammer_eta_scans if item.relation is not None
    ]
    lines.append(f"- {prefix}: `{scan.expression}`.")
    if scan.first_failure_power is None or scan.first_failure_coeff is None:
        lines.append(f"- {prefix}: matches `1` through the checked truncation.")
    else:
        lines.append(
            f"- {prefix}: `{scan.label} - 1` first fails at "
            f"`{series_symbol}^{scan.first_failure_power}` with coefficient "
            f"`{_format_expr(scan.first_failure_coeff)}`."
        )
    lines.append(
        f"- {prefix} self-polynomial uniqueness boxes: "
        f"`{len(scan.self_polynomial_scan.hits)}` / "
        f"`{len(scan.self_polynomial_scan.moduli_checked) * len(scan.self_polynomial_scan.fg_degree_values) * len(scan.self_polynomial_scan.t_degree_values)}` hit boxes."
    )
    lines.append(
        f"- {prefix} self-fractional-linear uniqueness boxes: "
        f"`{len(scan.self_fractional_linear_scan.hits)}` / "
        f"`{len(scan.self_fractional_linear_scan.moduli_checked) * len(scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
    )
    lines.append(
        f"- {prefix} self-quotient finite-product boxes: "
        f"`{len(self_product_hits)}` / `{len(scan.self_quotient_product_scans)}` hit boxes."
    )
    lines.append(
        f"- {prefix} eta templates: `{len(eta_hits)}` / `{len(scan.eta_scans)}` hit boxes."
    )
    lines.append(
        f"- {prefix} modular-unit / eta templates: "
        f"`{len(modular_hits)}` / `{len(scan.modular_unit_eta_scans)}` hit boxes."
    )
    lines.append(
        f"- {prefix} plus-Pochhammer templates: "
        f"`{len(plus_hits)}` / `{len(scan.self_plus_pochhammer_scans)}` hit boxes."
    )
    lines.append(
        f"- {prefix} plus-Pochhammer + eta templates: "
        f"`{len(plus_eta_hits)}` / `{len(scan.self_plus_pochhammer_eta_scans)}` hit boxes."
    )


def _append_named_gg_descendant_focus_lines(
    lines: list[str],
    *,
    prefix: str,
    descendant_scan: GGDescendantFocusedScan,
) -> None:
    direct_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                descendant_scan.direct_scans.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                descendant_scan.direct_scans.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                descendant_scan.direct_scans.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
        )
    )
    quotient_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                descendant_scan.quotient_scans.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                descendant_scan.quotient_scans.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                descendant_scan.quotient_scans.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
        )
    )
    lines.extend(
        [
            f"- {prefix} odd-prime descendant micro-scan: checked direct ladder "
            f"`{', '.join(descendant_scan.direct_labels)}` and quotient ladder "
            f"`{', '.join(descendant_scan.quotient_labels)}` through order "
            f"`{descendant_scan.order_checked}` with degree box `{descendant_scan.degree_values}` "
            f"and max exponent `{descendant_scan.max_abs_exponent}`.",
            f"- {prefix} odd-prime descendant direct micro-boxes: {direct_summary}.",
            f"- {prefix} odd-prime descendant quotient micro-boxes: {quotient_summary}.",
        ]
    )


def _append_named_coordinate_orbit_lines(
    lines: list[str],
    *,
    prefix: str,
    orbit_scan: NamedCoordinateOrbitScan,
) -> None:
    direct_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                orbit_scan.direct_scans.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.direct_scans.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.direct_scans.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
        )
    )
    quotient_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                orbit_scan.quotient_scans.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.quotient_scans.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.quotient_scans.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
        )
    )
    mixed_summary = "; ".join(
        (
            _format_tail_prefix_summary(
                orbit_scan.mixed_quotient_scans.polynomial_scans,
                label="polynomial",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.mixed_quotient_scans.multiplicative_scans,
                label="multiplicative",
                hit_predicate=lambda item: item.relation is not None,
            ),
            _format_tail_prefix_summary(
                orbit_scan.mixed_quotient_scans.fractional_linear_scans,
                label="fractional-linear",
                hit_predicate=lambda item: item.relation is not None,
            ),
        )
    )
    direct_labels = ", ".join(f"`{label}`" for label, _, _ in orbit_scan.ordered_basis_series)
    quotient_labels = ", ".join(f"`{label}`" for label, _, _ in orbit_scan.quotient_basis_series)
    mixed_labels = ", ".join(f"`{label}`" for label, _, _ in orbit_scan.mixed_quotient_basis_series)
    lines.extend(
        [
            f"- {prefix} focused `{orbit_scan.family_label}` direct ladder: {direct_labels}.",
            f"- {prefix} focused `{orbit_scan.family_label}` quotient ladder: {quotient_labels}.",
            f"- {prefix} focused `{orbit_scan.family_label}` mixed ladder: {mixed_labels}.",
            f"- {prefix} focused `{orbit_scan.family_label}` direct prefixes: {direct_summary}.",
            f"- {prefix} focused `{orbit_scan.family_label}` quotient prefixes: {quotient_summary}.",
            f"- {prefix} focused `{orbit_scan.family_label}` mixed prefixes: {mixed_summary}.",
        ]
    )


def _named_prefix_box_scan_hit_count(scans: NamedPrefixBoxScans) -> tuple[int, int]:
    total = (
        len(scans.polynomial_scans)
        + len(scans.multiplicative_scans)
        + len(scans.fractional_linear_scans)
        + len(scans.two_layer_fractional_linear_scans)
    )
    hits = (
        sum(1 for item in scans.polynomial_scans if item.relation is not None)
        + sum(1 for item in scans.multiplicative_scans if item.relation is not None)
        + sum(1 for item in scans.fractional_linear_scans if item.relation is not None)
        + sum(1 for item in scans.two_layer_fractional_linear_scans if item.total_hits > 0)
    )
    return hits, total


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


def _gg_descendant_preview(
    *,
    base_series: Series,
    order: int,
    supplemental_powers: tuple[int, ...] = (),
) -> GGDescendantPreview | None:
    descendant_powers = tuple(sorted({value for value in supplemental_powers if value >= 5}))
    if not descendant_powers:
        return None
    ordered_basis_entries = _gg_modular_equation_ordered_basis_series(
        base_series=base_series,
        order=order,
        supplemental_powers=descendant_powers,
    )
    quotient_basis_entries = _gg_modular_equation_quotient_basis_series(ordered_basis_entries)
    direct_labels = tuple(f"GG{power}" for power in descendant_powers)
    quotient_labels = tuple(f"Q_{power}" for power in descendant_powers)
    available_direct = {label for label, _, _ in ordered_basis_entries}
    available_quotient = {label for label, _, _ in quotient_basis_entries}
    return GGDescendantPreview(
        direct_labels=tuple(label for label in direct_labels if label in available_direct),
        quotient_labels=tuple(label for label in quotient_labels if label in available_quotient),
    )


def _scan_gg_descendant_focus_box(
    *,
    target_series: Series,
    gg_series: Series,
    order: int,
    supplemental_powers: tuple[int, ...] = (),
    degree_values: tuple[int, ...] = (1,),
    max_abs_exponent: int = 2,
) -> GGDescendantFocusedScan | None:
    descendant_powers = tuple(sorted({value for value in supplemental_powers if value >= 5}))
    if not descendant_powers:
        return None
    checked_order = min(order, 8)
    ordered_basis_entries = _gg_modular_equation_ordered_basis_series(
        base_series=gg_series[:checked_order],
        order=checked_order,
        supplemental_powers=descendant_powers,
    )
    direct_basis_series = tuple(
        (label, series)
        for label, _, series in ordered_basis_entries
        if label in {f"GG{power}" for power in descendant_powers}
    )
    quotient_basis_entries = _gg_modular_equation_quotient_basis_series(ordered_basis_entries)
    quotient_basis_series = tuple(
        (label, series)
        for label, _, series in quotient_basis_entries
        if label in {f"Q_{power}" for power in descendant_powers}
    )
    if not direct_basis_series or not quotient_basis_series:
        return None
    direct_scans = scan_named_prefix_boxes(
        target_series=target_series[:checked_order],
        ordered_basis_series=direct_basis_series,
        order=checked_order,
        degree_values=degree_values,
        max_abs_exponent=max_abs_exponent,
        solve_order=checked_order,
        include_two_layer=False,
    )
    quotient_scans = scan_named_prefix_boxes(
        target_series=target_series[:checked_order],
        ordered_basis_series=quotient_basis_series,
        order=checked_order,
        degree_values=degree_values,
        max_abs_exponent=max_abs_exponent,
        solve_order=checked_order,
        include_two_layer=False,
    )
    return GGDescendantFocusedScan(
        direct_labels=tuple(label for label, _ in direct_basis_series),
        quotient_labels=tuple(label for label, _ in quotient_basis_series),
        order_checked=checked_order,
        degree_values=tuple(sorted({degree for degree in degree_values if degree >= 1})),
        max_abs_exponent=max_abs_exponent,
        direct_scans=direct_scans,
        quotient_scans=quotient_scans,
    )


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


@lru_cache(maxsize=None)
def _morton_periodic_point_relations(order: int) -> tuple[tuple[str, PolynomialRelation], ...]:
    x, y = sp.symbols("F G")
    morton_f = _polynomial_relation_from_sympy_expr(
        expr=y**2 + (x**2 - 1) * y + x**2,
        variables=("F", "G"),
        order=order,
    )
    morton_g = _polynomial_relation_from_sympy_expr(
        expr=y**2 - (x**2 - 4 * x + 1) * y + x**2,
        variables=("F", "G"),
        order=order,
    )
    return (
        ("Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`", morton_f),
        ("Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`", morton_g),
        ("Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`", morton_f),
        ("Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`", morton_f),
    )


@lru_cache(maxsize=None)
def _morton_squared_coordinate_relations(order: int) -> tuple[tuple[str, PolynomialRelation], ...]:
    x, x2 = sp.symbols("X X_2")
    return (
        (
            "Morton Prop. 3.2 squared-coordinate template `X_2^2 - (X^2 - 4*X + 1)*X_2 + X^2`",
            _polynomial_relation_from_sympy_expr(
                expr=x2**2 - (x**2 - 4 * x + 1) * x2 + x**2,
                variables=("X", "X_2"),
                order=order,
            ),
        ),
    )


@lru_cache(maxsize=None)
def _morton_transformed_squared_coordinate_relations(
    order: int,
) -> tuple[tuple[str, PolynomialRelation], ...]:
    t_coord, t_coord2 = sp.symbols("T T_2")
    return (
        (
            "Morton Eq. (3.6) transformed squared-coordinate template `T^2 - (T_2^2 - 4*T_2 + 1)*T + T_2^2`",
            _polynomial_relation_from_sympy_expr(
                expr=t_coord**2 - (t_coord2**2 - 4 * t_coord2 + 1) * t_coord + t_coord2**2,
                variables=("T", "T_2"),
                order=order,
            ),
        ),
    )


@lru_cache(maxsize=None)
def _morton_weber_schlafli_relations(order: int) -> tuple[tuple[str, PolynomialRelation], ...]:
    p, p2 = sp.symbols("P P_2")
    return (
        (
            "Morton Weber-Schlafli template `P^2*P_2^2 + P^2 - 2*P_2`",
            _polynomial_relation_from_sympy_expr(
                expr=p**2 * p2**2 + p**2 - 2 * p2,
                variables=("P", "P_2"),
                order=order,
            ),
        ),
    )


@lru_cache(maxsize=None)
def _morton_weber_companion_relations(order: int) -> tuple[tuple[str, PolynomialRelation], ...]:
    b, b2, p = sp.symbols("B B_2 P")
    return (
        (
            "Morton Weber companion template `B^2 - B_2 - 4`",
            _polynomial_relation_from_sympy_expr(
                expr=b**2 - b2 - 4,
                variables=("B", "B_2"),
                order=order,
            ),
        ),
        (
            "Morton Weber companion template `B_2^4 - P^8 - 16*P^4`",
            _polynomial_relation_from_sympy_expr(
                expr=b2**4 - p**8 - 16 * p**4,
                variables=("B_2", "P"),
                order=order,
            ),
        ),
    )


def _first_nonzero_residual_term(
    residual: Series,
) -> tuple[int | None, sp.Expr | None]:
    for power, coeff in enumerate(residual):
        simplified = sp.simplify(coeff)
        if simplified != 0:
            return power, simplified
    return None, None


def _series_fractional_power_nonzero_constant(
    *,
    base: Series,
    exponent: sp.Rational,
) -> Series:
    order = len(base)
    if order < 1:
        raise ValueError("series must be non-empty")
    constant = sp.simplify(base[0])
    if constant == 0:
        raise ValueError("series constant term must be nonzero for fractional powers")

    delta = [sp.Integer(0) for _ in range(order)]
    for index in range(1, order):
        delta[index] = sp.simplify(base[index] / constant)

    result = [sp.Integer(0) for _ in range(order)]
    result[0] = sp.Integer(1)
    delta_power = delta[:]
    for degree in range(1, order):
        coefficient = sp.simplify(sp.binomial(exponent, degree))
        if coefficient != 0:
            result = series_add(
                result,
                [sp.simplify(coefficient * value) for value in delta_power],
            )
        delta_power = series_mul(delta_power, delta)

    scale = sp.simplify(constant**exponent)
    return [sp.simplify(scale * value) for value in result]


def _weber_schlafli_coordinate_series(
    *,
    target_series: Series,
    order: int,
) -> Series | None:
    if order < 1:
        return None
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0]) == 0:
        return None
    inverse_series = series_invert(target_series[:order])
    return [
        sp.simplify((inverse_series[index] - target_series[index]) / 2)
        for index in range(order)
    ]


def _weber_companion_coordinate_series(
    *,
    p_series: Series,
    order: int,
) -> Series | None:
    if order < 1:
        return None
    if len(p_series) < order:
        raise ValueError("p_series is shorter than requested order")

    p_trunc = p_series[:order]
    p_fourth = series_pow(p_trunc, 4)
    b2_source = p_fourth[:]
    b2_source[0] = sp.simplify(b2_source[0] + 16)
    b2_series = series_mul(
        p_trunc,
        _series_fractional_power_nonzero_constant(
            base=b2_source,
            exponent=sp.Rational(1, 4),
        ),
    )
    b_source = b2_series[:]
    b_source[0] = sp.simplify(b_source[0] + 4)
    return _series_fractional_power_nonzero_constant(
        base=b_source,
        exponent=sp.Rational(1, 2),
    )


@lru_cache(maxsize=None)
def _focused_weber_bridge_reference_ladder(
    *,
    order: int,
) -> WeberLeadingBridgeReferenceLadder | None:
    gg_series = list(
        _canonical_benchmark_series(
            "gollnitz_gordon_normalized",
            depth=order,
            order=order,
        )
    )
    p_series = _weber_schlafli_coordinate_series(
        target_series=gg_series,
        order=order,
    )
    if p_series is None:
        return None
    b_series = _weber_companion_coordinate_series(
        p_series=p_series,
        order=order,
    )
    if b_series is None:
        return None

    p_scan = _build_leading_normalized_coordinate_scan(
        source_label="P_ws",
        target_series=p_series,
        order=order,
        max_abs_exponent=4,
    )
    b_scan = _build_leading_normalized_coordinate_scan(
        source_label="B_ws",
        target_series=b_series,
        order=order,
        max_abs_exponent=4,
    )
    if p_scan is None or b_scan is None:
        return None

    q_pb_series = series_div(list(b_scan.normalized_series), list(p_scan.normalized_series))
    k_pb_series = _normalized_constant_one_followup_series(
        target_series=q_pb_series,
        order=order,
    )
    if k_pb_series is None:
        return None
    q_pk_series = series_div(k_pb_series, list(p_scan.normalized_series))
    l_pk_series = _normalized_constant_one_followup_series(
        target_series=q_pk_series,
        order=order,
    )
    if l_pk_series is None:
        return None

    return WeberLeadingBridgeReferenceLadder(
        q_pb_series=tuple(q_pb_series),
        k_pb_series=tuple(k_pb_series),
        q_pk_series=tuple(q_pk_series),
        l_pk_series=tuple(l_pk_series),
    )


def _morton_transformed_squared_coordinate_series(
    *,
    squared_coordinate_series: Series,
    order: int,
) -> Series | None:
    if order < 1:
        return None
    if len(squared_coordinate_series) < order:
        raise ValueError("squared_coordinate_series is shorter than requested order")

    sigma_squared = sp.expand((-1 + sp.sqrt(2)) ** 2)
    numerator = squared_coordinate_series[:order]
    numerator[0] = sp.simplify(numerator[0] - sigma_squared)
    denominator = [sp.simplify(sigma_squared * value) for value in squared_coordinate_series[:order]]
    denominator[0] = sp.simplify(denominator[0] - 1)
    if sp.simplify(denominator[0]) == 0:
        return None
    return series_div(numerator, denominator)


def _normalized_weber_g_coordinate_series(
    *,
    target_series: Series,
    order: int,
) -> Series | None:
    """Build the normalized Weber-g class-invariant coordinate on a GG-style series.

    Chan--Huang's `g_n^12` coordinate is naturally Laurent in `t`. Multiplying by
    `8 t` yields a constant-1 series which, on the true GG source, equals the
    eta product `(t^2; t^4)_inf^12`.
    """

    if order < 2:
        return None
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        return None

    one_series: Series = [sp.Integer(0) for _ in range(order)]
    one_series[0] = sp.Integer(1)
    t_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 1:
        t_series[1] = sp.Integer(1)
    t_squared_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 2:
        t_squared_series[2] = sp.Integer(1)

    target_squared = series_pow(target_series[:order], 2)
    t_target_squared = series_mul(t_series, target_squared)
    numerator = [sp.simplify(one_series[index] - t_target_squared[index]) for index in range(order)]
    numerator_squared = series_pow(numerator, 2)
    first_term = series_div(numerator_squared, target_squared)
    second_term = series_div(
        [sp.simplify(16 * value) for value in series_mul(t_squared_series, target_squared)],
        numerator_squared,
    )
    return [sp.simplify(first_term[index] - second_term[index]) for index in range(order)]


@lru_cache(maxsize=None)
def _normalized_weber_g_template_series_tuple(*, order: int) -> tuple[sp.Expr, ...]:
    if order < 2:
        raise ValueError("order must be at least 2")
    e2 = _eta_pochhammer_series(divisor=2, order=order)
    e4 = _eta_pochhammer_series(divisor=4, order=order)
    template = series_mul(series_pow(e2, 12), series_invert(series_pow(e4, 12)))
    return tuple(sp.simplify(value) for value in template)


def _normalized_weber_g_template_series(*, order: int) -> Series:
    return list(_normalized_weber_g_template_series_tuple(order=order))


def _normalized_weber_p_coordinate_series(
    *,
    target_series: Series,
    order: int,
) -> Series | None:
    """Build the normalized Weber-G class-invariant coordinate on a GG-style series."""

    if order < 2:
        return None
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")
    if sp.simplify(target_series[0] - 1) != 0:
        return None

    one_series: Series = [sp.Integer(0) for _ in range(order)]
    one_series[0] = sp.Integer(1)
    t_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 1:
        t_series[1] = sp.Integer(1)
    t_squared_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 2:
        t_squared_series[2] = sp.Integer(1)

    target_squared = series_pow(target_series[:order], 2)
    t_target_squared = series_mul(t_series, target_squared)
    numerator = [sp.simplify(one_series[index] - t_target_squared[index]) for index in range(order)]
    w_series = series_div(series_pow(numerator, 2), target_squared)
    w_squared = series_pow(w_series, 2)
    sqrt_argument = [
        sp.simplify(w_squared[index] - 16 * t_squared_series[index])
        for index in range(order)
    ]
    return series_div(
        w_squared,
        _series_fractional_power_nonzero_constant(
            base=sqrt_argument,
            exponent=sp.Rational(1, 2),
        ),
    )


@lru_cache(maxsize=None)
def _normalized_weber_p_template_series_tuple(*, order: int) -> tuple[sp.Expr, ...]:
    if order < 2:
        raise ValueError("order must be at least 2")
    pp2 = _plus_residue_pochhammer_series(residue=2, modulus=4, order=order)
    template = series_pow(pp2, 12)
    return tuple(sp.simplify(value) for value in template)


def _normalized_weber_p_template_series(*, order: int) -> Series:
    return list(_normalized_weber_p_template_series_tuple(order=order))


def scan_weber_class_invariant_box(
    *,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    max_abs_exponent: int = 8,
) -> WeberClassInvariantScan | None:
    coordinate_series = _normalized_weber_g_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if coordinate_series is None:
        return None

    template_series = _normalized_weber_g_template_series(order=order)
    correction_series = series_div(coordinate_series, template_series)
    correction_residual = [sp.simplify(value) for value in correction_series]
    correction_residual[0] = sp.simplify(correction_residual[0] - 1)
    first_failure_power, first_failure_coeff = _first_nonzero_residual_term(correction_residual)

    return WeberClassInvariantScan(
        label="g12_ws",
        expression=(
            "Z_g = ((1 - t*F^2)^2) / (4*t*F^2), "
            "g12_ws = 4*t*(Z_g - 1/Z_g)"
        ),
        template_label="Chan--Huang Weber g-coordinate template",
        template_expression="(t^2; t^4)_inf^12",
        template_hit=all(sp.simplify(value) == 0 for value in correction_residual),
        template_first_failure_power=first_failure_power,
        template_first_failure_coeff=first_failure_coeff,
        correction_label="G_g12_ws",
        correction_expression="G_g12_ws = g12_ws / (t^2; t^4)_inf^12",
        direct_eta_scans=tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=correction_series,
                levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        direct_modular_unit_eta_scans=tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=correction_series,
                moduli=moduli,
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        correction_self_plus_pochhammer_scans=tuple(
            scan_ratio_self_plus_pochhammer_relations(
                ratio_series=correction_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        correction_self_plus_pochhammer_eta_scans=tuple(
            scan_ratio_self_plus_pochhammer_eta_relations(
                ratio_series=correction_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
    )


def scan_weber_p_class_invariant_box(
    *,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    max_abs_exponent: int = 8,
) -> WeberClassInvariantScan | None:
    coordinate_series = _normalized_weber_p_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if coordinate_series is None:
        return None

    template_series = _normalized_weber_p_template_series(order=order)
    correction_series = series_div(coordinate_series, template_series)
    correction_residual = [sp.simplify(value) for value in correction_series]
    correction_residual[0] = sp.simplify(correction_residual[0] - 1)
    first_failure_power, first_failure_coeff = _first_nonzero_residual_term(correction_residual)

    return WeberClassInvariantScan(
        label="p12_ws",
        expression=(
            "Z_g = ((1 - t*F^2)^2) / (4*t*F^2), "
            "p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)"
        ),
        template_label="Chan--Huang Weber G-coordinate template",
        template_expression="(-t^2; t^4)_inf^12",
        template_hit=all(sp.simplify(value) == 0 for value in correction_residual),
        template_first_failure_power=first_failure_power,
        template_first_failure_coeff=first_failure_coeff,
        correction_label="G_p12_ws",
        correction_expression="G_p12_ws = p12_ws / (-t^2; t^4)_inf^12",
        direct_eta_scans=tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=correction_series,
                levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        direct_modular_unit_eta_scans=tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=correction_series,
                moduli=moduli,
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        correction_self_plus_pochhammer_scans=tuple(
            scan_ratio_self_plus_pochhammer_relations(
                ratio_series=correction_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        correction_self_plus_pochhammer_eta_scans=tuple(
            scan_ratio_self_plus_pochhammer_eta_relations(
                ratio_series=correction_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
    )


def scan_weber_class_invariant_bridge_box(
    *,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    max_abs_exponent: int = 8,
) -> WeberResidualBridgeScan | None:
    g_coordinate_series = _normalized_weber_g_coordinate_series(
        target_series=target_series,
        order=order,
    )
    p_coordinate_series = _normalized_weber_p_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if g_coordinate_series is None or p_coordinate_series is None:
        return None

    g_template_series = _normalized_weber_g_template_series(order=order)
    p_template_series = _normalized_weber_p_template_series(order=order)
    g_correction_series = series_div(g_coordinate_series, g_template_series)
    p_correction_series = series_div(p_coordinate_series, p_template_series)
    gg_series = _canonical_benchmark_series(
        "gollnitz_gordon_normalized",
        depth=order,
        order=order,
    )
    focused_named_gg_supplemental_powers = _supplemental_source_family_powers(
        smoke=max_abs_exponent <= 4
    ).get("GG", ())

    g_squared = series_pow(g_coordinate_series, 2)
    g_fourth = series_pow(g_squared, 2)
    p_squared = series_pow(p_coordinate_series, 2)
    p_fourth = series_pow(p_squared, 2)
    bridge_series = series_mul(g_fourth, p_squared)
    subtrahend = series_mul(g_squared, p_fourth)
    bridge_series = [
        sp.simplify(bridge_series[index] - subtrahend[index])
        for index in range(order)
    ]
    t_squared_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 2:
        t_squared_series[2] = sp.Integer(1)
    t_sixth_series: Series = [sp.Integer(0) for _ in range(order)]
    if order > 6:
        t_sixth_series[6] = sp.Integer(1)
    mixed_term = [
        sp.simplify(48 * value)
        for value in series_mul(t_squared_series, series_mul(g_squared, p_squared))
    ]
    bridge_series = [
        sp.simplify(bridge_series[index] + mixed_term[index] + 4096 * t_sixth_series[index])
        for index in range(order)
    ]
    bridge_first_failure_power, bridge_first_failure_coeff = _first_nonzero_residual_term(bridge_series)
    classical_product_coordinate_series = series_mul(g_correction_series, p_correction_series)
    classical_product_coordinate_scan = _scan_constant_one_series(
        label="G_f2_ws",
        expression=(
            "G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2)"
            " = G_g12_ws*G_p12_ws"
        ),
        target_series=classical_product_coordinate_series,
        order=order,
        eta_levels=eta_levels,
        moduli=moduli,
        max_abs_exponent=max_abs_exponent,
        followup_label="H_f2_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_solve_order=12,
        named_gg_include_weighted_coordinate_diagnostics=False,
        followup_named_gg_benchmark_name="gollnitz_gordon_normalized",
        followup_named_gg_series=gg_series,
        followup_named_gg_solve_order=12,
        followup_named_gg_include_weighted_coordinate_diagnostics=False,
    )
    canonical_j_numerator_series = series_pow(
        [
            sp.simplify((sp.Integer(16) if index == 0 else sp.Integer(0)) - value)
            for index, value in enumerate(classical_product_coordinate_series)
        ],
        3,
    )
    canonical_j_coordinate_series = series_div(
        canonical_j_numerator_series,
        [sp.simplify(3375 * value) for value in classical_product_coordinate_series],
    )
    canonical_j_coordinate_scan = _scan_constant_one_series(
        label="J_f2_ws",
        expression="J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)",
        target_series=canonical_j_coordinate_series,
        order=order,
        eta_levels=eta_levels,
        moduli=moduli,
        max_abs_exponent=max_abs_exponent,
        followup_label="H_J_f2_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_solve_order=12,
        named_gg_include_weighted_coordinate_diagnostics=False,
        followup_named_gg_benchmark_name="gollnitz_gordon_normalized",
        followup_named_gg_series=gg_series,
        followup_named_gg_solve_order=12,
        followup_named_gg_include_weighted_coordinate_diagnostics=False,
    )

    raw_quotient_series = series_div(p_coordinate_series, g_coordinate_series)
    raw_quotient_squared = series_pow(raw_quotient_series, 2)
    raw_quotient_fourth = series_pow(raw_quotient_squared, 2)
    quotient_coordinate_series = series_div(
        [sp.Integer(16) * value for value in t_squared_series],
        g_squared,
    )
    quotient_coordinate_cubic = series_pow(quotient_coordinate_series, 3)
    quotient_coordinate_bridge_series = [
        sp.simplify(raw_quotient_fourth[index] - raw_quotient_squared[index])
        for index in range(order)
    ]
    quotient_coordinate_mixed_term = series_mul(
        quotient_coordinate_series,
        raw_quotient_squared,
    )
    quotient_coordinate_bridge_series = [
        sp.simplify(
            quotient_coordinate_bridge_series[index]
            - 3 * quotient_coordinate_mixed_term[index]
            - quotient_coordinate_cubic[index]
        )
        for index in range(order)
    ]
    (
        quotient_coordinate_bridge_first_failure_power,
        quotient_coordinate_bridge_first_failure_coeff,
    ) = _first_nonzero_residual_term(quotient_coordinate_bridge_series)
    quotient_coordinate_template_series = series_div(
        [sp.Integer(1)] + [sp.Integer(0) for _ in range(order - 1)],
        series_pow(g_correction_series, 2),
    )
    quotient_coordinate_template_scan = _scan_constant_one_series(
        label="G_X_ws",
        expression=(
            "G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2)"
            " = 1 / G_g12_ws^2"
        ),
        target_series=quotient_coordinate_template_series,
        order=order,
        eta_levels=eta_levels,
        moduli=moduli,
        max_abs_exponent=max_abs_exponent,
        followup_label="H_X_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_solve_order=12,
        named_gg_include_weighted_coordinate_diagnostics=False,
        followup_named_gg_benchmark_name="gollnitz_gordon_normalized",
        followup_named_gg_series=gg_series,
        followup_named_gg_solve_order=12,
        followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
        followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
        followup_named_gg_descendant_focus_degree_values=(1,),
        followup_named_gg_descendant_focus_max_abs_exponent=2,
        followup_named_gg_include_weighted_coordinate_diagnostics=False,
    )
    anchor_canonical_j_coordinate_numerator_series = series_pow(
        [
            sp.simplify((sp.Integer(1) if index == 0 else sp.Integer(0)) + 16 * value)
            for index, value in enumerate(quotient_coordinate_template_series)
        ],
        3,
    )
    anchor_canonical_j_coordinate_series = series_div(
        anchor_canonical_j_coordinate_numerator_series,
        [sp.simplify(4913 * value) for value in series_pow(quotient_coordinate_template_series, 2)],
    )
    anchor_canonical_j_coordinate_scan = _scan_constant_one_series(
        label="J_X_ws",
        expression="J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)",
        target_series=anchor_canonical_j_coordinate_series,
        order=order,
        eta_levels=eta_levels,
        moduli=moduli,
        max_abs_exponent=max_abs_exponent,
        followup_label="H_J_X_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_solve_order=12,
        named_gg_include_weighted_coordinate_diagnostics=False,
        followup_named_gg_benchmark_name="gollnitz_gordon_normalized",
        followup_named_gg_series=gg_series,
        followup_named_gg_solve_order=12,
        followup_named_gg_include_weighted_coordinate_diagnostics=False,
    )

    alternate_anchor_canonical_j_coordinate_label = "J_X15_ws"
    alternate_anchor_canonical_j_coordinate_expression = (
        "J_X15_ws = (16*G_X_ws - 1)^3 / (3375*G_X_ws^2)"
    )
    alternate_anchor_canonical_j_coordinate_reason = (
        "The Weber cubic `j = (X - 16)^3 / X` can also be applied to the inverted "
        "branch `X = 1 / G_X_ws`, and then normalized by the CM-scale `-3375` so the "
        "true-source anchor again yields a constant-1 `j`-side coordinate."
    )
    alternate_anchor_canonical_j_coordinate_bridge_expression = (
        "3375*J_X15_ws*G_X_ws^2 - (16*G_X_ws - 1)^3 = 0"
    )
    alternate_anchor_canonical_j_coordinate_numerator_series = series_pow(
        [
            sp.simplify(
                (sp.Integer(16) * value)
                - (sp.Integer(1) if index == 0 else sp.Integer(0))
            )
            for index, value in enumerate(quotient_coordinate_template_series)
        ],
        3,
    )
    alternate_anchor_canonical_j_coordinate_series = series_div(
        alternate_anchor_canonical_j_coordinate_numerator_series,
        [sp.simplify(3375 * value) for value in series_pow(quotient_coordinate_template_series, 2)],
    )
    alternate_anchor_canonical_j_coordinate_residual = [
        sp.simplify(value - (sp.Integer(1) if index == 0 else sp.Integer(0)))
        for index, value in enumerate(alternate_anchor_canonical_j_coordinate_series)
    ]
    (
        alternate_anchor_canonical_j_coordinate_first_failure_power,
        alternate_anchor_canonical_j_coordinate_first_failure_coeff,
    ) = _first_nonzero_residual_term(alternate_anchor_canonical_j_coordinate_residual)
    canonical_j_anchor_bridge_scan = _scan_constant_one_series_pair_bridge(
        left_label=quotient_coordinate_template_scan.label,
        left_series=quotient_coordinate_template_series,
        right_label=canonical_j_coordinate_scan.label,
        right_series=canonical_j_coordinate_series,
        order=order,
        difference_label="D_XJ_ws",
        quotient_label="Q_XJ_ws",
        quotient_followup_label="K_XJ_ws",
        quotient_followup_bridge_left_label=quotient_coordinate_template_scan.label,
        quotient_followup_bridge_difference_label="D_XKJ_ws",
        quotient_followup_bridge_quotient_label="Q_XKJ_ws",
        quotient_followup_bridge_quotient_followup_label="L_XKJ_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_degree_values=(1, 2),
        named_gg_max_abs_exponent=max_abs_exponent,
        named_gg_solve_order=min(order, 24),
        named_gg_target_pairs=(("Q_XJ_ws", "K_XJ_ws"), ("Q_XKJ_ws", "L_XKJ_ws")),
        quotient_followup_named_gg_series=gg_series,
        quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
        quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
        quotient_followup_named_gg_descendant_focus_degree_values=(1,),
        quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
        nested_quotient_followup_named_gg_series=gg_series,
        nested_quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
        nested_quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
        nested_quotient_followup_named_gg_descendant_focus_degree_values=(1,),
        nested_quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
            named_gg_include_weighted_coordinate_diagnostics=False,
        )
    canonical_j_lift_bridge_scan: ConstantOnePairBridgeLiteScan | None = None
    canonical_j_alt_lift_bridge_scan: ConstantOnePairBridgeLiteScan | None = None
    if max_abs_exponent > 4:
        focused_named_gg_order = min(order, 16)
        focused_named_gg_solve_order = min(focused_named_gg_order, 12)
        lift_bridge_named_gg_supplemental_powers = focused_named_gg_supplemental_powers[:2]
        canonical_j_lift_bridge_scan = _scan_constant_one_series_pair_bridge_lite(
            left_label=anchor_canonical_j_coordinate_scan.label,
            left_series=anchor_canonical_j_coordinate_series,
            right_label=canonical_j_coordinate_scan.label,
            right_series=canonical_j_coordinate_series,
            order=order,
            difference_label="D_JX_ws",
            quotient_label="Q_JX_ws",
            solve_order=min(order, 24),
            eta_levels=eta_levels,
            moduli=moduli,
            named_gg_benchmark_name="gollnitz_gordon_normalized",
            named_gg_series=gg_series,
            named_gg_order=focused_named_gg_order,
            named_gg_degree_values=(1, 2),
            named_gg_max_abs_exponent=max_abs_exponent,
            named_gg_solve_order=focused_named_gg_solve_order,
            named_gg_supplemental_powers=(),
            named_gg_include_weighted_coordinate_diagnostics=False,
            quotient_followup_label="K_JX_ws",
        )
        canonical_j_alt_lift_bridge_scan = _scan_constant_one_series_pair_bridge_lite(
            left_label=alternate_anchor_canonical_j_coordinate_label,
            left_series=alternate_anchor_canonical_j_coordinate_series,
            right_label=canonical_j_coordinate_scan.label,
            right_series=canonical_j_coordinate_series,
            order=order,
            difference_label="D_JX15_ws",
            quotient_label="Q_JX15_ws",
            solve_order=min(order, 24),
            eta_levels=eta_levels,
            moduli=moduli,
            named_gg_benchmark_name="gollnitz_gordon_normalized",
            named_gg_series=gg_series,
            named_gg_order=focused_named_gg_order,
            named_gg_degree_values=(1, 2),
            named_gg_max_abs_exponent=max_abs_exponent,
            named_gg_solve_order=focused_named_gg_solve_order,
            named_gg_supplemental_powers=lift_bridge_named_gg_supplemental_powers,
            named_gg_include_weighted_coordinate_diagnostics=False,
            quotient_followup_label="K_JX15_ws",
        )
    primary_named_gg_modular_equation_scan: GGModularEquationScan | None = None
    if max_abs_exponent > 4:
        primary_named_gg_modular_equation_scan = scan_gg_modular_equation_box(
            target_series=g_correction_series[:focused_named_gg_order],
            benchmark_name="gollnitz_gordon_normalized",
            gg_series=gg_series[:focused_named_gg_order],
            order=focused_named_gg_order,
            degree_values=(1, 2),
            max_abs_exponent=max_abs_exponent,
            solve_order=focused_named_gg_solve_order,
            include_weighted_coordinate_diagnostics=False,
        )
    quotient_series = series_div(p_correction_series, g_correction_series)
    quotient_residual = [sp.simplify(value) for value in quotient_series]
    quotient_residual[0] = sp.simplify(quotient_residual[0] - 1)
    quotient_first_failure_power, quotient_first_failure_coeff = _first_nonzero_residual_term(
        quotient_residual
    )
    quotient_self_polynomial_scan = scan_self_polynomial_uniqueness_relations(
        target_series=quotient_series,
        moduli=moduli,
        order=order,
        fg_degree_values=(1, 2),
        t_degree_values=(1, 2, 3),
    )
    quotient_self_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
        target_series=quotient_series,
        moduli=moduli,
        order=order,
        t_degree_values=(1, 2, 3),
    )
    quotient_self_quotient_product_scans = tuple(
        scan_ratio_self_quotient_product_relations(
            ratio_series=quotient_series,
            moduli=moduli,
            order=order,
            max_abs_exponent=max_abs_exponent,
        )
    )
    normalized_followup: NormalizedResidualFollowupScan | None = None
    if quotient_first_failure_power is not None and quotient_first_failure_coeff is not None:
        followup_series: Series = [sp.Integer(0) for _ in range(order)]
        for index in range(order):
            source_index = index + quotient_first_failure_power
            if source_index >= order:
                break
            followup_series[index] = sp.simplify(
                quotient_residual[source_index] / quotient_first_failure_coeff
            )
        followup_residual = [sp.simplify(value) for value in followup_series]
        followup_residual[0] = sp.simplify(followup_residual[0] - 1)
        followup_first_failure_power, followup_first_failure_coeff = _first_nonzero_residual_term(
            followup_residual
        )
        normalized_followup_named_gg_descendant_preview: GGDescendantPreview | None = None
        normalized_followup_named_gg_descendant_focused_scan: GGDescendantFocusedScan | None = None
        if max_abs_exponent > 4:
            focused_named_gg_order = min(order, 16)
            normalized_followup_named_gg_descendant_preview = _gg_descendant_preview(
                base_series=gg_series[:focused_named_gg_order],
                order=focused_named_gg_order,
                supplemental_powers=focused_named_gg_supplemental_powers,
            )
            normalized_followup_named_gg_descendant_focused_scan = _scan_gg_descendant_focus_box(
                target_series=followup_series,
                gg_series=gg_series,
                order=order,
                supplemental_powers=focused_named_gg_supplemental_powers,
                degree_values=(1,),
                max_abs_exponent=2,
            )
        normalized_followup = NormalizedResidualFollowupScan(
            label="H_gp_ws",
            expression=(
                f"H_gp_ws = (R_gp_ws - 1) / "
                f"({_format_expr(quotient_first_failure_coeff)}*t^{quotient_first_failure_power})"
            ),
            first_failure_power=followup_first_failure_power,
            first_failure_coeff=followup_first_failure_coeff,
            self_polynomial_scan=scan_self_polynomial_uniqueness_relations(
                target_series=followup_series,
                moduli=moduli,
                order=order,
                fg_degree_values=(1, 2),
                t_degree_values=(1, 2, 3),
            ),
            self_fractional_linear_scan=scan_self_fractional_linear_uniqueness_relations(
                target_series=followup_series,
                moduli=moduli,
                order=order,
                t_degree_values=(1, 2, 3),
            ),
            self_quotient_product_scans=tuple(
                scan_ratio_self_quotient_product_relations(
                    ratio_series=followup_series,
                    moduli=moduli,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            eta_scans=tuple(
                scan_ratio_eta_quotient_relations(
                    ratio_series=followup_series,
                    levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            modular_unit_eta_scans=tuple(
                scan_ratio_modular_unit_eta_relations(
                    ratio_series=followup_series,
                    moduli=moduli,
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            self_plus_pochhammer_scans=tuple(
                scan_ratio_self_plus_pochhammer_relations(
                    ratio_series=followup_series,
                    moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            self_plus_pochhammer_eta_scans=tuple(
                scan_ratio_self_plus_pochhammer_eta_relations(
                    ratio_series=followup_series,
                    moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                    eta_levels=eta_levels,
                    order=order,
                    max_abs_exponent=max_abs_exponent,
                )
            ),
            named_gg_descendant_preview=normalized_followup_named_gg_descendant_preview,
            named_gg_descendant_focused_scan=normalized_followup_named_gg_descendant_focused_scan,
        )

    followup_bridge_scan: ConstantOnePairBridgeScan | None = None
    coordinate_followup = quotient_coordinate_template_scan.normalized_followup
    if coordinate_followup is not None and normalized_followup is not None:
        coordinate_followup_series = [sp.Integer(0) for _ in range(order)]
        coordinate_followup_series[0] = sp.Integer(1)
        if (
            quotient_coordinate_template_scan.first_failure_power is not None
            and quotient_coordinate_template_scan.first_failure_coeff is not None
        ):
            coordinate_template_residual = [
                sp.simplify(value - (sp.Integer(1) if index == 0 else sp.Integer(0)))
                for index, value in enumerate(quotient_coordinate_template_series)
            ]
            for index in range(order):
                source_index = index + quotient_coordinate_template_scan.first_failure_power
                if source_index >= order:
                    break
                coordinate_followup_series[index] = sp.simplify(
                    coordinate_template_residual[source_index]
                    / quotient_coordinate_template_scan.first_failure_coeff
                )

        residual_followup_series = [sp.Integer(0) for _ in range(order)]
        residual_followup_series[0] = sp.Integer(1)
        for index in range(order):
            source_index = index + quotient_first_failure_power
            if source_index >= order:
                break
            residual_followup_series[index] = sp.simplify(
                quotient_residual[source_index] / quotient_first_failure_coeff
            )

        followup_bridge_scan = _scan_constant_one_series_pair_bridge(
            left_label=coordinate_followup.label,
            left_series=coordinate_followup_series,
            right_label=normalized_followup.label,
            right_series=residual_followup_series,
            order=order,
            difference_label="D_XR_ws",
            quotient_label="Q_XR_ws",
            quotient_followup_label="K_XR_ws",
            quotient_followup_bridge_left_label=coordinate_followup.label,
            quotient_followup_bridge_difference_label="D_XK_ws",
            quotient_followup_bridge_quotient_label="Q_XK_ws",
            quotient_followup_bridge_quotient_followup_label="L_XK_ws",
            named_gg_benchmark_name="gollnitz_gordon_normalized",
            named_gg_series=gg_series,
            named_gg_degree_values=(1, 2),
            named_gg_max_abs_exponent=max_abs_exponent,
            named_gg_solve_order=min(order, 24),
            named_gg_target_pairs=(("Q_XR_ws", "K_XR_ws"), ("Q_XK_ws", "L_XK_ws")),
            quotient_followup_named_gg_series=gg_series,
            quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
            quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
            quotient_followup_named_gg_descendant_focus_degree_values=(1,),
            quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
            nested_quotient_followup_named_gg_series=gg_series,
            nested_quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
            nested_quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
            nested_quotient_followup_named_gg_descendant_focus_degree_values=(1,),
            nested_quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
            named_gg_include_weighted_coordinate_diagnostics=False,
        )

    return WeberResidualBridgeScan(
        primary_label="G_g12_ws",
        primary_expression="G_g12_ws = g12_ws / (t^2; t^4)_inf^12",
        primary_reason=(
            "The eta-side template `(t^2; t^4)_inf^12` is the simpler named source anchor, "
            "so `G_g12_ws` is the current primary residual and `G_p12_ws` is treated as its "
            "algebraically constrained companion."
        ),
        companion_label="G_p12_ws",
        companion_expression="G_p12_ws = p12_ws / (-t^2; t^4)_inf^12",
        exact_bridge_expression=(
            "g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0"
        ),
        exact_bridge_holds=bridge_first_failure_power is None,
        exact_bridge_first_failure_power=bridge_first_failure_power,
        exact_bridge_first_failure_coeff=bridge_first_failure_coeff,
        residual_bridge_expression=(
            "(t^2; t^4)_inf^48*(-t^2; t^4)_inf^24*G_g12_ws^4*G_p12_ws^2"
            " - (t^2; t^4)_inf^24*(-t^2; t^4)_inf^48*G_g12_ws^2*G_p12_ws^4"
            " + 48*t^2*(t^2; t^4)_inf^24*(-t^2; t^4)_inf^24*G_g12_ws^2*G_p12_ws^2"
            " + 4096*t^6 = 0"
        ),
        classical_product_coordinate_label="G_f2_ws",
        classical_product_coordinate_expression=(
            "G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2)"
            " = G_g12_ws*G_p12_ws"
        ),
        classical_product_coordinate_bridge_expression="G_f2_ws - G_g12_ws*G_p12_ws = 0",
        classical_product_coordinate_scan=classical_product_coordinate_scan,
        canonical_j_coordinate_label="J_f2_ws",
        canonical_j_coordinate_expression="J_f2_ws = (16 - G_f2_ws)^3 / (3375*G_f2_ws)",
        canonical_j_coordinate_reason=(
            "Yui--Zagier's Weber cubic `(X - 16)^3 = X*j` turns the "
            "classical Weber `f2` coordinate into a canonical `j`-side "
            "constant-1 object once the true-source normalization makes `G_f2_ws = 1`."
        ),
        canonical_j_coordinate_bridge_expression="3375*J_f2_ws*G_f2_ws - (16 - G_f2_ws)^3 = 0",
        canonical_j_coordinate_scan=canonical_j_coordinate_scan,
        canonical_j_anchor_bridge_scan=canonical_j_anchor_bridge_scan,
        canonical_j_lift_bridge_scan=canonical_j_lift_bridge_scan,
        alternate_anchor_canonical_j_coordinate_label=alternate_anchor_canonical_j_coordinate_label,
        alternate_anchor_canonical_j_coordinate_expression=alternate_anchor_canonical_j_coordinate_expression,
        alternate_anchor_canonical_j_coordinate_reason=alternate_anchor_canonical_j_coordinate_reason,
        alternate_anchor_canonical_j_coordinate_bridge_expression=alternate_anchor_canonical_j_coordinate_bridge_expression,
        alternate_anchor_canonical_j_coordinate_first_failure_power=alternate_anchor_canonical_j_coordinate_first_failure_power,
        alternate_anchor_canonical_j_coordinate_first_failure_coeff=alternate_anchor_canonical_j_coordinate_first_failure_coeff,
        canonical_j_alt_lift_bridge_scan=canonical_j_alt_lift_bridge_scan,
        quotient_coordinate_label="X_g_ws",
        quotient_coordinate_expression="X_g_ws = 16*t^2 / g12_ws^2",
        quotient_coordinate_bridge_expression=(
            "Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0, "
            "Q_gp_ws = p12_ws / g12_ws"
        ),
        quotient_coordinate_bridge_holds=quotient_coordinate_bridge_first_failure_power is None,
        quotient_coordinate_bridge_first_failure_power=quotient_coordinate_bridge_first_failure_power,
        quotient_coordinate_bridge_first_failure_coeff=quotient_coordinate_bridge_first_failure_coeff,
        quotient_coordinate_template_bridge_expression="G_X_ws*G_g12_ws^2 - 1 = 0",
        quotient_coordinate_template_scan=quotient_coordinate_template_scan,
        anchor_canonical_j_coordinate_label="J_X_ws",
        anchor_canonical_j_coordinate_expression="J_X_ws = (1 + 16*G_X_ws)^3 / (4913*G_X_ws^2)",
        anchor_canonical_j_coordinate_reason=(
            "The same signed Weber cubic also turns the squared `g`-side branch "
            "`1 / G_X_ws = G_g12_ws^2` into a second canonical `j`-side constant-1 "
            "object, so the template-normalized anchor now has its own direct `j` lift."
        ),
        anchor_canonical_j_coordinate_bridge_expression=(
            "4913*J_X_ws*G_X_ws^2 - (1 + 16*G_X_ws)^3 = 0"
        ),
        anchor_canonical_j_coordinate_scan=anchor_canonical_j_coordinate_scan,
        primary_named_gg_modular_equation_scan=primary_named_gg_modular_equation_scan,
        quotient_coordinate_template_named_gg_modular_equation_scan=quotient_coordinate_template_scan.named_gg_modular_equation_scan,
        quotient_label="R_gp_ws",
        quotient_expression="R_gp_ws = G_p12_ws / G_g12_ws",
        quotient_first_failure_power=quotient_first_failure_power,
        quotient_first_failure_coeff=quotient_first_failure_coeff,
        quotient_self_polynomial_scan=quotient_self_polynomial_scan,
        quotient_self_fractional_linear_scan=quotient_self_fractional_linear_scan,
        quotient_self_quotient_product_scans=quotient_self_quotient_product_scans,
        quotient_eta_scans=tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=quotient_series,
                levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        quotient_modular_unit_eta_scans=tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=quotient_series,
                moduli=moduli,
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        quotient_self_plus_pochhammer_scans=tuple(
            scan_ratio_self_plus_pochhammer_relations(
                ratio_series=quotient_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        quotient_self_plus_pochhammer_eta_scans=tuple(
            scan_ratio_self_plus_pochhammer_eta_relations(
                ratio_series=quotient_series,
                moduli=tuple(modulus for modulus in moduli if modulus <= 4),
                eta_levels=eta_levels,
                order=order,
                max_abs_exponent=max_abs_exponent,
            )
        ),
        normalized_followup=normalized_followup,
        followup_bridge_scan=followup_bridge_scan,
    )


def _build_weber_j_pb_bridge_scan(
    *,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    max_abs_exponent: int = 8,
) -> ConstantOnePairBridgeScan | None:
    p_coordinate_series = _weber_schlafli_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if p_coordinate_series is None:
        return None
    b_coordinate_series = _weber_companion_coordinate_series(
        p_series=p_coordinate_series,
        order=order,
    )
    if b_coordinate_series is None:
        return None
    p_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
        source_label="P_ws",
        target_series=p_coordinate_series,
        order=order,
        max_abs_exponent=min(max_abs_exponent, 4),
        eta_levels=eta_levels,
        moduli=moduli,
    )
    b_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
        source_label="B_ws",
        target_series=b_coordinate_series,
        order=order,
        max_abs_exponent=min(max_abs_exponent, 4),
        eta_levels=eta_levels,
        moduli=moduli,
    )
    if p_leading_normalized_scan is None or b_leading_normalized_scan is None:
        return None

    g_coordinate_series = _normalized_weber_g_coordinate_series(
        target_series=target_series,
        order=order,
    )
    p_named_coordinate_series = _normalized_weber_p_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if g_coordinate_series is None or p_named_coordinate_series is None:
        return None
    g_template_series = _normalized_weber_g_template_series(order=order)
    p_template_series = _normalized_weber_p_template_series(order=order)
    g_correction_series = series_div(g_coordinate_series, g_template_series)
    p_correction_series = series_div(p_named_coordinate_series, p_template_series)
    classical_product_coordinate_series = series_mul(g_correction_series, p_correction_series)
    canonical_j_numerator_series = series_pow(
        [
            sp.simplify((sp.Integer(16) if index == 0 else sp.Integer(0)) - value)
            for index, value in enumerate(classical_product_coordinate_series)
        ],
        3,
    )
    canonical_j_coordinate_series = series_div(
        canonical_j_numerator_series,
        [sp.simplify(3375 * value) for value in classical_product_coordinate_series],
    )
    q_pb_series = series_div(
        list(b_leading_normalized_scan.normalized_series),
        list(p_leading_normalized_scan.normalized_series),
    )
    gg_series = _canonical_benchmark_series(
        "gollnitz_gordon_normalized",
        depth=order,
        order=order,
    )
    focused_named_gg_supplemental_powers = _supplemental_source_family_powers(
        smoke=max_abs_exponent <= 4
    ).get("GG", ())
    return _scan_constant_one_series_pair_bridge(
        left_label="J_f2_ws",
        left_series=canonical_j_coordinate_series,
        right_label="Q_PB_ws",
        right_series=q_pb_series,
        order=order,
        difference_label="D_JPB_ws",
        quotient_label="Q_JPB_ws",
        quotient_followup_label="K_JPB_ws",
        quotient_followup_bridge_left_label="J_f2_ws",
        quotient_followup_bridge_difference_label="D_JKPB_ws",
        quotient_followup_bridge_quotient_label="Q_JKPB_ws",
        quotient_followup_bridge_quotient_followup_label="L_JKPB_ws",
        named_gg_benchmark_name="gollnitz_gordon_normalized",
        named_gg_series=gg_series,
        named_gg_degree_values=(1, 2),
        named_gg_max_abs_exponent=max_abs_exponent,
        named_gg_solve_order=min(order, 24),
        named_gg_target_pairs=(("Q_JPB_ws", "K_JPB_ws"), ("Q_JKPB_ws", "L_JKPB_ws")),
        quotient_followup_named_gg_series=gg_series,
        quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
        quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
        quotient_followup_named_gg_descendant_focus_degree_values=(1,),
        quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
        nested_quotient_followup_named_gg_series=gg_series,
        nested_quotient_followup_named_gg_descendant_preview_powers=focused_named_gg_supplemental_powers,
        nested_quotient_followup_named_gg_descendant_focus_powers=focused_named_gg_supplemental_powers,
        nested_quotient_followup_named_gg_descendant_focus_degree_values=(1,),
        nested_quotient_followup_named_gg_descendant_focus_max_abs_exponent=2,
        named_gg_include_weighted_coordinate_diagnostics=False,
    )


def _build_weber_j_lift_pivot_bridge_scans(
    *,
    target_series: Series,
    order: int,
    eta_levels: tuple[int, ...] = (1, 2, 4),
    moduli: tuple[int, ...] = (2, 3, 4),
    max_abs_exponent: int = 8,
) -> tuple[ConstantOnePairBridgeLiteScan, ...]:
    """Focused cross-rail comparisons using the `Q_JX15_ws` lift-bridge quotient as a pivot.

    Motivation (analogy): `Q_JX15_ws` is a tighter “lens” between two `j`-side lifts; this
    function checks whether that lens is just a thin coating on top of other already-tracked
    Weber seams (the `P/B` seam or the nested anchor seam) without widening anonymous boxes.
    """

    if order < 2:
        return ()
    if max_abs_exponent <= 4:
        return ()
    if len(target_series) < order:
        raise ValueError("target_series is shorter than requested order")

    g_coordinate_series = _normalized_weber_g_coordinate_series(
        target_series=target_series,
        order=order,
    )
    p_coordinate_series = _normalized_weber_p_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if g_coordinate_series is None or p_coordinate_series is None:
        return ()

    g_template_series = _normalized_weber_g_template_series(order=order)
    p_template_series = _normalized_weber_p_template_series(order=order)
    g_correction_series = series_div(g_coordinate_series, g_template_series)
    p_correction_series = series_div(p_coordinate_series, p_template_series)
    classical_product_coordinate_series = series_mul(g_correction_series, p_correction_series)

    canonical_j_numerator_series = series_pow(
        [
            sp.simplify((sp.Integer(16) if index == 0 else sp.Integer(0)) - value)
            for index, value in enumerate(classical_product_coordinate_series)
        ],
        3,
    )
    canonical_j_coordinate_series = series_div(
        canonical_j_numerator_series,
        [sp.simplify(3375 * value) for value in classical_product_coordinate_series],
    )

    quotient_coordinate_template_series = series_div(
        [sp.Integer(1)] + [sp.Integer(0) for _ in range(order - 1)],
        series_pow(g_correction_series, 2),
    )

    alternate_anchor_canonical_j_coordinate_series = series_div(
        series_pow(
            [
                sp.simplify(
                    (sp.Integer(16) * value)
                    - (sp.Integer(1) if index == 0 else sp.Integer(0))
                )
                for index, value in enumerate(quotient_coordinate_template_series)
            ],
            3,
        ),
        [sp.simplify(3375 * value) for value in series_pow(quotient_coordinate_template_series, 2)],
    )

    q_jx15_series = series_div(
        canonical_j_coordinate_series,
        alternate_anchor_canonical_j_coordinate_series,
    )
    _, q_jx15_first_failure_power, q_jx15_first_failure_coeff = _constant_one_residual_series(
        target_series=q_jx15_series,
        order=order,
    )
    if q_jx15_first_failure_power is None or q_jx15_first_failure_coeff is None:
        return ()

    gg_series = list(
        _canonical_benchmark_series(
            "gollnitz_gordon_normalized",
            depth=order,
            order=order,
        )
    )
    focused_named_gg_order = min(order, 16)
    focused_named_gg_solve_order = min(focused_named_gg_order, 12)

    scans: list[ConstantOnePairBridgeLiteScan] = []

    # (Step 3) Tie back to the Weber-Schlafli coordinate seam via the leading-normalized `Q_PB_ws`.
    p_ws_coordinate_series = _weber_schlafli_coordinate_series(
        target_series=target_series,
        order=order,
    )
    if p_ws_coordinate_series is not None:
        b_ws_coordinate_series = _weber_companion_coordinate_series(
            p_series=p_ws_coordinate_series,
            order=order,
        )
        if b_ws_coordinate_series is not None:
            p_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
                source_label="P_ws",
                target_series=p_ws_coordinate_series,
                order=order,
                max_abs_exponent=min(max_abs_exponent, 4),
                eta_levels=eta_levels,
                moduli=moduli,
            )
            b_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
                source_label="B_ws",
                target_series=b_ws_coordinate_series,
                order=order,
                max_abs_exponent=min(max_abs_exponent, 4),
                eta_levels=eta_levels,
                moduli=moduli,
            )
            if p_leading_normalized_scan is not None and b_leading_normalized_scan is not None:
                q_pb_series = series_div(
                    list(b_leading_normalized_scan.normalized_series),
                    list(p_leading_normalized_scan.normalized_series),
                )
                scans.append(
                    _scan_constant_one_series_pair_bridge_lite(
                        left_label="J_X15_ws",
                        left_series=alternate_anchor_canonical_j_coordinate_series,
                        right_label="Q_PB_ws",
                        right_series=q_pb_series,
                        order=order,
                        difference_label="D_JX15PB_ws",
                        quotient_label="Q_JX15PB_ws",
                        solve_order=min(order, 24),
                        eta_levels=eta_levels,
                        moduli=moduli,
                        named_gg_benchmark_name="gollnitz_gordon_normalized",
                        named_gg_series=gg_series,
                        named_gg_order=focused_named_gg_order,
                        named_gg_degree_values=(1, 2),
                        named_gg_max_abs_exponent=max_abs_exponent,
                        named_gg_solve_order=focused_named_gg_solve_order,
                        named_gg_supplemental_powers=(),
                        named_gg_include_weighted_coordinate_diagnostics=False,
                        quotient_followup_label="K_JX15PB_ws",
                    )
                )

                q_jpb_series = series_div(
                    q_pb_series,
                    canonical_j_coordinate_series,
                )
                scans.append(
                    _scan_constant_one_series_pair_bridge_lite(
                        left_label="Q_JX15_ws",
                        left_series=q_jx15_series,
                        right_label="Q_JPB_ws",
                        right_series=q_jpb_series,
                        order=order,
                        difference_label="D_JX15JPB_ws",
                        quotient_label="Q_JX15JPB_ws",
                        solve_order=min(order, 24),
                        eta_levels=eta_levels,
                        moduli=moduli,
                        named_gg_benchmark_name="gollnitz_gordon_normalized",
                        named_gg_series=gg_series,
                        named_gg_order=focused_named_gg_order,
                        named_gg_degree_values=(1, 2),
                        named_gg_max_abs_exponent=max_abs_exponent,
                        named_gg_solve_order=focused_named_gg_solve_order,
                        named_gg_supplemental_powers=(),
                        named_gg_include_weighted_coordinate_diagnostics=False,
                        quotient_followup_label="K_JX15JPB_ws",
                    )
                )

    # (Step 1) Compare against the nested anchor seam `Q_XKJ_ws` (defined by `K_XJ_ws / G_X_ws`).
    q_xj_series = series_div(
        canonical_j_coordinate_series,
        quotient_coordinate_template_series,
    )
    k_xj_series = _normalized_constant_one_followup_series(
        target_series=q_xj_series,
        order=order,
    )
    if k_xj_series is not None:
        q_xkj_series = series_div(
            k_xj_series,
            quotient_coordinate_template_series,
        )
        scans.append(
            _scan_constant_one_series_pair_bridge_lite(
                left_label="Q_JX15_ws",
                left_series=q_jx15_series,
                right_label="Q_XKJ_ws",
                right_series=q_xkj_series,
                order=order,
                difference_label="D_JX15XKJ_ws",
                quotient_label="Q_JX15XKJ_ws",
                solve_order=min(order, 24),
                eta_levels=eta_levels,
                moduli=moduli,
                named_gg_benchmark_name="gollnitz_gordon_normalized",
                named_gg_series=gg_series,
                named_gg_order=focused_named_gg_order,
                named_gg_degree_values=(1, 2),
                named_gg_max_abs_exponent=max_abs_exponent,
                named_gg_solve_order=focused_named_gg_solve_order,
                named_gg_supplemental_powers=(),
                named_gg_include_weighted_coordinate_diagnostics=False,
                quotient_followup_label="K_JX15XKJ_ws",
            )
        )

    return tuple(scans)


def scan_morton_periodic_point_box(
    *,
    target_series: Series,
    order: int,
    max_abs_exponent: int = 8,
) -> MortonPeriodicPointScan:
    if order < 2:
        return MortonPeriodicPointScan(template_results=())

    target_trunc = target_series[:order]
    target_q2 = benchmark_power_substitution_series(target_trunc, power=2, order=order)
    target_sq = series_pow(target_trunc, 2)
    target_q2_sq = series_pow(target_q2, 2)
    one_series = [sp.Integer(0) for _ in range(order)]
    one_series[0] = sp.Integer(1)
    neg_target_q2 = [sp.simplify(-value) for value in target_q2]
    positive_transform = series_div(
        series_add(one_series, neg_target_q2),
        series_add(one_series, target_q2),
    )
    negative_transform = series_div(
        series_add(target_q2, [sp.Integer(-1)] + [sp.Integer(0)] * (order - 1)),
        series_add(target_q2, one_series),
    )
    relation_map = dict(_morton_periodic_point_relations(order))
    template_specs = (
        (
            "Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`",
            relation_map["Morton Prop. 8.1 self-substitution template `f(Y, Y_2)`"],
            {"F": target_trunc, "G": target_q2},
        ),
        (
            "Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`",
            relation_map["Morton Prop. 8.2 squared template `g(Y^2, Y_2^2)`"],
            {"F": target_sq, "G": target_q2_sq},
        ),
        (
            "Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`",
            relation_map["Morton Cor. 7.4 transformed template `f(Y, (1-Y_2)/(1+Y_2))`"],
            {"F": target_trunc, "G": positive_transform},
        ),
        (
            "Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`",
            relation_map["Morton Cor. 7.5 transformed template `f(Y, (Y_2-1)/(Y_2+1))`"],
            {"F": target_trunc, "G": negative_transform},
        ),
    )
    template_results: list[MortonPeriodicPointTemplateResult] = []
    for label, relation, series_by_variable in template_specs:
        residual = _relation_residual_series(
            relation,
            series_by_variable=series_by_variable,
            order=order,
        )
        power, coeff = _first_nonzero_residual_term(residual)
        template_results.append(
            MortonPeriodicPointTemplateResult(
                label=label,
                first_failure_power=power,
                first_failure_coeff=coeff,
                hit=all(sp.simplify(value) == 0 for value in residual),
            )
        )
    named_coordinate_scans: list[MortonNamedCoordinateScan] = []
    leading_normalized_bridge_scans: list[ConstantOnePairBridgeScan] = []
    focused_weber_bridge_reference = (
        _focused_weber_bridge_reference_ladder(order=min(order, 16))
        if max_abs_exponent > 4
        else None
    )
    squared_coordinate_template_results: list[MortonPeriodicPointTemplateResult] = []
    for label, relation in _morton_squared_coordinate_relations(order):
        residual = _relation_residual_series(
            relation,
            series_by_variable={"X": target_sq, "X_2": target_q2_sq},
            order=order,
        )
        power, coeff = _first_nonzero_residual_term(residual)
        squared_coordinate_template_results.append(
            MortonPeriodicPointTemplateResult(
                label=label,
                first_failure_power=power,
                first_failure_coeff=coeff,
                hit=all(sp.simplify(value) == 0 for value in residual),
            )
        )
    named_coordinate_scans.append(
        MortonNamedCoordinateScan(
            family_label="squared",
            label="X_mt",
            expression="X_mt = F^2",
            template_results=tuple(squared_coordinate_template_results),
        )
    )
    transformed_squared_coordinate_series = _morton_transformed_squared_coordinate_series(
        squared_coordinate_series=target_sq,
        order=order,
    )
    if transformed_squared_coordinate_series is not None:
        transformed_squared_coordinate_q2 = benchmark_power_substitution_series(
            transformed_squared_coordinate_series,
            power=2,
            order=order,
        )
        transformed_squared_template_results: list[MortonPeriodicPointTemplateResult] = []
        for label, relation in _morton_transformed_squared_coordinate_relations(order):
            residual = _relation_residual_series(
                relation,
                series_by_variable={"T": transformed_squared_coordinate_series, "T_2": transformed_squared_coordinate_q2},
                order=order,
            )
            power, coeff = _first_nonzero_residual_term(residual)
            transformed_squared_template_results.append(
                MortonPeriodicPointTemplateResult(
                    label=label,
                    first_failure_power=power,
                    first_failure_coeff=coeff,
                    hit=all(sp.simplify(value) == 0 for value in residual),
                )
            )
        named_coordinate_scans.append(
            MortonNamedCoordinateScan(
                family_label="transformed squared",
                label="T_mt",
                expression="T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1), sigma = -1 + sqrt(2)",
                template_results=tuple(transformed_squared_template_results),
            )
        )
    weber_coordinate_series = _weber_schlafli_coordinate_series(
        target_series=target_trunc,
        order=order,
    )
    if weber_coordinate_series is not None:
        weber_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
            source_label="P_ws",
            target_series=weber_coordinate_series,
            order=order,
            max_abs_exponent=min(max_abs_exponent, 4),
        )
        weber_coordinate_q2 = benchmark_power_substitution_series(
            weber_coordinate_series,
            power=2,
            order=order,
        )
        weber_template_results: list[MortonPeriodicPointTemplateResult] = []
        for label, relation in _morton_weber_schlafli_relations(order):
            residual = _relation_residual_series(
                relation,
                series_by_variable={"P": weber_coordinate_series, "P_2": weber_coordinate_q2},
                order=order,
            )
            power, coeff = _first_nonzero_residual_term(residual)
            weber_template_results.append(
                MortonPeriodicPointTemplateResult(
                    label=label,
                    first_failure_power=power,
                    first_failure_coeff=coeff,
                    hit=all(sp.simplify(value) == 0 for value in residual),
                )
            )
        named_coordinate_scans.append(
            MortonNamedCoordinateScan(
                family_label="Weber-Schlafli",
                label="P_ws",
                expression="P_ws = (1/F - F) / 2",
                template_results=tuple(weber_template_results),
                leading_normalized_scan=weber_leading_normalized_scan,
            )
        )
        weber_companion_series = _weber_companion_coordinate_series(
            p_series=weber_coordinate_series,
            order=order,
        )
        if weber_companion_series is not None:
            weber_companion_leading_normalized_scan = _build_leading_normalized_coordinate_scan(
                source_label="B_ws",
                target_series=weber_companion_series,
                order=order,
                max_abs_exponent=min(max_abs_exponent, 4),
            )
            weber_companion_q2 = benchmark_power_substitution_series(
                weber_companion_series,
                power=2,
                order=order,
            )
            companion_template_results: list[MortonPeriodicPointTemplateResult] = []
            for label, relation in _morton_weber_companion_relations(order):
                if relation.variables == ("B", "B_2"):
                    variable_map = {
                        "B": weber_companion_series,
                        "B_2": weber_companion_q2,
                    }
                elif relation.variables == ("B_2", "P"):
                    variable_map = {
                        "B_2": weber_companion_q2,
                        "P": weber_coordinate_series,
                    }
                else:
                    raise ValueError(f"unexpected Weber companion variables: {relation.variables}")
                residual = _relation_residual_series(
                    relation,
                    series_by_variable=variable_map,
                    order=order,
                )
                power, coeff = _first_nonzero_residual_term(residual)
                companion_template_results.append(
                    MortonPeriodicPointTemplateResult(
                        label=label,
                        first_failure_power=power,
                        first_failure_coeff=coeff,
                        hit=all(sp.simplify(value) == 0 for value in residual),
                    )
                )
            named_coordinate_scans.append(
                MortonNamedCoordinateScan(
                    family_label="Weber-Schlafli",
                    label="B_ws",
                    expression="B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)",
                    template_results=tuple(companion_template_results),
                    leading_normalized_scan=weber_companion_leading_normalized_scan,
                )
            )
            if (
                weber_leading_normalized_scan is not None
                and weber_companion_leading_normalized_scan is not None
            ):
                bridge_scan = _scan_constant_one_series_pair_bridge(
                        left_label=weber_leading_normalized_scan.label,
                        left_series=list(weber_leading_normalized_scan.normalized_series),
                        right_label=weber_companion_leading_normalized_scan.label,
                        right_series=list(weber_companion_leading_normalized_scan.normalized_series),
                        order=order,
                        difference_label="D_PB_ws",
                        quotient_label="Q_PB_ws",
                        quotient_followup_label="K_PB_ws",
                        quotient_followup_bridge_left_label=weber_leading_normalized_scan.label,
                        quotient_followup_bridge_difference_label="D_PK_ws",
                        quotient_followup_bridge_quotient_label="Q_PK_ws",
                        quotient_followup_bridge_quotient_followup_label="L_PK_ws",
                    )
                q_pb_series = series_div(
                    list(weber_companion_leading_normalized_scan.normalized_series),
                    list(weber_leading_normalized_scan.normalized_series),
                )
                k_pb_series = _normalized_constant_one_followup_series(
                    target_series=q_pb_series,
                    order=order,
                )
                if (
                    focused_weber_bridge_reference is not None
                ):
                    focused_order = min(order, 16)
                    bridge_scan = replace(
                        bridge_scan,
                        quotient_named_coordinate_orbit_scan=scan_named_coordinate_orbit_box(
                            target_series=q_pb_series[:focused_order],
                            family_label="Weber bridge",
                            base_label="Q_PB_ref_ws",
                            base_expression="Q_PB_ref_ws = N_B_ref_ws / N_P_ref_ws",
                            base_series=list(focused_weber_bridge_reference.q_pb_series[:focused_order]),
                            order=focused_order,
                            powers=(2, 3, 4),
                            degree_values=(1, 2),
                            max_abs_exponent=2,
                            solve_order=min(focused_order, 12),
                        ),
                        quotient_followup_named_coordinate_orbit_scan=(
                            scan_named_coordinate_orbit_box(
                                target_series=k_pb_series[:focused_order],
                                family_label="Weber bridge",
                                base_label="K_PB_ref_ws",
                                base_expression="K_PB_ref_ws = (Q_PB_ref_ws - 1) / (leading term)",
                                base_series=list(focused_weber_bridge_reference.k_pb_series[:focused_order]),
                                order=focused_order,
                                powers=(2, 3, 4),
                                degree_values=(1, 2),
                                max_abs_exponent=2,
                                solve_order=min(focused_order, 12),
                            )
                            if k_pb_series is not None
                            else None
                        ),
                    )
                if bridge_scan.quotient_followup_bridge_scan is not None and k_pb_series is not None:
                    q_pk_series = series_div(
                        k_pb_series,
                        list(weber_leading_normalized_scan.normalized_series),
                    )
                    l_pk_series = _normalized_constant_one_followup_series(
                        target_series=q_pk_series,
                        order=order,
                    )
                    if focused_weber_bridge_reference is not None:
                        focused_order = min(order, 16)
                        bridge_scan = replace(
                            bridge_scan,
                            quotient_followup_bridge_scan=replace(
                                bridge_scan.quotient_followup_bridge_scan,
                                quotient_named_coordinate_orbit_scan=scan_named_coordinate_orbit_box(
                                    target_series=q_pk_series[:focused_order],
                                    family_label="Weber bridge",
                                    base_label="Q_PK_ref_ws",
                                    base_expression="Q_PK_ref_ws = K_PB_ref_ws / N_P_ref_ws",
                                    base_series=list(focused_weber_bridge_reference.q_pk_series[:focused_order]),
                                    order=focused_order,
                                    powers=(2, 3, 4),
                                    degree_values=(1, 2),
                                    max_abs_exponent=2,
                                    solve_order=min(focused_order, 12),
                                ),
                                quotient_followup_named_coordinate_orbit_scan=(
                                    scan_named_coordinate_orbit_box(
                                        target_series=l_pk_series[:focused_order],
                                        family_label="Weber bridge",
                                        base_label="L_PK_ref_ws",
                                        base_expression="L_PK_ref_ws = (Q_PK_ref_ws - 1) / (leading term)",
                                        base_series=list(focused_weber_bridge_reference.l_pk_series[:focused_order]),
                                        order=focused_order,
                                        powers=(2, 3, 4),
                                        degree_values=(1, 2),
                                        max_abs_exponent=2,
                                        solve_order=min(focused_order, 12),
                                    )
                                    if l_pk_series is not None
                                    else None
                                ),
                            ),
                        )
                leading_normalized_bridge_scans.append(bridge_scan)
    return MortonPeriodicPointScan(
        template_results=tuple(template_results),
        named_coordinate_scans=tuple(named_coordinate_scans),
        leading_normalized_bridge_scans=tuple(leading_normalized_bridge_scans),
    )


def _gg_exact_modular_template_hits(
    *,
    target_series: Series,
    ordered_basis_entries: tuple[tuple[str, str, Series], ...],
    quotient_basis_entries: tuple[tuple[str, str, Series], ...],
    order: int,
) -> tuple[
    tuple[str, ...],
    tuple[str, ...],
    tuple[tuple[str, int | None, sp.Expr | None], ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[tuple[str, int | None, sp.Expr | None], ...],
]:
    basis_by_label = {label: series for label, _, series in ordered_basis_entries}
    quotient_by_label = {label: series for label, _, series in quotient_basis_entries}
    relation_map = dict(_gg_exact_modular_relations(order))
    variable_series = [sp.Integer(0) for _ in range(order)]
    if order > 1:
        variable_series[1] = sp.Integer(1)

    exact_labels: list[str] = []
    exact_hits: list[str] = []
    exact_obstructions: list[tuple[str, int | None, sp.Expr | None]] = []
    if "GG3" in basis_by_label:
        label = "Chan--Huang Cor. 3.2(i) on (F, GG3)"
        exact_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(i) on (F, G3)"],
            series_by_variable={"F": target_series, "G": basis_by_label["GG3"], "T": variable_series},
            order=order,
        )
        exact_obstructions.append((label, *_first_nonzero_residual_term(residual)))
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
        exact_obstructions.append((label, *_first_nonzero_residual_term(residual)))
        if all(sp.simplify(value) == 0 for value in residual):
            exact_hits.append(label)

    quotient_labels: list[str] = []
    quotient_hits: list[str] = []
    quotient_obstructions: list[tuple[str, int | None, sp.Expr | None]] = []
    if "Q_3" in quotient_by_label:
        label = "Chan--Huang Cor. 3.2(i) on (F, Q_3)"
        quotient_labels.append(label)
        residual = _relation_residual_series(
            relation_map["Chan--Huang Cor. 3.2(i) on (F, Q_3)"],
            series_by_variable={"F": target_series, "Q": quotient_by_label["Q_3"], "T": variable_series},
            order=order,
        )
        quotient_obstructions.append((label, *_first_nonzero_residual_term(residual)))
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
        quotient_obstructions.append((label, *_first_nonzero_residual_term(residual)))
        if all(sp.simplify(value) == 0 for value in residual):
            quotient_hits.append(label)

    return (
        tuple(exact_labels),
        tuple(exact_hits),
        tuple(exact_obstructions),
        tuple(quotient_labels),
        tuple(quotient_hits),
        tuple(quotient_obstructions),
    )


def _gg_weighted_coordinate_diagnostics(
    *,
    target_series: Series,
    benchmark_name: str,
    gg_series: Series,
    quotient_basis_entries: tuple[tuple[str, str, Series], ...],
    order: int,
    coordinate_degree_values: tuple[int, ...] = (1, 2),
    coordinate_max_abs_exponent: int = 4,
    coordinate_solve_order: int | None = None,
    correction_eta_levels: tuple[int, ...] = (1, 2, 3, 4),
    correction_moduli: tuple[int, ...] = (2, 3, 4),
    correction_max_abs_exponent: int = 4,
    correction_source_families: tuple[tuple[str, str, Series], ...] = (),
    correction_source_powers: tuple[int, ...] = (2, 3, 4),
    correction_source_supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
) -> tuple[GGWeightedCoordinateDiagnostic, ...]:
    quotient_by_label = {label: series for label, _, series in quotient_basis_entries}
    q3_series = quotient_by_label.get("Q_3")
    q4_series = quotient_by_label.get("Q_4")
    if q3_series is None or q4_series is None:
        return ()

    quotient_ordered_basis_series = tuple(
        (label, series[:order]) for label, _, series in quotient_basis_entries
    )
    mixed_quotient_ordered_basis_series = (("GG", gg_series[:order]),) + quotient_ordered_basis_series

    weighted_series = series_div(series_pow(q3_series, 3), series_pow(q4_series, 2))
    weighted_log_series = [
        sp.simplify(3 * coeff_q3 - 2 * coeff_q4)
        for coeff_q3, coeff_q4 in zip(_series_log_coeffs(q3_series), _series_log_coeffs(q4_series))
    ]
    direct_difference = [
        sp.simplify(target_series[index] - weighted_series[index])
        for index in range(order)
    ]
    target_log_series = _series_log_coeffs(target_series)
    log_difference = [
        sp.simplify(target_log_series[index] - weighted_log_series[index])
        for index in range(order)
    ]
    weighted_correction_series = series_div(target_series, weighted_series)
    weighted_correction_gap = [
        sp.simplify(
            weighted_correction_series[index] - (sp.Integer(1) if index == 0 else sp.Integer(0))
        )
        for index in range(order)
    ]

    solve_order = min(order, 12)
    variable_series = [sp.Integer(0) for _ in range(solve_order)]
    if solve_order > 1:
        variable_series[1] = sp.Integer(1)

    polynomial_degree1_relation = None
    if solve_order >= 4:
        polynomial_degree1_relation = search_polynomial_relation(
            series_by_variable={
                "F": target_series[:solve_order],
                "W": weighted_series[:solve_order],
                "T": variable_series,
            },
            order=solve_order,
            max_total_degree=1,
            required_variable="F",
        )

    polynomial_degree2_relation = None
    if solve_order >= 10:
        polynomial_degree2_relation = search_polynomial_relation(
            series_by_variable={
                "F": target_series[:solve_order],
                "W": weighted_series[:solve_order],
                "T": variable_series,
            },
            order=solve_order,
            max_total_degree=2,
            required_variable="F",
        )

    fractional_linear_relation = None
    if solve_order >= 4:
        fractional_linear_relation = search_fractional_linear_relation(
            target_series=target_series[:solve_order],
            basis_series_by_variable={"W_34": weighted_series[:solve_order]},
            order=solve_order,
        )

    correction_eta_scans = tuple(
        scan_ratio_eta_quotient_relations(
            ratio_series=weighted_correction_series,
            levels=correction_eta_levels,
            order=solve_order,
            max_abs_exponent=correction_max_abs_exponent,
        )
    )
    correction_modular_unit_eta_scans = tuple(
        scan_ratio_modular_unit_eta_relations(
            ratio_series=weighted_correction_series,
            moduli=correction_moduli,
            eta_levels=correction_eta_levels,
            order=solve_order,
            max_abs_exponent=correction_max_abs_exponent,
        )
    )
    normalized_correction_gap = build_gap_normalized_series(target_series=weighted_correction_series)
    normalized_correction_label = None
    normalized_correction_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    normalized_correction_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    normalized_correction_source_family_eta_scans: tuple[SourceFamilyEtaCorrectionScan, ...] = ()
    second_normalized_correction_label = None
    second_normalized_correction_gap = None
    second_normalized_correction_eta_scans: tuple[EtaQuotientRelationScan, ...] = ()
    second_normalized_correction_modular_unit_eta_scans: tuple[ModularUnitEtaRelationScan, ...] = ()
    second_normalized_correction_source_family_eta_scans: tuple[SourceFamilyEtaCorrectionScan, ...] = ()
    second_normalized_correction_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...] = ()
    second_normalized_correction_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...] = ()
    second_normalized_correction_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_polynomial_scans: tuple[NamedPolynomialRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_multiplicative_scans: tuple[NamedMultiplicativeRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_fractional_linear_scans: tuple[NamedFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans: tuple[NamedTwoLayerFractionalLinearRelationScan, ...] = ()
    second_normalized_correction_explicit_transform_eta_scans: tuple[ExplicitSourceFamilyEtaCorrectionScan, ...] = ()
    if normalized_correction_gap is not None:
        normalized_correction_label = "G_W34"
        normalized_correction_series = list(normalized_correction_gap.normalized_series[:solve_order])
        normalized_correction_eta_scans = tuple(
            scan_ratio_eta_quotient_relations(
                ratio_series=normalized_correction_series,
                levels=correction_eta_levels,
                order=solve_order,
                max_abs_exponent=correction_max_abs_exponent,
            )
        )
        normalized_correction_modular_unit_eta_scans = tuple(
            scan_ratio_modular_unit_eta_relations(
                ratio_series=normalized_correction_series,
                moduli=correction_moduli,
                eta_levels=correction_eta_levels,
                order=solve_order,
                max_abs_exponent=correction_max_abs_exponent,
            )
        )
        normalized_correction_source_family_eta_scans = tuple(
            scan_source_family_eta_corrections(
                target_series=normalized_correction_series,
                ordered_base_families=tuple(
                    (label, benchmark_name, basis_series[:solve_order])
                    for label, benchmark_name, basis_series in correction_source_families
                ),
                powers=correction_source_powers,
                eta_levels=correction_eta_levels,
                order=solve_order,
                max_abs_exponent=correction_max_abs_exponent,
                supplemental_powers_by_family=correction_source_supplemental_powers_by_family,
            )
        )
        second_normalized_correction_gap = build_gap_normalized_series(
            target_series=list(normalized_correction_gap.normalized_series)
        )
        if second_normalized_correction_gap is not None:
            second_normalized_correction_label = "G2_W34"
            second_normalized_correction_series = list(
                second_normalized_correction_gap.normalized_series[:solve_order]
            )
            deeper_solve_order = coordinate_solve_order or solve_order
            second_normalized_correction_eta_scans = tuple(
                scan_ratio_eta_quotient_relations(
                    ratio_series=second_normalized_correction_series,
                    levels=correction_eta_levels,
                    order=solve_order,
                    max_abs_exponent=correction_max_abs_exponent,
                )
            )
            second_normalized_correction_modular_unit_eta_scans = tuple(
                scan_ratio_modular_unit_eta_relations(
                    ratio_series=second_normalized_correction_series,
                    moduli=correction_moduli,
                    eta_levels=correction_eta_levels,
                    order=solve_order,
                    max_abs_exponent=correction_max_abs_exponent,
                )
            )
            second_normalized_correction_source_family_eta_scans = tuple(
                scan_source_family_eta_corrections(
                    target_series=second_normalized_correction_series,
                    ordered_base_families=tuple(
                        (label, benchmark_name, basis_series[:solve_order])
                        for label, benchmark_name, basis_series in correction_source_families
                    ),
                    powers=correction_source_powers,
                    eta_levels=correction_eta_levels,
                    order=solve_order,
                    max_abs_exponent=correction_max_abs_exponent,
                    supplemental_powers_by_family=correction_source_supplemental_powers_by_family,
                )
            )
            second_normalized_correction_quotient_prefix_scans = scan_named_prefix_boxes(
                target_series=second_normalized_correction_series,
                ordered_basis_series=tuple(
                    (label, series[:solve_order]) for label, series in quotient_ordered_basis_series
                ),
                order=solve_order,
                degree_values=coordinate_degree_values,
                required_variable="F",
                max_abs_exponent=coordinate_max_abs_exponent,
                solve_order=deeper_solve_order,
            )
            second_normalized_correction_quotient_polynomial_scans = (
                second_normalized_correction_quotient_prefix_scans.polynomial_scans
            )
            second_normalized_correction_quotient_multiplicative_scans = (
                second_normalized_correction_quotient_prefix_scans.multiplicative_scans
            )
            second_normalized_correction_quotient_fractional_linear_scans = (
                second_normalized_correction_quotient_prefix_scans.fractional_linear_scans
            )
            second_normalized_correction_quotient_two_layer_fractional_linear_scans = (
                second_normalized_correction_quotient_prefix_scans.two_layer_fractional_linear_scans
            )
            second_normalized_correction_mixed_prefix_scans = scan_named_prefix_boxes(
                target_series=second_normalized_correction_series,
                ordered_basis_series=tuple(
                    (label, series[:solve_order]) for label, series in mixed_quotient_ordered_basis_series
                ),
                order=solve_order,
                degree_values=coordinate_degree_values,
                required_variable="F",
                max_abs_exponent=coordinate_max_abs_exponent,
                solve_order=deeper_solve_order,
            )
            second_normalized_correction_mixed_quotient_polynomial_scans = (
                second_normalized_correction_mixed_prefix_scans.polynomial_scans
            )
            second_normalized_correction_mixed_quotient_multiplicative_scans = (
                second_normalized_correction_mixed_prefix_scans.multiplicative_scans
            )
            second_normalized_correction_mixed_quotient_fractional_linear_scans = (
                second_normalized_correction_mixed_prefix_scans.fractional_linear_scans
            )
            second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans = (
                second_normalized_correction_mixed_prefix_scans.two_layer_fractional_linear_scans
            )
            second_normalized_correction_explicit_transform_eta_scans = tuple(
                scan_explicit_source_family_eta_correction_templates(
                    target_series=second_normalized_correction_series,
                    ordered_base_families=(("GG", benchmark_name, gg_series[:solve_order]),),
                    powers=correction_source_powers,
                    eta_levels=correction_eta_levels,
                    order=solve_order,
                    max_abs_exponent=correction_max_abs_exponent,
                    supplemental_powers_by_family=(
                        None
                        if correction_source_supplemental_powers_by_family is None
                        else {
                            "GG": correction_source_supplemental_powers_by_family.get("GG", ())
                        }
                    ),
                )
            )

    return (
        GGWeightedCoordinateDiagnostic(
            label="W_34",
            expression="Q_3^3 / Q_4^2",
            log_expression="3*log(Q_3) - 2*log(Q_4)",
            correction_expression="F / W_34",
            first_difference_power=_first_nonzero_residual_term(direct_difference)[0],
            first_difference_coeff=_first_nonzero_residual_term(direct_difference)[1],
            first_log_difference_power=_first_nonzero_residual_term(log_difference)[0],
            first_log_difference_coeff=_first_nonzero_residual_term(log_difference)[1],
            correction_first_gap_power=_first_nonzero_residual_term(weighted_correction_gap)[0],
            correction_first_gap_coeff=_first_nonzero_residual_term(weighted_correction_gap)[1],
            polynomial_degree1_relation=polynomial_degree1_relation,
            polynomial_degree2_relation=polynomial_degree2_relation,
            fractional_linear_relation=fractional_linear_relation,
            correction_eta_scans=correction_eta_scans,
            correction_modular_unit_eta_scans=correction_modular_unit_eta_scans,
            normalized_correction_label=normalized_correction_label,
            normalized_correction_gap=normalized_correction_gap,
            normalized_correction_eta_scans=normalized_correction_eta_scans,
            normalized_correction_modular_unit_eta_scans=normalized_correction_modular_unit_eta_scans,
            normalized_correction_source_family_eta_scans=normalized_correction_source_family_eta_scans,
            second_normalized_correction_label=second_normalized_correction_label,
            second_normalized_correction_gap=second_normalized_correction_gap,
            second_normalized_correction_eta_scans=second_normalized_correction_eta_scans,
            second_normalized_correction_modular_unit_eta_scans=second_normalized_correction_modular_unit_eta_scans,
            second_normalized_correction_source_family_eta_scans=second_normalized_correction_source_family_eta_scans,
            second_normalized_correction_quotient_polynomial_scans=second_normalized_correction_quotient_polynomial_scans,
            second_normalized_correction_quotient_multiplicative_scans=second_normalized_correction_quotient_multiplicative_scans,
            second_normalized_correction_quotient_fractional_linear_scans=second_normalized_correction_quotient_fractional_linear_scans,
            second_normalized_correction_quotient_two_layer_fractional_linear_scans=second_normalized_correction_quotient_two_layer_fractional_linear_scans,
            second_normalized_correction_mixed_quotient_polynomial_scans=second_normalized_correction_mixed_quotient_polynomial_scans,
            second_normalized_correction_mixed_quotient_multiplicative_scans=second_normalized_correction_mixed_quotient_multiplicative_scans,
            second_normalized_correction_mixed_quotient_fractional_linear_scans=second_normalized_correction_mixed_quotient_fractional_linear_scans,
            second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans=second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans,
            second_normalized_correction_explicit_transform_eta_scans=second_normalized_correction_explicit_transform_eta_scans,
        ),
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
    weighted_correction_eta_levels: tuple[int, ...] = (1, 2, 3, 4),
    weighted_correction_moduli: tuple[int, ...] = (2, 3, 4),
    weighted_correction_max_abs_exponent: int = 4,
    weighted_correction_source_families: tuple[tuple[str, str, Series], ...] = (),
    weighted_correction_source_powers: tuple[int, ...] = (2, 3, 4),
    weighted_correction_source_supplemental_powers_by_family: dict[str, tuple[int, ...]] | None = None,
    include_weighted_coordinate_diagnostics: bool = True,
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
        exact_polynomial_template_obstructions,
        quotient_exact_polynomial_template_labels,
        quotient_exact_polynomial_template_hits,
        quotient_exact_polynomial_template_obstructions,
    ) = _gg_exact_modular_template_hits(
        target_series=target_series,
        ordered_basis_entries=ordered_basis_entries,
        quotient_basis_entries=quotient_basis_series,
        order=order,
    )
    weighted_coordinate_diagnostics: tuple[GGWeightedCoordinateDiagnostic, ...] = ()
    if include_weighted_coordinate_diagnostics:
        weighted_coordinate_diagnostics = _gg_weighted_coordinate_diagnostics(
            target_series=target_series,
            benchmark_name=benchmark_name,
            gg_series=gg_series,
            quotient_basis_entries=quotient_basis_series,
            order=order,
            coordinate_degree_values=degree_values,
            coordinate_max_abs_exponent=max_abs_exponent,
            coordinate_solve_order=solve_order,
            correction_eta_levels=weighted_correction_eta_levels,
            correction_moduli=weighted_correction_moduli,
            correction_max_abs_exponent=weighted_correction_max_abs_exponent,
            correction_source_families=weighted_correction_source_families,
            correction_source_powers=weighted_correction_source_powers,
            correction_source_supplemental_powers_by_family=weighted_correction_source_supplemental_powers_by_family,
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
        exact_polynomial_template_obstructions=exact_polynomial_template_obstructions,
        polynomial_scans=direct_scans.polynomial_scans,
        multiplicative_scans=direct_scans.multiplicative_scans,
        fractional_linear_scans=direct_scans.fractional_linear_scans,
        two_layer_fractional_linear_scans=direct_scans.two_layer_fractional_linear_scans,
        quotient_basis_series=quotient_basis_series,
        quotient_exact_polynomial_template_labels=quotient_exact_polynomial_template_labels,
        quotient_exact_polynomial_template_hits=quotient_exact_polynomial_template_hits,
        quotient_exact_polynomial_template_obstructions=quotient_exact_polynomial_template_obstructions,
        quotient_polynomial_scans=quotient_scans.polynomial_scans,
        quotient_multiplicative_scans=quotient_scans.multiplicative_scans,
        quotient_fractional_linear_scans=quotient_scans.fractional_linear_scans,
        quotient_two_layer_fractional_linear_scans=quotient_scans.two_layer_fractional_linear_scans,
        weighted_coordinate_diagnostics=weighted_coordinate_diagnostics,
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
    reduced_tail_anchor: ReducedTailAnchor | None = None
    reduced_next_tail_anchor: ReducedTailAnchor | None = None
    reduced_tail_anchor_gap: GapNormalizedSeries | None = None
    reduced_next_tail_reciprocal_gap: GapNormalizedSeries | None = None
    reduced_tail_anchor_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_tail_anchor_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_tail_anchor_eta_scans: list[EtaQuotientRelationScan] = []
    reduced_tail_anchor_self_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_self_plus_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_self_signed_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_tail_anchor_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_tail_anchor_normalized_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_tail_anchor_normalized_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_tail_anchor_normalized_eta_scans: list[EtaQuotientRelationScan] = []
    reduced_tail_anchor_normalized_self_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_normalized_self_plus_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_normalized_self_signed_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_tail_anchor_normalized_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_next_tail_reciprocal_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_next_tail_reciprocal_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_next_tail_reciprocal_eta_scans: list[EtaQuotientRelationScan] = []
    reduced_next_tail_reciprocal_self_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_self_plus_scans: list[SelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_self_signed_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_tail_anchor_gap_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_tail_anchor_gap_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_tail_anchor_gap_eta_scans: list[EtaQuotientRelationScan] = []
    reduced_tail_anchor_gap_self_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_gap_self_plus_scans: list[SelfQuotientProductRelationScan] = []
    reduced_tail_anchor_gap_self_signed_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_tail_anchor_gap_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_tail_anchor_gap_source_family_eta_scans: list[SourceFamilyEtaCorrectionScan] = []
    reduced_tail_anchor_second_gap: GapNormalizedSeries | None = None
    reduced_tail_anchor_second_gap_source_family_eta_scans: list[SourceFamilyEtaCorrectionScan] = []
    reduced_next_tail_reciprocal_gap_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_next_tail_reciprocal_gap_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_next_tail_reciprocal_gap_eta_scans: list[EtaQuotientRelationScan] = []
    reduced_next_tail_reciprocal_gap_self_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_gap_self_plus_scans: list[SelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_gap_self_signed_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_next_tail_reciprocal_gap_self_signed_eta_scans: list[SelfSignedEtaRelationScan] = []
    reduced_next_tail_reciprocal_gap_source_family_eta_scans: list[SourceFamilyEtaCorrectionScan] = []
    reduced_next_tail_reciprocal_second_gap: GapNormalizedSeries | None = None
    reduced_next_tail_reciprocal_second_gap_source_family_eta_scans: list[SourceFamilyEtaCorrectionScan] = []
    reduced_object_self_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_object_self_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_ratio_self_polynomial_scan = SelfPolynomialUniquenessScan((), (), (), ())
    reduced_ratio_self_fractional_linear_scan = SelfFractionalLinearUniquenessScan((), (), ())
    reduced_object_mahler_scan = SelfMahlerLinearScan((), (), (), ())
    reduced_ratio_mahler_scan = SelfMahlerLinearScan((), (), (), ())
    reduced_ratio_self_quotient_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_ratio_self_plus_product_scans: list[SelfQuotientProductRelationScan] = []
    reduced_ratio_self_signed_product_scans: list[SignedSelfQuotientProductRelationScan] = []
    reduced_ratio_self_plus_pochhammer_scans: list[SelfPlusPochhammerRelationScan] = []
    reduced_ratio_self_plus_pochhammer_eta_scans: list[SelfPlusPochhammerEtaRelationScan] = []
    reduced_ratio_source_family_self_plus_pochhammer_eta_scans: list[SourceFamilySelfPlusPochhammerEtaCorrectionScan] = []
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
        reduced_tail_anchor = build_reduced_tail_anchor(
            reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
            symbol=sp.Symbol(series_symbol),
            start_stage=3,
            order=min(reduced_bridge_order, 24 if smoke else 24),
        )
        reduced_next_tail_anchor = build_reduced_tail_anchor(
            reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
            symbol=sp.Symbol(series_symbol),
            start_stage=4,
            order=min(reduced_bridge_order, 24 if smoke else 24),
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
        reduced_ratio_self_plus_pochhammer_scans = scan_ratio_self_plus_pochhammer_relations(
            ratio_series=reduced_ratio_series,
            moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 4),
            order=reduced_bridge_order,
            max_abs_exponent=6 if smoke else 8,
        )
        reduced_ratio_self_plus_pochhammer_eta_scans = scan_ratio_self_plus_pochhammer_eta_relations(
            ratio_series=reduced_ratio_series,
            moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 4),
            eta_levels=_eta_scan_levels((1, 2, 3, 4)),
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
        if reduced_tail_anchor is not None:
            tail_anchor_order = len(reduced_tail_anchor.tail_series)
            tail_anchor_moduli = tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 4)
            tail_anchor_eta_levels = _eta_scan_levels((1, 2, 3, 4, 6, 12))
            reduced_tail_anchor_gap = build_gap_normalized_series(target_series=list(reduced_tail_anchor.normalized_series))
            reduced_tail_anchor_polynomial_scan = scan_self_polynomial_uniqueness_relations(
                target_series=list(reduced_tail_anchor.tail_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                fg_degree_values=(1, 2),
                t_degree_values=(1, 2),
            )
            reduced_tail_anchor_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
                target_series=list(reduced_tail_anchor.tail_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                t_degree_values=(1, 2),
            )
            reduced_tail_anchor_eta_scans = scan_ratio_eta_quotient_relations(
                ratio_series=list(reduced_tail_anchor.tail_series),
                levels=tail_anchor_eta_levels,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_self_product_scans = scan_ratio_self_quotient_product_relations(
                ratio_series=list(reduced_tail_anchor.tail_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_self_plus_scans = scan_ratio_self_plus_product_relations(
                ratio_series=list(reduced_tail_anchor.tail_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_self_signed_scans = scan_ratio_self_signed_product_relations(
                ratio_series=list(reduced_tail_anchor.tail_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
                ratio_series=list(reduced_tail_anchor.tail_series),
                moduli=tuple(modulus for modulus in tail_anchor_moduli if modulus <= 3),
                eta_levels=_eta_scan_levels((1, 2, 3)),
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_normalized_polynomial_scan = scan_self_polynomial_uniqueness_relations(
                target_series=list(reduced_tail_anchor.normalized_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                fg_degree_values=(1, 2),
                t_degree_values=(1, 2),
            )
            reduced_tail_anchor_normalized_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
                target_series=list(reduced_tail_anchor.normalized_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                t_degree_values=(1, 2),
            )
            reduced_tail_anchor_normalized_eta_scans = scan_ratio_eta_quotient_relations(
                ratio_series=list(reduced_tail_anchor.normalized_series),
                levels=tail_anchor_eta_levels,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_normalized_self_product_scans = scan_ratio_self_quotient_product_relations(
                ratio_series=list(reduced_tail_anchor.normalized_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_normalized_self_plus_scans = scan_ratio_self_plus_product_relations(
                ratio_series=list(reduced_tail_anchor.normalized_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_normalized_self_signed_scans = scan_ratio_self_signed_product_relations(
                ratio_series=list(reduced_tail_anchor.normalized_series),
                moduli=tail_anchor_moduli,
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_tail_anchor_normalized_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
                ratio_series=list(reduced_tail_anchor.normalized_series),
                moduli=tuple(modulus for modulus in tail_anchor_moduli if modulus <= 3),
                eta_levels=_eta_scan_levels((1, 2, 3)),
                order=tail_anchor_order,
                max_abs_exponent=6 if smoke else 8,
            )
            if reduced_tail_anchor_gap is not None:
                reduced_tail_anchor_second_gap = build_gap_normalized_series(
                    target_series=list(reduced_tail_anchor_gap.normalized_series)
                )
                reduced_tail_anchor_gap_polynomial_scan = scan_self_polynomial_uniqueness_relations(
                    target_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tail_anchor_moduli,
                    order=tail_anchor_order,
                    fg_degree_values=(1, 2),
                    t_degree_values=(1, 2),
                )
                reduced_tail_anchor_gap_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
                    target_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tail_anchor_moduli,
                    order=tail_anchor_order,
                    t_degree_values=(1, 2),
                )
                reduced_tail_anchor_gap_eta_scans = scan_ratio_eta_quotient_relations(
                    ratio_series=list(reduced_tail_anchor_gap.normalized_series),
                    levels=tail_anchor_eta_levels,
                    order=tail_anchor_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_tail_anchor_gap_self_product_scans = scan_ratio_self_quotient_product_relations(
                    ratio_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tail_anchor_moduli,
                    order=tail_anchor_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_tail_anchor_gap_self_plus_scans = scan_ratio_self_plus_product_relations(
                    ratio_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tail_anchor_moduli,
                    order=tail_anchor_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_tail_anchor_gap_self_signed_scans = scan_ratio_self_signed_product_relations(
                    ratio_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tail_anchor_moduli,
                    order=tail_anchor_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_tail_anchor_gap_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
                    ratio_series=list(reduced_tail_anchor_gap.normalized_series),
                    moduli=tuple(modulus for modulus in tail_anchor_moduli if modulus <= 3),
                    eta_levels=_eta_scan_levels((1, 2, 3)),
                    order=tail_anchor_order,
                    max_abs_exponent=6 if smoke else 8,
                )
        if reduced_next_tail_anchor is not None:
            next_tail_reciprocal_series = series_invert(list(reduced_next_tail_anchor.normalized_series))
            next_tail_order = len(next_tail_reciprocal_series)
            next_tail_moduli = tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 4)
            next_tail_eta_levels = _eta_scan_levels((1, 2, 3, 4, 6, 12))
            reduced_next_tail_reciprocal_gap = build_gap_normalized_series(target_series=next_tail_reciprocal_series)
            reduced_next_tail_reciprocal_polynomial_scan = scan_self_polynomial_uniqueness_relations(
                target_series=next_tail_reciprocal_series,
                moduli=next_tail_moduli,
                order=next_tail_order,
                fg_degree_values=(1, 2),
                t_degree_values=(1, 2),
            )
            reduced_next_tail_reciprocal_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
                target_series=next_tail_reciprocal_series,
                moduli=next_tail_moduli,
                order=next_tail_order,
                t_degree_values=(1, 2),
            )
            reduced_next_tail_reciprocal_eta_scans = scan_ratio_eta_quotient_relations(
                ratio_series=next_tail_reciprocal_series,
                levels=next_tail_eta_levels,
                order=next_tail_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_next_tail_reciprocal_self_product_scans = scan_ratio_self_quotient_product_relations(
                ratio_series=next_tail_reciprocal_series,
                moduli=next_tail_moduli,
                order=next_tail_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_next_tail_reciprocal_self_plus_scans = scan_ratio_self_plus_product_relations(
                ratio_series=next_tail_reciprocal_series,
                moduli=next_tail_moduli,
                order=next_tail_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_next_tail_reciprocal_self_signed_scans = scan_ratio_self_signed_product_relations(
                ratio_series=next_tail_reciprocal_series,
                moduli=next_tail_moduli,
                order=next_tail_order,
                max_abs_exponent=6 if smoke else 8,
            )
            reduced_next_tail_reciprocal_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
                ratio_series=next_tail_reciprocal_series,
                moduli=tuple(modulus for modulus in next_tail_moduli if modulus <= 3),
                eta_levels=_eta_scan_levels((1, 2, 3)),
                order=next_tail_order,
                max_abs_exponent=6 if smoke else 8,
            )
            if reduced_next_tail_reciprocal_gap is not None:
                reduced_next_tail_reciprocal_second_gap = build_gap_normalized_series(
                    target_series=list(reduced_next_tail_reciprocal_gap.normalized_series)
                )
                reduced_next_tail_reciprocal_gap_polynomial_scan = scan_self_polynomial_uniqueness_relations(
                    target_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=next_tail_moduli,
                    order=next_tail_order,
                    fg_degree_values=(1, 2),
                    t_degree_values=(1, 2),
                )
                reduced_next_tail_reciprocal_gap_fractional_linear_scan = scan_self_fractional_linear_uniqueness_relations(
                    target_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=next_tail_moduli,
                    order=next_tail_order,
                    t_degree_values=(1, 2),
                )
                reduced_next_tail_reciprocal_gap_eta_scans = scan_ratio_eta_quotient_relations(
                    ratio_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    levels=next_tail_eta_levels,
                    order=next_tail_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_next_tail_reciprocal_gap_self_product_scans = scan_ratio_self_quotient_product_relations(
                    ratio_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=next_tail_moduli,
                    order=next_tail_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_next_tail_reciprocal_gap_self_plus_scans = scan_ratio_self_plus_product_relations(
                    ratio_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=next_tail_moduli,
                    order=next_tail_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_next_tail_reciprocal_gap_self_signed_scans = scan_ratio_self_signed_product_relations(
                    ratio_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=next_tail_moduli,
                    order=next_tail_order,
                    max_abs_exponent=6 if smoke else 8,
                )
                reduced_next_tail_reciprocal_gap_self_signed_eta_scans = scan_ratio_self_signed_eta_relations(
                    ratio_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
                    moduli=tuple(modulus for modulus in next_tail_moduli if modulus <= 3),
                    eta_levels=_eta_scan_levels((1, 2, 3)),
                    order=next_tail_order,
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
            _canonical_benchmark_series(
                benchmark_name,
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
    gap_tail_source_family_powers = (2, 3, 4)
    gap_tail_source_family_eta_levels = _eta_scan_levels((1, 2, 3, 4))
    if reduced_tail_anchor_gap is not None:
        gap_tail_order = len(reduced_tail_anchor_gap.normalized_series)
        reduced_tail_anchor_gap_source_family_eta_scans = scan_source_family_eta_corrections(
            target_series=list(reduced_tail_anchor_gap.normalized_series),
            ordered_base_families=tuple(
                (label, benchmark_name, basis_series[:gap_tail_order])
                for label, benchmark_name, basis_series in source_family_base_series
            ),
            powers=gap_tail_source_family_powers,
            eta_levels=gap_tail_source_family_eta_levels,
            order=gap_tail_order,
            max_abs_exponent=4 if smoke else 6,
        )
    if reduced_tail_anchor_second_gap is not None:
        gap_tail_order = len(reduced_tail_anchor_second_gap.normalized_series)
        reduced_tail_anchor_second_gap_source_family_eta_scans = scan_source_family_eta_corrections(
            target_series=list(reduced_tail_anchor_second_gap.normalized_series),
            ordered_base_families=tuple(
                (label, benchmark_name, basis_series[:gap_tail_order])
                for label, benchmark_name, basis_series in source_family_base_series
            ),
            powers=gap_tail_source_family_powers,
            eta_levels=gap_tail_source_family_eta_levels,
            order=gap_tail_order,
            max_abs_exponent=4 if smoke else 6,
        )
    if reduced_next_tail_reciprocal_gap is not None:
        gap_tail_order = len(reduced_next_tail_reciprocal_gap.normalized_series)
        reduced_next_tail_reciprocal_gap_source_family_eta_scans = scan_source_family_eta_corrections(
            target_series=list(reduced_next_tail_reciprocal_gap.normalized_series),
            ordered_base_families=tuple(
                (label, benchmark_name, basis_series[:gap_tail_order])
                for label, benchmark_name, basis_series in source_family_base_series
            ),
            powers=gap_tail_source_family_powers,
            eta_levels=gap_tail_source_family_eta_levels,
            order=gap_tail_order,
            max_abs_exponent=4 if smoke else 6,
        )
    if reduced_next_tail_reciprocal_second_gap is not None:
        gap_tail_order = len(reduced_next_tail_reciprocal_second_gap.normalized_series)
        reduced_next_tail_reciprocal_second_gap_source_family_eta_scans = scan_source_family_eta_corrections(
            target_series=list(reduced_next_tail_reciprocal_second_gap.normalized_series),
            ordered_base_families=tuple(
                (label, benchmark_name, basis_series[:gap_tail_order])
                for label, benchmark_name, basis_series in source_family_base_series
            ),
            powers=gap_tail_source_family_powers,
            eta_levels=gap_tail_source_family_eta_levels,
            order=gap_tail_order,
            max_abs_exponent=4 if smoke else 6,
        )
    if reduced_ratio_series is not None:
        reduced_ratio_source_family_self_plus_pochhammer_eta_scans = (
            scan_source_family_self_plus_pochhammer_eta_corrections(
                target_series=reduced_ratio_series,
                ordered_base_families=tuple(
                    (label, benchmark_name, basis_series[:reduced_bridge_order])
                    for label, benchmark_name, basis_series in source_family_base_series
                ),
                powers=(2, 3, 4),
                moduli=tuple(modulus for modulus in reduced_bridge_moduli if modulus <= 4),
                eta_levels=_eta_scan_levels((1, 2, 3, 4)),
                order=reduced_bridge_order,
                max_abs_exponent=4 if smoke else 6,
            )
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
    weighted_correction_source_families = tuple(
        entry for entry in source_family_base_series if entry[0] in {"RR", "GG"}
    )
    weighted_correction_source_supplemental_powers = {
        family_label: supplemental_source_family_powers.get(family_label, ())
        for family_label, _, _ in weighted_correction_source_families
        if family_label in supplemental_source_family_powers
    }
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
            weighted_correction_eta_levels=eta_scan_levels,
            weighted_correction_moduli=tuple(modulus for modulus in source_family_scan_powers if modulus <= 4) or (2, 3, 4),
            weighted_correction_max_abs_exponent=4 if smoke else 6,
            weighted_correction_source_families=weighted_correction_source_families,
            weighted_correction_source_powers=source_family_scan_powers,
            weighted_correction_source_supplemental_powers_by_family=weighted_correction_source_supplemental_powers,
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
    ratio_modular_unit_eta_scans: list[ModularUnitEtaRelationScan] = []
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
        ratio_modular_unit_eta_scans = scan_ratio_modular_unit_eta_relations(
            ratio_series=ratio_series,
            moduli=tuple(modulus for modulus in tuple(benchmark_power_series) if modulus <= 4) or (2, 3),
            eta_levels=_eta_scan_levels(tuple(benchmark_power_series)),
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
        if reduced_tail_anchor is not None:
            tail_state = _format_expr(reduced_tail_anchor.state_expr)
            tail_anchor_hits = any(
                (
                    reduced_tail_anchor_polynomial_scan.hits,
                    reduced_tail_anchor_fractional_linear_scan.hits,
                    [scan for scan in reduced_tail_anchor_eta_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_self_product_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_self_plus_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_self_signed_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_self_signed_eta_scans if scan.relation is not None],
                )
            )
            normalized_tail_anchor_hits = any(
                (
                    reduced_tail_anchor_normalized_polynomial_scan.hits,
                    reduced_tail_anchor_normalized_fractional_linear_scan.hits,
                    [scan for scan in reduced_tail_anchor_normalized_eta_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_normalized_self_product_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_normalized_self_plus_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_normalized_self_signed_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_normalized_self_signed_eta_scans if scan.relation is not None],
                )
            )
            lines.extend(
                [
                    "- We also anchored that tail law at the first stationary stage and scanned the concrete tail object itself.",
                    "",
                    "```text",
                    f"T_tail = T({tail_state})",
                    f"U_tail = T_tail / (1 + {tail_state})",
                    "```",
                    "",
                ]
            )
            if not tail_anchor_hits:
                lines.append(
                    "No anchored-tail hit was found in the scanned self-polynomial, self-fractional-linear, eta-quotient, finite-product, plus-product, signed-product, or signed-eta boxes."
                )
            else:
                lines.append("Anchored-tail hits were found in the scanned box family.")
            lines.append("")
            if not normalized_tail_anchor_hits:
                lines.append(
                    "No normalized anchored-tail hit was found in the same scanned box family after dividing by the visible factor `1 + x`."
                )
            else:
                lines.append("Normalized anchored-tail hits were found in the scanned box family.")
            lines.append("")
        if reduced_next_tail_anchor is not None:
            next_tail_state = _format_expr(reduced_next_tail_anchor.state_expr)
            next_tail_reciprocal_hits = any(
                (
                    reduced_next_tail_reciprocal_polynomial_scan.hits,
                    reduced_next_tail_reciprocal_fractional_linear_scan.hits,
                    [scan for scan in reduced_next_tail_reciprocal_eta_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_self_product_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_self_plus_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_self_signed_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_self_signed_eta_scans if scan.relation is not None],
                )
            )
            lines.extend(
                [
                    "- We also pushed one step deeper along the exact tail law and scanned the reciprocal-normalized next-tail object.",
                    "",
                    "```text",
                    f"R_tail = (1 + {next_tail_state}) / T({next_tail_state})",
                    "```",
                    "",
                ]
            )
            if not next_tail_reciprocal_hits:
                lines.append(
                    "No reciprocal-normalized next-tail hit was found in the scanned self-polynomial, self-fractional-linear, eta-quotient, finite-product, plus-product, signed-product, or signed-eta boxes."
                )
            else:
                lines.append("Reciprocal-normalized next-tail hits were found in the scanned box family.")
            lines.append("")
        if reduced_tail_anchor_gap is not None or reduced_next_tail_reciprocal_gap is not None:
            gap_tail_hits = any(
                (
                    reduced_tail_anchor_gap_polynomial_scan.hits,
                    reduced_tail_anchor_gap_fractional_linear_scan.hits,
                    [scan for scan in reduced_tail_anchor_gap_eta_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_gap_self_product_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_gap_self_plus_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_gap_self_signed_scans if scan.relation is not None],
                    [scan for scan in reduced_tail_anchor_gap_self_signed_eta_scans if scan.relation is not None],
                    reduced_next_tail_reciprocal_gap_polynomial_scan.hits,
                    reduced_next_tail_reciprocal_gap_fractional_linear_scan.hits,
                    [scan for scan in reduced_next_tail_reciprocal_gap_eta_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_gap_self_product_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_gap_self_plus_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_gap_self_signed_scans if scan.relation is not None],
                    [scan for scan in reduced_next_tail_reciprocal_gap_self_signed_eta_scans if scan.relation is not None],
                )
            )
            def _format_gap_formula(*, source_variable: str, target_variable: str, gap: GapNormalizedSeries) -> str:
                leading_coeff = _format_expr(gap.leading_coefficient)
                shift = gap.shift
                if leading_coeff == "1":
                    return f"{target_variable} = ({source_variable} - 1) / {series_symbol}^{shift}"
                if leading_coeff == "-1":
                    return f"{target_variable} = (1 - {source_variable}) / {series_symbol}^{shift}"
                return f"{target_variable} = ({source_variable} - 1) / ({leading_coeff}*{series_symbol}^{shift})"

            lines.extend(
                [
                    "- We also stripped off the first visible nonzero gap from the renormalized tail objects.",
                    "",
                    "```text",
                ]
            )
            if reduced_tail_anchor_gap is not None:
                lines.append(
                    _format_gap_formula(
                        source_variable="U_tail",
                        target_variable="G_tail",
                        gap=reduced_tail_anchor_gap,
                    )
                )
            if reduced_next_tail_reciprocal_gap is not None:
                lines.append(
                    _format_gap_formula(
                        source_variable="R_tail",
                        target_variable="H_tail",
                        gap=reduced_next_tail_reciprocal_gap,
                    )
                )
            lines.extend(["```", ""])
            if not gap_tail_hits:
                lines.append(
                    "No gap-normalized tail hit was found in the scanned self-polynomial, self-fractional-linear, eta-quotient, finite-product, plus-product, signed-product, or signed-eta boxes."
                )
            else:
                lines.append("Gap-normalized tail hits were found in the scanned box family.")
            lines.append("")
            lines.append(
                "Method compression: both exact tail residuals now fail twice in bounded boxes — first in their raw normalized form, and again after stripping the first visible nonzero gap (`t^3` for `U_tail - 1`, `t^4` for `1 - R_tail`)."
            )
            lines.append("")
            if reduced_tail_anchor_gap_source_family_eta_scans or reduced_next_tail_reciprocal_gap_source_family_eta_scans:
                lines.extend(
                    [
                        "We also checked whether the first-gap tail residuals look more like one nearby source-family core times a small eta tail than like a standalone tiny self-equation:",
                        "",
                        "```text",
                        "G_tail = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                        "H_tail = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                        "```",
                        "",
                    ]
                )

                def _flatten_source_family_eta_hits(
                    scans: list[SourceFamilyEtaCorrectionScan],
                ) -> list[tuple[str, str, str, int, MultiplicativeRelation]]:
                    hits: list[tuple[str, str, str, int, MultiplicativeRelation]] = []
                    for family_scan in scans:
                        for basis_scan in family_scan.direct_basis_scans:
                            for eta_scan in basis_scan.eta_scans:
                                if eta_scan.relation is not None:
                                    hits.append(
                                        (
                                            basis_scan.basis_label,
                                            basis_scan.basis_expression,
                                            "raw",
                                            eta_scan.level,
                                            eta_scan.relation,
                                        )
                                    )
                        for basis_scan in family_scan.quotient_basis_scans:
                            for eta_scan in basis_scan.eta_scans:
                                if eta_scan.relation is not None:
                                    hits.append(
                                        (
                                            basis_scan.basis_label,
                                            basis_scan.basis_expression,
                                            "quotient",
                                            eta_scan.level,
                                            eta_scan.relation,
                                        )
                                    )
                    return hits

                for target_label, target_variable, target_scans in (
                    ("`G_tail`", "G_tail", reduced_tail_anchor_gap_source_family_eta_scans),
                    ("`H_tail`", "H_tail", reduced_next_tail_reciprocal_gap_source_family_eta_scans),
                ):
                    if not target_scans:
                        continue
                    target_hits = _flatten_source_family_eta_hits(target_scans)
                    if not target_hits:
                        lines.append(
                            f"- {target_label}: no one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`."
                        )
                        continue
                    lines.append(f"- {target_label}: source-family eta-correction hits were found:")
                    for basis_label, basis_expression, basis_kind, level, relation in target_hits:
                        lines.append(
                            f"  - {basis_kind} basis `{basis_label}`, `N={level}`: `{_format_source_family_eta_correction(basis_expression=basis_expression, relation=relation, target_variable=target_variable, series_symbol=series_symbol)}`"
                        )
                lines.append("")

            if reduced_tail_anchor_second_gap is not None or reduced_next_tail_reciprocal_second_gap is not None:
                lines.extend(
                    [
                        "We then stripped one more visible nonzero gap from those first-gap residuals and kept only the more source-directed one-core eta-correction question:",
                        "",
                        "```text",
                    ]
                )
                if reduced_tail_anchor_second_gap is not None:
                    lines.append(
                        _format_gap_formula(
                            source_variable="G_tail",
                            target_variable="G2_tail",
                            gap=reduced_tail_anchor_second_gap,
                        )
                    )
                if reduced_next_tail_reciprocal_second_gap is not None:
                    lines.append(
                        _format_gap_formula(
                            source_variable="H_tail",
                            target_variable="H2_tail",
                            gap=reduced_next_tail_reciprocal_second_gap,
                        )
                    )
                lines.extend(
                    [
                        "```",
                        "",
                        "```text",
                        "G2_tail = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                        "H2_tail = T * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                        "```",
                        "",
                    ]
                )
                for target_label, target_variable, target_scans in (
                    ("`G2_tail`", "G2_tail", reduced_tail_anchor_second_gap_source_family_eta_scans),
                    ("`H2_tail`", "H2_tail", reduced_next_tail_reciprocal_second_gap_source_family_eta_scans),
                ):
                    if not target_scans:
                        continue
                    target_hits = _flatten_source_family_eta_hits(target_scans)
                    if not target_hits:
                        lines.append(
                            f"- {target_label}: no one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`."
                        )
                        continue
                    lines.append(f"- {target_label}: source-family eta-correction hits were found:")
                    for basis_label, basis_expression, basis_kind, level, relation in target_hits:
                        lines.append(
                            f"  - {basis_kind} basis `{basis_label}`, `N={level}`: `{_format_source_family_eta_correction(basis_expression=basis_expression, relation=relation, target_variable=target_variable, series_symbol=series_symbol)}`"
                        )
                lines.append("")
                lines.append(
                    "Second-gap compression: even after stripping a second visible gap from `G_tail` and `H_tail`, the one-core source-family eta-correction lane still stays empty in the scanned nearby basis ladders."
                )
                lines.append("")

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

        reduced_ratio_self_plus_pochhammer_hits = [
            scan for scan in reduced_ratio_self_plus_pochhammer_scans if scan.relation is not None
        ]
        if not reduced_ratio_self_plus_pochhammer_hits:
            lines.append("No reduced-ratio self plus-Pochhammer transfer hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self plus-Pochhammer transfer hits were found:")
            for scan in reduced_ratio_self_plus_pochhammer_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `m={scan.modulus}`: `{_format_self_plus_pochhammer_relation(scan.relation, modulus=scan.modulus, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        reduced_ratio_self_plus_pochhammer_eta_hits = [
            scan for scan in reduced_ratio_self_plus_pochhammer_eta_scans if scan.relation is not None
        ]
        if not reduced_ratio_self_plus_pochhammer_eta_hits:
            lines.append("No reduced-ratio self plus-Pochhammer + eta transfer hit was found in the scanned box.")
        else:
            lines.append("Reduced-ratio self plus-Pochhammer + eta transfer hits were found:")
            for scan in reduced_ratio_self_plus_pochhammer_eta_scans:
                if scan.relation is None:
                    continue
                lines.append(
                    f"- `m={scan.modulus}`, `N={scan.level}`: `{_format_self_plus_pochhammer_eta_relation(scan.relation, modulus=scan.modulus, target_variable='F_red', series_symbol=series_symbol)}`"
                )
        lines.append("")

        if reduced_ratio_source_family_self_plus_pochhammer_eta_scans:
            lines.extend(
                [
                    "We also checked whether removing one nearby source-family core reveals that same reverse-scale-aware mixed transfer box on the reduced ratio side:",
                    "",
                    "```text",
                    "G_red = F_red / T",
                    "G_red = G_red(t^m)^a * prod_r (-t^r; t^m)_inf^{e_r} * eta_tail",
                    "```",
                    "",
                ]
            )

            reduced_ratio_source_family_self_plus_pochhammer_eta_hits: list[
                tuple[str, str, str, int, int, MultiplicativeRelation]
            ] = []
            for family_scan in reduced_ratio_source_family_self_plus_pochhammer_eta_scans:
                for basis_scan in family_scan.direct_basis_scans:
                    for self_scan in basis_scan.self_scans:
                        if self_scan.relation is not None:
                            reduced_ratio_source_family_self_plus_pochhammer_eta_hits.append(
                                (
                                    basis_scan.basis_label,
                                    basis_scan.basis_expression,
                                    "raw",
                                    self_scan.modulus,
                                    self_scan.level,
                                    self_scan.relation,
                                )
                            )
                for basis_scan in family_scan.quotient_basis_scans:
                    for self_scan in basis_scan.self_scans:
                        if self_scan.relation is not None:
                            reduced_ratio_source_family_self_plus_pochhammer_eta_hits.append(
                                (
                                    basis_scan.basis_label,
                                    basis_scan.basis_expression,
                                    "quotient",
                                    self_scan.modulus,
                                    self_scan.level,
                                    self_scan.relation,
                                )
                            )

            if not reduced_ratio_source_family_self_plus_pochhammer_eta_hits:
                lines.append(
                    "- No one-core reduced-ratio self plus-Pochhammer + eta hit was found across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`."
                )
            else:
                lines.append("- One-core reduced-ratio self plus-Pochhammer + eta hits were found:")
                for basis_label, basis_expression, basis_kind, modulus, level, relation in (
                    reduced_ratio_source_family_self_plus_pochhammer_eta_hits
                ):
                    lines.append(
                        f"  - {basis_kind} basis `{basis_label}`, `m={modulus}`, `N={level}`:"
                        f" `G_red = F_red / ({basis_expression}); {_format_self_plus_pochhammer_eta_relation(relation, modulus=modulus, target_variable='G_red', series_symbol=series_symbol)}`"
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
            if gg_modular_equation_scan.weighted_coordinate_diagnostics:
                lines.append(
                    "- Weighted quotient-coordinate diagnostics now also test the first `3:2` Weber-style compression suggested by the `Q_3/Q_4` obstruction classes."
                )
                for diagnostic in gg_modular_equation_scan.weighted_coordinate_diagnostics:
                    eta_hit_levels = [
                        f"`N={scan.level}`"
                        for scan in diagnostic.correction_eta_scans
                        if scan.relation is not None
                    ]
                    modular_hit_boxes = [
                        f"`m={scan.modulus}, N={scan.level}`"
                        for scan in diagnostic.correction_modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    normalized_eta_hit_levels = [
                        f"`N={scan.level}`"
                        for scan in diagnostic.normalized_correction_eta_scans
                        if scan.relation is not None
                    ]
                    normalized_modular_hit_boxes = [
                        f"`m={scan.modulus}, N={scan.level}`"
                        for scan in diagnostic.normalized_correction_modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    normalized_source_hits = _flatten_source_family_eta_hits(
                        diagnostic.normalized_correction_source_family_eta_scans
                    )
                    normalized_source_family_labels = ", ".join(
                        f"`{family_scan.family_label}`"
                        for family_scan in diagnostic.normalized_correction_source_family_eta_scans
                    )
                    second_normalized_eta_hit_levels = [
                        f"`N={scan.level}`"
                        for scan in diagnostic.second_normalized_correction_eta_scans
                        if scan.relation is not None
                    ]
                    second_normalized_modular_hit_boxes = [
                        f"`m={scan.modulus}, N={scan.level}`"
                        for scan in diagnostic.second_normalized_correction_modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    second_normalized_source_hits = _flatten_source_family_eta_hits(
                        diagnostic.second_normalized_correction_source_family_eta_scans
                    )
                    second_normalized_source_family_labels = ", ".join(
                        f"`{family_scan.family_label}`"
                        for family_scan in diagnostic.second_normalized_correction_source_family_eta_scans
                    )
                    second_normalized_explicit_eta_hit_count = sum(
                        len(scan.hits)
                        for scan in diagnostic.second_normalized_correction_explicit_transform_eta_scans
                    )
                    second_normalized_quotient_polynomial_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_quotient_polynomial_scans,
                        label="polynomial",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_quotient_multiplicative_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_quotient_multiplicative_scans,
                        label="multiplicative",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_quotient_fractional_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_quotient_fractional_linear_scans,
                        label="fractional-linear",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_quotient_two_layer_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_quotient_two_layer_fractional_linear_scans,
                        label="two-layer fractional-linear",
                        hit_predicate=lambda item: item.total_hits > 0,
                    )
                    second_normalized_mixed_polynomial_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_mixed_quotient_polynomial_scans,
                        label="polynomial",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_mixed_multiplicative_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_mixed_quotient_multiplicative_scans,
                        label="multiplicative",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_mixed_fractional_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_mixed_quotient_fractional_linear_scans,
                        label="fractional-linear",
                        hit_predicate=lambda item: item.relation is not None,
                    )
                    second_normalized_mixed_two_layer_summary = _format_tail_prefix_summary(
                        diagnostic.second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans,
                        label="two-layer fractional-linear",
                        hit_predicate=lambda item: item.total_hits > 0,
                    )
                    lines.append(
                        f"- Weighted coordinate `{diagnostic.label} = {diagnostic.expression}`: "
                        + _format_weighted_coordinate_obstruction(
                            power=diagnostic.first_difference_power,
                            coeff=diagnostic.first_difference_coeff,
                            lhs=f"F - {diagnostic.label}",
                            series_symbol=series_symbol,
                        )
                        + "; "
                        + _format_weighted_coordinate_obstruction(
                            power=diagnostic.first_log_difference_power,
                            coeff=diagnostic.first_log_difference_coeff,
                            lhs=f"log(F) - ({diagnostic.log_expression})",
                            series_symbol=series_symbol,
                        )
                        + "."
                    )
                    if diagnostic.polynomial_degree1_relation is None:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: no candidate-dependent polynomial relation of total degree `<= 1` was found."
                        )
                    else:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: a candidate-dependent polynomial relation of total degree `<= 1` was found."
                        )
                    if diagnostic.polynomial_degree2_relation is None:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: no candidate-dependent polynomial relation of total degree `<= 2` was found."
                        )
                    else:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: a candidate-dependent polynomial relation of total degree `<= 2` was found."
                        )
                    if diagnostic.fractional_linear_relation is None:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: no one-coordinate fractional-linear closure was found."
                        )
                    else:
                        lines.append(
                            f"- Weighted coordinate `{diagnostic.label}`: a one-coordinate fractional-linear closure was found."
                        )
                    lines.append(
                        f"- Weighted correction `{diagnostic.correction_expression}`: "
                        + _format_weighted_coordinate_obstruction(
                            power=diagnostic.correction_first_gap_power,
                            coeff=diagnostic.correction_first_gap_coeff,
                            lhs=f"{diagnostic.correction_expression} - 1",
                            series_symbol=series_symbol,
                        )
                        + "."
                    )
                    if eta_hit_levels:
                        lines.append(
                            f"- Weighted correction `{diagnostic.correction_expression}`: eta-quotient hit(s) found at {', '.join(eta_hit_levels)}."
                        )
                    else:
                        lines.append(
                            f"- Weighted correction `{diagnostic.correction_expression}`: no eta-quotient hit was found in the checked small levels."
                        )
                    if modular_hit_boxes:
                        lines.append(
                            f"- Weighted correction `{diagnostic.correction_expression}`: modular-unit / eta hit(s) found at {', '.join(modular_hit_boxes)}."
                        )
                    else:
                        lines.append(
                            f"- Weighted correction `{diagnostic.correction_expression}`: no modular-unit / eta hit was found in the checked small boxes."
                        )
                    if diagnostic.normalized_correction_gap is not None and diagnostic.normalized_correction_label is not None:
                        lines.append(
                            f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: "
                            + _format_gap_normalization_formula(
                                source_variable=diagnostic.correction_expression,
                                target_variable=diagnostic.normalized_correction_label,
                                gap=diagnostic.normalized_correction_gap,
                                series_symbol=series_symbol,
                            )
                            + "."
                        )
                        if normalized_eta_hit_levels:
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: eta-quotient hit(s) found at {', '.join(normalized_eta_hit_levels)}."
                            )
                        else:
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: no eta-quotient hit was found in the checked small levels."
                            )
                        if normalized_modular_hit_boxes:
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: modular-unit / eta hit(s) found at {', '.join(normalized_modular_hit_boxes)}."
                            )
                        else:
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: no modular-unit / eta hit was found in the checked small boxes."
                            )
                        if normalized_source_hits:
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: source-family eta-correction hit(s) were found:"
                            )
                            for basis_label, basis_expression, basis_kind, level, relation in normalized_source_hits:
                                lines.append(
                                    "  - "
                                    + f"`{basis_kind}` basis `{basis_label}` at `N={level}`: "
                                    + _format_source_family_eta_correction(
                                        basis_expression=basis_expression,
                                        relation=relation,
                                        target_variable=diagnostic.normalized_correction_label,
                                        series_symbol=series_symbol,
                                    )
                                )
                        else:
                            family_text = normalized_source_family_labels or "the checked source-family"
                            lines.append(
                                f"- Normalized weighted correction `{diagnostic.normalized_correction_label}`: no one-core source-family eta-correction hit was found across {family_text} raw/quotient bases."
                            )
                    if (
                        diagnostic.second_normalized_correction_gap is not None
                        and diagnostic.second_normalized_correction_label is not None
                        and diagnostic.normalized_correction_label is not None
                    ):
                        lines.append(
                            f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: "
                            + _format_gap_normalization_formula(
                                source_variable=diagnostic.normalized_correction_label,
                                target_variable=diagnostic.second_normalized_correction_label,
                                gap=diagnostic.second_normalized_correction_gap,
                                series_symbol=series_symbol,
                            )
                            + "."
                        )
                        if second_normalized_eta_hit_levels:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: eta-quotient hit(s) found at {', '.join(second_normalized_eta_hit_levels)}."
                            )
                        else:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no eta-quotient hit was found in the checked small levels."
                            )
                        if second_normalized_modular_hit_boxes:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: modular-unit / eta hit(s) found at {', '.join(second_normalized_modular_hit_boxes)}."
                            )
                        else:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no modular-unit / eta hit was found in the checked small boxes."
                            )
                        if second_normalized_source_hits:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: source-family eta-correction hit(s) were found:"
                            )
                            for basis_label, basis_expression, basis_kind, level, relation in second_normalized_source_hits:
                                lines.append(
                                    "  - "
                                    + f"`{basis_kind}` basis `{basis_label}` at `N={level}`: "
                                    + _format_source_family_eta_correction(
                                        basis_expression=basis_expression,
                                        relation=relation,
                                        target_variable=diagnostic.second_normalized_correction_label,
                                        series_symbol=series_symbol,
                                    )
                                )
                        else:
                            family_text = second_normalized_source_family_labels or "the checked source-family"
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no one-core source-family eta-correction hit was found across {family_text} raw/quotient bases."
                            )
                        if second_normalized_explicit_eta_hit_count:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: explicit GG transform-template eta-correction hit(s) were found in `{second_normalized_explicit_eta_hit_count}` checked boxes."
                            )
                        else:
                            lines.append(
                                f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no explicit GG transform-template eta-correction hit was found in the checked small boxes."
                            )
                        lines.append(
                            f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}` quotient-coordinate prefixes: "
                            + "; ".join(
                                (
                                    second_normalized_quotient_polynomial_summary,
                                    second_normalized_quotient_multiplicative_summary,
                                    second_normalized_quotient_fractional_summary,
                                    second_normalized_quotient_two_layer_summary,
                                )
                            )
                            + "."
                        )
                        lines.append(
                            f"- Second normalized weighted correction `{diagnostic.second_normalized_correction_label}` mixed quotient-coordinate prefixes: "
                            + "; ".join(
                                (
                                    second_normalized_mixed_polynomial_summary,
                                    second_normalized_mixed_multiplicative_summary,
                                    second_normalized_mixed_fractional_summary,
                                    second_normalized_mixed_two_layer_summary,
                                )
                            )
                            + "."
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

    if ratio_modular_unit_eta_scans:
        lines.extend(
            [
                "",
                "## Ratio-Object Modular-Unit / Eta Scan",
                "",
                "We also checked whether the ratio object collapses into a small modular-unit / eta expression before introducing any self-copy lane:",
                "",
                "```text",
                "F = prod_r (1 - t^r)^{a_r} * prod_r (1 + t^r)^{b_r} * prod_{d|N} (t^d; t^d)_inf^{e_d}",
                "```",
                "",
                f"- `F = candidate / {record.closest_benchmark}`",
                "- This is a source-faithful modular-function recognition lane: a hit would be closer to an eta-quotient / modular-unit closed form than to another anonymous prefix fit.",
                "",
            ]
        )

        any_modular_unit_hit = any(scan.relation is not None for scan in ratio_modular_unit_eta_scans)
        if not any_modular_unit_hit:
            lines.append("No modular-unit / eta relation was found in any scanned modulus/level box.")
            lines.append("")

        no_hit_labels = [
            f"`m={scan.modulus}, N={scan.level}`"
            for scan in ratio_modular_unit_eta_scans
            if scan.error is None and scan.relation is None
        ]
        if not any_modular_unit_hit and no_hit_labels:
            lines.append(f"- No hit for modular-unit boxes {', '.join(no_hit_labels)}.")

        for scan in ratio_modular_unit_eta_scans:
            if scan.error is not None:
                lines.append(
                    f"- Modular-unit / eta box `m={scan.modulus}`, `N={scan.level}` skipped: {scan.error}"
                )
                continue
            if scan.relation is None:
                continue
            residual = _multiplicative_relation_residual_series(
                scan.relation,
                target_series=ratio_series,
                basis_series_by_variable=_modular_unit_eta_basis_series(
                    modulus=scan.modulus,
                    level=scan.level,
                    order=profile_order,
                ),
                order=profile_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Modular-unit / eta box `m={scan.modulus}`, `N={scan.level}` produced a candidate relation:",
                    "",
                    "```text",
                    _format_modular_unit_eta_relation(
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


def build_candidate_tail_family_note(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    depth: int = 40,
    series_order: int = 36,
    tail_stages: tuple[int, ...] = (3, 4, 5),
    max_gap_depth: int = 3,
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
    profile_depth = min(depth, 24 if smoke else depth)
    profile_order = min(series_order, 24 if smoke else series_order)
    active = _series_active_exponents(record.template) + _series_active_exponents(benchmark.canonical_template)
    step = 0
    for value in active:
        step = gcd(step, value)
    if step <= 0:
        step = 1

    reduced_candidate = record.template
    series_symbol = "q"
    variable_label = "q"
    if step > 1:
        maybe_candidate = reduce_template_by_step(record.template, step)
        maybe_benchmark = reduce_template_by_step(benchmark.canonical_template, step)
        if maybe_candidate is not None and maybe_benchmark is not None:
            reduced_candidate = maybe_candidate
            series_symbol = "t"
            variable_label = f"t = q^{step}"

    reduced_bridge_depth = min(profile_depth, 8 if smoke else 12)
    reduced_bridge_order = min(profile_order, 24 if smoke else 30)
    output_file = Path(output_path)
    build_started_at = perf_counter()
    lines = [
        f"# Tail-Family Note: `{record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{record.id}`",
        f"- Closest benchmark: `{record.closest_benchmark}`",
        f"- Variable view: `{variable_label}`",
        f"- Tail stages checked: `{', '.join(str(stage) for stage in tail_stages)}`",
        f"- Max gap depth checked: `{max_gap_depth}`",
        "",
    ]

    reduced_bridge_error: str | None = None
    reduced_reciprocal_witness = None
    reduced_tail_transfer_equation: ReducedTailTransferEquation | None = None
    try:
        reduced_reciprocal_witness, _ = _reduced_reciprocal_bridge(
            template=reduced_candidate,
            symbol=sp.Symbol(series_symbol),
            depth=reduced_bridge_depth,
            order=reduced_bridge_order,
        )
        reduced_tail_transfer_equation = detect_reduced_tail_transfer_equation(
            reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
            symbol=sp.Symbol(series_symbol),
        )
    except Exception as exc:
        reduced_bridge_error = str(exc)

    if reduced_bridge_error is not None or reduced_reciprocal_witness is None or reduced_tail_transfer_equation is None:
        lines.extend(
            [
                "Tail-family setup failed:",
                "",
                "```text",
                reduced_bridge_error or "unknown reduced-tail setup failure",
                "```",
                "",
            ]
        )
        output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    relation_lines = _format_reduced_tail_transfer_equation(reduced_tail_transfer_equation)
    lines.extend(
        [
            "## Exact Tail Family",
            "",
            "- From stage `3` onward, the reduced coefficients collapse into one stationary tail family.",
            "",
            "```text",
            *relation_lines,
            "```",
            "",
            "## Variable-Level Source-Core Recognition Lane",
            "",
            "We now treat the normalized tail family itself as the main intermediate object:",
            "",
            "```text",
            "U(x) = T(x) / (1 + x)",
            "```",
            "",
            "- For each sampled state `x = t^k`, we scan `U(x)` and its repeated gap-normalized residuals.",
            "- Unlike the broader `identify` note, this lane keeps only the source-driven one-core eta-correction question:",
            "",
            "```text",
            "Y = S * prod_{d|N} (t^d; t^d)_inf^{e_d}",
            "```",
            "",
            "- We also keep two more source-faithful direct recognition lanes on the same sampled objects:",
            "",
            "```text",
            "Y = prod_r (1 - t^r)^{a_r} * prod_r (1 + t^r)^{b_r} * prod_{d|N} (t^d; t^d)_inf^{e_d}",
            "f(Y, Y_2) = 0,   g(Y^2, Y_2^2) = 0,   f(Y, (1-Y_2)/(1+Y_2)) = 0",
            "```",
            "",
            "- The first is a direct modular-unit / eta lane.",
            "- The second is a Morton-2024-inspired periodic-point / algebraic-function lane built from the exact low-degree polynomials attached to the GG/Weber orbit.",
            "- Phase 2 now also isolates the source-faithful squared coordinate from Morton Prop. 3.2 / Theorem B:",
            "",
            "```text",
            "X_mt = Y^2",
            "X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 = 0",
            "T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1), sigma = -1 + sqrt(2)",
            "T_mt^2 - (T_mt,2^2 - 4*T_mt,2 + 1)*T_mt + T_mt,2^2 = 0",
            "```",
            "",
            "- Phase 2 now also opens the first deeper Weber-Schlafli coordinate lane on the same sampled objects:",
            "",
            "```text",
            "P_ws = (1/Y - Y) / 2",
            "B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)",
            "P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0",
            "B_ws^2 - B_ws,2 - 4 = 0",
            "B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0",
            "```",
            "",
            "- Phase 3 now adds the first Ramanujan-Weber class-invariant compression on the same normalized `GG` variable:",
            "",
            "```text",
            "Z_g = ((1 - t*Y^2)^2) / (4*t*Y^2)",
            "g12_ws = 4*t*(Z_g - 1/Z_g)",
            "p12_ws = (4*t*Z_g^2) / sqrt(Z_g^2 - 1)",
            "g12_ws ?= (t^2; t^4)_inf^12",
            "p12_ws ?= (-t^2; t^4)_inf^12",
            "R_gp_ws = G_p12_ws / G_g12_ws",
            "g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0",
            "```",
            "",
            "- Phase 4 now treats the eta-side residual `G_g12_ws` as the primary Weber hand-off, "
            "uses the exact algebraic bridge between `g12_ws` and `p12_ws` to keep `G_p12_ws` as a constrained companion, "
            "and then probes the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws` before reopening any broader search box.",
            "",
        ]
    )

    source_family_basis_catalog = _source_family_basis_catalog(record.closest_benchmark)
    source_family_base_series = tuple(
        (
            label,
            benchmark_name,
            _canonical_benchmark_series(
                benchmark_name,
                depth=profile_depth,
                order=reduced_bridge_order,
            ),
        )
        for label, benchmark_name in source_family_basis_catalog
    )
    gg_base_family_entry = next(
        (entry for entry in source_family_base_series if entry[0] == "GG"),
        None,
    )
    sample_scans = scan_tail_family_source_eta_ladder(
        reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
        symbol=sp.Symbol(series_symbol),
        ordered_base_families=source_family_base_series,
        start_stages=tail_stages,
        max_gap_depth=max_gap_depth,
        order=reduced_bridge_order,
        powers=(2, 3, 4),
        eta_levels=_eta_scan_levels((1, 2, 3, 4)),
        max_abs_exponent=4 if smoke else 6,
        gg_order=min(reduced_bridge_order, 20),
        gg_benchmark_name=None if gg_base_family_entry is None else gg_base_family_entry[1],
        gg_base_series=None if gg_base_family_entry is None else gg_base_family_entry[2],
        gg_degree_values=(1, 2),
        gg_max_abs_exponent=4 if smoke else 6,
        gg_solve_order=min(reduced_bridge_order, 14 if smoke else 18),
        gg_supplemental_powers=(),
        morton_order=min(reduced_bridge_order, 18 if smoke else 24),
    )

    if not sample_scans:
        lines.extend(
            [
                "No tail-family samples were available for the requested stage/gap profile.",
                "",
            ]
        )
    else:
        gg_lane_scan = next(
            (
                sample.gg_modular_equation_scan
                for sample in sample_scans
                if sample.gg_modular_equation_scan is not None
            ),
            None,
        )
        if gg_lane_scan is not None:
            gg_basis_descriptions = [
                f"`{label} = {expression}`"
                for label, expression, _ in gg_lane_scan.ordered_basis_series
            ]
            gg_quotient_descriptions = [
                f"`{label} = {expression}`"
                for label, expression, _ in gg_lane_scan.quotient_basis_series
            ]
            lines.extend(
                [
                    "- We also check a `GG/Weber modular-equation` lane on the same sampled `U(x)` objects and their gap residuals:",
                    "",
                    "```text",
                    "Y = T",
                    "Y = 1 / T",
                    "Y = T_i / T_j",
                    "P(Y, T_i) = 0",
                    "Y = prod_i T_i^e_i",
                    "Y = (1 + sum a_i*(T_i - 1)) / (1 + sum b_i*(T_i - 1))",
                    "```",
                    "",
                    f"- GG base benchmark: `{gg_lane_scan.benchmark_name}`",
                    f"- GG basis ladder: {', '.join(gg_basis_descriptions)}",
                    f"- Preferred quotient coordinates: {', '.join(gg_quotient_descriptions)}",
                    "- The narrowest exact quotient-coordinate lane keeps special attention on `Q_3 = GG(t^3)/GG(t)` and `Q_4 = GG(t^4)/GG(t)`, because those are the Chan--Huang exact modular-equation coordinates.",
                    "- Each sample below reports direct, quotient-coordinate, and mixed quotient-coordinate prefix summaries for this literature-driven lane.",
                    "",
                ]
            )
        for sample in sample_scans:
            hits = _flatten_source_family_eta_hits(sample.source_family_eta_scans)
            lines.extend(
                [
                    f"### `{sample.label}`",
                    "",
                    f"- Start stage: `{sample.start_stage}`",
                    f"- Gap depth: `{sample.gap_depth}`",
                    f"- State: `{_format_expr(sample.state_expr)}`",
                    "",
                    "```text",
                    sample.expression,
                    "```",
                    "",
                ]
            )
            if not hits:
                lines.append(
                    "- No one-core source-family eta-correction hit across raw/quotient basis choices from `RR`, `cubic`, `GG`, `S` through powers `2,3,4`."
                )
            else:
                lines.append("- One-core source-family eta-correction hits were found:")
                for basis_label, basis_expression, basis_kind, level, relation in hits:
                    lines.append(
                        f"  - {basis_kind} basis `{basis_label}`, `N={level}`: `{_format_source_family_eta_correction(basis_expression=basis_expression, relation=relation, target_variable=sample.label, series_symbol=series_symbol)}`"
                    )
            direct_eta_hits = [scan for scan in sample.direct_eta_scans if scan.relation is not None]
            direct_modular_unit_hits = [
                scan for scan in sample.direct_modular_unit_eta_scans if scan.relation is not None
            ]
            lines.append(
                f"- Direct eta-quotient templates: `{len(direct_eta_hits)}` / `{len(sample.direct_eta_scans)}` hit boxes."
            )
            if direct_eta_hits:
                for scan in direct_eta_hits:
                    lines.append(
                        f"  - eta level `N={scan.level}`: `{_format_eta_quotient_relation(scan.relation, target_variable=sample.label, series_symbol=series_symbol)}`"
                    )
            lines.append(
                f"- Direct modular-unit / eta templates: `{len(direct_modular_unit_hits)}` / `{len(sample.direct_modular_unit_eta_scans)}` hit boxes."
            )
            if direct_modular_unit_hits:
                for scan in direct_modular_unit_hits:
                    lines.append(
                        f"  - modular box `m={scan.modulus}`, `N={scan.level}`: `{_format_modular_unit_eta_relation(scan.relation, target_variable=sample.label, series_symbol=series_symbol)}`"
                    )
            if sample.morton_periodic_point_scan is not None:
                morton_scan = sample.morton_periodic_point_scan
                morton_hits = [item for item in morton_scan.template_results if item.hit]
                lines.append(
                    f"- Morton periodic-point / algebraic-function templates: `{len(morton_hits)}` / `{len(morton_scan.template_results)}` exact hits."
                )
                if morton_scan.template_results:
                    lines.append(
                        "- Morton obstruction witnesses: "
                        + "; ".join(
                            _format_exact_polynomial_obstruction(
                                (item.label, item.first_failure_power, item.first_failure_coeff),
                                series_symbol=series_symbol,
                            )
                            for item in morton_scan.template_results
                        )
                        + "."
                    )
                for coordinate_scan in morton_scan.named_coordinate_scans:
                    weber_hits = [item for item in coordinate_scan.template_results if item.hit]
                    lines.append(
                        f"- Morton {coordinate_scan.family_label} coordinate `{coordinate_scan.label}`: `{coordinate_scan.expression}`."
                    )
                    lines.append(
                        f"- Morton {coordinate_scan.family_label} coordinate templates on `{coordinate_scan.label}`: `{len(weber_hits)}` / `{len(coordinate_scan.template_results)}` exact hits."
                    )
                    if coordinate_scan.template_results:
                        lines.append(
                            f"- Morton {coordinate_scan.family_label} coordinate obstruction witnesses: "
                            + "; ".join(
                                _format_exact_polynomial_obstruction(
                                    (item.label, item.first_failure_power, item.first_failure_coeff),
                                    series_symbol=series_symbol,
                                )
                                for item in coordinate_scan.template_results
                            )
                            + "."
                        )
                    if coordinate_scan.leading_normalized_scan is not None:
                        _append_constant_one_scan_summary_lines(
                            lines,
                            prefix=(
                                f"Morton {coordinate_scan.family_label} leading-term-normalized coordinate "
                                f"`{coordinate_scan.leading_normalized_scan.label}`"
                            ),
                            scan=coordinate_scan.leading_normalized_scan.scan,
                            series_symbol=series_symbol,
                        )
                for bridge_scan in morton_scan.leading_normalized_bridge_scans:
                    bridge_polynomial_hits = [
                        item for item in bridge_scan.polynomial_scans if item.relation is not None
                    ]
                    lines.append(
                        f"- Morton Weber-Schlafli leading-normalized bridge difference `{bridge_scan.difference_label}`: "
                        f"`{bridge_scan.difference_expression}`."
                    )
                    if (
                        bridge_scan.difference_first_failure_power is None
                        or bridge_scan.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized bridge difference `{bridge_scan.difference_label}`: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized bridge difference `{bridge_scan.difference_label}`: "
                            f"first fails at `{series_symbol}^{bridge_scan.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(bridge_scan.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Morton Weber-Schlafli leading-normalized bridge quotient `{bridge_scan.quotient_label}`: "
                        f"`{bridge_scan.quotient_expression}`."
                    )
                    if (
                        bridge_scan.quotient_first_failure_power is None
                        or bridge_scan.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized bridge quotient `{bridge_scan.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized bridge quotient `{bridge_scan.quotient_label}`: "
                            f"`{bridge_scan.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{bridge_scan.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(bridge_scan.quotient_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Morton Weber-Schlafli leading-normalized bridge polynomial boxes: "
                        f"`{len(bridge_polynomial_hits)}` / `{len(bridge_scan.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        "- Morton Weber-Schlafli leading-normalized bridge fractional-linear box: "
                        + (
                            "`1` / `1` hit boxes."
                            if bridge_scan.fractional_linear_relation is not None
                            else "`0` / `1` hit boxes."
                        )
                    )
                    _append_constant_one_scan_summary_lines(
                        lines,
                        prefix=(
                            f"Morton Weber-Schlafli leading-normalized bridge quotient "
                            f"`{bridge_scan.quotient_scan.label}`"
                        ),
                        scan=bridge_scan.quotient_scan,
                        series_symbol=series_symbol,
                    )
                    if bridge_scan.quotient_named_coordinate_orbit_scan is not None:
                        _append_named_coordinate_orbit_lines(
                            lines,
                            prefix=(
                                f"Morton Weber-Schlafli leading-normalized bridge quotient "
                                f"`{bridge_scan.quotient_scan.label}`"
                            ),
                            orbit_scan=bridge_scan.quotient_named_coordinate_orbit_scan,
                        )
                    if bridge_scan.quotient_scan.normalized_followup is not None:
                        _append_constant_one_scan_summary_lines(
                            lines,
                            prefix=(
                                "Morton Weber-Schlafli leading-normalized bridge quotient normalized follow-up "
                                f"`{bridge_scan.quotient_scan.normalized_followup.label}`"
                            ),
                            scan=ConstantOneSeriesScan(
                                label=bridge_scan.quotient_scan.normalized_followup.label,
                                expression=bridge_scan.quotient_scan.normalized_followup.expression,
                                first_failure_power=bridge_scan.quotient_scan.normalized_followup.first_failure_power,
                                first_failure_coeff=bridge_scan.quotient_scan.normalized_followup.first_failure_coeff,
                                self_polynomial_scan=bridge_scan.quotient_scan.normalized_followup.self_polynomial_scan,
                                self_fractional_linear_scan=bridge_scan.quotient_scan.normalized_followup.self_fractional_linear_scan,
                                self_quotient_product_scans=bridge_scan.quotient_scan.normalized_followup.self_quotient_product_scans,
                                eta_scans=bridge_scan.quotient_scan.normalized_followup.eta_scans,
                                modular_unit_eta_scans=bridge_scan.quotient_scan.normalized_followup.modular_unit_eta_scans,
                                self_plus_pochhammer_scans=bridge_scan.quotient_scan.normalized_followup.self_plus_pochhammer_scans,
                                self_plus_pochhammer_eta_scans=bridge_scan.quotient_scan.normalized_followup.self_plus_pochhammer_eta_scans,
                                named_gg_modular_equation_scan=bridge_scan.quotient_scan.normalized_followup.named_gg_modular_equation_scan,
                                normalized_followup=None,
                                ),
                                series_symbol=series_symbol,
                            )
                        if bridge_scan.quotient_followup_named_coordinate_orbit_scan is not None:
                            _append_named_coordinate_orbit_lines(
                                lines,
                                prefix=(
                                    "Morton Weber-Schlafli leading-normalized bridge quotient normalized follow-up "
                                    f"`{bridge_scan.quotient_scan.normalized_followup.label}`"
                                ),
                                orbit_scan=bridge_scan.quotient_followup_named_coordinate_orbit_scan,
                            )
                    if bridge_scan.quotient_followup_bridge_scan is not None:
                        nested_bridge = bridge_scan.quotient_followup_bridge_scan
                        nested_bridge_polynomial_hits = [
                            item for item in nested_bridge.polynomial_scans if item.relation is not None
                        ]
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized nested bridge difference `{nested_bridge.difference_label}`: "
                            f"`{nested_bridge.difference_expression}`."
                        )
                        if (
                            nested_bridge.difference_first_failure_power is None
                            or nested_bridge.difference_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Morton Weber-Schlafli leading-normalized nested bridge difference `{nested_bridge.difference_label}`: "
                                "matches `0` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Morton Weber-Schlafli leading-normalized nested bridge difference `{nested_bridge.difference_label}`: "
                                f"first fails at `{series_symbol}^{nested_bridge.difference_first_failure_power}` "
                                f"with coefficient `{_format_expr(nested_bridge.difference_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized nested bridge quotient `{nested_bridge.quotient_label}`: "
                            f"`{nested_bridge.quotient_expression}`."
                        )
                        if (
                            nested_bridge.quotient_first_failure_power is None
                            or nested_bridge.quotient_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Morton Weber-Schlafli leading-normalized nested bridge quotient `{nested_bridge.quotient_label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Morton Weber-Schlafli leading-normalized nested bridge quotient `{nested_bridge.quotient_label}`: "
                                f"`{nested_bridge.quotient_label} - 1` first fails at "
                                f"`{series_symbol}^{nested_bridge.quotient_first_failure_power}` with coefficient "
                                f"`{_format_expr(nested_bridge.quotient_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Morton Weber-Schlafli leading-normalized nested bridge polynomial boxes: "
                            f"`{len(nested_bridge_polynomial_hits)}` / `{len(nested_bridge.polynomial_scans)}` hit boxes."
                        )
                        lines.append(
                            "- Morton Weber-Schlafli leading-normalized nested bridge fractional-linear box: "
                            + (
                                "`1` / `1` hit boxes."
                                if nested_bridge.fractional_linear_relation is not None
                                else "`0` / `1` hit boxes."
                            )
                        )
                        _append_constant_one_scan_summary_lines(
                            lines,
                            prefix=(
                                f"Morton Weber-Schlafli leading-normalized nested bridge quotient "
                                f"`{nested_bridge.quotient_scan.label}`"
                            ),
                            scan=nested_bridge.quotient_scan,
                            series_symbol=series_symbol,
                        )
                        if nested_bridge.quotient_named_coordinate_orbit_scan is not None:
                            _append_named_coordinate_orbit_lines(
                                lines,
                                prefix=(
                                    f"Morton Weber-Schlafli leading-normalized nested bridge quotient "
                                    f"`{nested_bridge.quotient_scan.label}`"
                                ),
                                orbit_scan=nested_bridge.quotient_named_coordinate_orbit_scan,
                            )
                        if nested_bridge.quotient_scan.normalized_followup is not None:
                            _append_constant_one_scan_summary_lines(
                                lines,
                                prefix=(
                                    "Morton Weber-Schlafli leading-normalized nested bridge quotient normalized follow-up "
                                    f"`{nested_bridge.quotient_scan.normalized_followup.label}`"
                                ),
                                scan=ConstantOneSeriesScan(
                                    label=nested_bridge.quotient_scan.normalized_followup.label,
                                    expression=nested_bridge.quotient_scan.normalized_followup.expression,
                                    first_failure_power=nested_bridge.quotient_scan.normalized_followup.first_failure_power,
                                    first_failure_coeff=nested_bridge.quotient_scan.normalized_followup.first_failure_coeff,
                                    self_polynomial_scan=nested_bridge.quotient_scan.normalized_followup.self_polynomial_scan,
                                    self_fractional_linear_scan=nested_bridge.quotient_scan.normalized_followup.self_fractional_linear_scan,
                                    self_quotient_product_scans=nested_bridge.quotient_scan.normalized_followup.self_quotient_product_scans,
                                    eta_scans=nested_bridge.quotient_scan.normalized_followup.eta_scans,
                                    modular_unit_eta_scans=nested_bridge.quotient_scan.normalized_followup.modular_unit_eta_scans,
                                    self_plus_pochhammer_scans=nested_bridge.quotient_scan.normalized_followup.self_plus_pochhammer_scans,
                                    self_plus_pochhammer_eta_scans=nested_bridge.quotient_scan.normalized_followup.self_plus_pochhammer_eta_scans,
                                    named_gg_modular_equation_scan=nested_bridge.quotient_scan.normalized_followup.named_gg_modular_equation_scan,
                                    normalized_followup=None,
                                    ),
                                    series_symbol=series_symbol,
                                )
                            if nested_bridge.quotient_followup_named_coordinate_orbit_scan is not None:
                                _append_named_coordinate_orbit_lines(
                                    lines,
                                    prefix=(
                                        "Morton Weber-Schlafli leading-normalized nested bridge quotient normalized follow-up "
                                        f"`{nested_bridge.quotient_scan.normalized_followup.label}`"
                                    ),
                                    orbit_scan=nested_bridge.quotient_followup_named_coordinate_orbit_scan,
                                )
            for weber_scan in (
                sample.weber_g_class_invariant_scan,
                sample.weber_p_class_invariant_scan,
            ):
                if weber_scan is None:
                    continue
                eta_hits = [scan for scan in weber_scan.direct_eta_scans if scan.relation is not None]
                modular_hits = [
                    scan for scan in weber_scan.direct_modular_unit_eta_scans if scan.relation is not None
                ]
                plus_hits = [
                    scan for scan in weber_scan.correction_self_plus_pochhammer_scans if scan.relation is not None
                ]
                plus_eta_hits = [
                    scan
                    for scan in weber_scan.correction_self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Weber class-invariant coordinate `{weber_scan.label}`: `{weber_scan.expression}`."
                )
                lines.append(
                    f"- Weber class-invariant template on `{weber_scan.label}`: "
                    f"`{weber_scan.template_expression}`."
                )
                lines.append(
                    "- Weber class-invariant obstruction witness: "
                    + _format_weber_class_invariant_obstruction(
                        weber_scan,
                        series_symbol=series_symbol,
                    )
                    + "."
                )
                lines.append(
                    f"- Weber class-invariant correction `{weber_scan.correction_label}`: "
                    f"`{weber_scan.correction_expression}`."
                )
                lines.append(
                    f"- Weber class-invariant correction eta templates: `{len(eta_hits)}` / "
                    f"`{len(weber_scan.direct_eta_scans)}` hit boxes."
                )
                if eta_hits:
                    for scan in eta_hits:
                        lines.append(
                            f"  - eta level `N={scan.level}`: "
                            f"`{_format_eta_quotient_relation(scan.relation, target_variable=weber_scan.correction_label, series_symbol=series_symbol)}`"
                        )
                lines.append(
                    f"- Weber class-invariant correction modular-unit / eta templates: "
                    f"`{len(modular_hits)}` / `{len(weber_scan.direct_modular_unit_eta_scans)}` hit boxes."
                )
                if modular_hits:
                    for scan in modular_hits:
                        lines.append(
                            f"  - modular box `m={scan.modulus}`, `N={scan.level}`: "
                            f"`{_format_modular_unit_eta_relation(scan.relation, target_variable=weber_scan.correction_label, series_symbol=series_symbol)}`"
                        )
                lines.append(
                    f"- Weber class-invariant correction plus-Pochhammer templates: "
                    f"`{len(plus_hits)}` / `{len(weber_scan.correction_self_plus_pochhammer_scans)}` hit boxes."
                )
                if plus_hits:
                    for scan in plus_hits:
                        lines.append(
                            f"  - plus box `m={scan.modulus}`: "
                            f"`{_format_self_plus_pochhammer_relation(scan.relation, target_variable=weber_scan.correction_label, series_symbol=series_symbol)}`"
                        )
                lines.append(
                    f"- Weber class-invariant correction plus-Pochhammer + eta templates: "
                    f"`{len(plus_eta_hits)}` / `{len(weber_scan.correction_self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if plus_eta_hits:
                    for scan in plus_eta_hits:
                        lines.append(
                            f"  - plus box `m={scan.modulus}`, `N={scan.level}`: "
                            f"`{_format_self_plus_pochhammer_eta_relation(scan.relation, modulus=scan.modulus, target_variable=weber_scan.correction_label, series_symbol=series_symbol)}`"
                        )
            if sample.weber_residual_bridge_scan is not None:
                bridge_scan = sample.weber_residual_bridge_scan
                quotient_eta_hits = [
                    scan for scan in bridge_scan.quotient_eta_scans if scan.relation is not None
                ]
                quotient_modular_hits = [
                    scan
                    for scan in bridge_scan.quotient_modular_unit_eta_scans
                    if scan.relation is not None
                ]
                quotient_plus_hits = [
                    scan
                    for scan in bridge_scan.quotient_self_plus_pochhammer_scans
                    if scan.relation is not None
                ]
                quotient_plus_eta_hits = [
                    scan
                    for scan in bridge_scan.quotient_self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                quotient_self_product_hits = [
                    scan
                    for scan in bridge_scan.quotient_self_quotient_product_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Weber residual bridge keeps `{bridge_scan.primary_label}` as the current primary residual: "
                    f"`{bridge_scan.primary_expression}`."
                )
                lines.append(f"- Weber residual bridge reason: {bridge_scan.primary_reason}")
                if bridge_scan.primary_named_gg_modular_equation_scan is not None:
                    lines.append(
                        "- Focused source-faithful residual pass: the primary Weber residual is now also "
                        "pushed directly through the named Chan--Huang / Cho--Koo--Park `GG` modular-equation "
                        "basis before moving on to the quotient residual."
                    )
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=f"Weber primary residual `{bridge_scan.primary_label}`",
                        gg_scan=bridge_scan.primary_named_gg_modular_equation_scan,
                        series_symbol=series_symbol,
                    )
                lines.append(
                    f"- Weber residual companion `{bridge_scan.companion_label}`: "
                    f"`{bridge_scan.companion_expression}`."
                )
                lines.append(
                    f"- Weber residual exact coordinate bridge: `{bridge_scan.exact_bridge_expression}`."
                )
                if bridge_scan.exact_bridge_holds:
                    lines.append(
                        "- Weber residual exact coordinate bridge verdict: matches through the checked truncation."
                    )
                else:
                    lines.append(
                        "- Weber residual exact coordinate bridge verdict: "
                        f"first fails at `{series_symbol}^{bridge_scan.exact_bridge_first_failure_power}` "
                        f"with coefficient `{_format_expr(bridge_scan.exact_bridge_first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Weber residual exact residual bridge: `{bridge_scan.residual_bridge_expression}`."
                )
                classical_product_scan = bridge_scan.classical_product_coordinate_scan
                classical_product_eta_hits = [
                    scan for scan in classical_product_scan.eta_scans if scan.relation is not None
                ]
                classical_product_modular_hits = [
                    scan
                    for scan in classical_product_scan.modular_unit_eta_scans
                    if scan.relation is not None
                ]
                classical_product_plus_hits = [
                    scan
                    for scan in classical_product_scan.self_plus_pochhammer_scans
                    if scan.relation is not None
                ]
                classical_product_plus_eta_hits = [
                    scan
                    for scan in classical_product_scan.self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                classical_product_self_product_hits = [
                    scan
                    for scan in classical_product_scan.self_quotient_product_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Classical Weber `f2` tri-product coordinate `{bridge_scan.classical_product_coordinate_label}`: "
                    f"`{bridge_scan.classical_product_coordinate_expression}`."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product bridge: "
                    f"`{bridge_scan.classical_product_coordinate_bridge_expression}`."
                )
                lines.append(
                    "- Classical Weber source reading: Berndt--Chan--Zhang identify "
                    "Ramanujan-Weber `G_n` / `g_n` with Weber `f` / `f1`, and "
                    "Yui--Zagier supplies the classical Weber `f`, `f1`, `f2` trio, "
                    "so the current `g12_ws` / `p12_ws` / `G_f2_ws` shell should be "
                    "read as that named Weber trio in the project's normalization "
                    "rather than as an anonymous product gadget."
                )
                if (
                    classical_product_scan.first_failure_power is None
                    or classical_product_scan.first_failure_coeff is None
                ):
                    lines.append(
                        f"- Classical Weber `f2` tri-product coordinate `{classical_product_scan.label}`: "
                        "matches `1` through the checked truncation."
                    )
                else:
                    lines.append(
                        f"- Classical Weber `f2` tri-product coordinate `{classical_product_scan.label}`: "
                        f"`{classical_product_scan.label} - 1` first fails at "
                        f"`{series_symbol}^{classical_product_scan.first_failure_power}` with coefficient "
                        f"`{_format_expr(classical_product_scan.first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Classical Weber `f2` tri-product self-polynomial uniqueness boxes: "
                    f"`{len(classical_product_scan.self_polynomial_scan.hits)}` / "
                    f"`{len(classical_product_scan.self_polynomial_scan.moduli_checked) * len(classical_product_scan.self_polynomial_scan.fg_degree_values) * len(classical_product_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product self-fractional-linear uniqueness boxes: "
                    f"`{len(classical_product_scan.self_fractional_linear_scan.hits)}` / "
                    f"`{len(classical_product_scan.self_fractional_linear_scan.moduli_checked) * len(classical_product_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product self-quotient finite-product boxes: "
                    f"`{len(classical_product_self_product_hits)}` / "
                    f"`{len(classical_product_scan.self_quotient_product_scans)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product eta templates: `{len(classical_product_eta_hits)}` / "
                    f"`{len(classical_product_scan.eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product modular-unit / eta templates: "
                    f"`{len(classical_product_modular_hits)}` / `{len(classical_product_scan.modular_unit_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product plus-Pochhammer templates: "
                    f"`{len(classical_product_plus_hits)}` / `{len(classical_product_scan.self_plus_pochhammer_scans)}` hit boxes."
                )
                lines.append(
                    f"- Classical Weber `f2` tri-product plus-Pochhammer + eta templates: "
                    f"`{len(classical_product_plus_eta_hits)}` / "
                    f"`{len(classical_product_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if classical_product_scan.named_gg_modular_equation_scan is not None:
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=(
                            f"Classical Weber `f2` tri-product coordinate "
                            f"`{classical_product_scan.label}`"
                        ),
                        gg_scan=classical_product_scan.named_gg_modular_equation_scan,
                        series_symbol=series_symbol,
                    )
                if classical_product_scan.normalized_followup is not None:
                    classical_product_followup = classical_product_scan.normalized_followup
                    classical_product_followup_self_product_hits = [
                        scan
                        for scan in classical_product_followup.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    classical_product_followup_eta_hits = [
                        scan for scan in classical_product_followup.eta_scans if scan.relation is not None
                    ]
                    classical_product_followup_modular_hits = [
                        scan
                        for scan in classical_product_followup.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    classical_product_followup_plus_hits = [
                        scan
                        for scan in classical_product_followup.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    classical_product_followup_plus_eta_hits = [
                        scan
                        for scan in classical_product_followup.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized follow-up `{classical_product_followup.label}`: "
                        f"`{classical_product_followup.expression}`."
                    )
                    if (
                        classical_product_followup.first_failure_power is None
                        or classical_product_followup.first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Classical Weber `f2` tri-product normalized follow-up `{classical_product_followup.label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Classical Weber `f2` tri-product normalized follow-up `{classical_product_followup.label}`: "
                            f"`{classical_product_followup.label} - 1` first fails at "
                            f"`{series_symbol}^{classical_product_followup.first_failure_power}` with coefficient "
                            f"`{_format_expr(classical_product_followup.first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized self-polynomial uniqueness boxes: "
                        f"`{len(classical_product_followup.self_polynomial_scan.hits)}` / "
                        f"`{len(classical_product_followup.self_polynomial_scan.moduli_checked) * len(classical_product_followup.self_polynomial_scan.fg_degree_values) * len(classical_product_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized self-fractional-linear uniqueness boxes: "
                        f"`{len(classical_product_followup.self_fractional_linear_scan.hits)}` / "
                        f"`{len(classical_product_followup.self_fractional_linear_scan.moduli_checked) * len(classical_product_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized self-quotient finite-product boxes: "
                        f"`{len(classical_product_followup_self_product_hits)}` / "
                        f"`{len(classical_product_followup.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized eta templates: `{len(classical_product_followup_eta_hits)}` / "
                        f"`{len(classical_product_followup.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized modular-unit / eta templates: "
                        f"`{len(classical_product_followup_modular_hits)}` / `{len(classical_product_followup.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized plus-Pochhammer templates: "
                        f"`{len(classical_product_followup_plus_hits)}` / `{len(classical_product_followup.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Classical Weber `f2` tri-product normalized plus-Pochhammer + eta templates: "
                        f"`{len(classical_product_followup_plus_eta_hits)}` / "
                        f"`{len(classical_product_followup.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                if classical_product_followup.named_gg_modular_equation_scan is not None:
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=(
                            f"Classical Weber `f2` tri-product normalized follow-up "
                                f"`{classical_product_followup.label}`"
                            ),
                            gg_scan=classical_product_followup.named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                canonical_j_scan = bridge_scan.canonical_j_coordinate_scan
                canonical_j_eta_hits = [
                    scan for scan in canonical_j_scan.eta_scans if scan.relation is not None
                ]
                canonical_j_modular_hits = [
                    scan
                    for scan in canonical_j_scan.modular_unit_eta_scans
                    if scan.relation is not None
                ]
                canonical_j_plus_hits = [
                    scan
                    for scan in canonical_j_scan.self_plus_pochhammer_scans
                    if scan.relation is not None
                ]
                canonical_j_plus_eta_hits = [
                    scan
                    for scan in canonical_j_scan.self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                canonical_j_self_product_hits = [
                    scan
                    for scan in canonical_j_scan.self_quotient_product_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Canonical Weber `j`-side coordinate `{bridge_scan.canonical_j_coordinate_label}`: "
                    f"`{bridge_scan.canonical_j_coordinate_expression}`."
                )
                lines.append(
                    f"- Canonical Weber `j`-side reason: {bridge_scan.canonical_j_coordinate_reason}"
                )
                lines.append(
                    f"- Canonical Weber `j`-side bridge: "
                    f"`{bridge_scan.canonical_j_coordinate_bridge_expression}`."
                )
                if (
                    canonical_j_scan.first_failure_power is None
                    or canonical_j_scan.first_failure_coeff is None
                ):
                    lines.append(
                        f"- Canonical Weber `j`-side coordinate `{canonical_j_scan.label}`: "
                        "matches `1` through the checked truncation."
                    )
                else:
                    lines.append(
                        f"- Canonical Weber `j`-side coordinate `{canonical_j_scan.label}`: "
                        f"`{canonical_j_scan.label} - 1` first fails at "
                        f"`{series_symbol}^{canonical_j_scan.first_failure_power}` with coefficient "
                        f"`{_format_expr(canonical_j_scan.first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Canonical Weber `j`-side self-polynomial uniqueness boxes: "
                    f"`{len(canonical_j_scan.self_polynomial_scan.hits)}` / "
                    f"`{len(canonical_j_scan.self_polynomial_scan.moduli_checked) * len(canonical_j_scan.self_polynomial_scan.fg_degree_values) * len(canonical_j_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side self-fractional-linear uniqueness boxes: "
                    f"`{len(canonical_j_scan.self_fractional_linear_scan.hits)}` / "
                    f"`{len(canonical_j_scan.self_fractional_linear_scan.moduli_checked) * len(canonical_j_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side self-quotient finite-product boxes: "
                    f"`{len(canonical_j_self_product_hits)}` / "
                    f"`{len(canonical_j_scan.self_quotient_product_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side eta templates: `{len(canonical_j_eta_hits)}` / "
                    f"`{len(canonical_j_scan.eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side modular-unit / eta templates: "
                    f"`{len(canonical_j_modular_hits)}` / `{len(canonical_j_scan.modular_unit_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side plus-Pochhammer templates: "
                    f"`{len(canonical_j_plus_hits)}` / `{len(canonical_j_scan.self_plus_pochhammer_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber `j`-side plus-Pochhammer + eta templates: "
                    f"`{len(canonical_j_plus_eta_hits)}` / "
                    f"`{len(canonical_j_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if canonical_j_scan.named_gg_modular_equation_scan is not None:
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=f"Canonical Weber `j`-side coordinate `{canonical_j_scan.label}`",
                        gg_scan=canonical_j_scan.named_gg_modular_equation_scan,
                        series_symbol=series_symbol,
                    )
                if canonical_j_scan.normalized_followup is not None:
                    canonical_j_followup = canonical_j_scan.normalized_followup
                    canonical_j_followup_self_product_hits = [
                        scan
                        for scan in canonical_j_followup.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    canonical_j_followup_eta_hits = [
                        scan for scan in canonical_j_followup.eta_scans if scan.relation is not None
                    ]
                    canonical_j_followup_modular_hits = [
                        scan
                        for scan in canonical_j_followup.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    canonical_j_followup_plus_hits = [
                        scan
                        for scan in canonical_j_followup.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    canonical_j_followup_plus_eta_hits = [
                        scan
                        for scan in canonical_j_followup.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Canonical Weber `j`-side normalized follow-up `{canonical_j_followup.label}`: "
                        f"`{canonical_j_followup.expression}`."
                    )
                    if (
                        canonical_j_followup.first_failure_power is None
                        or canonical_j_followup.first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side normalized follow-up `{canonical_j_followup.label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber `j`-side normalized follow-up `{canonical_j_followup.label}`: "
                            f"`{canonical_j_followup.label} - 1` first fails at "
                            f"`{series_symbol}^{canonical_j_followup.first_failure_power}` with coefficient "
                            f"`{_format_expr(canonical_j_followup.first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized self-polynomial uniqueness boxes: "
                        f"`{len(canonical_j_followup.self_polynomial_scan.hits)}` / "
                        f"`{len(canonical_j_followup.self_polynomial_scan.moduli_checked) * len(canonical_j_followup.self_polynomial_scan.fg_degree_values) * len(canonical_j_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized self-fractional-linear uniqueness boxes: "
                        f"`{len(canonical_j_followup.self_fractional_linear_scan.hits)}` / "
                        f"`{len(canonical_j_followup.self_fractional_linear_scan.moduli_checked) * len(canonical_j_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized self-quotient finite-product boxes: "
                        f"`{len(canonical_j_followup_self_product_hits)}` / "
                        f"`{len(canonical_j_followup.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized eta templates: `{len(canonical_j_followup_eta_hits)}` / "
                        f"`{len(canonical_j_followup.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized modular-unit / eta templates: "
                        f"`{len(canonical_j_followup_modular_hits)}` / `{len(canonical_j_followup.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized plus-Pochhammer templates: "
                        f"`{len(canonical_j_followup_plus_hits)}` / `{len(canonical_j_followup.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side normalized plus-Pochhammer + eta templates: "
                        f"`{len(canonical_j_followup_plus_eta_hits)}` / "
                        f"`{len(canonical_j_followup.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                    if canonical_j_followup.named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber `j`-side normalized follow-up "
                                f"`{canonical_j_followup.label}`"
                            ),
                            gg_scan=canonical_j_followup.named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                anchor_canonical_j_scan = bridge_scan.anchor_canonical_j_coordinate_scan
                anchor_canonical_j_eta_hits = [
                    scan for scan in anchor_canonical_j_scan.eta_scans if scan.relation is not None
                ]
                anchor_canonical_j_modular_hits = [
                    scan
                    for scan in anchor_canonical_j_scan.modular_unit_eta_scans
                    if scan.relation is not None
                ]
                anchor_canonical_j_plus_hits = [
                    scan
                    for scan in anchor_canonical_j_scan.self_plus_pochhammer_scans
                    if scan.relation is not None
                ]
                anchor_canonical_j_plus_eta_hits = [
                    scan
                    for scan in anchor_canonical_j_scan.self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                anchor_canonical_j_self_product_hits = [
                    scan
                    for scan in anchor_canonical_j_scan.self_quotient_product_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Canonical Weber anchor `j`-side coordinate `{bridge_scan.anchor_canonical_j_coordinate_label}`: "
                    f"`{bridge_scan.anchor_canonical_j_coordinate_expression}`."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side reason: {bridge_scan.anchor_canonical_j_coordinate_reason}"
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side bridge: "
                    f"`{bridge_scan.anchor_canonical_j_coordinate_bridge_expression}`."
                )
                if (
                    anchor_canonical_j_scan.first_failure_power is None
                    or anchor_canonical_j_scan.first_failure_coeff is None
                ):
                    lines.append(
                        f"- Canonical Weber anchor `j`-side coordinate `{anchor_canonical_j_scan.label}`: "
                        "matches `1` through the checked truncation."
                    )
                else:
                    lines.append(
                        f"- Canonical Weber anchor `j`-side coordinate `{anchor_canonical_j_scan.label}`: "
                        f"`{anchor_canonical_j_scan.label} - 1` first fails at "
                        f"`{series_symbol}^{anchor_canonical_j_scan.first_failure_power}` with coefficient "
                        f"`{_format_expr(anchor_canonical_j_scan.first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Canonical Weber anchor `j`-side self-polynomial uniqueness boxes: "
                    f"`{len(anchor_canonical_j_scan.self_polynomial_scan.hits)}` / "
                    f"`{len(anchor_canonical_j_scan.self_polynomial_scan.moduli_checked) * len(anchor_canonical_j_scan.self_polynomial_scan.fg_degree_values) * len(anchor_canonical_j_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side self-fractional-linear uniqueness boxes: "
                    f"`{len(anchor_canonical_j_scan.self_fractional_linear_scan.hits)}` / "
                    f"`{len(anchor_canonical_j_scan.self_fractional_linear_scan.moduli_checked) * len(anchor_canonical_j_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side self-quotient finite-product boxes: "
                    f"`{len(anchor_canonical_j_self_product_hits)}` / "
                    f"`{len(anchor_canonical_j_scan.self_quotient_product_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side eta templates: `{len(anchor_canonical_j_eta_hits)}` / "
                    f"`{len(anchor_canonical_j_scan.eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side modular-unit / eta templates: "
                    f"`{len(anchor_canonical_j_modular_hits)}` / `{len(anchor_canonical_j_scan.modular_unit_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side plus-Pochhammer templates: "
                    f"`{len(anchor_canonical_j_plus_hits)}` / `{len(anchor_canonical_j_scan.self_plus_pochhammer_scans)}` hit boxes."
                )
                lines.append(
                    f"- Canonical Weber anchor `j`-side plus-Pochhammer + eta templates: "
                    f"`{len(anchor_canonical_j_plus_eta_hits)}` / "
                    f"`{len(anchor_canonical_j_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if anchor_canonical_j_scan.named_gg_modular_equation_scan is not None:
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=(
                            f"Canonical Weber anchor `j`-side coordinate "
                            f"`{anchor_canonical_j_scan.label}`"
                        ),
                        gg_scan=anchor_canonical_j_scan.named_gg_modular_equation_scan,
                        series_symbol=series_symbol,
                    )
                if (
                    bridge_scan.alternate_anchor_canonical_j_coordinate_label is not None
                    and bridge_scan.alternate_anchor_canonical_j_coordinate_expression is not None
                    and bridge_scan.alternate_anchor_canonical_j_coordinate_reason is not None
                    and bridge_scan.alternate_anchor_canonical_j_coordinate_bridge_expression is not None
                ):
                    lines.append(
                        f"- Canonical Weber alternate anchor `j`-side coordinate `{bridge_scan.alternate_anchor_canonical_j_coordinate_label}`: "
                        f"`{bridge_scan.alternate_anchor_canonical_j_coordinate_expression}`."
                    )
                    lines.append(
                        f"- Canonical Weber alternate anchor `j`-side reason: {bridge_scan.alternate_anchor_canonical_j_coordinate_reason}"
                    )
                    lines.append(
                        f"- Canonical Weber alternate anchor `j`-side bridge: "
                        f"`{bridge_scan.alternate_anchor_canonical_j_coordinate_bridge_expression}`."
                    )
                    if (
                        bridge_scan.alternate_anchor_canonical_j_coordinate_first_failure_power is None
                        or bridge_scan.alternate_anchor_canonical_j_coordinate_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber alternate anchor `j`-side coordinate `{bridge_scan.alternate_anchor_canonical_j_coordinate_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber alternate anchor `j`-side coordinate `{bridge_scan.alternate_anchor_canonical_j_coordinate_label}`: "
                            f"`{bridge_scan.alternate_anchor_canonical_j_coordinate_label} - 1` first fails at "
                            f"`{series_symbol}^{bridge_scan.alternate_anchor_canonical_j_coordinate_first_failure_power}` with coefficient "
                            f"`{_format_expr(bridge_scan.alternate_anchor_canonical_j_coordinate_first_failure_coeff)}`."
                        )
                if anchor_canonical_j_scan.normalized_followup is not None:
                    anchor_canonical_j_followup = anchor_canonical_j_scan.normalized_followup
                    anchor_canonical_j_followup_self_product_hits = [
                        scan
                        for scan in anchor_canonical_j_followup.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    anchor_canonical_j_followup_eta_hits = [
                        scan
                        for scan in anchor_canonical_j_followup.eta_scans
                        if scan.relation is not None
                    ]
                    anchor_canonical_j_followup_modular_hits = [
                        scan
                        for scan in anchor_canonical_j_followup.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    anchor_canonical_j_followup_plus_hits = [
                        scan
                        for scan in anchor_canonical_j_followup.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    anchor_canonical_j_followup_plus_eta_hits = [
                        scan
                        for scan in anchor_canonical_j_followup.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized follow-up `{anchor_canonical_j_followup.label}`: "
                        f"`{anchor_canonical_j_followup.expression}`."
                    )
                    if (
                        anchor_canonical_j_followup.first_failure_power is None
                        or anchor_canonical_j_followup.first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber anchor `j`-side normalized follow-up `{anchor_canonical_j_followup.label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber anchor `j`-side normalized follow-up `{anchor_canonical_j_followup.label}`: "
                            f"`{anchor_canonical_j_followup.label} - 1` first fails at "
                            f"`{series_symbol}^{anchor_canonical_j_followup.first_failure_power}` with coefficient "
                            f"`{_format_expr(anchor_canonical_j_followup.first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized self-polynomial uniqueness boxes: "
                        f"`{len(anchor_canonical_j_followup.self_polynomial_scan.hits)}` / "
                        f"`{len(anchor_canonical_j_followup.self_polynomial_scan.moduli_checked) * len(anchor_canonical_j_followup.self_polynomial_scan.fg_degree_values) * len(anchor_canonical_j_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized self-fractional-linear uniqueness boxes: "
                        f"`{len(anchor_canonical_j_followup.self_fractional_linear_scan.hits)}` / "
                        f"`{len(anchor_canonical_j_followup.self_fractional_linear_scan.moduli_checked) * len(anchor_canonical_j_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized self-quotient finite-product boxes: "
                        f"`{len(anchor_canonical_j_followup_self_product_hits)}` / "
                        f"`{len(anchor_canonical_j_followup.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized eta templates: `{len(anchor_canonical_j_followup_eta_hits)}` / "
                        f"`{len(anchor_canonical_j_followup.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized modular-unit / eta templates: "
                        f"`{len(anchor_canonical_j_followup_modular_hits)}` / `{len(anchor_canonical_j_followup.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized plus-Pochhammer templates: "
                        f"`{len(anchor_canonical_j_followup_plus_hits)}` / `{len(anchor_canonical_j_followup.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber anchor `j`-side normalized plus-Pochhammer + eta templates: "
                        f"`{len(anchor_canonical_j_followup_plus_eta_hits)}` / "
                        f"`{len(anchor_canonical_j_followup.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                    if anchor_canonical_j_followup.named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber anchor `j`-side normalized follow-up "
                                f"`{anchor_canonical_j_followup.label}`"
                            ),
                            gg_scan=anchor_canonical_j_followup.named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                canonical_j_anchor_bridge = bridge_scan.canonical_j_anchor_bridge_scan
                if canonical_j_anchor_bridge is not None:
                    canonical_j_anchor_quotient = canonical_j_anchor_bridge.quotient_scan
                    lines.append(
                        f"- Canonical Weber `j`-side anchor bridge difference `{canonical_j_anchor_bridge.difference_label}`: "
                        f"`{canonical_j_anchor_bridge.difference_expression}`."
                    )
                    if (
                        canonical_j_anchor_bridge.difference_first_failure_power is None
                        or canonical_j_anchor_bridge.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            "- Canonical Weber `j`-side anchor bridge difference verdict: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            "- Canonical Weber `j`-side anchor bridge difference verdict: "
                            f"first fails at `{series_symbol}^{canonical_j_anchor_bridge.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(canonical_j_anchor_bridge.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side anchor bridge quotient `{canonical_j_anchor_bridge.quotient_label}`: "
                        f"`{canonical_j_anchor_bridge.quotient_expression}`."
                    )
                    if (
                        canonical_j_anchor_bridge.quotient_first_failure_power is None
                        or canonical_j_anchor_bridge.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side anchor bridge quotient `{canonical_j_anchor_bridge.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber `j`-side anchor bridge quotient `{canonical_j_anchor_bridge.quotient_label}`: "
                            f"`{canonical_j_anchor_bridge.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{canonical_j_anchor_bridge.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(canonical_j_anchor_bridge.quotient_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side anchor bridge polynomial boxes: "
                        f"`{sum(1 for scan in canonical_j_anchor_bridge.polynomial_scans if scan.relation is not None)}` / "
                        f"`{len(canonical_j_anchor_bridge.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side anchor bridge fractional-linear box: "
                        f"`{1 if canonical_j_anchor_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side anchor quotient self-polynomial uniqueness boxes: "
                        f"`{len(canonical_j_anchor_quotient.self_polynomial_scan.hits)}` / "
                        f"`{len(canonical_j_anchor_quotient.self_polynomial_scan.moduli_checked) * len(canonical_j_anchor_quotient.self_polynomial_scan.fg_degree_values) * len(canonical_j_anchor_quotient.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side anchor quotient self-fractional-linear uniqueness boxes: "
                        f"`{len(canonical_j_anchor_quotient.self_fractional_linear_scan.hits)}` / "
                        f"`{len(canonical_j_anchor_quotient.self_fractional_linear_scan.moduli_checked) * len(canonical_j_anchor_quotient.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    if canonical_j_anchor_bridge.quotient_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber `j`-side anchor bridge quotient "
                                f"`{canonical_j_anchor_bridge.quotient_label}`"
                            ),
                            gg_scan=canonical_j_anchor_bridge.quotient_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if canonical_j_anchor_quotient.normalized_followup is not None:
                        canonical_j_anchor_followup = canonical_j_anchor_quotient.normalized_followup
                        lines.append(
                            f"- Canonical Weber `j`-side anchor bridge quotient normalized follow-up `{canonical_j_anchor_followup.label}`: "
                            f"`{canonical_j_anchor_followup.expression}`."
                        )
                        if (
                            canonical_j_anchor_followup.first_failure_power is None
                            or canonical_j_anchor_followup.first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Canonical Weber `j`-side anchor bridge quotient normalized follow-up `{canonical_j_anchor_followup.label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Canonical Weber `j`-side anchor bridge quotient normalized follow-up `{canonical_j_anchor_followup.label}`: "
                                f"`{canonical_j_anchor_followup.label} - 1` first fails at "
                                f"`{series_symbol}^{canonical_j_anchor_followup.first_failure_power}` with coefficient "
                                f"`{_format_expr(canonical_j_anchor_followup.first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side anchor bridge quotient normalized self-polynomial uniqueness boxes: "
                            f"`{len(canonical_j_anchor_followup.self_polynomial_scan.hits)}` / "
                            f"`{len(canonical_j_anchor_followup.self_polynomial_scan.moduli_checked) * len(canonical_j_anchor_followup.self_polynomial_scan.fg_degree_values) * len(canonical_j_anchor_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side anchor bridge quotient normalized self-fractional-linear uniqueness boxes: "
                            f"`{len(canonical_j_anchor_followup.self_fractional_linear_scan.hits)}` / "
                            f"`{len(canonical_j_anchor_followup.self_fractional_linear_scan.moduli_checked) * len(canonical_j_anchor_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                        )
                    if canonical_j_anchor_bridge.quotient_followup_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                "Canonical Weber `j`-side anchor bridge quotient normalized follow-up "
                                f"`{canonical_j_anchor_quotient.normalized_followup.label}`"
                            ),
                            gg_scan=canonical_j_anchor_bridge.quotient_followup_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if canonical_j_anchor_bridge.quotient_followup_bridge_scan is not None:
                        canonical_j_anchor_nested_bridge = canonical_j_anchor_bridge.quotient_followup_bridge_scan
                        lines.append(
                            f"- Canonical Weber `j`-side nested anchor bridge difference `{canonical_j_anchor_nested_bridge.difference_label}`: "
                            f"`{canonical_j_anchor_nested_bridge.difference_expression}`."
                        )
                        if (
                            canonical_j_anchor_nested_bridge.difference_first_failure_power is None
                            or canonical_j_anchor_nested_bridge.difference_first_failure_coeff is None
                        ):
                            lines.append(
                                "- Canonical Weber `j`-side nested anchor bridge difference verdict: "
                                "matches `0` through the checked truncation."
                            )
                        else:
                            lines.append(
                                "- Canonical Weber `j`-side nested anchor bridge difference verdict: "
                                f"first fails at `{series_symbol}^{canonical_j_anchor_nested_bridge.difference_first_failure_power}` "
                                f"with coefficient `{_format_expr(canonical_j_anchor_nested_bridge.difference_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side nested anchor bridge quotient `{canonical_j_anchor_nested_bridge.quotient_label}`: "
                            f"`{canonical_j_anchor_nested_bridge.quotient_expression}`."
                        )
                        if (
                            canonical_j_anchor_nested_bridge.quotient_first_failure_power is None
                            or canonical_j_anchor_nested_bridge.quotient_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Canonical Weber `j`-side nested anchor bridge quotient `{canonical_j_anchor_nested_bridge.quotient_label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Canonical Weber `j`-side nested anchor bridge quotient `{canonical_j_anchor_nested_bridge.quotient_label}`: "
                                f"`{canonical_j_anchor_nested_bridge.quotient_label} - 1` first fails at "
                                f"`{series_symbol}^{canonical_j_anchor_nested_bridge.quotient_first_failure_power}` with coefficient "
                                f"`{_format_expr(canonical_j_anchor_nested_bridge.quotient_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side nested anchor bridge polynomial boxes: "
                            f"`{sum(1 for scan in canonical_j_anchor_nested_bridge.polynomial_scans if scan.relation is not None)}` / "
                            f"`{len(canonical_j_anchor_nested_bridge.polynomial_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side nested anchor bridge fractional-linear box: "
                            f"`{1 if canonical_j_anchor_nested_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                        )
                        if canonical_j_anchor_nested_bridge.quotient_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    f"Canonical Weber `j`-side nested anchor bridge quotient "
                                    f"`{canonical_j_anchor_nested_bridge.quotient_label}`"
                                ),
                                gg_scan=canonical_j_anchor_nested_bridge.quotient_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                canonical_j_lift_bridge = bridge_scan.canonical_j_lift_bridge_scan
                if canonical_j_lift_bridge is not None:
                    lines.append(
                        f"- Canonical Weber `j`-side lift-bridge difference `{canonical_j_lift_bridge.difference_label}`: "
                        f"`{canonical_j_lift_bridge.difference_expression}`."
                    )
                    if (
                        canonical_j_lift_bridge.difference_first_failure_power is None
                        or canonical_j_lift_bridge.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            "- Canonical Weber `j`-side lift-bridge difference verdict: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            "- Canonical Weber `j`-side lift-bridge difference verdict: "
                            f"first fails at `{series_symbol}^{canonical_j_lift_bridge.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(canonical_j_lift_bridge.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side lift-bridge quotient `{canonical_j_lift_bridge.quotient_label}`: "
                        f"`{canonical_j_lift_bridge.quotient_expression}`."
                    )
                    if (
                        canonical_j_lift_bridge.quotient_first_failure_power is None
                        or canonical_j_lift_bridge.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side lift-bridge quotient `{canonical_j_lift_bridge.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber `j`-side lift-bridge quotient `{canonical_j_lift_bridge.quotient_label}`: "
                            f"`{canonical_j_lift_bridge.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{canonical_j_lift_bridge.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(canonical_j_lift_bridge.quotient_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side lift-bridge polynomial boxes: "
                        f"`{sum(1 for scan in canonical_j_lift_bridge.polynomial_scans if scan.relation is not None)}` / "
                        f"`{len(canonical_j_lift_bridge.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side lift-bridge fractional-linear box: "
                        f"`{1 if canonical_j_lift_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                    )
                    if canonical_j_lift_bridge.quotient_eta_scans is not None:
                        lift_bridge_quotient_eta_hit_count = sum(
                            1 for scan in canonical_j_lift_bridge.quotient_eta_scans if scan.relation is not None
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side lift-bridge quotient eta templates: "
                            f"`{lift_bridge_quotient_eta_hit_count}` / `{len(canonical_j_lift_bridge.quotient_eta_scans)}` hit boxes."
                        )
                    if canonical_j_lift_bridge.quotient_modular_unit_eta_scans is not None:
                        lift_bridge_quotient_mu_eta_hit_count = sum(
                            1
                            for scan in canonical_j_lift_bridge.quotient_modular_unit_eta_scans
                            if scan.relation is not None
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side lift-bridge quotient modular-unit / eta templates: "
                            f"`{lift_bridge_quotient_mu_eta_hit_count}` / `{len(canonical_j_lift_bridge.quotient_modular_unit_eta_scans)}` hit boxes."
                        )
                    if canonical_j_lift_bridge.quotient_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber `j`-side lift-bridge quotient "
                                f"`{canonical_j_lift_bridge.quotient_label}`"
                            ),
                            gg_scan=canonical_j_lift_bridge.quotient_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if (
                        canonical_j_lift_bridge.quotient_followup_named_gg_modular_equation_scan is not None
                        and canonical_j_lift_bridge.quotient_followup_label is not None
                        and canonical_j_lift_bridge.quotient_followup_expression is not None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side lift-bridge quotient normalized follow-up `{canonical_j_lift_bridge.quotient_followup_label}`: "
                            f"`{canonical_j_lift_bridge.quotient_followup_expression}`."
                        )
                        if canonical_j_lift_bridge.quotient_followup_eta_scans is not None:
                            lift_bridge_followup_eta_hit_count = sum(
                                1
                                for scan in canonical_j_lift_bridge.quotient_followup_eta_scans
                                if scan.relation is not None
                            )
                            lines.append(
                                f"- Canonical Weber `j`-side lift-bridge quotient normalized follow-up `{canonical_j_lift_bridge.quotient_followup_label}` eta templates: "
                                f"`{lift_bridge_followup_eta_hit_count}` / `{len(canonical_j_lift_bridge.quotient_followup_eta_scans)}` hit boxes."
                            )
                        if canonical_j_lift_bridge.quotient_followup_modular_unit_eta_scans is not None:
                            lift_bridge_followup_mu_eta_hit_count = sum(
                                1
                                for scan in canonical_j_lift_bridge.quotient_followup_modular_unit_eta_scans
                                if scan.relation is not None
                            )
                            lines.append(
                                f"- Canonical Weber `j`-side lift-bridge quotient normalized follow-up `{canonical_j_lift_bridge.quotient_followup_label}` modular-unit / eta templates: "
                                f"`{lift_bridge_followup_mu_eta_hit_count}` / `{len(canonical_j_lift_bridge.quotient_followup_modular_unit_eta_scans)}` hit boxes."
                            )
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                "Canonical Weber `j`-side lift-bridge quotient normalized follow-up "
                                f"`{canonical_j_lift_bridge.quotient_followup_label}`"
                            ),
                            gg_scan=canonical_j_lift_bridge.quotient_followup_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                canonical_j_alt_lift_bridge = bridge_scan.canonical_j_alt_lift_bridge_scan
                if canonical_j_alt_lift_bridge is not None:
                    lines.append(
                        f"- Canonical Weber `j`-side alternate lift-bridge difference `{canonical_j_alt_lift_bridge.difference_label}`: "
                        f"`{canonical_j_alt_lift_bridge.difference_expression}`."
                    )
                    if (
                        canonical_j_alt_lift_bridge.difference_first_failure_power is None
                        or canonical_j_alt_lift_bridge.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            "- Canonical Weber `j`-side alternate lift-bridge difference verdict: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            "- Canonical Weber `j`-side alternate lift-bridge difference verdict: "
                            f"first fails at `{series_symbol}^{canonical_j_alt_lift_bridge.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(canonical_j_alt_lift_bridge.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side alternate lift-bridge quotient `{canonical_j_alt_lift_bridge.quotient_label}`: "
                        f"`{canonical_j_alt_lift_bridge.quotient_expression}`."
                    )
                    if (
                        canonical_j_alt_lift_bridge.quotient_first_failure_power is None
                        or canonical_j_alt_lift_bridge.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side alternate lift-bridge quotient `{canonical_j_alt_lift_bridge.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber `j`-side alternate lift-bridge quotient `{canonical_j_alt_lift_bridge.quotient_label}`: "
                            f"`{canonical_j_alt_lift_bridge.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{canonical_j_alt_lift_bridge.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(canonical_j_alt_lift_bridge.quotient_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side alternate lift-bridge polynomial boxes: "
                        f"`{sum(1 for scan in canonical_j_alt_lift_bridge.polynomial_scans if scan.relation is not None)}` / "
                        f"`{len(canonical_j_alt_lift_bridge.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side alternate lift-bridge fractional-linear box: "
                        f"`{1 if canonical_j_alt_lift_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                    )
                    if canonical_j_alt_lift_bridge.quotient_eta_scans is not None:
                        alt_lift_bridge_quotient_eta_hit_count = sum(
                            1 for scan in canonical_j_alt_lift_bridge.quotient_eta_scans if scan.relation is not None
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side alternate lift-bridge quotient eta templates: "
                            f"`{alt_lift_bridge_quotient_eta_hit_count}` / `{len(canonical_j_alt_lift_bridge.quotient_eta_scans)}` hit boxes."
                        )
                    if canonical_j_alt_lift_bridge.quotient_modular_unit_eta_scans is not None:
                        alt_lift_bridge_quotient_mu_eta_hit_count = sum(
                            1
                            for scan in canonical_j_alt_lift_bridge.quotient_modular_unit_eta_scans
                            if scan.relation is not None
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side alternate lift-bridge quotient modular-unit / eta templates: "
                            f"`{alt_lift_bridge_quotient_mu_eta_hit_count}` / `{len(canonical_j_alt_lift_bridge.quotient_modular_unit_eta_scans)}` hit boxes."
                        )
                    if canonical_j_alt_lift_bridge.quotient_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber `j`-side alternate lift-bridge quotient "
                                f"`{canonical_j_alt_lift_bridge.quotient_label}`"
                            ),
                            gg_scan=canonical_j_alt_lift_bridge.quotient_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if (
                        canonical_j_alt_lift_bridge.quotient_followup_named_gg_modular_equation_scan is not None
                        and canonical_j_alt_lift_bridge.quotient_followup_label is not None
                        and canonical_j_alt_lift_bridge.quotient_followup_expression is not None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side alternate lift-bridge quotient normalized follow-up `{canonical_j_alt_lift_bridge.quotient_followup_label}`: "
                            f"`{canonical_j_alt_lift_bridge.quotient_followup_expression}`."
                        )
                        if canonical_j_alt_lift_bridge.quotient_followup_eta_scans is not None:
                            alt_lift_bridge_followup_eta_hit_count = sum(
                                1
                                for scan in canonical_j_alt_lift_bridge.quotient_followup_eta_scans
                                if scan.relation is not None
                            )
                            lines.append(
                                f"- Canonical Weber `j`-side alternate lift-bridge quotient normalized follow-up `{canonical_j_alt_lift_bridge.quotient_followup_label}` eta templates: "
                                f"`{alt_lift_bridge_followup_eta_hit_count}` / `{len(canonical_j_alt_lift_bridge.quotient_followup_eta_scans)}` hit boxes."
                            )
                        if canonical_j_alt_lift_bridge.quotient_followup_modular_unit_eta_scans is not None:
                            alt_lift_bridge_followup_mu_eta_hit_count = sum(
                                1
                                for scan in canonical_j_alt_lift_bridge.quotient_followup_modular_unit_eta_scans
                                if scan.relation is not None
                            )
                            lines.append(
                                f"- Canonical Weber `j`-side alternate lift-bridge quotient normalized follow-up `{canonical_j_alt_lift_bridge.quotient_followup_label}` modular-unit / eta templates: "
                                f"`{alt_lift_bridge_followup_mu_eta_hit_count}` / `{len(canonical_j_alt_lift_bridge.quotient_followup_modular_unit_eta_scans)}` hit boxes."
                            )
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                "Canonical Weber `j`-side alternate lift-bridge quotient normalized follow-up "
                                f"`{canonical_j_alt_lift_bridge.quotient_followup_label}`"
                            ),
                            gg_scan=canonical_j_alt_lift_bridge.quotient_followup_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                if sample.weber_j_pb_bridge_scan is not None:
                    j_pb_bridge = sample.weber_j_pb_bridge_scan
                    j_pb_quotient = j_pb_bridge.quotient_scan
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge difference `{j_pb_bridge.difference_label}`: "
                        f"`{j_pb_bridge.difference_expression}`."
                    )
                    if (
                        j_pb_bridge.difference_first_failure_power is None
                        or j_pb_bridge.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            "- Canonical Weber `j`-side vs `P/B` bridge difference verdict: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            "- Canonical Weber `j`-side vs `P/B` bridge difference verdict: "
                            f"first fails at `{series_symbol}^{j_pb_bridge.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(j_pb_bridge.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge quotient `{j_pb_bridge.quotient_label}`: "
                        f"`{j_pb_bridge.quotient_expression}`."
                    )
                    if (
                        j_pb_bridge.quotient_first_failure_power is None
                        or j_pb_bridge.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` bridge quotient `{j_pb_bridge.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` bridge quotient `{j_pb_bridge.quotient_label}`: "
                            f"`{j_pb_bridge.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{j_pb_bridge.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(j_pb_bridge.quotient_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge polynomial boxes: "
                        f"`{sum(1 for scan in j_pb_bridge.polynomial_scans if scan.relation is not None)}` / "
                        f"`{len(j_pb_bridge.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge fractional-linear box: "
                        f"`{1 if j_pb_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge quotient self-polynomial uniqueness boxes: "
                        f"`{len(j_pb_quotient.self_polynomial_scan.hits)}` / "
                        f"`{len(j_pb_quotient.self_polynomial_scan.moduli_checked) * len(j_pb_quotient.self_polynomial_scan.fg_degree_values) * len(j_pb_quotient.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Canonical Weber `j`-side vs `P/B` bridge quotient self-fractional-linear uniqueness boxes: "
                        f"`{len(j_pb_quotient.self_fractional_linear_scan.hits)}` / "
                        f"`{len(j_pb_quotient.self_fractional_linear_scan.moduli_checked) * len(j_pb_quotient.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    if j_pb_bridge.quotient_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Canonical Weber `j`-side vs `P/B` bridge quotient "
                                f"`{j_pb_bridge.quotient_label}`"
                            ),
                            gg_scan=j_pb_bridge.quotient_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if j_pb_quotient.normalized_followup is not None:
                        j_pb_followup = j_pb_quotient.normalized_followup
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` bridge quotient normalized follow-up `{j_pb_followup.label}`: "
                            f"`{j_pb_followup.expression}`."
                        )
                        if (
                            j_pb_followup.first_failure_power is None
                            or j_pb_followup.first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Canonical Weber `j`-side vs `P/B` bridge quotient normalized follow-up `{j_pb_followup.label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Canonical Weber `j`-side vs `P/B` bridge quotient normalized follow-up `{j_pb_followup.label}`: "
                                f"`{j_pb_followup.label} - 1` first fails at "
                                f"`{series_symbol}^{j_pb_followup.first_failure_power}` with coefficient "
                                f"`{_format_expr(j_pb_followup.first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` bridge quotient normalized self-polynomial uniqueness boxes: "
                            f"`{len(j_pb_followup.self_polynomial_scan.hits)}` / "
                            f"`{len(j_pb_followup.self_polynomial_scan.moduli_checked) * len(j_pb_followup.self_polynomial_scan.fg_degree_values) * len(j_pb_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` bridge quotient normalized self-fractional-linear uniqueness boxes: "
                            f"`{len(j_pb_followup.self_fractional_linear_scan.hits)}` / "
                            f"`{len(j_pb_followup.self_fractional_linear_scan.moduli_checked) * len(j_pb_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                        )
                    if j_pb_bridge.quotient_followup_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                "Canonical Weber `j`-side vs `P/B` bridge quotient normalized follow-up "
                                f"`{j_pb_quotient.normalized_followup.label}`"
                            ),
                            gg_scan=j_pb_bridge.quotient_followup_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if j_pb_bridge.quotient_followup_bridge_scan is not None:
                        j_pb_nested_bridge = j_pb_bridge.quotient_followup_bridge_scan
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` nested bridge difference `{j_pb_nested_bridge.difference_label}`: "
                            f"`{j_pb_nested_bridge.difference_expression}`."
                        )
                        if (
                            j_pb_nested_bridge.difference_first_failure_power is None
                            or j_pb_nested_bridge.difference_first_failure_coeff is None
                        ):
                            lines.append(
                                "- Canonical Weber `j`-side vs `P/B` nested bridge difference verdict: "
                                "matches `0` through the checked truncation."
                            )
                        else:
                            lines.append(
                                "- Canonical Weber `j`-side vs `P/B` nested bridge difference verdict: "
                                f"first fails at `{series_symbol}^{j_pb_nested_bridge.difference_first_failure_power}` "
                                f"with coefficient `{_format_expr(j_pb_nested_bridge.difference_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` nested bridge quotient `{j_pb_nested_bridge.quotient_label}`: "
                            f"`{j_pb_nested_bridge.quotient_expression}`."
                        )
                        if (
                            j_pb_nested_bridge.quotient_first_failure_power is None
                            or j_pb_nested_bridge.quotient_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Canonical Weber `j`-side vs `P/B` nested bridge quotient `{j_pb_nested_bridge.quotient_label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Canonical Weber `j`-side vs `P/B` nested bridge quotient `{j_pb_nested_bridge.quotient_label}`: "
                                f"`{j_pb_nested_bridge.quotient_label} - 1` first fails at "
                                f"`{series_symbol}^{j_pb_nested_bridge.quotient_first_failure_power}` with coefficient "
                                f"`{_format_expr(j_pb_nested_bridge.quotient_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` nested bridge polynomial boxes: "
                            f"`{sum(1 for scan in j_pb_nested_bridge.polynomial_scans if scan.relation is not None)}` / "
                            f"`{len(j_pb_nested_bridge.polynomial_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Canonical Weber `j`-side vs `P/B` nested bridge fractional-linear box: "
                            f"`{1 if j_pb_nested_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                        )
                        if j_pb_nested_bridge.quotient_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    f"Canonical Weber `j`-side vs `P/B` nested bridge quotient "
                                    f"`{j_pb_nested_bridge.quotient_label}`"
                                ),
                                gg_scan=j_pb_nested_bridge.quotient_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                if sample.weber_j_lift_pivot_bridge_scans:
                    pivot_prefix = "Canonical Weber `j`-side pivot bridge"
                    for pivot_bridge in sample.weber_j_lift_pivot_bridge_scans:
                        lines.append(
                            f"- {pivot_prefix} difference `{pivot_bridge.difference_label}`: "
                            f"`{pivot_bridge.difference_expression}`."
                        )
                        if (
                            pivot_bridge.difference_first_failure_power is None
                            or pivot_bridge.difference_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- {pivot_prefix} difference verdict: matches `0` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- {pivot_prefix} difference verdict: first fails at "
                                f"`{series_symbol}^{pivot_bridge.difference_first_failure_power}` with coefficient "
                                f"`{_format_expr(pivot_bridge.difference_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- {pivot_prefix} quotient `{pivot_bridge.quotient_label}`: "
                            f"`{pivot_bridge.quotient_expression}`."
                        )
                        if (
                            pivot_bridge.quotient_first_failure_power is None
                            or pivot_bridge.quotient_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- {pivot_prefix} quotient `{pivot_bridge.quotient_label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- {pivot_prefix} quotient `{pivot_bridge.quotient_label}`: "
                                f"`{pivot_bridge.quotient_label} - 1` first fails at "
                                f"`{series_symbol}^{pivot_bridge.quotient_first_failure_power}` with coefficient "
                                f"`{_format_expr(pivot_bridge.quotient_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- {pivot_prefix} polynomial boxes: "
                            f"`{sum(1 for scan in pivot_bridge.polynomial_scans if scan.relation is not None)}` / "
                            f"`{len(pivot_bridge.polynomial_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- {pivot_prefix} fractional-linear box: "
                            f"`{1 if pivot_bridge.fractional_linear_relation is not None else 0}` / `1` hit boxes."
                        )
                        if pivot_bridge.quotient_eta_scans is not None:
                            pivot_quotient_eta_hit_count = sum(
                                1 for scan in pivot_bridge.quotient_eta_scans if scan.relation is not None
                            )
                            lines.append(
                                f"- {pivot_prefix} quotient eta templates: "
                                f"`{pivot_quotient_eta_hit_count}` / `{len(pivot_bridge.quotient_eta_scans)}` hit boxes."
                            )
                        if pivot_bridge.quotient_modular_unit_eta_scans is not None:
                            pivot_quotient_mu_eta_hit_count = sum(
                                1
                                for scan in pivot_bridge.quotient_modular_unit_eta_scans
                                if scan.relation is not None
                            )
                            lines.append(
                                f"- {pivot_prefix} quotient modular-unit / eta templates: "
                                f"`{pivot_quotient_mu_eta_hit_count}` / `{len(pivot_bridge.quotient_modular_unit_eta_scans)}` hit boxes."
                            )
                        if pivot_bridge.quotient_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    f"{pivot_prefix} quotient "
                                    f"`{pivot_bridge.quotient_label}`"
                                ),
                                gg_scan=pivot_bridge.quotient_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                        if (
                            pivot_bridge.quotient_followup_named_gg_modular_equation_scan is not None
                            and pivot_bridge.quotient_followup_label is not None
                            and pivot_bridge.quotient_followup_expression is not None
                        ):
                            lines.append(
                                f"- {pivot_prefix} quotient normalized follow-up `{pivot_bridge.quotient_followup_label}`: "
                                f"`{pivot_bridge.quotient_followup_expression}`."
                            )
                            if pivot_bridge.quotient_followup_eta_scans is not None:
                                pivot_followup_eta_hit_count = sum(
                                    1
                                    for scan in pivot_bridge.quotient_followup_eta_scans
                                    if scan.relation is not None
                                )
                                lines.append(
                                    f"- {pivot_prefix} quotient normalized follow-up `{pivot_bridge.quotient_followup_label}` eta templates: "
                                    f"`{pivot_followup_eta_hit_count}` / `{len(pivot_bridge.quotient_followup_eta_scans)}` hit boxes."
                                )
                            if pivot_bridge.quotient_followup_modular_unit_eta_scans is not None:
                                pivot_followup_mu_eta_hit_count = sum(
                                    1
                                    for scan in pivot_bridge.quotient_followup_modular_unit_eta_scans
                                    if scan.relation is not None
                                )
                                lines.append(
                                    f"- {pivot_prefix} quotient normalized follow-up `{pivot_bridge.quotient_followup_label}` modular-unit / eta templates: "
                                    f"`{pivot_followup_mu_eta_hit_count}` / `{len(pivot_bridge.quotient_followup_modular_unit_eta_scans)}` hit boxes."
                                )
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    f"{pivot_prefix} quotient normalized follow-up "
                                    f"`{pivot_bridge.quotient_followup_label}`"
                                ),
                                gg_scan=pivot_bridge.quotient_followup_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                lines.append(
                    f"- Weber residual quotient-coordinate `{bridge_scan.quotient_coordinate_label}`: "
                    f"`{bridge_scan.quotient_coordinate_expression}`."
                )
                lines.append(
                    f"- Weber residual exact quotient-coordinate bridge: "
                    f"`{bridge_scan.quotient_coordinate_bridge_expression}`."
                )
                if bridge_scan.quotient_coordinate_bridge_holds:
                    lines.append(
                        "- Weber residual exact quotient-coordinate bridge verdict: matches through the checked truncation."
                    )
                else:
                    lines.append(
                        "- Weber residual exact quotient-coordinate bridge verdict: "
                        f"first fails at `{series_symbol}^{bridge_scan.quotient_coordinate_bridge_first_failure_power}` "
                        f"with coefficient `{_format_expr(bridge_scan.quotient_coordinate_bridge_first_failure_coeff)}`."
                    )
                coordinate_scan = bridge_scan.quotient_coordinate_template_scan
                coordinate_eta_hits = [
                    scan for scan in coordinate_scan.eta_scans if scan.relation is not None
                ]
                coordinate_modular_hits = [
                    scan
                    for scan in coordinate_scan.modular_unit_eta_scans
                    if scan.relation is not None
                ]
                coordinate_plus_hits = [
                    scan
                    for scan in coordinate_scan.self_plus_pochhammer_scans
                    if scan.relation is not None
                ]
                coordinate_plus_eta_hits = [
                    scan
                    for scan in coordinate_scan.self_plus_pochhammer_eta_scans
                    if scan.relation is not None
                ]
                coordinate_self_product_hits = [
                    scan
                    for scan in coordinate_scan.self_quotient_product_scans
                    if scan.relation is not None
                ]
                lines.append(
                    f"- Weber quotient-coordinate template normalization `{coordinate_scan.label}`: "
                    f"`{coordinate_scan.expression}`."
                )
                lines.append(
                    f"- Weber quotient-coordinate template bridge: "
                    f"`{bridge_scan.quotient_coordinate_template_bridge_expression}`."
                )
                lines.append(
                    "- Weber quotient-coordinate route reading: with that classical "
                    "Weber source interpretation in place, the current "
                    "`Q_gp_ws` / `X_g_ws` / `G_X_ws` branch is the "
                    "source-faithful classical Weber quotient/template lane, so it "
                    "still outranks wider anonymous scans and the already-negative "
                    "product follow-up branch."
                )
                if (
                    coordinate_scan.first_failure_power is None
                    or coordinate_scan.first_failure_coeff is None
                ):
                    lines.append(
                        f"- Weber quotient-coordinate template normalization `{coordinate_scan.label}`: "
                        "matches `1` through the checked truncation."
                    )
                else:
                    lines.append(
                        f"- Weber quotient-coordinate template normalization `{coordinate_scan.label}`: "
                        f"`{coordinate_scan.label} - 1` first fails at "
                        f"`{series_symbol}^{coordinate_scan.first_failure_power}` with coefficient "
                        f"`{_format_expr(coordinate_scan.first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Weber quotient-coordinate template self-polynomial uniqueness boxes: "
                    f"`{len(coordinate_scan.self_polynomial_scan.hits)}` / "
                    f"`{len(coordinate_scan.self_polynomial_scan.moduli_checked) * len(coordinate_scan.self_polynomial_scan.fg_degree_values) * len(coordinate_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template self-fractional-linear uniqueness boxes: "
                    f"`{len(coordinate_scan.self_fractional_linear_scan.hits)}` / "
                    f"`{len(coordinate_scan.self_fractional_linear_scan.moduli_checked) * len(coordinate_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template self-quotient finite-product boxes: "
                    f"`{len(coordinate_self_product_hits)}` / "
                    f"`{len(coordinate_scan.self_quotient_product_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template eta templates: `{len(coordinate_eta_hits)}` / "
                    f"`{len(coordinate_scan.eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template modular-unit / eta templates: "
                    f"`{len(coordinate_modular_hits)}` / `{len(coordinate_scan.modular_unit_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template plus-Pochhammer templates: "
                    f"`{len(coordinate_plus_hits)}` / `{len(coordinate_scan.self_plus_pochhammer_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber quotient-coordinate template plus-Pochhammer + eta templates: "
                    f"`{len(coordinate_plus_eta_hits)}` / "
                    f"`{len(coordinate_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if coordinate_scan.named_gg_modular_equation_scan is not None:
                    lines.append(
                        "- Focused source-faithful quotient/template pass: the foreground `G_X_ws` lane is "
                        "now also checked directly against the named `GG` modular-equation coordinates, so "
                        "we test that template-normalized branch before falling back to the later residual quotient."
                    )
                    _append_named_gg_bridge_lines(
                        lines,
                        prefix=(
                            f"Weber quotient-coordinate template normalization "
                            f"`{coordinate_scan.label}`"
                        ),
                        gg_scan=coordinate_scan.named_gg_modular_equation_scan,
                        series_symbol=series_symbol,
                    )
                if coordinate_scan.normalized_followup is not None:
                    coordinate_followup_scan = coordinate_scan.normalized_followup
                    coordinate_followup_self_product_hits = [
                        scan
                        for scan in coordinate_followup_scan.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    coordinate_followup_eta_hits = [
                        scan
                        for scan in coordinate_followup_scan.eta_scans
                        if scan.relation is not None
                    ]
                    coordinate_followup_modular_hits = [
                        scan
                        for scan in coordinate_followup_scan.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    coordinate_followup_plus_hits = [
                        scan
                        for scan in coordinate_followup_scan.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    coordinate_followup_plus_eta_hits = [
                        scan
                        for scan in coordinate_followup_scan.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Weber quotient-coordinate normalized follow-up `{coordinate_followup_scan.label}`: "
                        f"`{coordinate_followup_scan.expression}`."
                    )
                    if (
                        coordinate_followup_scan.first_failure_power is None
                        or coordinate_followup_scan.first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Weber quotient-coordinate normalized follow-up `{coordinate_followup_scan.label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Weber quotient-coordinate normalized follow-up `{coordinate_followup_scan.label}`: "
                            f"`{coordinate_followup_scan.label} - 1` first fails at "
                            f"`{series_symbol}^{coordinate_followup_scan.first_failure_power}` with coefficient "
                            f"`{_format_expr(coordinate_followup_scan.first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Weber quotient-coordinate normalized self-polynomial uniqueness boxes: "
                        f"`{len(coordinate_followup_scan.self_polynomial_scan.hits)}` / "
                        f"`{len(coordinate_followup_scan.self_polynomial_scan.moduli_checked) * len(coordinate_followup_scan.self_polynomial_scan.fg_degree_values) * len(coordinate_followup_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized self-fractional-linear uniqueness boxes: "
                        f"`{len(coordinate_followup_scan.self_fractional_linear_scan.hits)}` / "
                        f"`{len(coordinate_followup_scan.self_fractional_linear_scan.moduli_checked) * len(coordinate_followup_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized self-quotient finite-product boxes: "
                        f"`{len(coordinate_followup_self_product_hits)}` / "
                        f"`{len(coordinate_followup_scan.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized eta templates: `{len(coordinate_followup_eta_hits)}` / "
                        f"`{len(coordinate_followup_scan.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized modular-unit / eta templates: "
                        f"`{len(coordinate_followup_modular_hits)}` / `{len(coordinate_followup_scan.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized plus-Pochhammer templates: "
                        f"`{len(coordinate_followup_plus_hits)}` / `{len(coordinate_followup_scan.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber quotient-coordinate normalized plus-Pochhammer + eta templates: "
                        f"`{len(coordinate_followup_plus_eta_hits)}` / "
                        f"`{len(coordinate_followup_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                    if coordinate_followup_scan.named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Weber quotient-coordinate normalized follow-up "
                                f"`{coordinate_followup_scan.label}`"
                            ),
                            gg_scan=coordinate_followup_scan.named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if coordinate_followup_scan.named_gg_descendant_preview is not None:
                        preview = coordinate_followup_scan.named_gg_descendant_preview
                        lines.append(
                            f"- Weber quotient-coordinate normalized follow-up `{coordinate_followup_scan.label}` "
                            f"odd-prime descendant preview: direct ladder `{', '.join(preview.direct_labels)}`; "
                            f"quotient ladder `{', '.join(preview.quotient_labels)}`."
                        )
                    if coordinate_followup_scan.named_gg_descendant_focused_scan is not None:
                        _append_named_gg_descendant_focus_lines(
                            lines,
                            prefix=(
                                f"Weber quotient-coordinate normalized follow-up "
                                f"`{coordinate_followup_scan.label}`"
                            ),
                            descendant_scan=coordinate_followup_scan.named_gg_descendant_focused_scan,
                        )
                if bridge_scan.followup_bridge_scan is not None:
                    followup_bridge = bridge_scan.followup_bridge_scan
                    polynomial_bridge_hits = [
                        scan for scan in followup_bridge.polynomial_scans if scan.relation is not None
                    ]
                    lines.append(
                        f"- Weber normalized follow-up bridge difference `{followup_bridge.difference_label}`: "
                        f"`{followup_bridge.difference_expression}`."
                    )
                    if (
                        followup_bridge.difference_first_failure_power is None
                        or followup_bridge.difference_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Weber normalized follow-up bridge difference `{followup_bridge.difference_label}`: "
                            "matches `0` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Weber normalized follow-up bridge difference `{followup_bridge.difference_label}`: "
                            f"first fails at `{series_symbol}^{followup_bridge.difference_first_failure_power}` "
                            f"with coefficient `{_format_expr(followup_bridge.difference_first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient `{followup_bridge.quotient_label}`: "
                        f"`{followup_bridge.quotient_expression}`."
                    )
                    if (
                        followup_bridge.quotient_first_failure_power is None
                        or followup_bridge.quotient_first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient `{followup_bridge.quotient_label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient `{followup_bridge.quotient_label}`: "
                            f"`{followup_bridge.quotient_label} - 1` first fails at "
                            f"`{series_symbol}^{followup_bridge.quotient_first_failure_power}` with coefficient "
                            f"`{_format_expr(followup_bridge.quotient_first_failure_coeff)}`."
                        )
                    quotient_bridge_scan = followup_bridge.quotient_scan
                    quotient_bridge_self_product_hits = [
                        scan
                        for scan in quotient_bridge_scan.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    quotient_bridge_eta_hits = [
                        scan for scan in quotient_bridge_scan.eta_scans if scan.relation is not None
                    ]
                    quotient_bridge_modular_hits = [
                        scan
                        for scan in quotient_bridge_scan.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    quotient_bridge_plus_hits = [
                        scan
                        for scan in quotient_bridge_scan.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    quotient_bridge_plus_eta_hits = [
                        scan
                        for scan in quotient_bridge_scan.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Weber normalized follow-up bridge polynomial boxes: "
                        f"`{len(polynomial_bridge_hits)}` / `{len(followup_bridge.polynomial_scans)}` hit boxes."
                    )
                    lines.append(
                        "- Weber normalized follow-up bridge fractional-linear box: "
                        + ("`1` / `1` hit boxes." if followup_bridge.fractional_linear_relation is not None else "`0` / `1` hit boxes.")
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient self-polynomial uniqueness boxes: "
                        f"`{len(quotient_bridge_scan.self_polynomial_scan.hits)}` / "
                        f"`{len(quotient_bridge_scan.self_polynomial_scan.moduli_checked) * len(quotient_bridge_scan.self_polynomial_scan.fg_degree_values) * len(quotient_bridge_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient self-fractional-linear uniqueness boxes: "
                        f"`{len(quotient_bridge_scan.self_fractional_linear_scan.hits)}` / "
                        f"`{len(quotient_bridge_scan.self_fractional_linear_scan.moduli_checked) * len(quotient_bridge_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient self-quotient finite-product boxes: "
                        f"`{len(quotient_bridge_self_product_hits)}` / "
                        f"`{len(quotient_bridge_scan.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient eta templates: `{len(quotient_bridge_eta_hits)}` / "
                        f"`{len(quotient_bridge_scan.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient modular-unit / eta templates: "
                        f"`{len(quotient_bridge_modular_hits)}` / `{len(quotient_bridge_scan.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient plus-Pochhammer templates: "
                        f"`{len(quotient_bridge_plus_hits)}` / `{len(quotient_bridge_scan.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber normalized follow-up bridge quotient plus-Pochhammer + eta templates: "
                        f"`{len(quotient_bridge_plus_eta_hits)}` / "
                        f"`{len(quotient_bridge_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                    if followup_bridge.quotient_named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=(
                                f"Weber normalized follow-up bridge quotient "
                                f"`{followup_bridge.quotient_label}`"
                            ),
                            gg_scan=followup_bridge.quotient_named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if quotient_bridge_scan.normalized_followup is not None:
                        quotient_bridge_followup = quotient_bridge_scan.normalized_followup
                        quotient_bridge_followup_self_product_hits = [
                            scan
                            for scan in quotient_bridge_followup.self_quotient_product_scans
                            if scan.relation is not None
                        ]
                        quotient_bridge_followup_eta_hits = [
                            scan for scan in quotient_bridge_followup.eta_scans if scan.relation is not None
                        ]
                        quotient_bridge_followup_modular_hits = [
                            scan
                            for scan in quotient_bridge_followup.modular_unit_eta_scans
                            if scan.relation is not None
                        ]
                        quotient_bridge_followup_plus_hits = [
                            scan
                            for scan in quotient_bridge_followup.self_plus_pochhammer_scans
                            if scan.relation is not None
                        ]
                        quotient_bridge_followup_plus_eta_hits = [
                            scan
                            for scan in quotient_bridge_followup.self_plus_pochhammer_eta_scans
                            if scan.relation is not None
                        ]
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized follow-up `{quotient_bridge_followup.label}`: "
                            f"`{quotient_bridge_followup.expression}`."
                        )
                        if (
                            quotient_bridge_followup.first_failure_power is None
                            or quotient_bridge_followup.first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Weber normalized follow-up bridge quotient normalized follow-up `{quotient_bridge_followup.label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Weber normalized follow-up bridge quotient normalized follow-up `{quotient_bridge_followup.label}`: "
                                f"`{quotient_bridge_followup.label} - 1` first fails at "
                                f"`{series_symbol}^{quotient_bridge_followup.first_failure_power}` with coefficient "
                                f"`{_format_expr(quotient_bridge_followup.first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized self-polynomial uniqueness boxes: "
                            f"`{len(quotient_bridge_followup.self_polynomial_scan.hits)}` / "
                            f"`{len(quotient_bridge_followup.self_polynomial_scan.moduli_checked) * len(quotient_bridge_followup.self_polynomial_scan.fg_degree_values) * len(quotient_bridge_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: "
                            f"`{len(quotient_bridge_followup.self_fractional_linear_scan.hits)}` / "
                            f"`{len(quotient_bridge_followup.self_fractional_linear_scan.moduli_checked) * len(quotient_bridge_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized self-quotient finite-product boxes: "
                            f"`{len(quotient_bridge_followup_self_product_hits)}` / "
                            f"`{len(quotient_bridge_followup.self_quotient_product_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized eta templates: `{len(quotient_bridge_followup_eta_hits)}` / "
                            f"`{len(quotient_bridge_followup.eta_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized modular-unit / eta templates: "
                            f"`{len(quotient_bridge_followup_modular_hits)}` / `{len(quotient_bridge_followup.modular_unit_eta_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized plus-Pochhammer templates: "
                            f"`{len(quotient_bridge_followup_plus_hits)}` / `{len(quotient_bridge_followup.self_plus_pochhammer_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber normalized follow-up bridge quotient normalized plus-Pochhammer + eta templates: "
                            f"`{len(quotient_bridge_followup_plus_eta_hits)}` / "
                            f"`{len(quotient_bridge_followup.self_plus_pochhammer_eta_scans)}` hit boxes."
                        )
                        if followup_bridge.quotient_followup_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    "Weber normalized follow-up bridge quotient normalized follow-up "
                                    f"`{quotient_bridge_followup.label}`"
                                ),
                                gg_scan=followup_bridge.quotient_followup_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                        if quotient_bridge_followup.named_gg_descendant_preview is not None:
                            preview = quotient_bridge_followup.named_gg_descendant_preview
                            lines.append(
                                f"- Weber normalized follow-up bridge quotient normalized follow-up `{quotient_bridge_followup.label}` "
                                f"odd-prime descendant preview: direct ladder `{', '.join(preview.direct_labels)}`; "
                                f"quotient ladder `{', '.join(preview.quotient_labels)}`."
                            )
                        if quotient_bridge_followup.named_gg_descendant_focused_scan is not None:
                            _append_named_gg_descendant_focus_lines(
                                lines,
                                prefix=(
                                    f"Weber normalized follow-up bridge quotient normalized follow-up "
                                    f"`{quotient_bridge_followup.label}`"
                                ),
                                descendant_scan=quotient_bridge_followup.named_gg_descendant_focused_scan,
                            )
                    if followup_bridge.quotient_followup_bridge_scan is not None:
                        quotient_followup_bridge = followup_bridge.quotient_followup_bridge_scan
                        quotient_followup_polynomial_hits = [
                            scan
                            for scan in quotient_followup_bridge.polynomial_scans
                            if scan.relation is not None
                        ]
                        lines.append(
                            f"- Weber quotient-follow-up bridge difference `{quotient_followup_bridge.difference_label}`: "
                            f"`{quotient_followup_bridge.difference_expression}`."
                        )
                        if (
                            quotient_followup_bridge.difference_first_failure_power is None
                            or quotient_followup_bridge.difference_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Weber quotient-follow-up bridge difference `{quotient_followup_bridge.difference_label}`: "
                                "matches `0` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Weber quotient-follow-up bridge difference `{quotient_followup_bridge.difference_label}`: "
                                f"first fails at `{series_symbol}^{quotient_followup_bridge.difference_first_failure_power}` "
                                f"with coefficient `{_format_expr(quotient_followup_bridge.difference_first_failure_coeff)}`."
                            )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient `{quotient_followup_bridge.quotient_label}`: "
                            f"`{quotient_followup_bridge.quotient_expression}`."
                        )
                        if (
                            quotient_followup_bridge.quotient_first_failure_power is None
                            or quotient_followup_bridge.quotient_first_failure_coeff is None
                        ):
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient `{quotient_followup_bridge.quotient_label}`: "
                                "matches `1` through the checked truncation."
                            )
                        else:
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient `{quotient_followup_bridge.quotient_label}`: "
                                f"`{quotient_followup_bridge.quotient_label} - 1` first fails at "
                                f"`{series_symbol}^{quotient_followup_bridge.quotient_first_failure_power}` with coefficient "
                                f"`{_format_expr(quotient_followup_bridge.quotient_first_failure_coeff)}`."
                            )
                        quotient_followup_bridge_scan = quotient_followup_bridge.quotient_scan
                        quotient_followup_bridge_self_product_hits = [
                            scan
                            for scan in quotient_followup_bridge_scan.self_quotient_product_scans
                            if scan.relation is not None
                        ]
                        quotient_followup_bridge_eta_hits = [
                            scan
                            for scan in quotient_followup_bridge_scan.eta_scans
                            if scan.relation is not None
                        ]
                        quotient_followup_bridge_modular_hits = [
                            scan
                            for scan in quotient_followup_bridge_scan.modular_unit_eta_scans
                            if scan.relation is not None
                        ]
                        quotient_followup_bridge_plus_hits = [
                            scan
                            for scan in quotient_followup_bridge_scan.self_plus_pochhammer_scans
                            if scan.relation is not None
                        ]
                        quotient_followup_bridge_plus_eta_hits = [
                            scan
                            for scan in quotient_followup_bridge_scan.self_plus_pochhammer_eta_scans
                            if scan.relation is not None
                        ]
                        lines.append(
                            f"- Weber quotient-follow-up bridge polynomial boxes: "
                            f"`{len(quotient_followup_polynomial_hits)}` / `{len(quotient_followup_bridge.polynomial_scans)}` hit boxes."
                        )
                        lines.append(
                            "- Weber quotient-follow-up bridge fractional-linear box: "
                            + (
                                "`1` / `1` hit boxes."
                                if quotient_followup_bridge.fractional_linear_relation is not None
                                else "`0` / `1` hit boxes."
                            )
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient self-polynomial uniqueness boxes: "
                            f"`{len(quotient_followup_bridge_scan.self_polynomial_scan.hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.self_polynomial_scan.moduli_checked) * len(quotient_followup_bridge_scan.self_polynomial_scan.fg_degree_values) * len(quotient_followup_bridge_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient self-fractional-linear uniqueness boxes: "
                            f"`{len(quotient_followup_bridge_scan.self_fractional_linear_scan.hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.self_fractional_linear_scan.moduli_checked) * len(quotient_followup_bridge_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient self-quotient finite-product boxes: "
                            f"`{len(quotient_followup_bridge_self_product_hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.self_quotient_product_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient eta templates: "
                            f"`{len(quotient_followup_bridge_eta_hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.eta_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient modular-unit / eta templates: "
                            f"`{len(quotient_followup_bridge_modular_hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.modular_unit_eta_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient plus-Pochhammer templates: "
                            f"`{len(quotient_followup_bridge_plus_hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.self_plus_pochhammer_scans)}` hit boxes."
                        )
                        lines.append(
                            f"- Weber quotient-follow-up bridge quotient plus-Pochhammer + eta templates: "
                            f"`{len(quotient_followup_bridge_plus_eta_hits)}` / "
                            f"`{len(quotient_followup_bridge_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                        )
                        if quotient_followup_bridge.quotient_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    f"Weber quotient-follow-up bridge quotient `{quotient_followup_bridge.quotient_label}`"
                                ),
                                gg_scan=quotient_followup_bridge.quotient_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                        if quotient_followup_bridge_scan.normalized_followup is not None:
                            quotient_followup_bridge_followup = quotient_followup_bridge_scan.normalized_followup
                            quotient_followup_bridge_followup_self_product_hits = [
                                scan
                                for scan in quotient_followup_bridge_followup.self_quotient_product_scans
                                if scan.relation is not None
                            ]
                            quotient_followup_bridge_followup_eta_hits = [
                                scan
                                for scan in quotient_followup_bridge_followup.eta_scans
                                if scan.relation is not None
                            ]
                            quotient_followup_bridge_followup_modular_hits = [
                                scan
                                for scan in quotient_followup_bridge_followup.modular_unit_eta_scans
                                if scan.relation is not None
                            ]
                            quotient_followup_bridge_followup_plus_hits = [
                                scan
                                for scan in quotient_followup_bridge_followup.self_plus_pochhammer_scans
                                if scan.relation is not None
                            ]
                            quotient_followup_bridge_followup_plus_eta_hits = [
                                scan
                                for scan in quotient_followup_bridge_followup.self_plus_pochhammer_eta_scans
                                if scan.relation is not None
                            ]
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized follow-up `{quotient_followup_bridge_followup.label}`: "
                                f"`{quotient_followup_bridge_followup.expression}`."
                            )
                            if (
                                quotient_followup_bridge_followup.first_failure_power is None
                                or quotient_followup_bridge_followup.first_failure_coeff is None
                            ):
                                lines.append(
                                    f"- Weber quotient-follow-up bridge quotient normalized follow-up `{quotient_followup_bridge_followup.label}`: "
                                    "matches `1` through the checked truncation."
                                )
                            else:
                                lines.append(
                                    f"- Weber quotient-follow-up bridge quotient normalized follow-up `{quotient_followup_bridge_followup.label}`: "
                                    f"`{quotient_followup_bridge_followup.label} - 1` first fails at "
                                    f"`{series_symbol}^{quotient_followup_bridge_followup.first_failure_power}` with coefficient "
                                    f"`{_format_expr(quotient_followup_bridge_followup.first_failure_coeff)}`."
                                )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized self-polynomial uniqueness boxes: "
                                f"`{len(quotient_followup_bridge_followup.self_polynomial_scan.hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.self_polynomial_scan.moduli_checked) * len(quotient_followup_bridge_followup.self_polynomial_scan.fg_degree_values) * len(quotient_followup_bridge_followup.self_polynomial_scan.t_degree_values)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized self-fractional-linear uniqueness boxes: "
                                f"`{len(quotient_followup_bridge_followup.self_fractional_linear_scan.hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.self_fractional_linear_scan.moduli_checked) * len(quotient_followup_bridge_followup.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized self-quotient finite-product boxes: "
                                f"`{len(quotient_followup_bridge_followup_self_product_hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.self_quotient_product_scans)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized eta templates: "
                                f"`{len(quotient_followup_bridge_followup_eta_hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.eta_scans)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized modular-unit / eta templates: "
                                f"`{len(quotient_followup_bridge_followup_modular_hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.modular_unit_eta_scans)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized plus-Pochhammer templates: "
                                f"`{len(quotient_followup_bridge_followup_plus_hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.self_plus_pochhammer_scans)}` hit boxes."
                            )
                            lines.append(
                                f"- Weber quotient-follow-up bridge quotient normalized plus-Pochhammer + eta templates: "
                                f"`{len(quotient_followup_bridge_followup_plus_eta_hits)}` / "
                                f"`{len(quotient_followup_bridge_followup.self_plus_pochhammer_eta_scans)}` hit boxes."
                            )
                            if quotient_followup_bridge_followup.named_gg_descendant_preview is not None:
                                preview = quotient_followup_bridge_followup.named_gg_descendant_preview
                                lines.append(
                                    f"- Weber quotient-follow-up bridge quotient normalized follow-up `{quotient_followup_bridge_followup.label}` "
                                    f"odd-prime descendant preview: direct ladder `{', '.join(preview.direct_labels)}`; "
                                    f"quotient ladder `{', '.join(preview.quotient_labels)}`."
                                )
                            if quotient_followup_bridge_followup.named_gg_descendant_focused_scan is not None:
                                _append_named_gg_descendant_focus_lines(
                                    lines,
                                    prefix=(
                                        "Weber quotient-follow-up bridge quotient normalized follow-up "
                                        f"`{quotient_followup_bridge_followup.label}`"
                                    ),
                                    descendant_scan=quotient_followup_bridge_followup.named_gg_descendant_focused_scan,
                                )
                        if quotient_followup_bridge.quotient_followup_named_gg_modular_equation_scan is not None:
                            _append_named_gg_bridge_lines(
                                lines,
                                prefix=(
                                    "Weber quotient-follow-up bridge quotient normalized follow-up "
                                    f"`{quotient_followup_bridge.quotient_scan.normalized_followup.label}`"
                                ),
                                gg_scan=quotient_followup_bridge.quotient_followup_named_gg_modular_equation_scan,
                                series_symbol=series_symbol,
                            )
                lines.append(
                    f"- Weber residual quotient diagnostic `{bridge_scan.quotient_label}`: "
                    f"`{bridge_scan.quotient_expression}`."
                )
                if (
                    bridge_scan.quotient_first_failure_power is None
                    or bridge_scan.quotient_first_failure_coeff is None
                ):
                    lines.append(
                        f"- Weber residual quotient diagnostic `{bridge_scan.quotient_label}`: "
                        "matches `1` through the checked truncation."
                    )
                else:
                    lines.append(
                        f"- Weber residual quotient diagnostic `{bridge_scan.quotient_label}`: "
                        f"`{bridge_scan.quotient_expression} - 1` first fails at "
                        f"`{series_symbol}^{bridge_scan.quotient_first_failure_power}` with coefficient "
                        f"`{_format_expr(bridge_scan.quotient_first_failure_coeff)}`."
                    )
                lines.append(
                    f"- Weber residual quotient self-polynomial uniqueness boxes: "
                    f"`{len(bridge_scan.quotient_self_polynomial_scan.hits)}` / "
                    f"`{len(bridge_scan.quotient_self_polynomial_scan.moduli_checked) * len(bridge_scan.quotient_self_polynomial_scan.fg_degree_values) * len(bridge_scan.quotient_self_polynomial_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient self-fractional-linear uniqueness boxes: "
                    f"`{len(bridge_scan.quotient_self_fractional_linear_scan.hits)}` / "
                    f"`{len(bridge_scan.quotient_self_fractional_linear_scan.moduli_checked) * len(bridge_scan.quotient_self_fractional_linear_scan.t_degree_values)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient self-quotient finite-product boxes: "
                    f"`{len(quotient_self_product_hits)}` / "
                    f"`{len(bridge_scan.quotient_self_quotient_product_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient eta templates: `{len(quotient_eta_hits)}` / "
                    f"`{len(bridge_scan.quotient_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient modular-unit / eta templates: "
                    f"`{len(quotient_modular_hits)}` / `{len(bridge_scan.quotient_modular_unit_eta_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient plus-Pochhammer templates: "
                    f"`{len(quotient_plus_hits)}` / `{len(bridge_scan.quotient_self_plus_pochhammer_scans)}` hit boxes."
                )
                lines.append(
                    f"- Weber residual quotient plus-Pochhammer + eta templates: "
                    f"`{len(quotient_plus_eta_hits)}` / "
                    f"`{len(bridge_scan.quotient_self_plus_pochhammer_eta_scans)}` hit boxes."
                )
                if bridge_scan.normalized_followup is not None:
                    followup_scan = bridge_scan.normalized_followup
                    followup_self_product_hits = [
                        scan
                        for scan in followup_scan.self_quotient_product_scans
                        if scan.relation is not None
                    ]
                    followup_eta_hits = [
                        scan for scan in followup_scan.eta_scans if scan.relation is not None
                    ]
                    followup_modular_hits = [
                        scan
                        for scan in followup_scan.modular_unit_eta_scans
                        if scan.relation is not None
                    ]
                    followup_plus_hits = [
                        scan
                        for scan in followup_scan.self_plus_pochhammer_scans
                        if scan.relation is not None
                    ]
                    followup_plus_eta_hits = [
                        scan
                        for scan in followup_scan.self_plus_pochhammer_eta_scans
                        if scan.relation is not None
                    ]
                    lines.append(
                        f"- Weber residual normalized follow-up `{followup_scan.label}`: "
                        f"`{followup_scan.expression}`."
                    )
                    if (
                        followup_scan.first_failure_power is None
                        or followup_scan.first_failure_coeff is None
                    ):
                        lines.append(
                            f"- Weber residual normalized follow-up `{followup_scan.label}`: "
                            "matches `1` through the checked truncation."
                        )
                    else:
                        lines.append(
                            f"- Weber residual normalized follow-up `{followup_scan.label}`: "
                            f"`{followup_scan.label} - 1` first fails at "
                            f"`{series_symbol}^{followup_scan.first_failure_power}` with coefficient "
                            f"`{_format_expr(followup_scan.first_failure_coeff)}`."
                        )
                    lines.append(
                        f"- Weber residual normalized self-polynomial uniqueness boxes: "
                        f"`{len(followup_scan.self_polynomial_scan.hits)}` / "
                        f"`{len(followup_scan.self_polynomial_scan.moduli_checked) * len(followup_scan.self_polynomial_scan.fg_degree_values) * len(followup_scan.self_polynomial_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized self-fractional-linear uniqueness boxes: "
                        f"`{len(followup_scan.self_fractional_linear_scan.hits)}` / "
                        f"`{len(followup_scan.self_fractional_linear_scan.moduli_checked) * len(followup_scan.self_fractional_linear_scan.t_degree_values)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized self-quotient finite-product boxes: "
                        f"`{len(followup_self_product_hits)}` / "
                        f"`{len(followup_scan.self_quotient_product_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized eta templates: `{len(followup_eta_hits)}` / "
                        f"`{len(followup_scan.eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized modular-unit / eta templates: "
                        f"`{len(followup_modular_hits)}` / `{len(followup_scan.modular_unit_eta_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized plus-Pochhammer templates: "
                        f"`{len(followup_plus_hits)}` / `{len(followup_scan.self_plus_pochhammer_scans)}` hit boxes."
                    )
                    lines.append(
                        f"- Weber residual normalized plus-Pochhammer + eta templates: "
                        f"`{len(followup_plus_eta_hits)}` / `{len(followup_scan.self_plus_pochhammer_eta_scans)}` hit boxes."
                    )
                    if followup_scan.named_gg_modular_equation_scan is not None:
                        _append_named_gg_bridge_lines(
                            lines,
                            prefix=f"Weber residual normalized follow-up `{followup_scan.label}`",
                            gg_scan=followup_scan.named_gg_modular_equation_scan,
                            series_symbol=series_symbol,
                        )
                    if followup_scan.named_gg_descendant_preview is not None:
                        preview = followup_scan.named_gg_descendant_preview
                        lines.append(
                            f"- Weber residual normalized follow-up `{followup_scan.label}` "
                            f"odd-prime descendant preview: direct ladder `{', '.join(preview.direct_labels)}`; "
                            f"quotient ladder `{', '.join(preview.quotient_labels)}`."
                        )
                    if followup_scan.named_gg_descendant_focused_scan is not None:
                        _append_named_gg_descendant_focus_lines(
                            lines,
                            prefix=f"Weber residual normalized follow-up `{followup_scan.label}`",
                            descendant_scan=followup_scan.named_gg_descendant_focused_scan,
                        )
                    coordinate_followup_scan = bridge_scan.quotient_coordinate_template_scan.normalized_followup
                    if (
                        coordinate_followup_scan is not None
                        and coordinate_followup_scan.named_gg_descendant_focused_scan is not None
                        and followup_scan.named_gg_descendant_focused_scan is not None
                    ):
                        left_direct_hits, left_direct_total = _named_prefix_box_scan_hit_count(
                            coordinate_followup_scan.named_gg_descendant_focused_scan.direct_scans
                        )
                        right_direct_hits, right_direct_total = _named_prefix_box_scan_hit_count(
                            followup_scan.named_gg_descendant_focused_scan.direct_scans
                        )
                        left_quotient_hits, left_quotient_total = _named_prefix_box_scan_hit_count(
                            coordinate_followup_scan.named_gg_descendant_focused_scan.quotient_scans
                        )
                        right_quotient_hits, right_quotient_total = _named_prefix_box_scan_hit_count(
                            followup_scan.named_gg_descendant_focused_scan.quotient_scans
                        )
                        lines.append(
                            f"- Weber odd-prime descendant comparison `{coordinate_followup_scan.label}` vs `{followup_scan.label}`: "
                            f"shared direct ladder `{', '.join(coordinate_followup_scan.named_gg_descendant_focused_scan.direct_labels)}` "
                            f"and shared quotient ladder `{', '.join(coordinate_followup_scan.named_gg_descendant_focused_scan.quotient_labels)}`."
                        )
                        lines.append(
                            f"- Weber odd-prime descendant direct micro-box comparison: "
                            f"`{coordinate_followup_scan.label}` has `{left_direct_hits}` / `{left_direct_total}` hits, "
                            f"`{followup_scan.label}` has `{right_direct_hits}` / `{right_direct_total}` hits."
                        )
                        lines.append(
                            f"- Weber odd-prime descendant quotient micro-box comparison: "
                            f"`{coordinate_followup_scan.label}` has `{left_quotient_hits}` / `{left_quotient_total}` hits, "
                            f"`{followup_scan.label}` has `{right_quotient_hits}` / `{right_quotient_total}` hits."
                        )
            if sample.gg_modular_equation_scan is not None:
                gg_scan = sample.gg_modular_equation_scan
                direct_prefix_summary = "; ".join(
                    (
                        _format_tail_prefix_summary(
                            gg_scan.polynomial_scans,
                            label="polynomial",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.multiplicative_scans,
                            label="multiplicative",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.fractional_linear_scans,
                            label="fractional-linear",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.two_layer_fractional_linear_scans,
                            label="two-layer fractional-linear",
                            hit_predicate=lambda item: item.total_hits > 0,
                        ),
                    )
                )
                quotient_prefix_summary = "; ".join(
                    (
                        _format_tail_prefix_summary(
                            gg_scan.quotient_polynomial_scans,
                            label="polynomial",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.quotient_multiplicative_scans,
                            label="multiplicative",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.quotient_fractional_linear_scans,
                            label="fractional-linear",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.quotient_two_layer_fractional_linear_scans,
                            label="two-layer fractional-linear",
                            hit_predicate=lambda item: item.total_hits > 0,
                        ),
                    )
                )
                mixed_prefix_summary = "; ".join(
                    (
                        _format_tail_prefix_summary(
                            gg_scan.mixed_quotient_polynomial_scans,
                            label="polynomial",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.mixed_quotient_multiplicative_scans,
                            label="multiplicative",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.mixed_quotient_fractional_linear_scans,
                            label="fractional-linear",
                            hit_predicate=lambda item: item.relation is not None,
                        ),
                        _format_tail_prefix_summary(
                            gg_scan.mixed_quotient_two_layer_fractional_linear_scans,
                            label="two-layer fractional-linear",
                            hit_predicate=lambda item: item.total_hits > 0,
                        ),
                    )
                )
                exact_template_summary = (
                    f"`{len(gg_scan.hit_templates)}` / `{len(gg_scan.checked_templates)}` exact template hits"
                )
                if gg_scan.hit_templates:
                    exact_template_summary += (
                        f" ({', '.join(f'`{label}`' for label in gg_scan.hit_templates)})"
                    )
                direct_exact_summary = (
                    f"`{len(gg_scan.exact_polynomial_template_hits)}` / `{len(gg_scan.exact_polynomial_template_labels)}` exact Chan--Huang direct hits"
                )
                if gg_scan.exact_polynomial_template_hits:
                    direct_exact_summary += (
                        f" ({', '.join(f'`{label}`' for label in gg_scan.exact_polynomial_template_hits)})"
                    )
                quotient_exact_summary = (
                    f"`{len(gg_scan.quotient_exact_polynomial_template_hits)}` / `{len(gg_scan.quotient_exact_polynomial_template_labels)}` exact Chan--Huang quotient-coordinate hits"
                )
                if gg_scan.quotient_exact_polynomial_template_hits:
                    quotient_exact_summary += (
                        f" ({', '.join(f'`{label}`' for label in gg_scan.quotient_exact_polynomial_template_hits)})"
                    )
                lines.extend(
                    [
                        f"- GG direct / reciprocal / quotient templates: {exact_template_summary}.",
                        f"- GG direct exact modular-equation templates: {direct_exact_summary}.",
                        f"- GG quotient exact modular-equation templates: {quotient_exact_summary}.",
                        "- GG narrow quotient-coordinate exact lane focuses on `Q_3` and `Q_4` before any broader quotient-prefix interpretation.",
                        "- GG exact quotient-coordinate obstruction witnesses: "
                        + "; ".join(
                            _format_exact_polynomial_obstruction(
                                obstruction,
                                series_symbol=series_symbol,
                            )
                            for obstruction in gg_scan.quotient_exact_polynomial_template_obstructions
                        )
                        + ".",
                        *[
                            "- GG weighted quotient-coordinate diagnostic "
                            + f"`{diagnostic.label} = {diagnostic.expression}`: "
                            + _format_weighted_coordinate_obstruction(
                                power=diagnostic.first_difference_power,
                                coeff=diagnostic.first_difference_coeff,
                                lhs=f"F - {diagnostic.label}",
                                series_symbol=series_symbol,
                            )
                            + "; "
                            + _format_weighted_coordinate_obstruction(
                                power=diagnostic.first_log_difference_power,
                                coeff=diagnostic.first_log_difference_coeff,
                                lhs=f"log(F) - ({diagnostic.log_expression})",
                                series_symbol=series_symbol,
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                        ],
                        *[
                            "- GG weighted correction "
                            + f"`{diagnostic.correction_expression}`: "
                            + _format_weighted_coordinate_obstruction(
                                power=diagnostic.correction_first_gap_power,
                                coeff=diagnostic.correction_first_gap_coeff,
                                lhs=f"{diagnostic.correction_expression} - 1",
                                series_symbol=series_symbol,
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                        ],
                        *[
                            "- GG normalized weighted correction "
                            + f"`{diagnostic.normalized_correction_label}`: "
                            + _format_gap_normalization_formula(
                                source_variable=diagnostic.correction_expression,
                                target_variable=diagnostic.normalized_correction_label,
                                gap=diagnostic.normalized_correction_gap,
                                series_symbol=series_symbol,
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.normalized_correction_label is not None
                            and diagnostic.normalized_correction_gap is not None
                        ],
                        *[
                            "- GG second normalized weighted correction "
                            + f"`{diagnostic.second_normalized_correction_label}`: "
                            + _format_gap_normalization_formula(
                                source_variable=diagnostic.normalized_correction_label,
                                target_variable=diagnostic.second_normalized_correction_label,
                                gap=diagnostic.second_normalized_correction_gap,
                                series_symbol=series_symbol,
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.normalized_correction_label is not None
                            and diagnostic.second_normalized_correction_label is not None
                            and diagnostic.second_normalized_correction_gap is not None
                        ],
                        *[
                            f"- GG weighted quotient-coordinate `{diagnostic.label}`: no degree-`<= 2` polynomial or one-coordinate fractional-linear closure was found."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.polynomial_degree1_relation is None
                            and diagnostic.polynomial_degree2_relation is None
                            and diagnostic.fractional_linear_relation is None
                        ],
                        *[
                            f"- GG weighted correction `{diagnostic.correction_expression}`: no eta-quotient or modular-unit / eta hit was found in the checked small correction boxes."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if not any(scan.relation is not None for scan in diagnostic.correction_eta_scans)
                            and not any(scan.relation is not None for scan in diagnostic.correction_modular_unit_eta_scans)
                        ],
                        *[
                            f"- GG normalized weighted correction `{diagnostic.normalized_correction_label}`: no eta-quotient or modular-unit / eta hit was found in the checked small normalized boxes."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.normalized_correction_label is not None
                            and not any(scan.relation is not None for scan in diagnostic.normalized_correction_eta_scans)
                            and not any(scan.relation is not None for scan in diagnostic.normalized_correction_modular_unit_eta_scans)
                        ],
                        *[
                            f"- GG normalized weighted correction `{diagnostic.normalized_correction_label}`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.normalized_correction_label is not None
                            and not _flatten_source_family_eta_hits(diagnostic.normalized_correction_source_family_eta_scans)
                        ],
                        *[
                            f"- GG second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no eta-quotient or modular-unit / eta hit was found in the checked small second-normalized boxes."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.second_normalized_correction_label is not None
                            and not any(scan.relation is not None for scan in diagnostic.second_normalized_correction_eta_scans)
                            and not any(scan.relation is not None for scan in diagnostic.second_normalized_correction_modular_unit_eta_scans)
                        ],
                        *[
                            f"- GG second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no one-core source-family eta-correction hit was found across `RR`, `GG` raw/quotient bases."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.second_normalized_correction_label is not None
                            and not _flatten_source_family_eta_hits(diagnostic.second_normalized_correction_source_family_eta_scans)
                        ],
                        *[
                            f"- GG second normalized weighted correction `{diagnostic.second_normalized_correction_label}`: no explicit GG transform-template eta-correction hit was found in the checked small boxes."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.second_normalized_correction_label is not None
                            and not any(scan.hits for scan in diagnostic.second_normalized_correction_explicit_transform_eta_scans)
                        ],
                        *[
                            "- GG second normalized weighted correction "
                            + f"`{diagnostic.second_normalized_correction_label}` quotient-coordinate prefixes: "
                            + "; ".join(
                                (
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_quotient_polynomial_scans,
                                        label="polynomial",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_quotient_multiplicative_scans,
                                        label="multiplicative",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_quotient_fractional_linear_scans,
                                        label="fractional-linear",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_quotient_two_layer_fractional_linear_scans,
                                        label="two-layer fractional-linear",
                                        hit_predicate=lambda item: item.total_hits > 0,
                                    ),
                                )
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.second_normalized_correction_label is not None
                        ],
                        *[
                            "- GG second normalized weighted correction "
                            + f"`{diagnostic.second_normalized_correction_label}` mixed quotient-coordinate prefixes: "
                            + "; ".join(
                                (
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_mixed_quotient_polynomial_scans,
                                        label="polynomial",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_mixed_quotient_multiplicative_scans,
                                        label="multiplicative",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_mixed_quotient_fractional_linear_scans,
                                        label="fractional-linear",
                                        hit_predicate=lambda item: item.relation is not None,
                                    ),
                                    _format_tail_prefix_summary(
                                        diagnostic.second_normalized_correction_mixed_quotient_two_layer_fractional_linear_scans,
                                        label="two-layer fractional-linear",
                                        hit_predicate=lambda item: item.total_hits > 0,
                                    ),
                                )
                            )
                            + "."
                            for diagnostic in gg_scan.weighted_coordinate_diagnostics
                            if diagnostic.second_normalized_correction_label is not None
                        ],
                        f"- GG direct prefixes: {direct_prefix_summary}.",
                        f"- GG quotient-coordinate prefixes: {quotient_prefix_summary}.",
                        f"- GG mixed quotient-coordinate prefixes: {mixed_prefix_summary}.",
                    ]
                )
            lines.append("")

        hit_count = sum(1 for sample in sample_scans if _flatten_source_family_eta_hits(sample.source_family_eta_scans))
        direct_eta_hit_count = sum(1 for sample in sample_scans if any(scan.relation is not None for scan in sample.direct_eta_scans))
        direct_modular_unit_hit_count = sum(
            1
            for sample in sample_scans
            if any(scan.relation is not None for scan in sample.direct_modular_unit_eta_scans)
        )
        gg_hit_count = sum(
            1
            for sample in sample_scans
            if sample.gg_modular_equation_scan is not None
            and _gg_modular_equation_scan_has_hit(sample.gg_modular_equation_scan)
        )
        gg_exact_quotient_hit_count = sum(
            1
            for sample in sample_scans
            if sample.gg_modular_equation_scan is not None
            and _gg_modular_equation_scan_has_exact_quotient_hit(sample.gg_modular_equation_scan)
        )
        morton_hit_count = sum(
            1
            for sample in sample_scans
            if sample.morton_periodic_point_scan is not None
            and _morton_periodic_point_scan_has_hit(sample.morton_periodic_point_scan)
        )
        morton_weber_hit_count = sum(
            1
            for sample in sample_scans
            if sample.morton_periodic_point_scan is not None
            and any(
                item.hit
                for coordinate_scan in sample.morton_periodic_point_scan.named_coordinate_scans
                for item in coordinate_scan.template_results
            )
        )
        weber_g_class_invariant_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_g_class_invariant_scan is not None
            and _weber_class_invariant_scan_has_hit(sample.weber_g_class_invariant_scan)
        )
        weber_p_class_invariant_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_p_class_invariant_scan is not None
            and _weber_class_invariant_scan_has_hit(sample.weber_p_class_invariant_scan)
        )
        weber_classical_product_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_residual_bridge_scan is not None
            and _constant_one_series_scan_has_hit(
                sample.weber_residual_bridge_scan.classical_product_coordinate_scan
            )
        )
        weber_anchor_canonical_j_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_residual_bridge_scan is not None
            and _constant_one_series_scan_has_hit(
                sample.weber_residual_bridge_scan.anchor_canonical_j_coordinate_scan
            )
        )
        weber_residual_quotient_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_residual_bridge_scan is not None
            and (
                any(scan.relation is not None for scan in sample.weber_residual_bridge_scan.quotient_eta_scans)
                or bool(sample.weber_residual_bridge_scan.quotient_self_polynomial_scan.hits)
                or bool(sample.weber_residual_bridge_scan.quotient_self_fractional_linear_scan.hits)
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.quotient_self_quotient_product_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.quotient_modular_unit_eta_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.quotient_self_plus_pochhammer_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.quotient_self_plus_pochhammer_eta_scans
                )
            )
        )
        weber_residual_followup_hit_count = sum(
            1
            for sample in sample_scans
            if sample.weber_residual_bridge_scan is not None
            and sample.weber_residual_bridge_scan.normalized_followup is not None
            and (
                bool(sample.weber_residual_bridge_scan.normalized_followup.self_polynomial_scan.hits)
                or bool(sample.weber_residual_bridge_scan.normalized_followup.self_fractional_linear_scan.hits)
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.normalized_followup.self_quotient_product_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.normalized_followup.eta_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.normalized_followup.modular_unit_eta_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.normalized_followup.self_plus_pochhammer_scans
                )
                or any(
                    scan.relation is not None
                    for scan in sample.weber_residual_bridge_scan.normalized_followup.self_plus_pochhammer_eta_scans
                )
            )
        )
        lines.extend(
            [
                "## Tail Verdict",
                "",
                f"- Samples checked: `{len(sample_scans)}`",
                f"- Source-core eta hits found: `{hit_count}`",
                f"- Direct eta-quotient sample hits found: `{direct_eta_hit_count}`",
                f"- Direct modular-unit / eta sample hits found: `{direct_modular_unit_hit_count}`",
                f"- GG/Weber modular-equation sample hits found: `{gg_hit_count}`",
                f"- GG exact quotient-coordinate sample hits found: `{gg_exact_quotient_hit_count}`",
                f"- Morton periodic-point / algebraic-function sample hits found: `{morton_hit_count}`",
                f"- Morton named-coordinate sample hits found: `{morton_weber_hit_count}`",
                f"- Weber g-class-invariant sample hits found: `{weber_g_class_invariant_hit_count}`",
                f"- Weber G-class-invariant sample hits found: `{weber_p_class_invariant_hit_count}`",
                f"- Classical Weber `f2` tri-product sample hits found: `{weber_classical_product_hit_count}`",
                f"- Canonical Weber anchor `j`-side sample hits found: `{weber_anchor_canonical_j_hit_count}`",
                f"- Weber residual-quotient sample hits found: `{weber_residual_quotient_hit_count}`",
                f"- Weber residual-follow-up sample hits found: `{weber_residual_followup_hit_count}`",
                "- Current reading: the tail-family ladder remains structurally informative, but the sampled `U(x)` objects and their deeper gap residuals still do not collapse into the first direct eta / modular-unit boxes, the first nearby one-core eta-correction boxes, the direct Morton algebraic-function templates, the first Weber-Schlafli coordinate / companion templates, the first Ramanujan-Weber class-invariant compression boxes, the classical Weber `f2` tri-product box, the anchor-derived canonical Weber `j` box, the focused Weber residual-quotient box, the normalized Weber residual-follow-up box, or the first literature-driven GG/Weber modular-equation boxes.",
                "",
                f"- Build elapsed seconds: `{perf_counter() - build_started_at:.2f}`",
                "",
            ]
        )

    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_candidate_tail_operator_note(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    depth: int = 40,
    series_order: int = 36,
    tail_stages: tuple[int, ...] = (3, 4, 5),
    max_gap_depth: int = 3,
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
    profile_depth = min(depth, 24 if smoke else depth)
    profile_order = min(series_order, 24 if smoke else series_order)
    active = _series_active_exponents(record.template) + _series_active_exponents(benchmark.canonical_template)
    step = 0
    for value in active:
        step = gcd(step, value)
    if step <= 0:
        step = 1

    reduced_candidate = record.template
    series_symbol = "q"
    variable_label = "q"
    if step > 1:
        maybe_candidate = reduce_template_by_step(record.template, step)
        maybe_benchmark = reduce_template_by_step(benchmark.canonical_template, step)
        if maybe_candidate is not None and maybe_benchmark is not None:
            reduced_candidate = maybe_candidate
            series_symbol = "t"
            variable_label = f"t = q^{step}"

    reduced_bridge_depth = min(profile_depth, 8 if smoke else 12)
    reduced_bridge_order = min(profile_order, 24 if smoke else 30)
    output_file = Path(output_path)
    build_started_at = perf_counter()
    lines = [
        f"# Tail-Operator Note: `{record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{record.id}`",
        f"- Closest benchmark: `{record.closest_benchmark}`",
        f"- Variable view: `{variable_label}`",
        f"- Tail stages checked: `{', '.join(str(stage) for stage in tail_stages)}`",
        f"- Max gap depth checked: `{max_gap_depth}`",
        "",
    ]

    reduced_bridge_error: str | None = None
    reduced_reciprocal_witness = None
    reduced_tail_transfer_equation: ReducedTailTransferEquation | None = None
    try:
        reduced_reciprocal_witness, _ = _reduced_reciprocal_bridge(
            template=reduced_candidate,
            symbol=sp.Symbol(series_symbol),
            depth=reduced_bridge_depth,
            order=reduced_bridge_order,
        )
        reduced_tail_transfer_equation = detect_reduced_tail_transfer_equation(
            reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
            symbol=sp.Symbol(series_symbol),
        )
    except Exception as exc:
        reduced_bridge_error = str(exc)

    if reduced_bridge_error is not None or reduced_reciprocal_witness is None or reduced_tail_transfer_equation is None:
        lines.extend(
            [
                "Tail-operator setup failed:",
                "",
                "```text",
                reduced_bridge_error or "unknown reduced-tail setup failure",
                "```",
                "",
            ]
        )
        output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    relation_lines = _format_reduced_tail_transfer_equation(reduced_tail_transfer_equation)
    lines.extend(
        [
            "## Exact Tail Family",
            "",
            "- The operator lane starts from the same exact stationary tail law rather than from a fresh anonymous ansatz.",
            "",
            "```text",
            *relation_lines,
            "```",
            "",
            "## Operator Lane",
            "",
            "We now ask a more recurrence-first question on the sampled tail ladder:",
            "",
            "```text",
            "A_0(t) + A_1(t)*Y(t) + A_2(t)*Y(t^m) + A_3(t)*Y(t^(m^2)) = 0",
            "```",
            "",
            "- This is the current affine q-difference / Mahler-style operator box.",
            "- The goal is not another family label; it is a compact operator statement that could later support an operator-factorization or uniqueness proof.",
            "",
        ]
    )

    sample_scans = scan_tail_family_source_eta_ladder(
        reduced_coeffs=reduced_reciprocal_witness.reduction.reduced_coeffs,
        symbol=sp.Symbol(series_symbol),
        ordered_base_families=(),
        start_stages=tail_stages,
        max_gap_depth=max_gap_depth,
        order=reduced_bridge_order,
        powers=(2, 3),
        eta_levels=(),
        max_abs_exponent=4 if smoke else 6,
        gg_benchmark_name=None,
        gg_base_series=None,
        morton_order=min(reduced_bridge_order, 18 if smoke else 24),
    )

    total_hits = 0
    for sample in sample_scans:
        sample_order = min(len(sample.series), 18 if smoke else 24)
        mahler_scan = scan_self_mahler_linear_relations(
            target_series=list(sample.series[:sample_order]),
            moduli=(2,) if smoke else (2, 3),
            levels_checked=(2,) if smoke else (2, 3),
            order=sample_order,
            t_degree_values=(1, 2) if smoke else (1, 2, 3),
        )
        total_hits += len(mahler_scan.hits)
        lines.extend(
            [
                f"### `{sample.label}`",
                "",
                f"- Start stage: `{sample.start_stage}`",
                f"- Gap depth: `{sample.gap_depth}`",
                f"- State: `{_format_expr(sample.state_expr)}`",
                "",
                "```text",
                sample.expression,
                "```",
                "",
                f"- Moduli checked: {', '.join(f'`m={value}`' for value in mahler_scan.moduli_checked) if mahler_scan.moduli_checked else '`none`'}",
                f"- Recurrence depths checked: {', '.join(f'`levels={value}`' for value in mahler_scan.levels_checked) if mahler_scan.levels_checked else '`none`'}",
                f"- Polynomial t-degrees checked: {', '.join(f'`deg_t={value}`' for value in mahler_scan.t_degree_values) if mahler_scan.t_degree_values else '`none`'}",
            ]
        )
        if not mahler_scan.hits:
            lines.append("- No affine q-difference / Mahler operator hit was found in the scanned box.")
            lines.append("")
            continue
        lines.append(f"- Affine q-difference / Mahler operator hits found: `{len(mahler_scan.hits)}`")
        lines.append("")
        for hit in mahler_scan.hits:
            relation_series_by_variable = {
                "T": _t_series(order=sample_order),
                "F": list(sample.series[:sample_order]),
            }
            for variable in hit.relation.variables:
                if variable.startswith("G"):
                    relation_series_by_variable[variable] = benchmark_power_substitution_series(
                        list(sample.series[:sample_order]),
                        power=int(variable.removeprefix("G")),
                        order=sample_order,
                    )
            residual = _relation_residual_series(
                hit.relation,
                series_by_variable=relation_series_by_variable,
                order=sample_order,
            )
            residual_ok = all(sp.simplify(value) == 0 for value in residual)
            lines.extend(
                [
                    f"- Modulus `m={hit.modulus}`, recurrence depth `{hit.levels}`, `deg_t={hit.max_t_degree}`:",
                    "",
                    "```text",
                    _format_self_mahler_linear_relation(
                        hit.relation,
                        target_variable=sample.label,
                        series_symbol=series_symbol,
                    ),
                    "```",
                    "",
                    f"  Verified by exact series re-expansion modulo `{series_symbol}^{sample_order}`: `{residual_ok}`",
                    "",
                ]
            )

    lines.extend(
        [
            "## Operator Verdict",
            "",
            f"- Samples checked: `{len(sample_scans)}`",
            f"- Total affine q-difference / Mahler hits found: `{total_hits}`",
            "- Current reading: the exact tail law is strong enough to justify an operator-first endgame, but the first low-degree affine q-difference box is still mostly a diagnostic lane rather than a final theorem.",
            "",
            f"- Build elapsed seconds: `{perf_counter() - build_started_at:.2f}`",
            "",
        ]
    )

    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
