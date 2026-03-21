from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import sympy as sp

from ramanujan_discovery.benchmarks import BenchmarkDefinition, get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.research import (
    BauerMuirDirectObstruction,
    ContinuedFractionCoeffs,
    ConvergentFactorEquivalenceWitness,
    HeineCor2CFContractionObstruction,
    Page43F2EquivalenceObstruction,
    Page43F2UnitLambdaShiftEquivalenceObstruction,
    Page43F4EquivalenceObstruction,
    Page43F4UnitLambdaShiftEquivalenceObstruction,
    Page43ParameterHit,
    Page43PolynomialSinglePrefactorObstruction,
    SubsequenceContractionHit,
    _template_reciprocal_coeffs,
    arithmetic_subsequence_contraction_search,
    convergent_factor_equivalence_witness,
    continued_fraction_convergents,
    direct_bauer_muir_obstruction,
    heine_cor2cf_a_zero_contraction_obstruction,
    page43_f2_zero_shift_equivalence_obstruction,
    page43_f2_unit_a_shift_equivalence_obstruction,
    page43_f2_unit_ab_shift_equivalence_obstruction,
    page43_f2_unit_ab_lambda_shift_equivalence_obstruction,
    page43_f2_unit_a_lambda_shift_equivalence_obstruction,
    page43_f2_unit_b_shift_equivalence_obstruction,
    page43_f2_unit_b_lambda_shift_equivalence_obstruction,
    page43_f2_unit_lambda_shift_equivalence_obstruction,
    page43_f4_zero_shift_equivalence_obstruction,
    page43_f4_unit_a_shift_equivalence_obstruction,
    page43_f4_unit_ab_shift_equivalence_obstruction,
    page43_f4_unit_ab_lambda_shift_equivalence_obstruction,
    page43_f4_unit_a_lambda_shift_equivalence_obstruction,
    page43_f4_unit_b_shift_equivalence_obstruction,
    page43_f4_unit_b_lambda_shift_equivalence_obstruction,
    page43_f4_unit_lambda_shift_equivalence_obstruction,
    page43_monomial_parameter_search,
    page43_rational_parameter_search,
    page43_zero_shift_polynomial_single_prefactor_obstructions,
    page43_zero_shift_reciprocal_single_prefactor_obstructions,
    parity_contraction_coeffs,
    reduce_template_by_step,
)
from ramanujan_discovery.storage import read_candidates


@dataclass(frozen=True)
class FormalizationBuildProfile:
    label: str
    factor_depth: int
    subsequence_stages: int
    max_page43_shift: int
    page43_stages: int
    page43_rational_max_nontrivial_profiles: int


@dataclass(frozen=True)
class FormalizationContext:
    record: CandidateRecord
    benchmark: BenchmarkDefinition
    step: int
    build_profile: str
    max_stride: int
    subsequence_stages: int
    page43_max_shift: int
    page43_stages: int
    page43_rational_max_shift: int
    page43_rational_stages: int
    page43_rational_max_nontrivial_profiles: int
    reduced_candidate: QCFTemplate | None
    reduced_benchmark: QCFTemplate | None
    target_b0: sp.Expr | None
    target_a_terms: list[sp.Expr]
    target_b_terms: list[sp.Expr]
    factor_depth: int
    factor_witness: ConvergentFactorEquivalenceWitness | None
    rr_direct: BauerMuirDirectObstruction | None
    cubic_direct: BauerMuirDirectObstruction | None
    cubic_odd: ContinuedFractionCoeffs | None
    cubic_even: ContinuedFractionCoeffs | None
    heine_cor2cf: HeineCor2CFContractionObstruction | None
    f2_equivalence: Page43F2EquivalenceObstruction | None
    f4_equivalence: Page43F4EquivalenceObstruction | None
    f2_unit_a_shift_equivalence: Page43F2EquivalenceObstruction | None
    f4_unit_a_shift_equivalence: Page43F4EquivalenceObstruction | None
    f2_unit_b_shift_equivalence: Page43F2EquivalenceObstruction | None
    f4_unit_b_shift_equivalence: Page43F4EquivalenceObstruction | None
    f2_unit_ab_shift_equivalence: Page43F2EquivalenceObstruction | None
    f4_unit_ab_shift_equivalence: Page43F4EquivalenceObstruction | None
    f2_unit_ab_lambda_shift_equivalence: Page43F2UnitLambdaShiftEquivalenceObstruction | None
    f4_unit_ab_lambda_shift_equivalence: Page43F4UnitLambdaShiftEquivalenceObstruction | None
    f2_unit_a_lambda_shift_equivalence: Page43F2UnitLambdaShiftEquivalenceObstruction | None
    f4_unit_a_lambda_shift_equivalence: Page43F4UnitLambdaShiftEquivalenceObstruction | None
    f2_unit_b_lambda_shift_equivalence: Page43F2UnitLambdaShiftEquivalenceObstruction | None
    f4_unit_b_lambda_shift_equivalence: Page43F4UnitLambdaShiftEquivalenceObstruction | None
    f2_unit_lambda_shift_equivalence: Page43F2UnitLambdaShiftEquivalenceObstruction | None
    f4_unit_lambda_shift_equivalence: Page43F4UnitLambdaShiftEquivalenceObstruction | None
    f2_polynomial_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction]
    f4_polynomial_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction]
    f2_reciprocal_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction]
    f4_reciprocal_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction]
    rr_subsequence_hits: list[SubsequenceContractionHit]
    cubic_subsequence_hits: list[SubsequenceContractionHit]
    f2_hits: list[Page43ParameterHit]
    f4_hits: list[Page43ParameterHit]
    f2_rational_hits: list[Page43ParameterHit]
    f4_rational_hits: list[Page43ParameterHit]
    looks_like_hero: bool


def _find_candidate(records: list[CandidateRecord], candidate_id: str) -> CandidateRecord:
    for record in records:
        if record.id == candidate_id:
            return record
    raise KeyError(f"unknown candidate id: {candidate_id}")


def _format_expr(expr) -> str:
    return str(sp.expand(expr)).replace("**", "^")


def _format_fraction_expr(expr) -> str:
    return str(sp.cancel(expr)).replace("**", "^")


def _format_symbolic_solution(
    solution: dict[sp.Symbol, sp.Expr] | None,
    *,
    symbols: tuple[sp.Symbol, ...],
) -> str:
    if not solution:
        return "no constant solution"
    return ", ".join(f"{symbol} = {_format_expr(solution[symbol])}" for symbol in symbols if symbol in solution)


def _format_prefactor_case_label(obstruction: Page43PolynomialSinglePrefactorObstruction) -> str:
    return (
        f"`phi_a = {obstruction.a_profile.replace(' / ', '/').replace(' + ', '+')}`, "
        f"`phi_b = {obstruction.b_profile.replace(' / ', '/').replace(' + ', '+')}`, "
        f"`phi_lambda = {obstruction.lambda_profile.replace(' / ', '/').replace(' + ', '+')}`"
    )


def _formalization_build_profile(*, smoke: bool) -> FormalizationBuildProfile:
    if smoke:
        return FormalizationBuildProfile(
            label="smoke",
            factor_depth=4,
            subsequence_stages=2,
            max_page43_shift=1,
            page43_stages=2,
            page43_rational_max_nontrivial_profiles=1,
        )
    return FormalizationBuildProfile(
        label="full",
        factor_depth=8,
        subsequence_stages=3,
        max_page43_shift=3,
        page43_stages=3,
        page43_rational_max_nontrivial_profiles=3,
    )


def _page43_rational_prefactor_max_shift(*, smoke: bool) -> int:
    return 1


def _power_term(scale: int, shift: int, step: int, variable: str) -> str:
    if shift == step:
        exponent = "n" if step == 1 else f"{step}n"
    else:
        exponent = f"{shift} + {step}(n-1)"

    if scale == 1:
        return f"{variable}^{exponent}"
    return f"{scale}*{variable}^{exponent}"


def _reciprocal_rule_lines(template: QCFTemplate, variable: str) -> list[str]:
    normalized = template.normalized()
    numerator_terms = [_power_term(normalized.numerator_scale, normalized.numerator_q_shift, normalized.numerator_q_step, variable)]
    if normalized.numerator_extra_scale != 0:
        numerator_terms.append(
            _power_term(
                normalized.numerator_extra_scale,
                normalized.numerator_extra_q_shift,
                normalized.numerator_extra_q_step,
                variable,
            )
        )
    numerator = " + ".join(numerator_terms)

    denominator = str(normalized.denominator_constant)
    if normalized.denominator_scale != 0:
        denominator += " + " + _power_term(
            normalized.denominator_scale,
            normalized.denominator_q_shift,
            normalized.denominator_q_step,
            variable,
        )

    return [
        f"b0 = {sp.Rational(normalized.base_denominator, normalized.top_constant)}",
        f"a_n = {numerator}",
        f"b_n = {denominator}",
    ]


def _build_formalization_context(
    *,
    input_path: str,
    candidate_id: str,
    max_stride: int,
    smoke: bool,
) -> FormalizationContext:
    records = read_candidates(input_path)
    record = _find_candidate(records, candidate_id)
    benchmark = get_benchmark(record.closest_benchmark)
    step = benchmark.canonical_template.numerator_q_step or 1
    profile = _formalization_build_profile(smoke=smoke)
    effective_max_stride = min(max_stride, 2) if smoke else max_stride

    reduced_candidate = reduce_template_by_step(record.template.normalized(), step=step)
    reduced_benchmark = reduce_template_by_step(benchmark.canonical_template.normalized(), step=step)
    page43_rational_max_shift = _page43_rational_prefactor_max_shift(smoke=smoke)

    if reduced_candidate is None or reduced_benchmark is None:
        return FormalizationContext(
            record=record,
            benchmark=benchmark,
            step=step,
            build_profile=profile.label,
            max_stride=effective_max_stride,
            subsequence_stages=profile.subsequence_stages,
            page43_max_shift=profile.max_page43_shift,
            page43_stages=profile.page43_stages,
            page43_rational_max_shift=page43_rational_max_shift,
            page43_rational_stages=profile.page43_stages,
            page43_rational_max_nontrivial_profiles=profile.page43_rational_max_nontrivial_profiles,
            reduced_candidate=reduced_candidate,
            reduced_benchmark=reduced_benchmark,
            target_b0=None,
            target_a_terms=[],
            target_b_terms=[],
            factor_depth=0,
            factor_witness=None,
            rr_direct=None,
            cubic_direct=None,
            cubic_odd=None,
            cubic_even=None,
            heine_cor2cf=None,
            f2_equivalence=None,
            f4_equivalence=None,
            f2_unit_a_shift_equivalence=None,
            f4_unit_a_shift_equivalence=None,
            f2_unit_b_shift_equivalence=None,
            f4_unit_b_shift_equivalence=None,
            f2_unit_ab_shift_equivalence=None,
            f4_unit_ab_shift_equivalence=None,
            f2_unit_ab_lambda_shift_equivalence=None,
            f4_unit_ab_lambda_shift_equivalence=None,
            f2_unit_a_lambda_shift_equivalence=None,
            f4_unit_a_lambda_shift_equivalence=None,
            f2_unit_b_lambda_shift_equivalence=None,
            f4_unit_b_lambda_shift_equivalence=None,
            f2_unit_lambda_shift_equivalence=None,
            f4_unit_lambda_shift_equivalence=None,
            f2_polynomial_prefactor_obstructions=[],
            f4_polynomial_prefactor_obstructions=[],
            f2_reciprocal_prefactor_obstructions=[],
            f4_reciprocal_prefactor_obstructions=[],
            rr_subsequence_hits=[],
            cubic_subsequence_hits=[],
            f2_hits=[],
            f4_hits=[],
            f2_rational_hits=[],
            f4_rational_hits=[],
            looks_like_hero=False,
        )

    t = sp.Symbol("t")
    target_b0, target_a_terms, target_b_terms = _template_reciprocal_coeffs(reduced_candidate, q=t, depth=4)
    factor_depth = profile.factor_depth
    factor_b0, factor_a_terms, factor_b_terms = _template_reciprocal_coeffs(reduced_candidate, q=t, depth=factor_depth)
    factor_witness = convergent_factor_equivalence_witness(
        b0=factor_b0,
        a_terms=factor_a_terms,
        b_terms=factor_b_terms,
    )
    looks_like_hero = (
        sp.simplify(target_b0 - 1) == 0
        and sp.simplify(target_a_terms[1] - (t + t**2)) == 0
        and sp.simplify(target_b_terms[1] - (1 + t)) == 0
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
    heine_cor2cf = None
    f2_equivalence = None
    f4_equivalence = None
    f2_unit_a_shift_equivalence = None
    f4_unit_a_shift_equivalence = None
    f2_unit_b_shift_equivalence = None
    f4_unit_b_shift_equivalence = None
    f2_unit_ab_shift_equivalence = None
    f4_unit_ab_shift_equivalence = None
    f2_unit_ab_lambda_shift_equivalence = None
    f4_unit_ab_lambda_shift_equivalence = None
    f2_unit_a_lambda_shift_equivalence = None
    f4_unit_a_lambda_shift_equivalence = None
    f2_unit_b_lambda_shift_equivalence = None
    f4_unit_b_lambda_shift_equivalence = None
    f2_unit_lambda_shift_equivalence = None
    f4_unit_lambda_shift_equivalence = None
    f2_polynomial_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction] = []
    f4_polynomial_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction] = []
    f2_reciprocal_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction] = []
    f4_reciprocal_prefactor_obstructions: list[Page43PolynomialSinglePrefactorObstruction] = []
    if looks_like_hero:
        b_cor, lam_cor = sp.symbols("b lambda")
        heine_cor2cf = heine_cor2cf_a_zero_contraction_obstruction(
            b=b_cor,
            lam=lam_cor,
            q=t,
            depth=12,
        )
        f2_equivalence = page43_f2_zero_shift_equivalence_obstruction(q=t)
        f4_equivalence = page43_f4_zero_shift_equivalence_obstruction(q=t)
        f2_unit_a_shift_equivalence = page43_f2_unit_a_shift_equivalence_obstruction(q=t)
        f4_unit_a_shift_equivalence = page43_f4_unit_a_shift_equivalence_obstruction(q=t)
        f2_unit_b_shift_equivalence = page43_f2_unit_b_shift_equivalence_obstruction(q=t)
        f4_unit_b_shift_equivalence = page43_f4_unit_b_shift_equivalence_obstruction(q=t)
        f2_unit_ab_shift_equivalence = page43_f2_unit_ab_shift_equivalence_obstruction(q=t)
        f4_unit_ab_shift_equivalence = page43_f4_unit_ab_shift_equivalence_obstruction(q=t)
        f2_unit_ab_lambda_shift_equivalence = page43_f2_unit_ab_lambda_shift_equivalence_obstruction(q=t)
        f4_unit_ab_lambda_shift_equivalence = page43_f4_unit_ab_lambda_shift_equivalence_obstruction(q=t)
        f2_unit_a_lambda_shift_equivalence = page43_f2_unit_a_lambda_shift_equivalence_obstruction(q=t)
        f4_unit_a_lambda_shift_equivalence = page43_f4_unit_a_lambda_shift_equivalence_obstruction(q=t)
        f2_unit_b_lambda_shift_equivalence = page43_f2_unit_b_lambda_shift_equivalence_obstruction(q=t)
        f4_unit_b_lambda_shift_equivalence = page43_f4_unit_b_lambda_shift_equivalence_obstruction(q=t)
        f2_unit_lambda_shift_equivalence = page43_f2_unit_lambda_shift_equivalence_obstruction(q=t)
        f4_unit_lambda_shift_equivalence = page43_f4_unit_lambda_shift_equivalence_obstruction(q=t)
        f2_polynomial_prefactor_obstructions = page43_zero_shift_polynomial_single_prefactor_obstructions(
            family="f2",
            q=t,
        )
        f4_polynomial_prefactor_obstructions = page43_zero_shift_polynomial_single_prefactor_obstructions(
            family="f4",
            q=t,
        )
        f2_reciprocal_prefactor_obstructions = page43_zero_shift_reciprocal_single_prefactor_obstructions(
            family="f2",
            q=t,
        )
        f4_reciprocal_prefactor_obstructions = page43_zero_shift_reciprocal_single_prefactor_obstructions(
            family="f4",
            q=t,
        )

    rr_subsequence_hits = arithmetic_subsequence_contraction_search(
        source_label="RR reciprocal",
        source_template=get_benchmark("rogers_ramanujan_normalized").canonical_template,
        target_template=reduced_candidate,
        q=t,
        max_stride=effective_max_stride,
        stages=profile.subsequence_stages,
    )
    cubic_subsequence_hits = arithmetic_subsequence_contraction_search(
        source_label="cubic reciprocal",
        source_template=get_benchmark("ramanujan_cubic_normalized").canonical_template,
        target_template=reduced_candidate,
        q=t,
        max_stride=effective_max_stride,
        stages=profile.subsequence_stages,
    )
    f2_hits = page43_monomial_parameter_search(
        family="f2",
        target_template=reduced_candidate,
        q=t,
        max_shift=profile.max_page43_shift,
        stages=profile.page43_stages,
    )
    f4_hits = page43_monomial_parameter_search(
        family="f4",
        target_template=reduced_candidate,
        q=t,
        max_shift=profile.max_page43_shift,
        stages=profile.page43_stages,
    )
    page43_rational_stages = profile.page43_stages
    f2_rational_hits = page43_rational_parameter_search(
        family="f2",
        target_template=reduced_candidate,
        q=t,
        max_shift=page43_rational_max_shift,
        stages=page43_rational_stages,
        max_nontrivial_profiles=profile.page43_rational_max_nontrivial_profiles,
    )
    f4_rational_hits = page43_rational_parameter_search(
        family="f4",
        target_template=reduced_candidate,
        q=t,
        max_shift=page43_rational_max_shift,
        stages=page43_rational_stages,
        max_nontrivial_profiles=profile.page43_rational_max_nontrivial_profiles,
    )

    return FormalizationContext(
        record=record,
        benchmark=benchmark,
        step=step,
        build_profile=profile.label,
        max_stride=effective_max_stride,
        subsequence_stages=profile.subsequence_stages,
        page43_max_shift=profile.max_page43_shift,
        page43_stages=profile.page43_stages,
        page43_rational_max_shift=page43_rational_max_shift,
        page43_rational_stages=page43_rational_stages,
        page43_rational_max_nontrivial_profiles=profile.page43_rational_max_nontrivial_profiles,
        reduced_candidate=reduced_candidate,
        reduced_benchmark=reduced_benchmark,
        target_b0=target_b0,
        target_a_terms=target_a_terms,
        target_b_terms=target_b_terms,
        factor_depth=factor_depth,
        factor_witness=factor_witness,
        rr_direct=rr_direct,
        cubic_direct=cubic_direct,
        cubic_odd=cubic_odd,
        cubic_even=cubic_even,
        heine_cor2cf=heine_cor2cf,
        f2_equivalence=f2_equivalence,
        f4_equivalence=f4_equivalence,
        f2_unit_a_shift_equivalence=f2_unit_a_shift_equivalence,
        f4_unit_a_shift_equivalence=f4_unit_a_shift_equivalence,
        f2_unit_b_shift_equivalence=f2_unit_b_shift_equivalence,
        f4_unit_b_shift_equivalence=f4_unit_b_shift_equivalence,
        f2_unit_ab_shift_equivalence=f2_unit_ab_shift_equivalence,
        f4_unit_ab_shift_equivalence=f4_unit_ab_shift_equivalence,
        f2_unit_ab_lambda_shift_equivalence=f2_unit_ab_lambda_shift_equivalence,
        f4_unit_ab_lambda_shift_equivalence=f4_unit_ab_lambda_shift_equivalence,
        f2_unit_a_lambda_shift_equivalence=f2_unit_a_lambda_shift_equivalence,
        f4_unit_a_lambda_shift_equivalence=f4_unit_a_lambda_shift_equivalence,
        f2_unit_b_lambda_shift_equivalence=f2_unit_b_lambda_shift_equivalence,
        f4_unit_b_lambda_shift_equivalence=f4_unit_b_lambda_shift_equivalence,
        f2_unit_lambda_shift_equivalence=f2_unit_lambda_shift_equivalence,
        f4_unit_lambda_shift_equivalence=f4_unit_lambda_shift_equivalence,
        f2_polynomial_prefactor_obstructions=f2_polynomial_prefactor_obstructions,
        f4_polynomial_prefactor_obstructions=f4_polynomial_prefactor_obstructions,
        f2_reciprocal_prefactor_obstructions=f2_reciprocal_prefactor_obstructions,
        f4_reciprocal_prefactor_obstructions=f4_reciprocal_prefactor_obstructions,
        rr_subsequence_hits=rr_subsequence_hits,
        cubic_subsequence_hits=cubic_subsequence_hits,
        f2_hits=f2_hits,
        f4_hits=f4_hits,
        f2_rational_hits=f2_rational_hits,
        f4_rational_hits=f4_rational_hits,
        looks_like_hero=looks_like_hero,
    )


def _write_formalization_note(context: FormalizationContext, output_path: str) -> None:
    lines: list[str] = [
        f"# Formalization Prep: `{context.record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{context.record.id}`",
        f"- Closest benchmark: `{context.record.closest_benchmark}` ({context.record.closest_benchmark_digits} shared digits)",
        f"- Candidate template: `{context.record.template.signature()}`",
        f"- Benchmark template: `{context.benchmark.canonical_template.signature()}`",
        f"- Build profile: `{context.build_profile}`",
        "",
        "## Current Theorem Status",
        "",
        "- No complete source theorem is identified yet.",
        "- This candidate is therefore **not ready** for a full Lean/Coq formalization of a final identity.",
        "- The award-track Lean scaffold now packages an exact intermediate waypoint: finite convergents agree with the reduced-by-factor model over the rational-function reverse-equivalence layer, the page-43 nearest-shift cube is ruled out as a source-family-specific exact lane family, the full zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active is also ruled out exactly in both page-43 families, the nearest RR/cubic arithmetic-subsequence source lanes are excluded exactly, and the direct RR/cubic plus Heine-`cor2cf` low-stage mismatch layers are packaged as named local obstructions.",
        "- The correct near-term target is to formalize exact local lemmas and keep bounded search evidence clearly separated from theorem-grade statements.",
    ]

    if context.reduced_candidate is None or context.reduced_benchmark is None:
        lines.extend(
            [
                "",
                "## Blocker",
                "",
                f"- The candidate does not cleanly reduce in the benchmark step `{context.step}`, so there is no obvious one-variable `t = q^{context.step}` theorem statement to prepare yet.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "## Exact Objects To Formalize",
                "",
                f"- Reduced variable: `t = q^{context.step}`",
                "- Target reciprocal object:",
                "",
                "```text",
                f"C(t) = {sp.Rational(context.reduced_candidate.base_denominator, context.reduced_candidate.top_constant)} + K_(n>=1) a_n / b_n",
                *_reciprocal_rule_lines(context.reduced_candidate, "t"),
                "```",
                "",
                "- Closest benchmark reciprocal object:",
                "",
                "```text",
                f"R(t) = {sp.Rational(context.reduced_benchmark.base_denominator, context.reduced_benchmark.top_constant)} + K_(n>=1) a_n / b_n",
                *_reciprocal_rule_lines(context.reduced_benchmark, "t"),
                "```",
                "",
                "## Exact Reduction And Equivalence Witness",
                "",
                f"- Exact convergent gcd factors were checked through stage `{context.factor_depth}`.",
                "- First common factors:",
                "",
                "```text",
                *[
                    f"g{n} = {_format_expr(context.factor_witness.reduction.gcd_factors[n])}"
                    for n in range(1, min(5, len(context.factor_witness.reduction.gcd_factors)))
                ],
                "```",
                "",
                "- After cancellation, the induced reduced-by-factor object begins:",
                "",
                "```text",
                f"b0_red = {_format_expr(context.factor_witness.reduction.reduced_coeffs.b0)}",
                *[
                    item
                    for n in range(
                        1,
                        min(5, len(context.factor_witness.reduction.reduced_coeffs.a_terms)),
                    )
                    for item in (
                        f"a{n}_red = {_format_expr(context.factor_witness.reduction.reduced_coeffs.a_terms[n])}",
                        f"b{n}_red = {_format_expr(context.factor_witness.reduction.reduced_coeffs.b_terms[n])}",
                    )
                ],
                "```",
                "",
                "- Reverse equivalence transform stage scales:",
                "",
                "```text",
                *[
                    f"r{n} = {_format_fraction_expr(context.factor_witness.scale_terms[n])}"
                    for n in range(1, min(5, len(context.factor_witness.scale_terms)))
                ],
                "```",
                "",
                "- These reverse scales are rational functions in `t`; the finite-stage reverse step now lives in Lean over `RatFunc Rat`, but a final identity still needs an infinite-object bridge and a fuller fraction-field coefficient layer beyond the current exact waypoint.",
                "",
                "## Exact Lemma Candidates",
                "",
                "### Direct 1-Step Bauer-Muir Obstructions",
                "",
                (
                    f"- RR source: `w0 = {_format_expr(context.rr_direct.forced_w0)}`, "
                    f"`w1 = {_format_expr(context.rr_direct.forced_w1)}`, so the first transformed numerator is "
                    f"`{_format_expr(context.rr_direct.transformed_a1)}` instead of `{_format_expr(context.rr_direct.target_a1)}`."
                ),
                (
                    f"- Cubic source: `w0 = {_format_expr(context.cubic_direct.forced_w0)}`, "
                    f"`w1 = {_format_expr(context.cubic_direct.forced_w1)}`, "
                    f"`w2 = {_format_expr(context.cubic_direct.forced_w2)}`, so the second transformed numerator is "
                    f"`{_format_expr(context.cubic_direct.transformed_a2)}` instead of `{_format_expr(context.cubic_direct.target_a2)}`."
                ),
                "",
                "### Simple Cubic Contraction Obstructions",
                "",
                (
                    f"- Odd contraction: initial term is `{_format_expr(context.cubic_odd.b0)}` instead of the target "
                    f"`{_format_expr(context.target_b0)}`."
                ),
                (
                    f"- Even contraction: first numerator does match as `{_format_expr(context.cubic_even.a_terms[1])}`, "
                    f"but the first denominator is `{_format_expr(context.cubic_even.b_terms[1])}` instead of "
                    f"`{_format_expr(context.target_b_terms[1])}`."
                ),
                "",
            ]
        )
        if context.heine_cor2cf is not None:
            lines.extend(
                [
                    "### Heine `cor2cf` Odd/Even Branch Obstructions",
                    "",
                    "- Lean mirror module: `proofs/Proofs/HeroCaseHeineCor2cf.lean`.",
                    (
                        f"- In the relevant `a = 0` lane, the odd part already has initial term "
                        f"`{_format_expr(context.heine_cor2cf.odd_part.b0)}` instead of "
                        f"`{_format_expr(context.target_b0)}`."
                    ),
                    (
                        f"- The even part keeps initial term `{_format_expr(context.heine_cor2cf.even_part.b0)}`, "
                        f"but its first numerator is `{_format_expr(context.heine_cor2cf.even_part.a_terms[1])}` "
                        f"instead of `{_format_expr(context.target_a_terms[1])}`."
                    ),
                    (
                        "- The odd-of-even branch changes the initial term to "
                        f"`{_format_fraction_expr(context.heine_cor2cf.even_odd_part.b0)}`, "
                        "so it fails before the first nontrivial numerator."
                    ),
                    (
                        "- The even-of-even branch keeps initial term `1`, but its first numerator is "
                        f"`{_format_expr(context.heine_cor2cf.even_even_part.a_terms[1])}`. "
                        "That numerator has no `t^2` term, so it cannot equal the target `t + t^2`."
                    ),
                    "",
                ]
            )
        if context.f2_equivalence is not None:
            lines.extend(
                [
                    "### Exact `f2` / `gcf3` `n`-Dependent Equivalence Lane",
                    "",
                    "- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.",
                    "- Prioritized source-family-specific lane: zero-shift `f2` / `gcf3` under arbitrary `n`-dependent equivalence factors.",
                    "- Write `m = t^(n-1)` and enforce the exact necessary identity",
                    "",
                    "```text",
                    "alpha_n * (1 + t^(n-1)) = t^n * beta_(n-1) * beta_n",
                    "```",
                    "",
                    "- In the zero-shift `f2` / `gcf3` lane this becomes",
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
                    _format_expr(context.f2_equivalence.residual_polynomial),
                    "```",
                    "",
                    (
                        "- `m^3` coefficient is "
                        f"`{_format_expr(context.f2_equivalence.m_coefficients[3])}`; "
                        "exact vanishing forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient is "
                        f"`{_format_expr(context.f2_equivalence.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes "
                        f"`{_format_expr(context.f2_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    "- So no arbitrary `n`-dependent equivalence transformation sends the hero case into this zero-shift `f2` / `gcf3` lane.",
                    "",
                ]
            )
        if context.f4_equivalence is not None:
            lines.extend(
                [
                    "### Exact `f4` / `gcf2` `n`-Dependent Equivalence Lane",
                    "",
                    "- The same Lean mirror module also now covers the zero-shift `f4` / `gcf2` lane.",
                    "- In that lane the necessary identity becomes",
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
                    _format_expr(context.f4_equivalence.residual_polynomial),
                    "```",
                    "",
                    (
                        "- `m^0` coefficient is "
                        f"`{_format_expr(context.f4_equivalence.m_coefficients[0])}`; "
                        "exact vanishing forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient is "
                        f"`{_format_expr(sp.simplify(context.f4_equivalence.m_coefficients[3].subs(context.f4_equivalence.forced_a_solution)))}`; "
                        "exact vanishing forces `b = 0`."
                    ),
                    (
                        "- Then the `m^1` coefficient is "
                        f"`{_format_expr(context.f4_equivalence.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- With `a = b = 0`, `lambda = 1`, the `m^2` coefficient becomes "
                        f"`{_format_expr(context.f4_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    "- So no arbitrary `n`-dependent equivalence transformation sends the hero case into this zero-shift `f4` / `gcf2` lane either.",
                    "",
                ]
            )
        if (
            context.f2_unit_lambda_shift_equivalence is not None
            and context.f4_unit_lambda_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Unit-Shift `lambda` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the next nearby shift choice `lambda -> lambda*t` with zero `a`/`b` shifts.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient is unchanged from the zero-shift lane, so exact vanishing still forces "
                        "`a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_lambda_shift_equivalence.impossible_m1_coefficient)}`; "
                        "no constant `lambda` can make that polynomial vanish identically."
                    ),
                    (
                        "- For `f4/gcf2`, exact vanishing again forces "
                        "`a = 0`, then `b = 0`, and the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f4_unit_lambda_shift_equivalence.impossible_m1_coefficient)}`."
                    ),
                    "- So the first nonzero `lambda`-shift nearest lanes already fail before any final `m^2` obstruction is needed.",
                    "",
                ]
            )
        if (
            context.f2_unit_a_lambda_shift_equivalence is not None
            and context.f4_unit_a_lambda_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Mixed Unit-Shift `a`/`lambda` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the first mixed nearby shift choice `a -> a*t`, `lambda -> lambda*t` with zero `b` shift.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_a_lambda_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_a_lambda_shift_equivalence.impossible_m1_coefficient)}`; "
                        "no constant `lambda` can make that polynomial vanish identically."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_a_lambda_shift_equivalence.m_coefficients[0])}` forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_a_lambda_shift_equivalence.m_coefficients[3].subs(context.f4_unit_a_lambda_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f4_unit_a_lambda_shift_equivalence.impossible_m1_coefficient)}`."
                    ),
                    "- So the first mixed unit-`a` / unit-`lambda` lanes already fail before any final `m^2` obstruction is needed.",
                    "",
                ]
            )
        if (
            context.f2_unit_b_lambda_shift_equivalence is not None
            and context.f4_unit_b_lambda_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Mixed Unit-Shift `b`/`lambda` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the first mixed nearby shift choice `b -> b*t`, `lambda -> lambda*t` with zero `a` shift.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_b_lambda_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_b_lambda_shift_equivalence.impossible_m1_coefficient)}`; "
                        "no constant `lambda` can make that polynomial vanish identically."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_b_lambda_shift_equivalence.m_coefficients[0])}` still forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_b_lambda_shift_equivalence.m_coefficients[3].subs(context.f4_unit_b_lambda_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f4_unit_b_lambda_shift_equivalence.impossible_m1_coefficient)}`."
                    ),
                    "- So the first mixed unit-`b` / unit-`lambda` lanes also fail before any final `m^2` obstruction is needed.",
                    "",
                ]
            )
        if (
            context.f2_unit_ab_lambda_shift_equivalence is not None
            and context.f4_unit_ab_lambda_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Mixed Unit-Shift `a`/`b`/`lambda` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the first full nearby shift choice `a -> a*t`, `b -> b*t`, `lambda -> lambda*t`.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_ab_lambda_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_ab_lambda_shift_equivalence.impossible_m1_coefficient)}`; "
                        "no constant `lambda` can make that polynomial vanish identically."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_ab_lambda_shift_equivalence.m_coefficients[0])}` forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_ab_lambda_shift_equivalence.m_coefficients[3].subs(context.f4_unit_ab_lambda_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and the surviving `m^1` coefficient becomes "
                        f"`{_format_expr(context.f4_unit_ab_lambda_shift_equivalence.impossible_m1_coefficient)}`."
                    ),
                    "- So the first full three-parameter nearest lane also fails before any final `m^2` obstruction is needed.",
                    "- Together with the zero-shift and the other seven nearest-shift cases, this closes the full `{0,1}^3` nearest-shift cube at the current page-43 exact-equivalence level.",
                    "- Lean now exposes that cube not only as summary theorems, but also as Bool-parameterized theorems over the shift bits.",
                    "",
                ]
            )
        if (
            context.f2_unit_a_shift_equivalence is not None
            and context.f4_unit_a_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Unit-Shift `a` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the next nearby shift choice `a -> a*t` with zero `b`/`lambda` shifts.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_a_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_a_shift_equivalence.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f2_unit_a_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_a_shift_equivalence.m_coefficients[0])}` already forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_a_shift_equivalence.m_coefficients[3].subs(context.f4_unit_a_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and then the `m^1` coefficient "
                        f"`{_format_expr(context.f4_unit_a_shift_equivalence.reduced_m1_coefficient)}` forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f4_unit_a_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    "- So the first nonzero `a`-shift nearest lanes also fail by an exact final obstruction.",
                    "",
                ]
            )
        if (
            context.f2_unit_b_shift_equivalence is not None
            and context.f4_unit_b_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Unit-Shift `b` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the next nearby shift choice `b -> b*t` with zero `a`/`lambda` shifts.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_b_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_b_shift_equivalence.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f2_unit_b_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_b_shift_equivalence.m_coefficients[0])}` still forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_b_shift_equivalence.m_coefficients[3].subs(context.f4_unit_b_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and then the `m^1` coefficient "
                        f"`{_format_expr(context.f4_unit_b_shift_equivalence.reduced_m1_coefficient)}` forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f4_unit_b_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    "- So the first nonzero `b`-shift nearest lanes also fail by an exact final obstruction.",
                    "",
                ]
            )
        if (
            context.f2_unit_ab_shift_equivalence is not None
            and context.f4_unit_ab_shift_equivalence is not None
        ):
            lines.extend(
                [
                    "### Exact Mixed Unit-Shift `a`/`b` Page-43 Lanes",
                    "",
                    "- The same exact-equivalence layer now also covers the first mixed nearby shift choice `a -> a*t`, `b -> b*t` with zero `lambda` shift.",
                    (
                        "- For `f2/gcf3`, the `m^3` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_ab_shift_equivalence.m_coefficients[3])}`, "
                        "and exact vanishing still forces `a = 0`, `b = 0`."
                    ),
                    (
                        "- After that specialization, the `m^1` coefficient becomes "
                        f"`{_format_expr(context.f2_unit_ab_shift_equivalence.reduced_m1_coefficient)}`; "
                        "exact vanishing forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f2_unit_ab_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    (
                        "- For `f4/gcf2`, the constant coefficient "
                        f"`{_format_expr(context.f4_unit_ab_shift_equivalence.m_coefficients[0])}` already forces `a = 0`."
                    ),
                    (
                        "- After that, the `m^3` coefficient "
                        f"`{_format_expr(sp.simplify(context.f4_unit_ab_shift_equivalence.m_coefficients[3].subs(context.f4_unit_ab_shift_equivalence.forced_a_solution)))}` "
                        "forces `b = 0`, and then the `m^1` coefficient "
                        f"`{_format_expr(context.f4_unit_ab_shift_equivalence.reduced_m1_coefficient)}` forces `lambda = 1`."
                    ),
                    (
                        "- The surviving `m^2` coefficient is then "
                        f"`{_format_expr(context.f4_unit_ab_shift_equivalence.final_m2_coefficient)}`, still nonzero."
                    ),
                    "- So the first mixed unit-`a` / unit-`b` lane also fails by an exact final obstruction.",
                    "",
                ]
            )
        if (
            context.f2_polynomial_prefactor_obstructions
            and context.f4_polynomial_prefactor_obstructions
        ):
            a_sym, b_sym, lam_sym = sp.symbols("a b lambda")
            lines.extend(
                [
                    "### Exact Zero-Shift Polynomial Single-Prefactor Page-43 Lanes",
                    "",
                    "- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.",
                    "- This exact layer now covers the whole polynomial sub-box `phi in {1, 1+t}` with zero shifts and at most one non-plain prefactor active.",
                ]
            )
            for obstruction in context.f2_polynomial_prefactor_obstructions:
                label = _format_prefactor_case_label(obstruction)
                if obstruction.failure_stage == "denominator":
                    lines.append(f"- `f2/gcf3`, {label}: denominator matching is already incompatible.")
                elif obstruction.failure_stage == "stage1_numerator":
                    lines.append(
                        f"- `f2/gcf3`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "but the remaining stage-1 numerator difference is "
                        f"`{_format_expr(obstruction.reduced_a1_difference)}`."
                    )
                else:
                    lines.append(
                        f"- `f2/gcf3`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "stage-1 numerator matching then forces "
                        f"`{_format_symbolic_solution(obstruction.forced_lambda_solution, symbols=(lam_sym,))}`, "
                        "but the stage-2 numerator becomes "
                        f"`{_format_expr(obstruction.final_stage2_a)}` instead of "
                        f"`{_format_expr(obstruction.target_stage2_a)}`."
                    )
            for obstruction in context.f4_polynomial_prefactor_obstructions:
                label = _format_prefactor_case_label(obstruction)
                if obstruction.failure_stage == "denominator":
                    lines.append(f"- `f4/gcf2`, {label}: denominator matching is already incompatible.")
                elif obstruction.failure_stage == "stage1_numerator":
                    lines.append(
                        f"- `f4/gcf2`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "but the remaining stage-1 numerator difference is "
                        f"`{_format_expr(obstruction.reduced_a1_difference)}`."
                    )
                else:
                    lines.append(
                        f"- `f4/gcf2`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "stage-1 numerator matching then forces "
                        f"`{_format_symbolic_solution(obstruction.forced_lambda_solution, symbols=(lam_sym,))}`, "
                        "but the stage-2 numerator becomes "
                        f"`{_format_expr(obstruction.final_stage2_a)}` instead of "
                        f"`{_format_expr(obstruction.target_stage2_a)}`."
                    )
            lines.extend(
                [
                    "- So the whole zero-shift polynomial single-prefactor sub-box is excluded exactly in both page-43 families.",
                    "",
                ]
            )
        if (
            context.f2_reciprocal_prefactor_obstructions
            and context.f4_reciprocal_prefactor_obstructions
        ):
            a_sym, b_sym = sp.symbols("a b")
            lines.extend(
                [
                    "### Exact Zero-Shift Reciprocal Single-Prefactor Page-43 Lanes",
                    "",
                    "- Lean mirror module: `proofs/Proofs/HeroCasePage43Equivalence.lean`.",
                    "- This exact layer now also covers the reciprocal sub-box `phi in {1/(1+t)}` via cross-multiplied coefficient identities.",
                ]
            )
            for obstruction in context.f2_reciprocal_prefactor_obstructions:
                label = _format_prefactor_case_label(obstruction)
                if obstruction.failure_stage == "denominator":
                    lines.append(f"- `f2/gcf3`, {label}: cross-multiplied denominator matching is already incompatible.")
                else:
                    lines.append(
                        f"- `f2/gcf3`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "but the remaining cross-multiplied stage-1 numerator difference has numerator "
                        f"`{_format_expr(sp.together(obstruction.reduced_a1_difference).as_numer_denom()[0])}`."
                    )
            for obstruction in context.f4_reciprocal_prefactor_obstructions:
                label = _format_prefactor_case_label(obstruction)
                if obstruction.failure_stage == "denominator":
                    lines.append(f"- `f4/gcf2`, {label}: cross-multiplied denominator matching is already incompatible.")
                else:
                    lines.append(
                        f"- `f4/gcf2`, {label}: denominator matching forces "
                        f"`{_format_symbolic_solution(obstruction.forced_ab_solution, symbols=(a_sym, b_sym))}`, "
                        "but the remaining cross-multiplied stage-1 numerator difference has numerator "
                        f"`{_format_expr(sp.together(obstruction.reduced_a1_difference).as_numer_denom()[0])}`."
                    )
            lines.extend(
                [
                    "- So the whole zero-shift reciprocal single-prefactor sub-box is excluded exactly in both page-43 families.",
                    "- Taken together, the current zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` is now fully exactified at theorem grade.",
                    "",
                ]
            )
        lines.extend(
            [
                "## Bounded Exact Exclusion Results",
                "",
                (
                    f"- Arithmetic subsequence contractions up to stride `{context.max_stride}` "
                    f"with stage comparison depth `{context.subsequence_stages}`: "
                    f"RR hits `{len(context.rr_subsequence_hits)}`, cubic hits `{len(context.cubic_subsequence_hits)}`."
                ),
                "- These are exact statements for the bounded class being checked, but they do not identify a final source theorem.",
                (
                    f"- Page-43 monomial substitutions in the current "
                    f"`[-{context.page43_max_shift},{context.page43_max_shift}]` shift box "
                    f"with `{context.page43_stages}` matched stages: "
                    f"`f2/gcf3` hits `{len(context.f2_hits)}`, `f4/gcf2` hits `{len(context.f4_hits)}`."
                ),
                (
                    "- Page-43 low-complexity rational-prefactor box: "
                    f"`phi in {{1, 1+t, 1/(1+t)}}` with at most "
                    f"`{context.page43_rational_max_nontrivial_profiles}` non-plain prefactor"
                    f"{'' if context.page43_rational_max_nontrivial_profiles == 1 else 's'} active, "
                    f"shift box `[-{context.page43_rational_max_shift},{context.page43_rational_max_shift}]`, "
                    f"and `{context.page43_rational_stages}` matched stages: "
                    f"`f2/gcf3` hits `{len(context.f2_rational_hits)}`, "
                    f"`f4/gcf2` hits `{len(context.f4_rational_hits)}`."
                ),
                "- These are bounded symbolic searches, useful for narrowing the theorem statement but not substitutes for a full origin proof.",
            ]
        )
        if (
            context.f2_polynomial_prefactor_obstructions
            and context.f4_polynomial_prefactor_obstructions
            and context.f2_reciprocal_prefactor_obstructions
            and context.f4_reciprocal_prefactor_obstructions
        ):
            lines.append(
                "- Within that bounded prefactor box, the full zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active is now theorem-grade exact."
            )
        if (
            context.rr_subsequence_hits
            or context.cubic_subsequence_hits
            or context.f2_hits
            or context.f4_hits
            or context.f2_rational_hits
            or context.f4_rational_hits
        ):
            lines.extend(["", "```text"])
            for hit in context.rr_subsequence_hits + context.cubic_subsequence_hits:
                lines.append(f"{hit.source_label}: stride={hit.stride}, offset={hit.offset}")
            for hit in context.f2_hits + context.f4_hits:
                lines.append(
                    f"{hit.family}: A={hit.a_shift}, B={hit.b_shift}, L={hit.lambda_shift}, "
                    f"alpha={_format_expr(hit.a_coeff)}, beta={_format_expr(hit.b_coeff)}, "
                    f"lambda={_format_expr(hit.lambda_coeff)}"
                )
            for hit in context.f2_rational_hits + context.f4_rational_hits:
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
                "## Formalization Order",
                "",
                "1. Formalize generalized continued fractions and convergent recurrence for finite truncations.",
                "2. Reuse the exact convergent-factor reduction theorem for the candidate-side local model.",
                "3. Bridge the current rational-function reverse-equivalence layer from finite convergents to the infinite object.",
                "4. Extend `Proofs/HeroCasePage43Equivalence.lean` beyond the currently formalized nearest-shift cube plus the exact zero-shift single-prefactor box, especially if shifted or multi-prefactor rational lanes become theorem-relevant.",
                "5. Formalize the direct 1-step Bauer-Muir obstruction lemmas against the reduced target.",
                "6. Formalize odd/even contraction reconstruction together with the cubic and Heine-`cor2cf` low-stage mismatch lemmas.",
                "7. Defer the bounded search exclusions until a final theorem statement makes them clearly necessary.",
                "8. Do not start a full Lean/Coq origin theorem until a unique source family or exact identity is identified.",
            ]
        )

        if context.looks_like_hero:
            lines.extend(
                [
                    "",
                    "## Why This Is Still Not A Full Theorem",
                    "",
                    "- The current exact lemmas only rule out nearby transforms and simple contraction sources.",
                    "- They do not prove what the candidate *is*.",
                    "- A full formal proof needs a final theorem statement of the form `C(t) = known_object(t)` or a uniquely characterizing theorem that has not been found yet.",
                ]
            )
            award_target = "proofs/Proofs/HeroCaseFinalIdentity.lean"
            lines.extend(
                [
                    "",
                    "## Award-Track Endgame Hook",
                    "",
                    "- This candidate matches the current hero-case structural signature in reduced variable `t`.",
                    f"- Award-track target module (Lean scaffold): `{award_target}`",
                    "- Current state: the module compiles, carries status marker `exclusion_waypoint`, and now exposes the current exact waypoint as a certificate object `currentExactWaypointCertificate`, whose two fields are proved by `finiteConvergentReductionWaypoint_true` and `knownSourceOrbitExclusionWaypoint_true` and then repackaged by `exactWaypointStatement_true`. The lower-level exact ingredients now also explicitly include `mortonSquaredCoordinateExcluded`, `weberSchlafliCoordinateExcluded`, and the older `directLocalObstructions` / `simpleCor2cfBranchesExcluded` layer, while the research-side hand-off now threads through `mortonNamedCoordinateResearchWaypoint_true`, `weberSchlafliBridgeResearchWaypoint_true`, `namedWeberCoordinateBridgeResearchWaypoint_true`, `currentNamedWeberCoordinateBridgeResearchCertificate`, `ggWeightedOrbitResearchWaypoint_true`, `currentGGWeightedOrbitResearchCertificate`, `currentNamedWeberOrbitResearchCertificate`, and `namedWeberOrbitResearchWaypoint_true`; the frontier is now compressed into two dedicated research-side blocks, one coordinate-bridge certificate for `Morton named stack -> P_ws -> Weber bridge` and one weighted-`GG` orbit certificate for `GGQ34 -> GG weighted correction`, before the combined exact-plus-research frontier is exposed as `currentRecognitionFrontierCertificate` / `currentRecognitionFrontierWaypoint_true`.",
                    "- Current Weber-Schlafli exact shell: `proofs/Proofs/HeroCaseWeberSchlafliCoordinateObstruction.lean` packages `P_ws = (1/F - F) / 2` and the first Weber-Schlafli template `P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0` as a seven-class repeated witness table `(2,-1)`, `(2,3)`, `(2,8)`, `(4,3)`, `(6,3)`, `(8,3)`, `(10,3)` rather than a universal one-class obstruction; that table is now also compressed one step further into coefficient support `{-1, 3, 8}`, the coefficient-`3` even ladder `2,4,6,8,10`, and the two exceptional power-`2` classes for coefficients `-1` and `8`.",
                    "- Current Morton named-coordinate waypoint: `proofs/Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean` packages the unified source-faithful stack `X_mt = F^2`, `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`, `P_ws = (1/F - F) / 2`, and `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`, with universal obstruction classes already visible in the `X_mt`, `T_mt`, and `B_ws` lanes and the standalone exact `P_ws` shell reused inside the larger stack.",
                    "- Current Morton square-coordinate exact shell: `proofs/Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean` still packages the explicit source-faithful lane `X_mt = F^2`, the Proposition 3.2 template `X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 = 0`, and the current universal constant-term obstruction `(t^0, 4)` across all sampled tail-family objects.",
                    "- Current modular-coordinate waypoint: `proofs/Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean` records the first failure of `F / W_34`, the normalized probe `G_W34`, the deeper follow-up `G2_W34`, and the current empty small correction boxes at both normalized layers.",
                    "- Current Weber residual waypoint: `proofs/Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean` records the primary residual choice `G_g12_ws`, the exact coordinate bridge tying `G_p12_ws` back to it, the classical Weber `f2` tri-product coordinate `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2) = G_g12_ws*G_p12_ws` together with its normalized follow-up `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)`, and now also records the source-backed reading that Berndt--Chan--Zhang identify Ramanujan-Weber `G_n` / `g_n` with classical Weber `f` / `f1` while Yui--Zagier supplies the classical Weber `f`, `f1`, `f2` trio, so this `g12_ws` / `p12_ws` / `G_f2_ws` shell should be read as that named Weber trio in the project's normalization. The same waypoint then keeps the derived quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2`, the exact quotient-coordinate bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0` with `Q_gp_ws = p12_ws / g12_ws`, the template-normalized coordinate `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`, its normalized follow-up `H_X_ws = (G_X_ws - 1) / (4*t^1)`, the direct follow-up bridge `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws`, the quotient-bridge follow-up `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, the next stripped comparison back to the template-normalized lane `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws`, its quotient follow-up `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)`, the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws`, the normalized follow-up `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`, the theorem-shaped return-bridge closure shell saying `Q_XK_ws` and `L_XK_ws` still have `0` hits in the checked direct / quotient / mixed named `GG` modular-equation boxes, and now also packages the exact Chan--Huang obstruction quartets on the later return bridge: for `Q_XK_ws`, direct failures `(-9/2, -6)` against `GG3/GG4` and quotient failures `(-3, -9/2)` against `Q_3/Q_4`; for `L_XK_ws`, direct failures `(593/10, 1186/15)` and quotient failures `(593/15, 593/10)`, alongside the current empty theorem-shaped uniqueness / closure boxes on the hero-side follow-up objects.",
                    "- `finalIdentityStatement` is still a placeholder because no final closed form is identified yet.",
                    "- Only replace the placeholder after a concrete closed form is identified.",
                ]
            )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _lean_namespace_component(raw: str) -> str:
    pieces = re.findall(r"[A-Za-z0-9]+", raw)
    if not pieces:
        return "GeneratedCandidate"
    normalized = "".join(piece[:1].upper() + piece[1:] for piece in pieces)
    if normalized[0].isdigit():
        normalized = f"Candidate{normalized}"
    return normalized


def _lean_rat(value: sp.Expr | int) -> str:
    rational = sp.Rational(value)
    if rational.q == 1:
        return f"({rational.p} : Rat)"
    return f"(({rational.p} : Rat) / ({rational.q} : Rat))"


def _lean_nat_expr(shift: int, step: int, index_var: str = "k") -> str:
    if step == 0:
        return str(shift)
    if shift == 0:
        return f"{step} * {index_var}"
    return f"{shift} + {step} * {index_var}"


def _lean_poly_term(scale: int, shift: int, step: int, index_var: str = "k") -> str:
    exponent = _lean_nat_expr(shift, step, index_var=index_var)
    return f"C {_lean_rat(scale)} * X ^ ({exponent})"


def _lean_poly_term_at_stage(scale: int, exponent: int) -> str:
    if exponent == 0:
        return f"C {_lean_rat(scale)}"
    base = "X" if exponent == 1 else f"X ^ {exponent}"
    if scale == 1:
        return base
    if scale == -1:
        return f"-{base}"
    return f"C {_lean_rat(scale)} * {base}"


def _lean_template_a_expr(template: QCFTemplate, index_var: str = "k") -> str:
    terms = [_lean_poly_term(template.numerator_scale, template.numerator_q_shift, template.numerator_q_step, index_var=index_var)]
    if template.numerator_extra_scale != 0:
        terms.append(
            _lean_poly_term(
                template.numerator_extra_scale,
                template.numerator_extra_q_shift,
                template.numerator_extra_q_step,
                index_var=index_var,
            )
        )
    return " + ".join(terms)


def _lean_template_b_expr(template: QCFTemplate, index_var: str = "k") -> str:
    terms = [f"C {_lean_rat(template.denominator_constant)}"]
    if template.denominator_scale != 0:
        terms.append(
            _lean_poly_term(
                template.denominator_scale,
                template.denominator_q_shift,
                template.denominator_q_step,
                index_var=index_var,
            )
        )
    return " + ".join(terms)


def _lean_template_b_var(template: QCFTemplate) -> str:
    if template.denominator_scale == 0 or template.denominator_q_step == 0:
        return "_"
    return "k"


def _lean_template_a_stage_expr(template: QCFTemplate, stage: int) -> str:
    exponent = template.numerator_q_shift + template.numerator_q_step * stage
    terms = [_lean_poly_term_at_stage(template.numerator_scale, exponent)]
    if template.numerator_extra_scale != 0:
        extra_exponent = template.numerator_extra_q_shift + template.numerator_extra_q_step * stage
        terms.append(_lean_poly_term_at_stage(template.numerator_extra_scale, extra_exponent))
    return " + ".join(terms)


def _lean_template_b_stage_expr(template: QCFTemplate, stage: int) -> str:
    terms = [f"C {_lean_rat(template.denominator_constant)}"]
    if template.denominator_scale != 0:
        exponent = template.denominator_q_shift + template.denominator_q_step * stage
        terms.append(_lean_poly_term_at_stage(template.denominator_scale, exponent))
    return " + ".join(terms)


def _sympy_poly_to_lean(expr, variable: str = "X") -> str:
    poly = sp.Poly(sp.expand(expr), sp.Symbol("t"), domain="QQ")
    pieces: list[str] = []
    for exponent_tuple, coeff in poly.terms():
        exponent = exponent_tuple[0]
        coeff = sp.Rational(coeff)
        integer_coeff = coeff.q == 1
        if exponent == 0:
            term = str(coeff.p) if integer_coeff else f"C {_lean_rat(coeff)}"
        elif coeff == 1:
            term = variable if exponent == 1 else f"{variable} ^ {exponent}"
        elif coeff == -1:
            base = variable if exponent == 1 else f"{variable} ^ {exponent}"
            term = f"-{base}"
        else:
            base = variable if exponent == 1 else f"{variable} ^ {exponent}"
            if integer_coeff:
                term = f"{coeff.p} * {base}"
            else:
                term = f"C {_lean_rat(coeff)} * {base}"
        pieces.append(f"({term})")
    if not pieces:
        return "(0 : QPoly)"
    return " + ".join(pieces)


def _write_lean_skeleton(context: FormalizationContext, output_path: str) -> None:
    namespace_name = _lean_namespace_component(context.record.id)
    candidate_convergents: list[tuple[sp.Expr, sp.Expr]] = []
    if context.target_b0 is not None and context.target_a_terms and context.target_b_terms:
        candidate_convergents = continued_fraction_convergents(
            b0=context.target_b0,
            a_terms=context.target_a_terms[:3],
            b_terms=context.target_b_terms[:3],
        )
    lines: list[str] = [
        "import Mathlib",
        "import Proofs.GeneralizedCF",
        "",
        "open Polynomial",
        "",
        "namespace Proofs",
        "namespace Generated",
        f"namespace {namespace_name}",
        "",
        "/-!",
        "Auto-generated by `python -m ramanujan_discovery formalize`.",
        f"Candidate id: {context.record.id}",
        f"Closest benchmark: {context.record.closest_benchmark}",
        f"Reduced variable: t = q^{context.step}",
        f"Candidate signature: {context.record.template.signature()}",
        f"Benchmark signature: {context.benchmark.canonical_template.signature()}",
        "-/",
        "",
        "abbrev QPoly := Polynomial Rat",
        "",
        "noncomputable section",
        "",
    ]

    if context.reduced_candidate is None or context.reduced_benchmark is None:
        lines.extend(
            [
                "/-!",
                "Blocker:",
                f"The candidate does not cleanly reduce in step q^{context.step}, so there is",
                "no one-variable theorem skeleton to generate yet.",
                "-/",
                "",
                "end",
                "",
                f"end {namespace_name}",
                "end Generated",
                "end Proofs",
            ]
        )
    else:
        lines.extend(
            [
                f"def candidateData : Proofs.GCFData QPoly := {{",
                f"  b0 := C {_lean_rat(sp.Rational(context.reduced_candidate.base_denominator, context.reduced_candidate.top_constant))}",
                f"  a := fun k => {_lean_template_a_expr(context.reduced_candidate)}",
                f"  b := fun {_lean_template_b_var(context.reduced_candidate)} => {_lean_template_b_expr(context.reduced_candidate)}",
                "}",
                "",
                f"def benchmarkData : Proofs.GCFData QPoly := {{",
                f"  b0 := C {_lean_rat(sp.Rational(context.reduced_benchmark.base_denominator, context.reduced_benchmark.top_constant))}",
                f"  a := fun k => {_lean_template_a_expr(context.reduced_benchmark)}",
                f"  b := fun {_lean_template_b_var(context.reduced_benchmark)} => {_lean_template_b_expr(context.reduced_benchmark)}",
                "}",
                "",
                f"def rrForcedStage1Numerator : QPoly := {_sympy_poly_to_lean(context.rr_direct.transformed_a1)}",
                f"def cubicForcedStage2Numerator : QPoly := {_sympy_poly_to_lean(context.cubic_direct.transformed_a2)}",
                f"def cubicOddInitialDenominator : QPoly := {_sympy_poly_to_lean(context.cubic_odd.b0)}",
                f"def cubicEvenFirstDenominator : QPoly := {_sympy_poly_to_lean(context.cubic_even.b_terms[1])}",
                "",
                "theorem candidate_stage0_a :",
                f"    candidateData.a 0 = {_lean_template_a_stage_expr(context.reduced_candidate, 0)} := by",
                "  simp [candidateData]",
                "",
                "theorem candidate_stage1_a :",
                f"    candidateData.a 1 = {_lean_template_a_stage_expr(context.reduced_candidate, 1)} := by",
                "  simp [candidateData]",
                "",
                "theorem candidate_stage0_b :",
                f"    candidateData.b 0 = {_lean_template_b_stage_expr(context.reduced_candidate, 0)} := by",
                "  simp [candidateData]",
                "",
                "theorem candidate_first_convergent_num :",
                f"    Proofs.continuantNum candidateData 1 = {_sympy_poly_to_lean(candidate_convergents[1][0])} := by",
                "  simp [Proofs.continuantNum, candidateData]",
                "  ring_nf",
                "",
                "theorem candidate_first_convergent_den :",
                f"    Proofs.continuantDen candidateData 1 = {_sympy_poly_to_lean(candidate_convergents[1][1])} := by",
                "  simp [Proofs.continuantDen, candidateData]",
                "  ring_nf",
                "",
                "theorem candidate_second_convergent_num :",
                f"    Proofs.continuantNum candidateData 2 = {_sympy_poly_to_lean(candidate_convergents[2][0])} := by",
                "  simp [Proofs.continuantNum, candidateData]",
                "  ring_nf",
                "",
                "theorem candidate_second_convergent_den :",
                f"    Proofs.continuantDen candidateData 2 = {_sympy_poly_to_lean(candidate_convergents[2][1])} := by",
                "  simp [Proofs.continuantDen, candidateData]",
                "  ring_nf",
                "",
                "/-!",
                "Theorems below are generated from exact symbolic witnesses recovered by the",
                "Python formalization pipeline.",
                "",
                f"RR direct witness: w0 = {_format_expr(context.rr_direct.forced_w0)}, w1 = {_format_expr(context.rr_direct.forced_w1)}, transformed a1 = {_format_expr(context.rr_direct.transformed_a1)}, target a1 = {_format_expr(context.rr_direct.target_a1)}",
                f"Cubic direct witness: w0 = {_format_expr(context.cubic_direct.forced_w0)}, w1 = {_format_expr(context.cubic_direct.forced_w1)}, w2 = {_format_expr(context.cubic_direct.forced_w2)}, transformed a2 = {_format_expr(context.cubic_direct.transformed_a2)}, target a2 = {_format_expr(context.cubic_direct.target_a2)}",
                f"First reverse-equivalence scales: {', '.join(f'r{n} = {_format_fraction_expr(context.factor_witness.scale_terms[n])}' for n in range(1, min(5, len(context.factor_witness.scale_terms))))}",
                "These scales are rational functions, and Proofs/RationalEquivalence.lean already formalizes the finite-stage reverse equivalence transform over RatFunc Rat; the remaining gap is a fuller fraction-field coefficient layer for the infinite-object bridge.",
                "Proofs/HeroCaseFinalIdentity.lean packages that exact waypoint via the certificate object `currentExactWaypointCertificate`, whose two fields are established by `finiteConvergentReductionWaypoint_true` and `knownSourceOrbitExclusionWaypoint_true` and then bundled by `exactWaypointStatement_true`; the lower-level ingredients still include `heroConvergentRatFunc_eq_reducedHeroConvergentRatFunc`, `reverseEquivalenceRecoversHeroData`, `page43PolynomialPrefactorExcluded`, `page43ReciprocalPrefactorExcluded`, `nearestArithmeticSubsequenceSourcesExcluded`, `mortonSquaredCoordinateExcluded`, `weberSchlafliCoordinateExcluded`, `directLocalObstructions`, and `simpleCor2cfBranchesExcluded`. The same module also now threads the research-side hand-off through `mortonNamedCoordinateResearchWaypoint_true`, `weberSchlafliBridgeResearchWaypoint_true`, `namedWeberCoordinateBridgeResearchWaypoint_true`, `currentNamedWeberCoordinateBridgeResearchCertificate`, `ggWeightedOrbitResearchWaypoint_true`, `currentGGWeightedOrbitResearchCertificate`, `currentNamedWeberOrbitResearchCertificate`, and `namedWeberOrbitResearchWaypoint_true`; the frontier is now compressed into two dedicated research-side blocks, one coordinate-bridge certificate for `Morton named stack -> P_ws -> Weber bridge` and one weighted-`GG` orbit certificate for `GGQ34 -> GG weighted correction`, and the module exposes the combined exact-plus-research frontier as `currentRecognitionFrontierCertificate` / `currentRecognitionFrontierWaypoint_true`.",
                "The same proof workspace also now records the standalone Weber-Schlafli exact shell in `Proofs/HeroCaseWeberSchlafliCoordinateObstruction.lean`: `P_ws = (1/F - F) / 2`, the first template `P_ws^2 * P_ws,2^2 + P_ws^2 - 2*P_ws,2 = 0`, the current seven-class repeated witness table `(2,-1)`, `(2,3)`, `(2,8)`, `(4,3)`, `(6,3)`, `(8,3)`, `(10,3)`, and the compressed profile saying the coefficient support is exactly `{-1, 3, 8}`, the coefficient `3` sits on the even ladder `2,4,6,8,10`, and the exceptional coefficients `-1` and `8` only occur at power `2`.",
                "The same proof workspace also now records the unified Morton named-coordinate shell in `Proofs/HeroCaseMortonNamedCoordinateWaypoint.lean`: `X_mt = F^2`, `T_mt = (X_mt - sigma^2) / (sigma^2*X_mt - 1)`, `P_ws = (1/F - F) / 2`, and `B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)`, with universal obstruction classes already visible in the `X_mt`, `T_mt`, and `B_ws` lanes and the standalone exact `P_ws` shell reused inside the larger stack.",
                "The same proof workspace also now records the explicit Morton square-coordinate shell in `Proofs/HeroCaseMortonSquaredCoordinateObstruction.lean`: named coordinate `X_mt = F^2`, the source-faithful Proposition 3.2 template `X_mt,2^2 - (X_mt^2 - 4*X_mt + 1)*X_mt,2 + X_mt^2 = 0`, and the current universal constant-term obstruction `(t^0, 4)` across all sampled tail-family objects.",
                "The same proof workspace also now records the current weighted `GG` correction ladder in `Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`: first failure of `F / W_34`, normalized follow-up `G_W34`, deeper follow-up `G2_W34`, and the current no-small-correction verdicts at both normalized layers.",
                "The same proof workspace also now records the focused Weber residual bridge in `Proofs/HeroCaseWeberClassInvariantBridgeWaypoint.lean`: primary residual `G_g12_ws`, companion `G_p12_ws`, exact coordinate bridge `g12_ws^4*p12_ws^2 - g12_ws^2*p12_ws^4 + 48*t^2*g12_ws^2*p12_ws^2 + 4096*t^6 = 0`, the classical Weber `f2` tri-product coordinate `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2) = G_g12_ws*G_p12_ws`, its normalized follow-up `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)`, and the source-backed reading that Berndt--Chan--Zhang identify Ramanujan-Weber `G_n` / `g_n` with classical Weber `f` / `f1` while Yui--Zagier supplies the classical Weber `f`, `f1`, `f2` trio, so this `g12_ws` / `p12_ws` / `G_f2_ws` shell should now be read as that named Weber trio in the project's normalization. The same waypoint then keeps the derived quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2`, the exact quotient-coordinate bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0` with `Q_gp_ws = p12_ws / g12_ws`, the template-normalized coordinate `G_X_ws = X_g_ws*(t^2; t^4)_inf^24 / (16*t^2) = 1 / G_g12_ws^2`, the normalized follow-up `H_X_ws = (G_X_ws - 1) / (4*t^1)`, the direct follow-up bridge `D_XR_ws = H_gp_ws - H_X_ws` and `Q_XR_ws = H_gp_ws / H_X_ws`, the quotient-bridge follow-up `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, the next bridge back to the template-normalized lane `D_XK_ws = K_XR_ws - H_X_ws` and `Q_XK_ws = K_XR_ws / H_X_ws`, the quotient-follow-up `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)`, the focused quotient `R_gp_ws = G_p12_ws / G_g12_ws`, the normalized follow-up `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`, the theorem-shaped return-bridge closure shell saying `Q_XK_ws` and `L_XK_ws` still have `0` hits in the checked direct / quotient / mixed named `GG` modular-equation boxes, now also exposes the exact Chan--Huang obstruction quartets for those two return-bridge objects, and still records the current no-hit verdicts in the first self-polynomial, self-fractional-linear, self-quotient-product, eta, modular-unit, and plus-Pochhammer boxes on the hero-side follow-up objects.",
                (
                    "Current source-family-specific exact lanes: zero-shift f2/gcf3 and f4/gcf2"
                    " n-dependent equivalence, plus the nearest unit-a-shift, unit-b-shift, mixed unit-a/unit-b-shift, mixed unit-a/unit-lambda-shift, mixed unit-b/unit-lambda-shift, mixed unit-a/unit-b/unit-lambda-shift, and unit-lambda-shift lanes,"
                    " formalized in Proofs/HeroCasePage43Equivalence.lean,"
                    f" with zero-shift final surviving obstruction coefficients {_format_expr(context.f2_equivalence.final_m2_coefficient)}"
                    f" and {_format_expr(context.f4_equivalence.final_m2_coefficient)}, unit-a-shift final"
                    f" surviving coefficients {_format_expr(context.f2_unit_a_shift_equivalence.final_m2_coefficient)}"
                    f" and {_format_expr(context.f4_unit_a_shift_equivalence.final_m2_coefficient)}, unit-b-shift final"
                    f" surviving coefficients {_format_expr(context.f2_unit_b_shift_equivalence.final_m2_coefficient)}"
                    f" and {_format_expr(context.f4_unit_b_shift_equivalence.final_m2_coefficient)}, mixed unit-a/unit-b-shift final"
                    f" surviving coefficients {_format_expr(context.f2_unit_ab_shift_equivalence.final_m2_coefficient)}"
                    f" and {_format_expr(context.f4_unit_ab_shift_equivalence.final_m2_coefficient)}, mixed unit-a/unit-lambda-shift impossible"
                    f" m1 coefficients {_format_expr(context.f2_unit_a_lambda_shift_equivalence.impossible_m1_coefficient)}"
                    f" and {_format_expr(context.f4_unit_a_lambda_shift_equivalence.impossible_m1_coefficient)}, mixed unit-b/unit-lambda-shift impossible"
                    f" m1 coefficients {_format_expr(context.f2_unit_b_lambda_shift_equivalence.impossible_m1_coefficient)}"
                    f" and {_format_expr(context.f4_unit_b_lambda_shift_equivalence.impossible_m1_coefficient)}, mixed unit-a/unit-b/unit-lambda-shift impossible"
                    f" m1 coefficients {_format_expr(context.f2_unit_ab_lambda_shift_equivalence.impossible_m1_coefficient)}"
                    f" and {_format_expr(context.f4_unit_ab_lambda_shift_equivalence.impossible_m1_coefficient)}, and shifted impossible"
                    f" m1 coefficients {_format_expr(context.f2_unit_lambda_shift_equivalence.impossible_m1_coefficient)}"
                    f" and {_format_expr(context.f4_unit_lambda_shift_equivalence.impossible_m1_coefficient)}."
                    " These lanes are also packaged there as the nearest-shift cube summary theorems"
                    " `noNearestShiftCubeF2ExactEquivalence`, `noNearestShiftCubeF4ExactEquivalence`, and"
                    " `noNearestShiftCubeExactEquivalence`, together with the Bool-parameterized"
                    " variants `noNearestShiftCubeF2ExactEquivalenceFor`,"
                    " `noNearestShiftCubeF4ExactEquivalenceFor`, and"
                    " `noNearestShiftCubeExactEquivalenceFor`."
                    if (
                        context.f2_equivalence is not None
                        and context.f4_equivalence is not None
                        and context.f2_unit_a_shift_equivalence is not None
                        and context.f4_unit_a_shift_equivalence is not None
                        and context.f2_unit_b_shift_equivalence is not None
                        and context.f4_unit_b_shift_equivalence is not None
                        and context.f2_unit_ab_shift_equivalence is not None
                        and context.f4_unit_ab_shift_equivalence is not None
                        and context.f2_unit_a_lambda_shift_equivalence is not None
                        and context.f4_unit_a_lambda_shift_equivalence is not None
                        and context.f2_unit_b_lambda_shift_equivalence is not None
                        and context.f4_unit_b_lambda_shift_equivalence is not None
                        and context.f2_unit_ab_lambda_shift_equivalence is not None
                        and context.f4_unit_ab_lambda_shift_equivalence is not None
                        and context.f2_unit_lambda_shift_equivalence is not None
                        and context.f4_unit_lambda_shift_equivalence is not None
                    )
                    else (
                        "Current source-family-specific exact lane: zero-shift f2/gcf3 n-dependent equivalence,"
                        f" whose final surviving obstruction coefficient is {_format_expr(context.f2_equivalence.final_m2_coefficient)}."
                        if context.f2_equivalence is not None
                        else "No source-family-specific exact lane has been attached yet."
                    )
                ),
                (
                    "The same module also now excludes the full zero-shift single-prefactor box `phi in {1, 1+t, 1/(1+t)}` with at most one non-plain prefactor active in both page-43 families, via the polynomial theorems `noZeroShiftPolynomialSinglePrefactorF2DirectMatches`, `noZeroShiftPolynomialSinglePrefactorF4DirectMatches`, `noZeroShiftPolynomialSinglePrefactorDirectMatches` and the reciprocal cross-multiplied theorems `noZeroShiftReciprocalSinglePrefactorF2CrossMatches`, `noZeroShiftReciprocalSinglePrefactorF4CrossMatches`, `noZeroShiftReciprocalSinglePrefactorCrossMatches`."
                    if (
                        context.f2_polynomial_prefactor_obstructions
                        and context.f4_polynomial_prefactor_obstructions
                        and context.f2_reciprocal_prefactor_obstructions
                        and context.f4_reciprocal_prefactor_obstructions
                    )
                    else ""
                ),
                "-/",
                "",
                "",
                "theorem rr_direct_obstruction :",
                "    rrForcedStage1Numerator ≠ candidateData.a 0 := by",
                "  intro h",
                "  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h",
                "  norm_num [rrForcedStage1Numerator, candidateData] at hEval",
                "",
                "theorem cubic_direct_obstruction :",
                "    cubicForcedStage2Numerator ≠ candidateData.a 1 := by",
                "  intro h",
                "  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h",
                "  norm_num [cubicForcedStage2Numerator, candidateData] at hEval",
                "",
                "theorem odd_contraction_obstruction :",
                "    cubicOddInitialDenominator ≠ candidateData.b0 := by",
                "  intro h",
                "  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h",
                "  norm_num [cubicOddInitialDenominator, candidateData] at hEval",
                "",
                "theorem even_contraction_obstruction :",
                "    cubicEvenFirstDenominator ≠ candidateData.b 0 := by",
                "  intro h",
                "  have hEval := congrArg (fun p : QPoly => Polynomial.eval (1 : Rat) p) h",
                "  norm_num [cubicEvenFirstDenominator, candidateData] at hEval",
                "",
                "/-!",
                "Suggested next theorem extensions:",
                "",
                "1. Bridge the current rational-function reverse-equivalence layer from",
                "   finite convergents to the infinite object.",
                "2. Upgrade the focused Weber residual bridge around `G_g12_ws`,",
                "   `G_p12_ws`, the classical Weber `f2` tri-product coordinate",
                "   `G_f2_ws = (g12_ws*p12_ws*(-t^4; t^4)_inf^12) / (64*t^2)`,",
                "   its normalized follow-up `H_f2_ws = (G_f2_ws - 1) / (-4*t^1)`,",
                "   the source-backed classical Weber reading of that shell, the derived",
                "   quotient-coordinate `X_g_ws = 16*t^2 / g12_ws^2`,",
                "   the exact raw quotient bridge `Q_gp_ws^4 - (1 + 3*X_g_ws)*Q_gp_ws^2 - X_g_ws^3 = 0`,",
                "   the template-normalized coordinate `G_X_ws = 1 / G_g12_ws^2`,",
                "   the direct normalized-follow-up bridge `D_XR_ws = H_gp_ws - H_X_ws` /",
                "   `Q_XR_ws = H_gp_ws / H_X_ws`, its quotient follow-up",
                "   `K_XR_ws = (Q_XR_ws - 1) / (-24*t^2)`, the return bridge",
                "   `D_XK_ws = K_XR_ws - H_X_ws` / `Q_XK_ws = K_XR_ws / H_X_ws`,",
                "   its quotient follow-up `L_XK_ws = (Q_XK_ws - 1) / (-5/2*t^1)`, the named-GG return-bridge closure shell",
                "   saying those two objects still have `0` hits in the checked direct / quotient / mixed named `GG` boxes, `R_gp_ws = G_p12_ws / G_g12_ws`, and the normalized",
                "   follow-ups `H_X_ws = (G_X_ws - 1) / (4*t^1)` and `H_gp_ws = (R_gp_ws - 1) / (96*t^3)`",
                "   from a research waypoint shell into an exact uniqueness / closure theorem.",
                "3. Extend the page-43 exact layer beyond the currently formalized",
                "   nearest-shift cube plus the exact zero-shift single-prefactor box,",
                "   especially if shifted or multi-prefactor rational lanes matter later.",
                "4. Compare candidate convergents against nearby benchmark convergents.",
                "5. Formalize the Bauer-Muir transform algebra itself instead of injecting only",
                "   the recovered witnesses.",
                "6. Attempt a final source theorem only after a unique identity is known.",
                "-/",
                "",
                "end",
                "",
                f"end {namespace_name}",
                "end Generated",
                "end Proofs",
            ]
        )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_candidate_formalization_assets(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    max_stride: int = 4,
    lean_output_path: str | None = None,
    smoke: bool = False,
) -> None:
    context = _build_formalization_context(
        input_path=input_path,
        candidate_id=candidate_id,
        max_stride=max_stride,
        smoke=smoke,
    )
    _write_formalization_note(context, output_path)
    if lean_output_path:
        _write_lean_skeleton(context, lean_output_path)


def build_candidate_formalization_note(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    max_stride: int = 4,
    smoke: bool = False,
) -> None:
    build_candidate_formalization_assets(
        input_path=input_path,
        candidate_id=candidate_id,
        output_path=output_path,
        max_stride=max_stride,
        lean_output_path=None,
        smoke=smoke,
    )


def build_candidate_lean_skeleton(
    *,
    input_path: str,
    candidate_id: str,
    output_path: str,
    max_stride: int = 4,
    smoke: bool = False,
) -> None:
    context = _build_formalization_context(
        input_path=input_path,
        candidate_id=candidate_id,
        max_stride=max_stride,
        smoke=smoke,
    )
    _write_lean_skeleton(context, output_path)
