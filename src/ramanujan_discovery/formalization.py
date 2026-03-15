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
    Page43MonomialHit,
    SubsequenceContractionHit,
    _template_reciprocal_coeffs,
    arithmetic_subsequence_contraction_search,
    convergent_factor_equivalence_witness,
    continued_fraction_convergents,
    direct_bauer_muir_obstruction,
    page43_monomial_parameter_search,
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
    rr_subsequence_hits: list[SubsequenceContractionHit]
    cubic_subsequence_hits: list[SubsequenceContractionHit]
    f2_hits: list[Page43MonomialHit]
    f4_hits: list[Page43MonomialHit]
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


def _formalization_build_profile(*, smoke: bool) -> FormalizationBuildProfile:
    if smoke:
        return FormalizationBuildProfile(
            label="smoke",
            factor_depth=4,
            subsequence_stages=2,
            max_page43_shift=1,
            page43_stages=2,
        )
    return FormalizationBuildProfile(
        label="full",
        factor_depth=8,
        subsequence_stages=3,
        max_page43_shift=3,
        page43_stages=3,
    )


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
            rr_subsequence_hits=[],
            cubic_subsequence_hits=[],
            f2_hits=[],
            f4_hits=[],
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

    return FormalizationContext(
        record=record,
        benchmark=benchmark,
        step=step,
        build_profile=profile.label,
        max_stride=effective_max_stride,
        subsequence_stages=profile.subsequence_stages,
        page43_max_shift=profile.max_page43_shift,
        page43_stages=profile.page43_stages,
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
        rr_subsequence_hits=rr_subsequence_hits,
        cubic_subsequence_hits=cubic_subsequence_hits,
        f2_hits=f2_hits,
        f4_hits=f4_hits,
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
                "- These reverse scales are rational functions in `t`, so a full formalization of this step likely needs a fraction-field coefficient layer in addition to the current polynomial one.",
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
                "- These are bounded symbolic searches, useful for narrowing the theorem statement but not substitutes for a full origin proof.",
            ]
        )
        if context.rr_subsequence_hits or context.cubic_subsequence_hits or context.f2_hits or context.f4_hits:
            lines.extend(["", "```text"])
            for hit in context.rr_subsequence_hits + context.cubic_subsequence_hits:
                lines.append(f"{hit.source_label}: stride={hit.stride}, offset={hit.offset}")
            for hit in context.f2_hits + context.f4_hits:
                lines.append(
                    f"{hit.family}: A={hit.a_shift}, B={hit.b_shift}, L={hit.lambda_shift}, "
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
                "3. Add a rational-function or fraction-field coefficient layer for the reverse equivalence transform.",
                "4. Formalize the direct 1-step Bauer-Muir obstruction lemmas against the reduced target.",
                "5. Formalize odd/even contraction reconstruction and the cubic denominator mismatch lemma.",
                "6. Defer the bounded search exclusions until a final theorem statement makes them clearly necessary.",
                "7. Do not start a full Lean/Coq origin theorem until a unique source family or exact identity is identified.",
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
                "These scales are rational functions, so formalizing this reverse step will likely require a fraction-field coefficient layer.",
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
                "1. Lift the coefficient domain from polynomials to rational functions and",
                "   formalize the reverse equivalence transform.",
                "2. Compare candidate convergents against nearby benchmark convergents.",
                "3. Formalize the Bauer-Muir transform algebra itself instead of injecting only",
                "   the recovered witnesses.",
                "4. Attempt a final source theorem only after a unique identity is known.",
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
