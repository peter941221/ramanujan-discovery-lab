from __future__ import annotations

from pathlib import Path

import sympy as sp

from ramanujan_discovery.benchmarks import BENCHMARKS, get_benchmark
from ramanujan_discovery.models import CandidateRecord, QCFTemplate
from ramanujan_discovery.storage import read_candidates


def _series_expr(template: QCFTemplate, depth: int, order: int, q_symbol: sp.Symbol | None = None):
    if depth < 1:
        raise ValueError("depth must be at least 1")
    if order < 2:
        raise ValueError("order must be at least 2")

    q = q_symbol if q_symbol is not None else sp.Symbol("q")
    tail = None

    for index in range(depth, 0, -1):
        position = index - 1
        numerator_exponent = template.numerator_q_shift + template.numerator_q_step * position
        numerator = template.numerator_scale * q**numerator_exponent
        if template.numerator_extra_scale != 0:
            extra_exponent = template.numerator_extra_q_shift + template.numerator_extra_q_step * position
            numerator += template.numerator_extra_scale * q**extra_exponent

        denominator = sp.Integer(template.denominator_constant)
        if template.denominator_scale != 0:
            denominator_exponent = template.denominator_q_shift + template.denominator_q_step * position
            denominator += template.denominator_scale * q**denominator_exponent

        if tail is None:
            tail = sp.series(numerator / denominator, q, 0, order).removeO()
        else:
            tail = sp.series(numerator / (denominator + tail), q, 0, order).removeO()

    return sp.expand(
        sp.series(
            sp.Integer(template.top_constant) / (sp.Integer(template.base_denominator) + tail),
            q,
            0,
            order,
        ).removeO()
    )


def _first_difference_order(expr) -> int | None:
    exponents = _difference_exponents(expr)
    if not exponents:
        return None
    return min(exponents)


def _difference_exponents(expr) -> list[int]:
    expanded = sp.expand(expr)
    if expanded == 0:
        return []

    q = sp.Symbol("q")
    terms = expanded.as_ordered_terms()
    exponents: list[int] = []
    for term in terms:
        coefficient, exponent = term.as_coeff_exponent(q)
        if coefficient != 0 and exponent.is_integer:
            exponents.append(int(exponent))

    return sorted(set(exponents))


def _low_order_multiplier(candidate_series, benchmark_series, terms: int = 3):
    exponents = _difference_exponents(candidate_series - benchmark_series)[:terms]
    if not exponents:
        return None

    q = sp.Symbol("q")
    unknowns = sp.symbols(f"m0:{len(exponents)}")
    modifier = sp.Integer(1)
    for unknown, exponent in zip(unknowns, exponents):
        modifier += unknown * q**exponent

    max_order = max(exponents) + 1
    residual = sp.expand(sp.series(benchmark_series * modifier - candidate_series, q, 0, max_order).removeO())
    equations = [sp.Eq(residual.coeff(q, exponent), 0) for exponent in exponents]
    solutions = sp.solve(equations, unknowns, dict=True)
    if not solutions:
        return None

    fitted_modifier = sp.expand(modifier.subs(solutions[0]))
    fitted_residual = sp.expand(
        sp.series(benchmark_series * fitted_modifier - candidate_series, q, 0, max_order).removeO()
    )
    return fitted_modifier, fitted_residual


def _format_expr(expr) -> str:
    return str(sp.expand(expr)).replace("**", "^")


def _find_candidate(records: list[CandidateRecord], candidate_id: str) -> CandidateRecord:
    for record in records:
        if record.id == candidate_id:
            return record
    raise KeyError(f"unknown candidate id: {candidate_id}")


def _same_step_alternative(record: CandidateRecord) -> tuple[str, str] | None:
    closest = record.closest_benchmark
    if closest not in BENCHMARKS:
        return None

    step = BENCHMARKS[closest].canonical_template.numerator_q_step
    if closest.startswith("rogers_ramanujan"):
        prefix = "ramanujan_cubic"
    elif closest.startswith("ramanujan_cubic"):
        prefix = "rogers_ramanujan"
    else:
        return None

    for name, benchmark in BENCHMARKS.items():
        template = benchmark.canonical_template
        if name.startswith(prefix) and template.numerator_q_step == step:
            return name, benchmark.description
    return None


def analyze_candidate(
    input_path: str,
    candidate_id: str,
    depth: int = 12,
    series_order: int = 31,
) -> dict[str, object]:
    records = read_candidates(input_path)
    record = _find_candidate(records, candidate_id)
    benchmark = get_benchmark(record.closest_benchmark)

    candidate_series = _series_expr(record.template, depth=depth, order=series_order)
    benchmark_series = _series_expr(benchmark.canonical_template, depth=depth, order=series_order)
    delta_series = sp.expand(candidate_series - benchmark_series)
    first_difference = _first_difference_order(delta_series)
    multiplier_fit = _low_order_multiplier(candidate_series, benchmark_series)

    alternative_name = None
    alternative_description = None
    alternative_series = None
    alternative_delta = None
    alternative_first_difference = None
    alternative_match = _same_step_alternative(record)
    if alternative_match is not None:
        alternative_name, alternative_description = alternative_match
        alternative_series = _series_expr(
            get_benchmark(alternative_name).canonical_template,
            depth=depth,
            order=series_order,
        )
        alternative_delta = sp.expand(candidate_series - alternative_series)
        alternative_first_difference = _first_difference_order(alternative_delta)

    return {
        "record": record,
        "benchmark": benchmark,
        "candidate_series": candidate_series,
        "benchmark_series": benchmark_series,
        "delta_series": delta_series,
        "first_difference": first_difference,
        "multiplier_fit": multiplier_fit,
        "alternative_name": alternative_name,
        "alternative_description": alternative_description,
        "alternative_series": alternative_series,
        "alternative_delta": alternative_delta,
        "alternative_first_difference": alternative_first_difference,
        "benchmark_template": benchmark.canonical_template.normalized(),
        "candidate_template": record.template.normalized(),
        "depth": depth,
        "series_order": series_order,
    }


def _render_expr(expr, math_format: str) -> str:
    expanded = sp.expand(expr)
    if math_format == "unicode":
        return sp.pretty(expanded, use_unicode=True)
    if math_format == "latex":
        return sp.latex(expanded)
    return _format_expr(expanded)


def build_candidate_terminal_summary(
    input_path: str,
    candidate_id: str,
    depth: int = 12,
    series_order: int = 31,
    math_format: str = "unicode",
) -> str:
    context = analyze_candidate(
        input_path=input_path,
        candidate_id=candidate_id,
        depth=depth,
        series_order=series_order,
    )
    record = context["record"]
    assert isinstance(record, CandidateRecord)

    first_difference = context["first_difference"]
    alternative_name = context["alternative_name"]
    alternative_first_difference = context["alternative_first_difference"]
    multiplier_fit = context["multiplier_fit"]

    lines = [
        f"[{record.id}]",
        f"closest benchmark: {record.closest_benchmark} ({record.closest_benchmark_digits} shared digits)",
        (
            f"first divergence vs closest: q^{first_difference}"
            if first_difference is not None
            else "first divergence vs closest: exact at this order"
        ),
    ]

    if alternative_name is not None:
        lines.append(
            (
                f"first divergence vs {alternative_name}: q^{alternative_first_difference}"
                if alternative_first_difference is not None
                else f"first divergence vs {alternative_name}: exact at this order"
            )
        )

    if multiplier_fit is not None:
        fitted_modifier, _ = multiplier_fit
        lines.extend(
            [
                "low-order multiplicative fit:",
                f"candidate / {record.closest_benchmark} =",
                _render_expr(fitted_modifier, math_format),
            ]
        )

    if math_format == "latex":
        lines.append("copyable LaTeX above; plain terminals do not natively render TeX.")

    return "\n".join(lines)


def build_candidate_analysis_note(
    input_path: str,
    candidate_id: str,
    output_path: str,
    depth: int = 12,
    series_order: int = 31,
) -> None:
    records = read_candidates(input_path)
    record = _find_candidate(records, candidate_id)
    benchmark = get_benchmark(record.closest_benchmark)

    candidate_series = _series_expr(record.template, depth=depth, order=series_order)
    benchmark_series = _series_expr(benchmark.canonical_template, depth=depth, order=series_order)
    delta_series = sp.expand(candidate_series - benchmark_series)
    first_difference = _first_difference_order(delta_series)
    multiplier_fit = _low_order_multiplier(candidate_series, benchmark_series)

    alternative_name = None
    alternative_description = None
    alternative_series = None
    alternative_delta = None
    alternative_first_difference = None
    alternative_match = _same_step_alternative(record)
    if alternative_match is not None:
        alternative_name, alternative_description = alternative_match
        alternative_series = _series_expr(
            get_benchmark(alternative_name).canonical_template,
            depth=depth,
            order=series_order,
        )
        alternative_delta = sp.expand(candidate_series - alternative_series)
        alternative_first_difference = _first_difference_order(alternative_delta)

    benchmark_template = benchmark.canonical_template.normalized()
    candidate_template = record.template.normalized()

    lines = [
        f"# Hero Case Analysis: `{record.id}`",
        "",
        "## Snapshot",
        "",
        f"- Candidate id: `{record.id}`",
        f"- Closest benchmark: `{record.closest_benchmark}` ({record.closest_benchmark_digits} shared digits)",
        f"- Benchmark kind: `{benchmark.kind}`",
        f"- Novelty status: `{record.novelty_status}`",
        f"- Family bucket: `{record.family_bucket}`",
        f"- Equivalence key: `{record.equivalence_key}`",
        "",
        "## Structural Delta vs Closest Benchmark",
        "",
        f"- Candidate template: `{record.template.signature()}`",
        f"- Benchmark template: `{benchmark.canonical_template.signature()}`",
        f"- Main numerator shift delta: `{candidate_template.numerator_q_shift - benchmark_template.numerator_q_shift}`",
        f"- Main numerator step delta: `{candidate_template.numerator_q_step - benchmark_template.numerator_q_step}`",
        f"- Extra numerator present: `{candidate_template.numerator_extra_scale != 0}`",
        f"- Denominator perturbation present: `{candidate_template.denominator_scale != 0}`",
        "",
        "## Symbolic q-Series",
        "",
        f"- Candidate series (depth `{depth}`, order `{series_order}`):",
        "```text",
        _format_expr(candidate_series),
        "```",
        f"- Benchmark series `{record.closest_benchmark}`:",
        "```text",
        _format_expr(benchmark_series),
        "```",
        "- Candidate minus benchmark:",
        "```text",
        _format_expr(delta_series),
        "```",
        f"- First divergence order: `{first_difference if first_difference is not None else 'exact at this order'}`",
    ]

    if multiplier_fit is not None:
        fitted_modifier, fitted_residual = multiplier_fit
        lines.extend(
            [
                "",
                "## Low-Order Multiplicative Fit",
                "",
                "- Fitted modifier multiplying the closest benchmark:",
                "```text",
                _format_expr(fitted_modifier),
                "```",
                "- Truncated residual after this fit:",
                "```text",
                _format_expr(fitted_residual),
                "```",
            ]
        )

    if alternative_name is not None and alternative_series is not None and alternative_delta is not None:
        lines.extend(
            [
                "",
                "## Same-Step Cross-Family Comparison",
                "",
                f"- Alternative benchmark: `{alternative_name}`",
                f"- Description: {alternative_description}",
                "```text",
                _format_expr(alternative_series),
                "```",
                "- Candidate minus alternative benchmark:",
                "```text",
                _format_expr(alternative_delta),
                "```",
                (
                    f"- First divergence order vs `{alternative_name}`: "
                    f"`{alternative_first_difference if alternative_first_difference is not None else 'exact at this order'}`"
                ),
            ]
        )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
