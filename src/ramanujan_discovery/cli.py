from __future__ import annotations

import argparse

from ramanujan_discovery.analysis import build_candidate_analysis_note, build_candidate_terminal_summary
from ramanujan_discovery.config import SearchConfig, VerificationConfig
from ramanujan_discovery.discovery import discover_candidates
from ramanujan_discovery.formalization import build_candidate_formalization_assets
from ramanujan_discovery.identification import (
    build_candidate_identification_note,
    build_candidate_tail_operator_note,
    build_candidate_tail_family_note,
)
from ramanujan_discovery.research import build_candidate_research_note
from ramanujan_discovery.reporting import build_report, build_site
from ramanujan_discovery.storage import write_candidates
from ramanujan_discovery.verification import verify_candidates


def _parse_q_values(raw: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in raw.split(",") if item.strip())

def _parse_int_list(raw: str) -> tuple[int, ...]:
    return tuple(int(item.strip()) for item in raw.split(",") if item.strip())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ramanujan-discovery")
    subparsers = parser.add_subparsers(dest="command", required=True)

    discover = subparsers.add_parser("discover", help="Search the built-in q-continued-fraction template grid.")
    discover.add_argument("--depth", type=int, default=36)
    discover.add_argument("--precision", type=int, default=80)
    discover.add_argument("--budget-hours", type=float, default=0.1)
    discover.add_argument("--q-values", type=str, default="0.05,0.09,0.13")
    discover.add_argument("--min-digits", type=int, default=18)
    discover.add_argument("--min-stability", type=int, default=18)
    discover.add_argument("--min-review-stability", type=int, default=28)
    discover.add_argument("--max-per-target", type=int, default=4)
    discover.add_argument("--max-review-candidates", type=int, default=6)
    discover.add_argument("--out", type=str, required=True)

    verify = subparsers.add_parser("verify", help="Recompute candidates at higher precision.")
    verify.add_argument("--in", dest="input_path", type=str, required=True)
    verify.add_argument("--depth", type=int, default=52)
    verify.add_argument("--precision", type=int, default=160)
    verify.add_argument("--q-values", type=str, default="0.04,0.11,0.17,0.23")
    verify.add_argument("--min-digits", type=int, default=30)
    verify.add_argument("--min-review-stability", type=int, default=42)
    verify.add_argument("--max-review-candidates", type=int, default=6)
    verify.add_argument("--out", type=str, required=True)

    report = subparsers.add_parser("report", help="Render a Markdown report from verified candidates.")
    report.add_argument("--in", dest="input_path", type=str, required=True)
    report.add_argument("--out", type=str, required=True)

    analyze = subparsers.add_parser("analyze", help="Write a focused Markdown analysis for one candidate.")
    analyze.add_argument("--in", dest="input_path", type=str, required=True)
    analyze.add_argument("--candidate-id", type=str, required=True)
    analyze.add_argument("--depth", type=int, default=12)
    analyze.add_argument("--series-order", type=int, default=31)
    analyze.add_argument("--stdout-format", choices=("unicode", "plain", "latex", "none"), default="unicode")
    analyze.add_argument("--out", type=str, required=True)

    research = subparsers.add_parser(
        "research",
        help="Run heavier candidate research helpers (ratio-series fits, literature coefficient checks).",
    )
    research.add_argument("--in", dest="input_path", type=str, required=True)
    research.add_argument("--candidate-id", type=str, required=True)
    research.add_argument("--depth", type=int, default=40)
    research.add_argument("--series-order", type=int, default=151)
    research.add_argument("--smoke", action="store_true", help="Use a reduced-cost profile for fast local validation.")
    research.add_argument("--out", type=str, required=True)

    formalize = subparsers.add_parser(
        "formalize",
        help="Write a formalization-prep note for one candidate without invoking a proof assistant.",
    )
    formalize.add_argument("--in", dest="input_path", type=str, required=True)
    formalize.add_argument("--candidate-id", type=str, required=True)
    formalize.add_argument("--max-stride", type=int, default=4)
    formalize.add_argument("--smoke", action="store_true", help="Use a reduced-cost profile for fast local validation.")
    formalize.add_argument("--out", type=str, required=True)
    formalize.add_argument("--lean-out", type=str, default=None)

    identify = subparsers.add_parser(
        "identify",
        help="Attempt to relate a candidate to its nearest benchmark (small algebraic relation guesses).",
    )
    identify.add_argument("--in", dest="input_path", type=str, required=True)
    identify.add_argument("--candidate-id", type=str, required=True)
    identify.add_argument("--depth", type=int, default=40)
    identify.add_argument("--series-order", type=int, default=90)
    identify.add_argument("--max-degree", type=int, default=4)
    identify.add_argument(
        "--benchmark-powers",
        type=str,
        default="",
        help="Comma-separated k>=2 values to include Bk = benchmark_recip(t^k) in a multivariate relation search.",
    )
    identify.add_argument("--smoke", action="store_true", help="Use a reduced-cost profile for fast local validation.")
    identify.add_argument("--out", type=str, required=True)

    tail_note = subparsers.add_parser(
        "tail-note",
        help="Write a tail-family-first source-recognition note for one candidate.",
    )
    tail_note.add_argument("--in", dest="input_path", type=str, required=True)
    tail_note.add_argument("--candidate-id", type=str, required=True)
    tail_note.add_argument("--depth", type=int, default=40)
    tail_note.add_argument("--series-order", type=int, default=36)
    tail_note.add_argument(
        "--tail-stages",
        type=str,
        default="3,4,5",
        help="Comma-separated start stages whose normalized tail samples should be scanned.",
    )
    tail_note.add_argument("--max-gap-depth", type=int, default=3)
    tail_note.add_argument("--smoke", action="store_true", help="Use a reduced-cost profile for fast local validation.")
    tail_note.add_argument("--out", type=str, required=True)

    tail_operator_note = subparsers.add_parser(
        "tail-operator-note",
        help="Write a tail-family-first operator / q-difference note for one candidate.",
    )
    tail_operator_note.add_argument("--in", dest="input_path", type=str, required=True)
    tail_operator_note.add_argument("--candidate-id", type=str, required=True)
    tail_operator_note.add_argument("--depth", type=int, default=40)
    tail_operator_note.add_argument("--series-order", type=int, default=36)
    tail_operator_note.add_argument(
        "--tail-stages",
        type=str,
        default="3,4,5",
        help="Comma-separated start stages whose normalized tail samples should be scanned.",
    )
    tail_operator_note.add_argument("--max-gap-depth", type=int, default=3)
    tail_operator_note.add_argument("--smoke", action="store_true", help="Use a reduced-cost profile for fast local validation.")
    tail_operator_note.add_argument("--out", type=str, required=True)

    site = subparsers.add_parser("site", help="Render a GitHub Pages-friendly static site.")
    site.add_argument("--in", dest="input_path", type=str, required=True)
    site.add_argument("--out-dir", type=str, required=True)
    site.add_argument("--title", type=str, default="Ramanujan Discovery Lab")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "discover":
        config = SearchConfig(
            depth=args.depth,
            precision=args.precision,
            budget_hours=args.budget_hours,
            q_values=_parse_q_values(args.q_values),
            min_discovery_digits=args.min_digits,
            min_stability_digits=args.min_stability,
            min_review_stability=args.min_review_stability,
            max_per_target=args.max_per_target,
            max_review_candidates=args.max_review_candidates,
        )
        records = discover_candidates(config)
        write_candidates(args.out, records)
        return 0

    if args.command == "verify":
        config = VerificationConfig(
            depth=args.depth,
            precision=args.precision,
            q_values=_parse_q_values(args.q_values),
            min_verified_digits=args.min_digits,
            min_review_stability=args.min_review_stability,
            max_review_candidates=args.max_review_candidates,
        )
        records = verify_candidates(args.input_path, config)
        write_candidates(args.out, records)
        return 0

    if args.command == "report":
        build_report(args.input_path, args.out)
        return 0

    if args.command == "analyze":
        build_candidate_analysis_note(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            depth=args.depth,
            series_order=args.series_order,
        )
        if args.stdout_format != "none":
            print(
                build_candidate_terminal_summary(
                    input_path=args.input_path,
                    candidate_id=args.candidate_id,
                    depth=args.depth,
                    series_order=args.series_order,
                    math_format=args.stdout_format,
                )
            )
        return 0

    if args.command == "research":
        build_candidate_research_note(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            depth=args.depth,
            series_order=args.series_order,
            smoke=args.smoke,
        )
        return 0

    if args.command == "formalize":
        build_candidate_formalization_assets(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            max_stride=args.max_stride,
            lean_output_path=args.lean_out,
            smoke=args.smoke,
        )
        return 0

    if args.command == "identify":
        build_candidate_identification_note(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            depth=args.depth,
            series_order=args.series_order,
            max_degree=args.max_degree,
            benchmark_powers=_parse_int_list(args.benchmark_powers),
            smoke=args.smoke,
        )
        return 0

    if args.command == "tail-note":
        build_candidate_tail_family_note(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            depth=args.depth,
            series_order=args.series_order,
            tail_stages=_parse_int_list(args.tail_stages),
            max_gap_depth=args.max_gap_depth,
            smoke=args.smoke,
        )
        return 0

    if args.command == "tail-operator-note":
        build_candidate_tail_operator_note(
            input_path=args.input_path,
            candidate_id=args.candidate_id,
            output_path=args.out,
            depth=args.depth,
            series_order=args.series_order,
            tail_stages=_parse_int_list(args.tail_stages),
            max_gap_depth=args.max_gap_depth,
            smoke=args.smoke,
        )
        return 0

    if args.command == "site":
        build_site(args.input_path, args.out_dir, title=args.title)
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2
