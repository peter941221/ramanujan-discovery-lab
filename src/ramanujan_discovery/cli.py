from __future__ import annotations

import argparse

from ramanujan_discovery.config import SearchConfig, VerificationConfig
from ramanujan_discovery.discovery import discover_candidates
from ramanujan_discovery.reporting import build_report, build_site
from ramanujan_discovery.storage import write_candidates
from ramanujan_discovery.verification import verify_candidates


def _parse_q_values(raw: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in raw.split(",") if item.strip())


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

    if args.command == "site":
        build_site(args.input_path, args.out_dir, title=args.title)
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2
