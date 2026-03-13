from dataclasses import dataclass


@dataclass(frozen=True)
class SearchConfig:
    depth: int = 36
    precision: int = 80
    budget_hours: float = 0.1
    q_values: tuple[float, ...] = (0.05, 0.09, 0.13)
    min_discovery_digits: int = 18
    min_stability_digits: int = 18
    min_review_stability: int = 28
    max_per_target: int = 4
    max_review_candidates: int = 6


@dataclass(frozen=True)
class VerificationConfig:
    depth: int = 52
    precision: int = 160
    q_values: tuple[float, ...] = (0.04, 0.11, 0.17, 0.23)
    min_verified_digits: int = 30
    min_review_stability: int = 42
    max_review_candidates: int = 6
