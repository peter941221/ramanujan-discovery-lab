from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QCFTemplate:
    top_constant: int = 1
    base_denominator: int = 1
    base_denominator_q_scale: int = 0
    base_denominator_q_shift: int = 1
    numerator_scale: int = 1
    numerator_q_shift: int = 1
    numerator_q_step: int = 1
    numerator_extra_scale: int = 0
    numerator_extra_q_shift: int = 1
    numerator_extra_q_step: int = 1
    denominator_constant: int = 1
    denominator_scale: int = 0
    denominator_q_shift: int = 1
    denominator_q_step: int = 1

    def _normalized_parts(self) -> tuple[int, int, int, int, int, int, int, int, int, int, int, int, int, int]:
        base_denominator_q_shift = self.base_denominator_q_shift
        if self.base_denominator_q_scale == 0:
            base_denominator_q_shift = 0

        numerator_extra_q_shift = self.numerator_extra_q_shift
        numerator_extra_q_step = self.numerator_extra_q_step
        if self.numerator_extra_scale == 0:
            numerator_extra_q_shift = 0
            numerator_extra_q_step = 0

        denominator_q_shift = self.denominator_q_shift
        denominator_q_step = self.denominator_q_step
        if self.denominator_scale == 0:
            denominator_q_shift = 0
            denominator_q_step = 0

        return (
            self.top_constant,
            self.base_denominator,
            self.base_denominator_q_scale,
            base_denominator_q_shift,
            self.numerator_scale,
            self.numerator_q_shift,
            self.numerator_q_step,
            self.numerator_extra_scale,
            numerator_extra_q_shift,
            numerator_extra_q_step,
            self.denominator_constant,
            self.denominator_scale,
            denominator_q_shift,
            denominator_q_step,
        )

    def normalized(self) -> "QCFTemplate":
        return QCFTemplate(*self._normalized_parts())

    def signature(self) -> str:
        (
            top_constant,
            base_denominator,
            base_denominator_q_scale,
            base_denominator_q_shift,
            numerator_scale,
            numerator_q_shift,
            numerator_q_step,
            numerator_extra_scale,
            numerator_extra_q_shift,
            numerator_extra_q_step,
            denominator_constant,
            denominator_scale,
            denominator_q_shift,
            denominator_q_step,
        ) = self._normalized_parts()
        base_part = f"top={top_constant};base={base_denominator}"
        if base_denominator_q_scale != 0:
            base_part += f";bqs={base_denominator_q_scale};bqr={base_denominator_q_shift}"
        return (
            "qcf:"
            f"{base_part};"
            f"ns={numerator_scale};nr={numerator_q_shift};nk={numerator_q_step};"
            f"nes={numerator_extra_scale};ner={numerator_extra_q_shift};nek={numerator_extra_q_step};"
            f"dc={denominator_constant};ds={denominator_scale};"
            f"dr={denominator_q_shift};dk={denominator_q_step}"
        )

    def latex(self) -> str:
        base_denominator = str(self.base_denominator)
        if self.base_denominator_q_scale != 0:
            base_denominator += (
                " + "
                f"{self.base_denominator_q_scale} q^{{{self.base_denominator_q_shift}}}"
            )

        numerator = (
            f"{self.numerator_scale} q^{{{self.numerator_q_shift} + {self.numerator_q_step}(n-1)}}"
        )
        if self.numerator_extra_scale != 0:
            numerator += (
                " + "
                f"{self.numerator_extra_scale} q^{{{self.numerator_extra_q_shift} + "
                f"{self.numerator_extra_q_step}(n-1)}}"
            )

        return (
            r"\cfrac{1}{"
            f"{base_denominator}"
            r" + \cfrac{"
            f"{numerator}"
            r"}{"
            f"{self.denominator_constant} + "
            f"{self.denominator_scale} q^{{{self.denominator_q_shift} + {self.denominator_q_step}(n-1)}}"
            r" + \ddots}}"
        )

    def to_dict(self) -> dict[str, int]:
        return {
            "top_constant": self.top_constant,
            "base_denominator": self.base_denominator,
            "base_denominator_q_scale": self.base_denominator_q_scale,
            "base_denominator_q_shift": self.base_denominator_q_shift,
            "numerator_scale": self.numerator_scale,
            "numerator_q_shift": self.numerator_q_shift,
            "numerator_q_step": self.numerator_q_step,
            "numerator_extra_scale": self.numerator_extra_scale,
            "numerator_extra_q_shift": self.numerator_extra_q_shift,
            "numerator_extra_q_step": self.numerator_extra_q_step,
            "denominator_constant": self.denominator_constant,
            "denominator_scale": self.denominator_scale,
            "denominator_q_shift": self.denominator_q_shift,
            "denominator_q_step": self.denominator_q_step,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, int]) -> "QCFTemplate":
        return cls(**payload)

    def complexity_score(self) -> int:
        score = 0
        if self.base_denominator_q_scale != 0:
            score += 2
        if self.numerator_scale != 1:
            score += 1
        if self.numerator_q_shift != 1:
            score += 1
        if self.numerator_q_step != 1:
            score += 1
        if self.numerator_extra_scale != 0:
            score += 2
        if self.denominator_scale != 0:
            score += 2
        if self.denominator_constant != 1:
            score += 1
        return score

    def exploratory_family_key(self, closest_benchmark: str) -> str:
        if self.base_denominator_q_scale != 0 and self.numerator_extra_scale == 0 and self.denominator_scale == 0:
            family_kind = "base_perturbed_family"
        elif self.numerator_extra_scale == 0 and self.denominator_scale == 0:
            family_kind = "single_numerator_family"
        elif self.numerator_extra_scale != 0 and self.denominator_scale == 0:
            family_kind = "double_numerator_family"
        elif self.numerator_extra_scale == 0 and self.denominator_scale != 0:
            family_kind = "denominator_perturbed_family"
        else:
            family_kind = "hybrid_perturbed_family"

        extra_ratio = 0
        if self.numerator_extra_scale != 0 and self.numerator_q_step != 0:
            extra_ratio = self.numerator_extra_q_step // self.numerator_q_step

        return (
            f"{family_kind}::"
            f"num_scale={self.numerator_scale}::"
            f"extra_ratio={extra_ratio}::"
            f"den_scale={self.denominator_scale}"
        )


@dataclass(frozen=True)
class BenchmarkDefinition:
    name: str
    kind: str
    description: str
    canonical_template: QCFTemplate


@dataclass
class CandidateRecord:
    id: str
    template: QCFTemplate
    q_values: list[float]
    value_estimates: list[str]
    matched_target: str
    closest_benchmark: str
    closest_benchmark_digits: int
    family_bucket: str
    equivalence_key: str
    benchmark_kind: str
    digits_agree: int
    stability_score: int
    novelty_status: str
    notes: str

    def to_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "template": self.template.to_dict(),
            "q_values": self.q_values,
            "value_estimates": self.value_estimates,
            "matched_target": self.matched_target,
            "closest_benchmark": self.closest_benchmark,
            "closest_benchmark_digits": self.closest_benchmark_digits,
            "family_bucket": self.family_bucket,
            "equivalence_key": self.equivalence_key,
            "benchmark_kind": self.benchmark_kind,
            "digits_agree": self.digits_agree,
            "stability_score": self.stability_score,
            "novelty_status": self.novelty_status,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "CandidateRecord":
        template = QCFTemplate.from_dict(dict(payload["template"]))
        family_bucket = str(payload.get("family_bucket", "legacy"))
        return cls(
            id=str(payload["id"]),
            template=template,
            q_values=[float(value) for value in payload["q_values"]],
            value_estimates=[str(value) for value in payload["value_estimates"]],
            matched_target=str(payload["matched_target"]),
            closest_benchmark=str(payload.get("closest_benchmark", payload["matched_target"])),
            closest_benchmark_digits=int(payload.get("closest_benchmark_digits", payload["digits_agree"])),
            family_bucket=family_bucket,
            equivalence_key=str(payload.get("equivalence_key", family_bucket or template.signature())),
            benchmark_kind=str(payload["benchmark_kind"]),
            digits_agree=int(payload["digits_agree"]),
            stability_score=int(payload["stability_score"]),
            novelty_status=str(payload["novelty_status"]),
            notes=str(payload["notes"]),
        )
