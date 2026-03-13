from __future__ import annotations

from mpmath import mp

from ramanujan_discovery.continued_fraction import evaluate_qcf, q_pochhammer
from ramanujan_discovery.models import BenchmarkDefinition, QCFTemplate


def _rr_template(step: int) -> QCFTemplate:
    return QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=step,
        numerator_q_step=step,
        denominator_constant=1,
        denominator_scale=0,
        denominator_q_shift=1,
        denominator_q_step=1,
    )


RR_TEMPLATE = _rr_template(1)
RR_Q2_TEMPLATE = _rr_template(2)
RR_Q3_TEMPLATE = _rr_template(3)
RR_Q4_TEMPLATE = _rr_template(4)

def _cubic_template(step: int) -> QCFTemplate:
    return QCFTemplate(
        numerator_scale=1,
        numerator_q_shift=step,
        numerator_q_step=step,
        numerator_extra_scale=1,
        numerator_extra_q_shift=2 * step,
        numerator_extra_q_step=2 * step,
        denominator_constant=1,
        denominator_scale=0,
        denominator_q_shift=1,
        denominator_q_step=1,
    )


CUBIC_TEMPLATE = _cubic_template(1)
CUBIC_Q2_TEMPLATE = _cubic_template(2)
CUBIC_Q3_TEMPLATE = _cubic_template(3)

SHIFTED_FIXTURE_TEMPLATE = QCFTemplate(
    numerator_scale=1,
    numerator_q_shift=2,
    numerator_q_step=1,
    denominator_constant=1,
    denominator_scale=0,
    denominator_q_shift=1,
    denominator_q_step=1,
)

DENOMINATOR_PERTURBED_TEMPLATE = QCFTemplate(
    numerator_scale=1,
    numerator_q_shift=1,
    numerator_q_step=1,
    denominator_constant=1,
    denominator_scale=1,
    denominator_q_shift=1,
    denominator_q_step=1,
)

BENCHMARKS: dict[str, BenchmarkDefinition] = {
    "rogers_ramanujan_normalized": BenchmarkDefinition(
        name="rogers_ramanujan_normalized",
        kind="classical",
        description="Normalized Rogers-Ramanujan continued fraction verified by its product formula.",
        canonical_template=RR_TEMPLATE,
    ),
    "rogers_ramanujan_q2_normalized": BenchmarkDefinition(
        name="rogers_ramanujan_q2_normalized",
        kind="classical_family",
        description="Rogers-Ramanujan family benchmark evaluated at q^2.",
        canonical_template=RR_Q2_TEMPLATE,
    ),
    "rogers_ramanujan_q3_normalized": BenchmarkDefinition(
        name="rogers_ramanujan_q3_normalized",
        kind="classical_family",
        description="Rogers-Ramanujan family benchmark evaluated at q^3.",
        canonical_template=RR_Q3_TEMPLATE,
    ),
    "rogers_ramanujan_q4_normalized": BenchmarkDefinition(
        name="rogers_ramanujan_q4_normalized",
        kind="classical_family",
        description="Rogers-Ramanujan family benchmark evaluated at q^4.",
        canonical_template=RR_Q4_TEMPLATE,
    ),
    "ramanujan_cubic_normalized": BenchmarkDefinition(
        name="ramanujan_cubic_normalized",
        kind="classical",
        description="Normalized Ramanujan cubic continued fraction verified by an infinite-product formula.",
        canonical_template=CUBIC_TEMPLATE,
    ),
    "ramanujan_cubic_q2_normalized": BenchmarkDefinition(
        name="ramanujan_cubic_q2_normalized",
        kind="classical_family",
        description="Ramanujan cubic continued fraction benchmark evaluated at q^2.",
        canonical_template=CUBIC_Q2_TEMPLATE,
    ),
    "ramanujan_cubic_q3_normalized": BenchmarkDefinition(
        name="ramanujan_cubic_q3_normalized",
        kind="classical_family",
        description="Ramanujan cubic continued fraction benchmark evaluated at q^3.",
        canonical_template=CUBIC_Q3_TEMPLATE,
    ),
    "shifted_rr_fixture": BenchmarkDefinition(
        name="shifted_rr_fixture",
        kind="internal_fixture",
        description="Internal regression fixture with numerator shift 2.",
        canonical_template=SHIFTED_FIXTURE_TEMPLATE,
    ),
    "denominator_perturbed_fixture": BenchmarkDefinition(
        name="denominator_perturbed_fixture",
        kind="internal_fixture",
        description="Internal regression fixture with a q-dependent denominator term.",
        canonical_template=DENOMINATOR_PERTURBED_TEMPLATE,
    ),
}


def benchmark_names() -> tuple[str, ...]:
    return tuple(BENCHMARKS.keys())


def get_benchmark(name: str) -> BenchmarkDefinition:
    return BENCHMARKS[name]


def _rr_product(q, precision: int, step: int):
    q_value = mp.mpf(str(q))
    q_power = q_value**step
    q_step = q_power**5
    numerator = q_pochhammer(q_power, q_step, precision) * q_pochhammer(q_power**4, q_step, precision)
    denominator = q_pochhammer(q_power**2, q_step, precision) * q_pochhammer(q_power**3, q_step, precision)
    return numerator / denominator


def _cubic_product(q, precision: int, step: int):
    q_value = mp.mpf(str(q))
    q_power = q_value**step
    q_step = q_power**6
    numerator = q_pochhammer(q_power, q_step, precision) * q_pochhammer(q_power**5, q_step, precision)
    denominator = q_pochhammer(q_power**3, q_step, precision) ** 2
    return numerator / denominator


def target_value(name: str, q: float, precision: int, depth: int):
    with mp.workdps(precision):
        if name == "rogers_ramanujan_normalized":
            return _rr_product(q, precision, step=1)

        if name == "rogers_ramanujan_q2_normalized":
            return _rr_product(q, precision, step=2)

        if name == "rogers_ramanujan_q3_normalized":
            return _rr_product(q, precision, step=3)

        if name == "rogers_ramanujan_q4_normalized":
            return _rr_product(q, precision, step=4)

        if name == "ramanujan_cubic_normalized":
            return _cubic_product(q, precision, step=1)

        if name == "ramanujan_cubic_q2_normalized":
            return _cubic_product(q, precision, step=2)

        if name == "ramanujan_cubic_q3_normalized":
            return _cubic_product(q, precision, step=3)

    if name == "shifted_rr_fixture":
        return evaluate_qcf(SHIFTED_FIXTURE_TEMPLATE, q=q, depth=depth, precision=precision)

    if name == "denominator_perturbed_fixture":
        return evaluate_qcf(DENOMINATOR_PERTURBED_TEMPLATE, q=q, depth=depth, precision=precision)

    raise KeyError(f"unknown benchmark: {name}")
