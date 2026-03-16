# `cb60fd71d1d7`: An Audit-First Case Study in Experimental Mathematics

## One-Sentence Version

We built a local-first discovery pipeline for Ramanujan-style `q`-continued fractions, found a strong Rogers-Ramanujan-adjacent candidate, and then used a mix of symbolic elimination and Lean-checked local lemmas to show that the obvious explanations are not enough.

## Suggested Title Variants

1. `cb60fd71d1d7`: An Audit-First Experimental Mathematics Case Study
2. Searching for Ramanujan-Style Continued Fractions with Exact Local Proofs
3. A Rogers-Ramanujan-Adjacent Hero Case from a CLI-First Discovery Pipeline

## Abstract

This note presents a conservative case study from the `Ramanujan Discovery Lab`.
The project searches a bounded family of `q`-continued-fraction templates,
verifies stable candidates numerically against a benchmark catalog, and then
audits the strongest unexplained examples using symbolic calculations and Lean.

Our current hero case is the candidate `cb60fd71d1d7`, whose nearest benchmark
is the Rogers-Ramanujan continued fraction at `q^3`. Symbolically, the candidate
agrees with `RR(q^3)` through `q^9` and first diverges at `q^12`. In the reduced
variable `t = q^3`, its reciprocal takes the clean form

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n).
```

We do not claim a new continued-fraction identity. Instead, we report a stronger
and more honest result: the candidate survives several rings of nearby
explanations. Current exact or bounded eliminations include direct
constant-parameter matches to nearby published families, simple odd/even
contractions, short constrained Bauer-Muir chains, and low-degree algebraic
relations against the nearest Rogers-Ramanujan benchmark object. On the proof
side, Lean formalizes finite convergent recurrences, exact convergent-factor
reduction, and local obstruction lemmas around the hero case.

The main contribution is therefore methodological: a reproducible workflow for
moving from heuristic discovery to conservative public claims, with a clear line
between theorem-grade statements and bounded-search evidence.

## Why This Is Worth Sharing

- The candidate is structured, not random numerical noise.
- The project separates `known`, `fixture`, and `review` cases conservatively.
- The public claim is modest but strong: we have a good unexplained object plus
  exact local structure and multiple eliminated nearby origins.
- This is a useful pattern for AI-assisted mathematics: search broadly, filter
  aggressively, and publish only what survives exact audit.

## The Core Object

The current hero-case template is

```text
1 / (1 + (q^3 + q^6)/(1 + q^3 + (q^6 + q^12)/(1 + q^6 + ...))).
```

Its nearest built-in benchmark is `rogers_ramanujan_q3_normalized`.

Low-order symbolic comparison:

- matches `RR(q^3)` through `q^9`
- first divergence from `RR(q^3)` occurs at `q^12`
- remains materially closer to `RR(q^3)` than to the cubic benchmark in the
  first visible symbolic orders

## What We Can Say Exactly

In reduced variable `t = q^3`, the reciprocal object is

```text
C(t) = 1 + K_{n>=1} (t^n + t^(2n)) / (1 + t^n).
```

Current exact machine-checkable structure includes:

1. finite-truncation convergent recurrences for generalized continued fractions
2. exact common-factor reduction of the hero-case convergents
3. an equivalence witness reconstructing the original reciprocal from the
   reduced-by-factor object
4. coefficient-level local obstruction lemmas against nearby source families

## What We Are Not Claiming

- not a proof of novelty
- not a final source identity
- not a publishable "new formula" announcement yet

The honest current status is:

```text
high-value unexplained case
+ exact local structure
+ stronger-than-usual elimination evidence
```

## Links

- Project site: `https://peter941221.github.io/ramanujan-discovery-lab/`
- Repo: `https://github.com/peter941221/ramanujan-discovery-lab`
- Hero summary: `CB60FD71D1D7_PUBLIC_SUMMARY.md`
- Hero transform audit: `CB60FD71D1D7_TRANSFORM_AUDIT.md`
- Hero identification note: `CB60FD71D1D7_IDENTIFICATION_NOTE.md`
- Hero case study: `CB60FD71D1D7_CASE_STUDY.md`

## Suggested CTA

If you work on continued fractions, modular objects, or formalized mathematics,
the most useful feedback is not "is this novel?" in the abstract, but one of:

1. a source family that naturally produces `1 + q^n` partial denominators with a
   two-term numerator pattern
2. a functional equation or modular parametrization that this object should satisfy
3. a better transform/contraction lens than the current nearby-family audit

