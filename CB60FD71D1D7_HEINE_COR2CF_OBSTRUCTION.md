# `cb60fd71d1d7` Heine `cor2cf` Obstruction

Status date: `2026-03-16`

## Goal

Record an exact low-stage obstruction for the Heine-derived `cor2cf` lane that
shows up naturally in the Rogers-Ramanujan neighborhood audit.

This note does **not** claim literature closure.

It only isolates one nearby source mechanism and shows that its relevant
1-step / 2-step odd-even contraction branches do not produce the hero-case
target.

## Relevant Specialization

In the transform audit, matching the hero target's initial term `1` forces the
nearby `cor2cf` branch to the `a = 0` specialization.

In reduced variable `t = q^3`, that source begins:

```text
1 + (lambda*t)/1 + (b*t + lambda*t^2)/1 + (lambda*t^3)/1 + (b*t^2 + lambda*t^4)/1 + ...
```

Equivalently, as coefficient data:

```text
b0 = 1
a1 = lambda*t
b1 = 1
a2 = b*t + lambda*t^2
b2 = 1
a3 = lambda*t^3
b3 = 1
a4 = b*t^2 + lambda*t^4
b4 = 1
```

The hero-case target reciprocal is

```text
H(t) = 1 + K_(n>=1) (t^n + t^(2n)) / (1 + t^n).
```

So the first target data are:

```text
b0 = 1
a1 = t + t^2
b1 = 1 + t
```

## One-Step Odd / Even Parts

### Odd part

The odd part of the specialized `cor2cf` source has initial term

```text
d0 = 1 + lambda*t.
```

That is already incompatible with the target initial term `1`.

### Even part

The even part keeps initial term `1`, but its first numerator is

```text
a1_even = lambda*t
```

instead of the target

```text
t + t^2.
```

So the source is not a 1-step odd or even contraction of the hero target.

## Two-Step Branches

The only natural two-step follow-up comes from contracting the even part again.

### Odd of even

The odd part of the even part has initial term

```text
(1 + b*t + lambda*t + lambda*t^2) / (1 + b*t + lambda*t^2).
```

This is not identically `1`, so that branch fails before the first nontrivial
numerator.

### Even of even

The even part of the even part keeps initial term `1`, but its first numerator
is

```text
lambda*t*(1 + b*t^2 + lambda*t^3 + lambda*t^4).
```

Its `t^2` coefficient is `0`, so it cannot equal the target first numerator

```text
t + t^2.
```

## Conclusion

The relevant Heine-derived `cor2cf` specialization is **not**

- a 1-step odd contraction
- a 1-step even contraction
- the odd-of-even two-step branch
- the even-of-even two-step branch

for the hero-case target.

This strengthens the earlier transform audit by turning that local argument into
an exact symbolic artifact rather than a prose-only elimination.

## Validation Hook

These formulas are now reproduced by the code path:

- `src/ramanujan_discovery/research.py`

and regression-checked in:

- `tests/test_research.py`

Useful commands:

```powershell
$env:PYTHONPATH='src'; pytest tests/test_research.py -q
$env:PYTHONPATH='src'; python -m ramanujan_discovery research --in results/verified.jsonl --candidate-id cb60fd71d1d7 --depth 40 --series-order 151 --out CB60FD71D1D7_RESEARCH_NOTE.md
```
