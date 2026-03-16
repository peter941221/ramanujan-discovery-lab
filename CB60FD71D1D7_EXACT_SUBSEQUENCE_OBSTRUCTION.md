# `cb60fd71d1d7` Exact Subsequence Obstruction

Status date: `2026-03-16`

## Goal

Upgrade the old bounded arithmetic-subsequence scan

```text
stride <= 4, exact rational sample points, stages <= 3
```

to a stronger exact argument for the two nearby source families that matter most:

- the reduced Rogers-Ramanujan reciprocal
- the reduced cubic reciprocal

The target reciprocal in `t = q^3` is

```text
H(t) = 1 + K_(n>=1) (t^n + t^(2n)) / (1 + t^n).
```

Its first two convergents are:

```text
H_0(t) = 1
H_1(t) = 1 + t
```

So any arithmetic subsequence contraction of a source continued fraction that
matches `H(t)` must satisfy:

1. the subsequence start value equals `1`
2. the next subsequence value equals `1 + t`

That is already enough to rule out the nearby RR and cubic sources exactly.

## RR Reciprocal

Let `P_n / Q_n` be the convergents of the reduced RR reciprocal

```text
R(t) = 1 + K_(n>=1) t^n / 1.
```

The continuants satisfy

```text
P_(n+2) = P_(n+1) + t^(n+2) P_n
Q_(n+2) = Q_(n+1) + t^(n+2) Q_n.
```

Define

```text
Delta_n = P_n - Q_n.
```

Then `Delta_(n+2) = Delta_(n+1) + t^(n+2) Delta_n`, with

```text
Delta_0 = 0
Delta_1 = t.
```

Because the extra term has degree at least `2`, it never changes the linear
coefficient. Therefore:

```text
[t] Delta_n = 1    for every n >= 1.
```

Hence `P_n != Q_n` for every `n >= 1`, so no positive-offset subsequence can
start at `H_0(t) = 1`.

Now define

```text
E_n = P_n - (1 + t) Q_n.
```

The same continuant recurrence gives

```text
E_(n+2) = E_(n+1) + t^(n+2) E_n.
```

Direct calculation at the first relevant stages gives

```text
[t^3] E_2 = -1
[t^3] E_3 = -1.
```

For all later stages, the added term has degree at least `4`, so it cannot
change the `t^3` coefficient. Therefore:

```text
[t^3] E_n = -1    for every n >= 2.
```

Hence `P_n != (1 + t) Q_n` for every `n >= 2`, so no zero-offset subsequence
with stride `>= 2` can hit `H_1(t) = 1 + t`.

Conclusion:

```text
No arithmetic subsequence contraction of the reduced RR reciprocal
can equal the hero-case target.
```

## Cubic Reciprocal

Let `P_n / Q_n` be the convergents of the reduced cubic reciprocal

```text
C(t) = 1 + K_(n>=1) (t^n + t^(2n)) / 1.
```

Its continuants satisfy

```text
P_(n+2) = P_(n+1) + (t^(n+2) + t^(2n+4)) P_n
Q_(n+2) = Q_(n+1) + (t^(n+2) + t^(2n+4)) Q_n.
```

Set again

```text
Delta_n = P_n - Q_n.
```

Then

```text
Delta_(n+2) = Delta_(n+1) + (t^(n+2) + t^(2n+4)) Delta_n
```

with

```text
Delta_0 = 0
Delta_1 = t + t^2.
```

The added term has degree at least `2`, so it never changes the linear term.
Therefore:

```text
[t] Delta_n = 1    for every n >= 1.
```

So again, no positive-offset subsequence can start at `1`.

For the first-step obstruction, define

```text
E_n = P_n - (1 + t) Q_n.
```

Then

```text
E_(n+2) = E_(n+1) + (t^(n+2) + t^(2n+4)) E_n.
```

Direct calculation gives

```text
[t^3] E_2 = -1
[t^3] E_3 = -1.
```

All later correction terms have degree at least `4`, so the `t^3` coefficient
stays fixed:

```text
[t^3] E_n = -1    for every n >= 2.
```

Hence no zero-offset subsequence with stride `>= 2` can hit `1 + t`.

Conclusion:

```text
No arithmetic subsequence contraction of the reduced cubic reciprocal
can equal the hero-case target.
```

## What This Strengthens

This is stronger than the old bounded scan because it removes the

- stride bound
- sample-point bound
- stage-depth bound

for these two specific nearby source families.

## Validation Hook

This note is backed by symbolic regression checks in:

- `tests/test_research.py`

Current regression commands:

```powershell
$env:PYTHONPATH='src'; pytest tests/test_research.py -q
Set-Location proofs; lake build
```
