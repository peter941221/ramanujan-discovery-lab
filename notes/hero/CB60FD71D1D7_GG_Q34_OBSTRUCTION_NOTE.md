# `cb60fd71d1d7` GG Quotient-Coordinate Obstruction Note

Status date: `2026-03-19`

## Purpose

Compress the current tail-family-first `GG` quotient-coordinate evidence into a
small reusable obstruction layer.

The exact local question is:

```text
can a sampled tail-family object Y
match the Chan--Huang quotient-coordinate lanes
Q_3 = GG(t^3) / GG(t)
Q_4 = GG(t^4) / GG(t)
```

where `Y` ranges over the exact tail-family ladder

```text
U_t2, U_t2_g1, U_t2_g2, U_t2_g3,
U_t3, U_t3_g1, U_t3_g2, U_t3_g3,
U_t4, U_t4_g1, U_t4_g2, U_t4_g3
```

from `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`.

## Exact Outcome

- Sample count: `12`
- Exact Chan--Huang quotient-coordinate hits on `Q_3`: `0`
- Exact Chan--Huang quotient-coordinate hits on `Q_4`: `0`

Instead of only recording "no hit", the current pipeline now records the first
visible residual term where each exact quotient-coordinate template fails.

## Obstruction Witness Table

| Tail object | `Q_3` first failure | `Q_4` first failure |
| --- | --- | --- |
| `U_t2` | `t^1` with coefficient `2` | `t^1` with coefficient `3` |
| `U_t2_g1` | `t^1` with coefficient `4` | `t^1` with coefficient `6` |
| `U_t2_g2` | `t^2` with coefficient `-4` | `t^2` with coefficient `-6` |
| `U_t2_g3` | `t^1` with coefficient `6` | `t^1` with coefficient `9` |
| `U_t3` | `t^1` with coefficient `2` | `t^1` with coefficient `3` |
| `U_t3_g1` | `t^1` with coefficient `2` | `t^1` with coefficient `3` |
| `U_t3_g2` | `t^2` with coefficient `-2` | `t^2` with coefficient `-3` |
| `U_t3_g3` | `t^1` with coefficient `4` | `t^1` with coefficient `6` |
| `U_t4` | `t^1` with coefficient `2` | `t^1` with coefficient `3` |
| `U_t4_g1` | `t^1` with coefficient `2` | `t^1` with coefficient `3` |
| `U_t4_g2` | `t^2` with coefficient `-2` | `t^2` with coefficient `-3` |
| `U_t4_g3` | `t^1` with coefficient `4` | `t^1` with coefficient `6` |

## Pattern Compression

The witness table is not random noise.  It compresses into a few small
obstruction classes:

```text
Class A:  Q_3 -> (t^1,  2),  Q_4 -> (t^1,  3)
Class B:  Q_3 -> (t^1,  4),  Q_4 -> (t^1,  6)
Class C:  Q_3 -> (t^2, -2),  Q_4 -> (t^2, -3)
Class D:  Q_3 -> (t^2, -4),  Q_4 -> (t^2, -6)
Class E:  Q_3 -> (t^1,  6),  Q_4 -> (t^1,  9)
```

So the current local evidence says:

- the tail-family ladder does not merely miss in an unstructured way
- it misses in a small number of repeated quotient-coordinate obstruction
  shapes

## Renormalized Compression

The five raw classes compress one step further.

If `R_3(Y)` and `R_4(Y)` denote the exact residual series for the `Q_3` and
`Q_4` Chan--Huang quotient-coordinate templates on a sampled tail object `Y`,
then every current class fits the same leading-order shape

```text
R_3(Y) = 2*rho(Y)*t^m(Y) + ...
R_4(Y) = 3*rho(Y)*t^m(Y) + ...
```

with a shared failure order `m(Y)` and a shared obstruction scalar `rho(Y)`.

For the current `12` samples, this tail-derived coordinate

```text
Theta_tail(Y) = (m(Y), rho(Y))
```

only takes the five values

```text
(1,  1), (1,  2), (2, -1), (2, -2), (1,  3)
```

Analogy:

- the old table listed the jammed teeth on two separate keys
- this renormalized view says both keys are really jamming on one shared tooth
  height `rho`, seen through two fixed magnifications `2` and `3`

## Tail-Derived Quotient Coordinate Verdict

This means the obstruction classes **do** unify at the first visible
residual layer.

So there is a meaningful prototype coordinate:

```text
Theta_tail(Y) = (leading failure order, shared obstruction scalar)
```

But the current evidence still does **not** justify treating `Theta_tail`
itself as the final named source coordinate.

Working reason:

- the leading `3:2` ratio is stable
- but the remaining residual after canceling that shared leading term is still
  sample-dependent rather than collapsing to one obvious universal next object

So Step 3 succeeds only as a **diagnostic compression**, not yet as a final
recognition coordinate.

## Weighted Coordinate Follow-Up

We now tested the first obvious weighted quotient coordinate suggested by the
shared leading `3:2` obstruction ratio:

```text
W_34 = Q_3^3 / Q_4^2
log W_34 = 3*log(Q_3) - 2*log(Q_4)
```

Current outcome on the full `12`-sample tail ladder:

- no sampled tail-family object matches `W_34` exactly
- no sampled tail-family object matches `log W_34` at the checked truncation
- no degree-`<= 2` polynomial relation in `(F, W_34, t)` was found
- no one-coordinate fractional-linear closure in `W_34` was found
- the first multiplicative correction `F / W_34` also keeps a visible first
  gap on every sampled object, so the weighted lane still does not collapse at
  the first modular-tail correction
- after normalizing that first gap into the sample-level object `G_W34`, the
  checked small eta-quotient, modular-unit / eta, and one-core `RR/GG`
  source-family eta-correction boxes still report `0` hits on the current
  sampled ladder
- after stripping one more visible gap into the sample-level object `G2_W34`,
  the same checked small eta-quotient, modular-unit / eta, and one-core
  `RR/GG` source-family eta-correction boxes still report `0` hits on the
  current sampled ladder

So the first weighted quotient coordinate and its first normalized correction
are useful as **probes**, but they still do not collapse the tail-family lane
into a named source object.

## Direction Consequence

The next source-driven move should therefore be:

- keep the exact `Q_3 / Q_4` obstruction layer
- keep `Theta_tail` as the tail-family diagnostic coordinate
- record that the first weighted follow-up `W_34 = Q_3^3 / Q_4^2`, its first
  correction `F / W_34`, the normalized residual `G_W34`, and the deeper
  residual `G2_W34` all still miss in the current small recognition boxes
- shift the next positive-recognition attempt to a deeper `Weber /
  modular-curve` parameter that respects the same `3:2` weighting, instead of
  widening anonymous `GG` prefix scans again

Natural next parameter shapes to inspect are weighted quotient combinations
such as a `Q_3^3 / Q_4^2`-type modular coordinate or an equivalent
modular-function ratio from the Chan--Huang / Weber orbit.

## Why This Matters

Analogy:

- the old output was "the key does not open the lock"
- the new output is "the key always jams at one of five teeth-patterns"

That is a better theorem target because it suggests the next proof object
should package exact residual-coefficient obstructions, not just another
bounded prefix scan.

## Immediate Theorem Target

The next exact lane should aim to prove a statement of the form:

```text
for each sampled tail-family object Y in the current ladder,
the Chan--Huang exact quotient-coordinate identities for Q_3 and Q_4
already fail at the first listed residual coefficient
```

This is still weaker than a final source-identity theorem, but it is stronger
and more reusable than the current plain `0`-hit summaries.

For the weighted follow-up, the current theorem-shaped companion target is
slightly different:

```text
record the first correction label F / W_34,
record the normalized follow-up label G_W34,
record the deeper follow-up label G2_W34,
and package the current no-hit correction verdicts
before committing to a final positive modular coordinate
```

## Lean Waypoint

The proof workspace now has a first landing spot for this obstruction layer:

- file: `proofs/Proofs/HeroCaseGGQuotientCoordinateObstruction.lean`
- current payload:
  - `TailSample`
  - `witnessFor`
  - `classFor`
  - `sampleWitness_true`
  - `q3FirstFailure_true`
  - `q4FirstFailure_true`
  - `pairedFirstFailure_true`
  - `tailDerivedCoordinate_true`
  - `sharedLeadingObstruction_true`
- `everyTailSampleHasFirstFailure_true`
- `everyTailSampleHasSharedLeadingObstruction_true`

The proof workspace now also has a second landing spot for the weighted
follow-up layer:

- file: `proofs/Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`
- current payload:
  - `weightedCoordinateLabel`
  - `weightedCorrectionExpression`
  - `normalizedWeightedCorrectionLabel`
  - `normalizedWeightedCorrectionExpression`
  - `secondNormalizedWeightedCorrectionLabel`
  - `secondNormalizedWeightedCorrectionExpression`
  - `currentWaypoint`
  - `currentWaypoint_true`

Interpretation:

- this is not yet the final modular-equation proof
- it is the first sample-indexed theorem layer saying:

```text
for every currently sampled tail-family object,
the exact Chan--Huang Q_3 / Q_4 lane already fails at the recorded first
residual term, and that leading failure compresses into one shared
tail-derived obstruction scalar
```

## Companion Files

- `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`
- `notes/hero/CB60FD71D1D7_IDENTIFICATION_NOTE.md`
- `notes/hero/CB60FD71D1D7_DIRECTION_REVIEW.md`
- `notes/hero/CB60FD71D1D7_AWARD_TRACK_EXECUTION_PLAN.md`
- `proofs/Proofs/HeroCaseGGQuotientCoordinateObstruction.lean`
- `proofs/Proofs/HeroCaseGGWeightedCorrectionWaypoint.lean`
