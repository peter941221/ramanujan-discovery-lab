# `cb60fd71d1d7` Weber Companion Obstruction Note

Status date: `2026-03-20`

## Purpose

Compress the first companion-coordinate follow-up inside the
Weber-Schlafli / Morton lane into a small reusable obstruction layer.

The exact local question is:

```text
if Y is one sampled tail-family object,
and P_ws = (1/Y - Y) / 2 is the first Weber-Schlafli coordinate,
does the companion coordinate B = b(4τ)
explain the next source-faithful layer?
```

The `Y` objects are the same exact tail-family ladder:

```text
U_t2, U_t2_g1, U_t2_g2, U_t2_g3,
U_t3, U_t3_g1, U_t3_g2, U_t3_g3,
U_t4, U_t4_g1, U_t4_g2, U_t4_g3
```

from `notes/hero/CB60FD71D1D7_TAIL_FAMILY_NOTE.md`.

## Source-Faithful Coordinate Definition

This note uses the Akkarapakam--Morton low-degree Weber relations:

```text
2*p(8τ) = 1/v(τ) - v(τ)
b(τ)^4 = p(τ)^8 + 16 p(τ)^4
b(τ/2)^2 = b(τ) + 4
```

The coded local coordinates are therefore:

```text
P_ws = (1/Y - Y) / 2
B_ws = sqrt(root4(P_ws^8 + 16*P_ws^4) + 4)
```

and the first exact companion templates are:

```text
B_ws^2 - B_ws,2 - 4 = 0
B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0
```

## Exact Outcome

- Sample count: `12`
- Exact hits for `B_ws^2 - B_ws,2 - 4 = 0`: `0`
- Exact hits for `B_ws,2^4 - P_ws^8 - 16*P_ws^4 = 0`: `0`

Unlike the earlier `P_ws` lane, the `B_ws` misses do not spread across a few
different failure classes.

They collapse immediately to one universal constant-term obstruction.

## Obstruction Witness Table

| Tail object | `B_ws^2 - B_ws,2 - 4` first failure | `B_ws,2^4 - P_ws^8 - 16*P_ws^4` first failure |
| --- | --- | --- |
| `U_t2` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t2_g1` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t2_g2` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t2_g3` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t3` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t3_g1` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t3_g2` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t3_g3` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t4` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t4_g1` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t4_g2` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |
| `U_t4_g3` | `t^0` with coefficient `-2` | `t^0` with coefficient `16` |

## Pattern Compression

The witness table collapses to one universal class:

```text
Companion Class Omega:
  B_ws^2 - B_ws,2 - 4 -> (t^0, -2)
  B_ws,2^4 - P_ws^8 - 16*P_ws^4 -> (t^0, 16)
```

This means the companion lane is not merely "still 0-hit."

It fails at the threshold, before any higher-order coefficient matching has a
chance to matter.

## Interpretation

Analogy:

- the `P_ws` lane was like opening the gearbox and checking whether the first
  moving parts line up
- the `B_ws` lane is like trying the companion gear and discovering that the
  teeth already clash at the axle center

So the local takeaway is:

- the direct `B = b(4τ)` companion closure does not explain the sampled
  tail-family ladder
- the failure is exact and uniform across all current samples
- if the project stays in the same named `GG/Weber` orbit, the next move
  should be another named Weber coordinate or a theorem-grade packaging of
  this universal constant-term obstruction, not a broader anonymous scan

## Current Handoff

This note upgrades the Phase-2 Weber story from:

```text
P_ws also misses
```

to:

```text
P_ws misses in a few repeated low-order classes
B_ws misses in one universal constant-term class
```

So the orbit is now narrower, even without a positive recognition hit.
