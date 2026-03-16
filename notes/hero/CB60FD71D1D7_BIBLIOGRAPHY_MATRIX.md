# `cb60fd71d1d7` Bibliography Matrix

Status date: `2026-03-16`

## Purpose

This matrix is not a claim that the literature is closed.

It is a working coverage map for the hero case:

- what source families have been checked
- which primary sources anchor those checks
- whether the current exact pattern

```text
(q^(3n) + q^(6n)) / (1 + q^(3n))
```

has been located there

## Matrix

| Family / lane | Primary source | Link | What it covers for this project | Current outcome |
| --- | --- | --- | --- | --- |
| Page-43 ratio family, `H`, `H_1`, nearby RR/cubic neighborhoods | Bowman, Mc Laughlin, Wyshinski, "Generalized Rogers-Ramanujan continued fractions, Ramanujan's continued fraction, and Hirschhorn's continued fraction" | https://arxiv.org/abs/1901.00584 | Main nearby source families, Ramanujan notebook page-43 material, denominator-perturbed RR and cubic-adjacent lanes | No exact hit for the hero pattern |
| Hirschhorn special-issue exposition, `S(q)`, `GG(q)` neighborhood | Bowman, Mc Laughlin, Wyshinski, 2005 special-issue PDF | https://www.wcupa.edu/sciences-mathematics/mathematics/documents/hirschhornnov08.pdf | Older `1 + q^n` denominator examples and an alternate presentation of the page-43 neighborhood; both `S(q)` and `GG(q)` are now encoded as normalized project benchmarks with verified product formulas | No exact stage-pattern hit found; the named source-family multiplicative scan over `RR`, `cubic`, `GG`, and `S` still has no hero-case hit |
| Bauer-Muir transforms, `gcf2`, `gcf3`, `gcf4`, transformed RR-type fractions | Lee, Mc Laughlin, Sohn, "Bauer-Muir transformations and generalizations of Rogers-Ramanujan type continued fractions" | https://arxiv.org/abs/1906.11991 | Main transform lane currently audited against the hero case | No direct constant-parameter or tiny-chain hit found |
| Orthogonal-polynomial / J-fraction viewpoint | Ismail, Stanton, "Orthogonal Polynomials and q-Continued Fractions" | http://www-users.math.umn.edu/~isman/papers/ramanujan.pdf | Structural recognition lane for RR-type fractions as continued fractions from recurrences | Relevant framework, but no direct match extracted yet |
| Weber / Gollnitz-Gordon neighborhood | Adiga, Kim, et al., "Ramanujan-Weber Class Invariants G_n and g_n" | https://www.filomat.org/index.php/filomat/article/download/4026/4026/31-13-2-4026.pdf | Class-invariant and Gollnitz-Gordon-adjacent continued fractions | No direct coefficient-level match extracted in current pass |
| Recent RR modular-product lane | Aricheta, Guadalupe, "A Remark on Modularity of Certain Products of the Rogers-Ramanujan Continued Fraction" | https://arxiv.org/abs/2405.06678 | Modern modularity constraints for products built from the RR continued fraction | No exact hero-pattern source identity extracted; widens the active RR-product landscape |
| Recent RR order-10 modular-function lane | Aricheta, Guadalupe, "On Certain Order 10 Modular Functions Involving the Rogers-Ramanujan Continued Fraction" | https://arxiv.org/abs/2404.05756 | Recent modular-function identities around the RR orbit | No exact hit recovered in the scanned formulas; confirms nearby RR structure remains active |
| Recent RR modular-equation lane | Akkarapakam, Morton, "A remark on modular equations involving the Rogers-Ramanujan continued fraction via 5-dissections" | https://arxiv.org/abs/2410.14149 | Modern modular-equation identities for the RR continued fraction | No direct coefficient-level match extracted in this pass |
| Recent machine-discovery RR lane | Yamamoto, "Proof and Generalization of Conjectures of Ramanujan Machine" | https://arxiv.org/abs/2403.09729 | Modern discovery/proof workflow on RR-adjacent conjectural identities | Relevant for novelty burden and recognition posture, but no exact hero source identity |
| Recent RR-adjacent sign / product behavior | Baruah, Sarma, "Sign Patterns and Congruences for Certain Infinite Products involving the Rogers-Ramanujan Continued Fraction" | https://arxiv.org/abs/2503.08517 | Modern RR-adjacent infinite-product behavior around the RR continued fraction | Not an exact source match, but another reminder that RR-adjacent work is current |
| Recent RR-adjacent conjecture follow-up | Ghoshal, Jana, "Affirmation of Certain Conjectures on the Rogers-Ramanujan Continued Fraction" | https://arxiv.org/abs/2503.17950 | Recent follow-up on RR continued-fraction conjectures | Not a source match; used only as evidence that the area remains active |
| Survey of RR continued fraction properties, values, identities, generalizations | Berndt, Rebaka, "The Rogers-Ramanujan continued fraction" | https://arxiv.org/abs/2512.19952 | Recent survey-level snapshot of the RR continued-fraction literature | No exact hit seen, but this widens the novelty burden because the literature is still active |

## Coverage Read

What is reasonably covered now:

- the immediate page-43 / RR / cubic / Bauer-Muir neighborhood
- several older `1 + q^n` denominator examples
- a first named source-family correction scan over `RR`, `cubic`, `GG`, and `S`
- a verified post-2020 RR bibliography slice covering modular products, modular functions, modular equations, and conjecture/proof papers

What is **not** covered well enough yet:

- a broad systematic bibliography across all modern Ramanujan-style `q`-continued fractions
- deeper nontrivial substitution / contraction / transform chains outside the current nearby families
- parameterized Gordon/Hirschhorn transform families beyond the base benchmark objects
- exact title-normalized coverage for every recent GG / Weber / RR orbit paper that might still sit just outside the current verified list

## Decision Use

This matrix is meant to support one narrow decision:

```text
Can the project already say "new discovery"?
```

Current answer:

```text
No.
```

Reason:

- the nearby-family coverage is now meaningful
- the post-2020 RR slice is better grounded because the rows above use verified titles/links
- but the global literature closure is still not strong enough

## Companion Notes

- `CB60FD71D1D7_LITERATURE_LOG.md`
- `CB60FD71D1D7_EXACT_SUBSEQUENCE_OBSTRUCTION.md`
- `CB60FD71D1D7_NOVELTY_GATE.md`
