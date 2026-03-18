# `cb60fd71d1d7` Bibliography Matrix

Status date: `2026-03-18`

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
| Hirschhorn special-issue exposition, `S(q)`, `GG(q)` neighborhood | Bowman, Mc Laughlin, Wyshinski, 2005 special-issue PDF | https://www.wcupa.edu/sciences-mathematics/mathematics/documents/hirschhornnov08.pdf | Older `1 + q^n` denominator examples and an alternate presentation of the page-43 neighborhood; both `S(q)` and `GG(q)` are now encoded as normalized project benchmarks with verified product formulas | No exact stage-pattern hit found; neither the named mixed-family boxes, the per-family powered ladders over `GG/GG2/GG3/GG4` and `S/S2/S3/S4`, nor the explicit direct / reciprocal / quotient template box in those same ladders explains the hero case |
| Classical Gordon / Hirschhorn orbit | Basil Gordon, *Some continued fractions of the Rogers-Ramanujan type* (Duke Math. J. 32 (1965), 741-748); M. D. Hirschhorn, *A continued fraction* (Duke Math. J. 41 (1974), 27-33); M. D. Hirschhorn, *A continued fraction of Ramanujan* (J. Austral. Math. Soc. Ser. A 29 (1980), 80-86) | Gordon 1965 DOI: https://doi.org/10.1215/S0012-7094-65-03278-3 ; Hirschhorn 1974 DOI: https://doi.org/10.1215/S0012-7094-74-04104-0 ; Hirschhorn 1980 DOI/Cambridge: https://doi.org/10.1017/S1446788700020954 | The older exact citation spine behind the `GG` / `S` neighborhood and related `1 + q^n` denominator continued fractions | Primary-link hygiene is now stronger, but no exact source hit or transform identity has been extracted yet |
| Modern Hirschhorn structural lane | Bhatnagar, Ismail, "Orthogonal polynomials associated with a continued fraction of Hirschhorn" | https://arxiv.org/abs/1901.09985 | Direct structural analysis of Hirschhorn's continued fraction, including convergents, generating functions, orthogonality, and Stieltjes-transform data; explicitly notes RR and Ramanujan-generalization special cases | No exact hero-pattern source identity extracted; widens the Hirschhorn-side recognition burden beyond base product formulas |
| GG power / substitution identity lane | Vasuki, Srivatsa Kumar, "Certain identities for Ramanujan-Göllnitz-Gordon continued fraction" | https://doi.org/10.1016/j.cam.2005.03.038 | Explicit identities relating `H(q)` to `H(q^3)`, `H(q^5)`, `H(q^7)`, and `H(q^11)` inside the GG orbit | No exact coefficient-level hero match extracted; confirms the GG orbit has nontrivial power/substitution relations beyond the current simple basis scans |
| GG base modular-equation lane | Chan, Huang, "On the Ramanujan-Göllnitz-Gordon continued fraction" | https://mrc.sdu.edu.cn/ziliao/8.pdf | Treats `H(q)` as a modular object, derives identities linking `H(q)` with `H(-q)` and `H(q^2)`, and explains how modular equations produce relations between `H(q)` and `H(q^n)`, with explicit emphasis on the `q^3` and `q^4` lanes | Strong evidence that the next recognition pass should encode structured `GG(q), GG(-q), GG(q^2), GG(q^3), GG(q^4)` templates rather than only wider anonymous basis boxes |
| GG modular-unit / odd-prime modular-equation lane | Cho, Koo, Park, "Arithmetic of the Ramanujan-Göllnitz-Gordon continued fraction" | https://doi.org/10.1016/j.jnt.2008.09.018 | Extends Chan--Huang / Vasuki--Srivatsa Kumar style modular equations to all odd primes and treats the GG object as a modular unit with arithmetic structure | Strong evidence that the current local `GG` recognition story is still structurally incomplete even after the powered-ladder and quotient-template passes |
| GG coefficient / evaluation property lane | Yuttanan, "New properties for the Ramanujan-Göllnitz-Gordon continued fraction" | https://eudml.org/doc/278973 | Additional exact properties of the GG continued fraction, including coefficient and evaluation behavior beyond the earlier product-identity lane | Widens the burden against declaring the GG orbit exhausted; supports a literature-driven modular / structural recognition pass before widening generic search boxes again |
| Bauer-Muir transforms, `gcf2`, `gcf3`, `gcf4`, transformed RR-type fractions | Lee, Mc Laughlin, Sohn, "Bauer-Muir transformations and generalizations of Rogers-Ramanujan type continued fractions" | https://arxiv.org/abs/1906.11991 | Main transform lane currently audited against the hero case | No direct constant-parameter or tiny-chain hit found |
| Orthogonal-polynomial / J-fraction viewpoint | Ismail, Stanton, "Orthogonal Polynomials and q-Continued Fractions" | http://www-users.math.umn.edu/~isman/papers/ramanujan.pdf | Structural recognition lane for RR-type fractions as continued fractions from recurrences | Relevant framework, but no direct match extracted yet |
| Weber / Gollnitz-Gordon neighborhood | Adiga, Kim, et al., "Ramanujan-Weber Class Invariants G_n and g_n" | https://www.filomat.org/index.php/filomat/article/download/4026/4026/31-13-2-4026.pdf | Class-invariant and Gollnitz-Gordon-adjacent continued fractions | No direct coefficient-level match extracted in current pass |
| Notebook-wide `q`-continued-fraction survey lane | Bhatnagar, "Ramanujan's `q`-continued fractions" | https://arxiv.org/abs/2208.12656 | Reorganizes Ramanujan's notebook `q`-continued fractions and summarizes later work across multiple continued-fraction orders and families | No exact hero-pattern hit extracted; widens the notebook-level novelty burden |
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
- a modern structural Hirschhorn paper and a direct GG power-relation paper
- a base GG modular-equation paper, an odd-prime modular-equation / modular-unit paper, and a later GG property paper
- a first named source-family correction scan over `RR`, `cubic`, `GG`, and `S`
- a per-family powered source-family pass over `RR/RR2/RR3/RR4`,
  `cubic/cubic2/cubic3/cubic4`, `GG/GG2/GG3/GG4/GG5/GG7/GG11`, and `S/S2/S3/S4`
- a per-family low-degree polynomial scan over those same powered family ladders
- a per-family two-layer fractional-linear scan over those same powered family ladders
- a within-family quotient-ladder scan `Qk = Tk / T1` over those same powered family ladders
- a quotient-ladder two-layer fractional-linear scan over those same powered family ladders
- a mixed quotient-basis scan combining `T1` with `Qk = Tk / T1` over those same powered family ladders
- a mixed quotient-basis two-layer fractional-linear scan over those same powered family ladders
- an explicit `GG` / `S` transform-template pass over direct objects,
  reciprocals, and pairwise quotients inside those literature-family ladders
- a notebook-level survey source plus a verified post-2020 RR bibliography slice covering modular products, modular functions, modular equations, and conjecture/proof papers

What is **not** covered well enough yet:

- a broad systematic bibliography across all modern Ramanujan-style `q`-continued fractions
- deeper nontrivial substitution / contraction / transform chains outside the current nearby families
- reciprocation / quotient / transform templates in the Gordon/Hirschhorn orbit beyond the current powered family ladders, quotient ladders, quotient-ladder two-layer boxes, mixed quotient-basis boxes, mixed quotient-basis two-layer boxes, and two-layer within-family corrections
- an explicit modular-equation-style recognition lane for the `GG` orbit that treats `GG(q)`, `GG(-q)`, `GG(q^2)`, `GG(q^3)`, and `GG(q^4)` as first-class literature objects
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
