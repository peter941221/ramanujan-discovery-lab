# `cb60fd71d1d7` Bibliography Matrix

Status date: `2026-03-20`

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
| GG base modular-equation lane | Chan, Huang, "On the Ramanujan-Göllnitz-Gordon continued fraction" | https://mrc.sdu.edu.cn/ziliao/8.pdf | Treats `H(q)` as a modular object, derives identities linking `H(q)` with `H(-q)` and `H(q^2)`, and explains how modular equations produce relations between `H(q)` and `H(q^n)`, with explicit emphasis on the `q^3` and `q^4` lanes | The project now has a tail-family-first `GG` lane with explicit `GG(t), GG(-t), GG(t^2), GG(t^3), GG(t^4)` coordinates plus exact `Q_3 = GG(t^3)/GG(t)` and `Q_4 = GG(t^4)/GG(t)` checks; after those misses, the first weighted follow-up `Q_3^3 / Q_4^2` also stays negative, and even the later Weber return-bridge objects `Q_XK_ws` and `L_XK_ws` still give `0` hits when pushed back through the same named `GG` modular-equation box. The current exact witness quartets `(-9/2, -6 ; -3, -9/2)` for `Q_XK_ws` and `(593/10, 1186/15 ; 593/15, 593/10)` for `L_XK_ws` make that negative reading sharper: the return bridge is still not a hidden rename of the present Chan--Huang basis, so the next move should narrow toward a deeper modular-curve / Weber parameter rather than broader anonymous boxes |
| Ramanujan-Weber / classical-Weber naming bridge | Berndt, Chan, Zhang, "Ramanujan's class invariants, Kronecker's limit formula, and modular equations" | https://mrc.sdu.edu.cn/ziliao/10.pdf | Explicitly identifies Ramanujan-Weber `G_n` / `g_n` with classical Weber `f` / `f1`, giving the dictionary that lets the project read its current `g12_ws` / `p12_ws` shell in named classical-Weber terms when combined with Yui--Zagier's `f`, `f1`, `f2` trio | Does not identify the hero object directly, but it upgrades the current `g12_ws` / `p12_ws` / `G_f2_ws` shell from an anonymous algebra gadget to a source-backed classical Weber trio in project normalization |
| GG modular-unit / odd-prime modular-equation lane | Cho, Koo, Park, "Arithmetic of the Ramanujan-Göllnitz-Gordon continued fraction" | https://doi.org/10.1016/j.jnt.2008.09.018 | Extends Chan--Huang / Vasuki--Srivatsa Kumar style modular equations to all odd primes and treats the GG object as a modular unit with arithmetic structure | Still supports staying in the `GG` / Weber neighborhood, but after the `Q_3` / `Q_4` and `Q_3^3 / Q_4^2` tail-lane misses the best next step is now a smaller quotient-coordinate / modular-curve style lane, not a larger odd-prime prefix catalog |
| Weber singular-value anchor | Yui, Zagier, "On the singular values of Weber modular functions" | https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1090/S0025-5718-97-00854-5/fulltext.pdf | Classical Weber modular-function source spine behind later modular-curve and class-invariant coordinates | Does not identify the hero object directly, but strengthens the choice to move from raw GG quotient coordinates to Weber-flavored modular-function coordinates |
| GG coefficient / evaluation property lane | Yuttanan, "New properties for the Ramanujan-Göllnitz-Gordon continued fraction" | https://eudml.org/doc/278973 | Additional exact properties of the GG continued fraction, including coefficient and evaluation behavior beyond the earlier product-identity lane | Widens the burden against declaring the GG orbit exhausted; supports a literature-driven modular / structural recognition pass before widening generic search boxes again |
| Bauer-Muir transforms, `gcf2`, `gcf3`, `gcf4`, transformed RR-type fractions | Lee, Mc Laughlin, Sohn, "Bauer-Muir transformations and generalizations of Rogers-Ramanujan type continued fractions" | https://arxiv.org/abs/1906.11991 | Main transform lane currently audited against the hero case | No direct constant-parameter or tiny-chain hit found |
| Orthogonal-polynomial / J-fraction viewpoint | Ismail, Stanton, "Orthogonal Polynomials and q-Continued Fractions" | http://www-users.math.umn.edu/~isman/papers/ramanujan.pdf | Structural recognition lane for RR-type fractions as continued fractions from recurrences | Relevant framework, but no direct match extracted yet |
| Weber / Gollnitz-Gordon neighborhood | Adiga, Kim, et al., "Ramanujan-Weber Class Invariants G_n and g_n" | https://www.filomat.org/index.php/filomat/article/download/4026/4026/31-13-2-4026.pdf | Class-invariant and Gollnitz-Gordon-adjacent continued fractions | This lane now supports a source-faithful Ramanujan-Weber class-invariant compression: on true normalized `GG`, the derived coordinates recover the named objects `(t^2; t^4)_inf^12` and `(-t^2; t^4)_inf^12`; however, on the hero side the later return-bridge objects still do not collapse back to a named literature coordinate, so the pass currently yields a constructive hand-off rather than a final identification |
| Notebook-wide `q`-continued-fraction survey lane | Bhatnagar, "Ramanujan's `q`-continued fractions" | https://arxiv.org/abs/2208.12656 | Reorganizes Ramanujan's notebook `q`-continued fractions and summarizes later work across multiple continued-fraction orders and families | No exact hero-pattern hit extracted; widens the notebook-level novelty burden |
| Recent RR modular-product lane | Guadalupe, "Modularity of certain products of the Rogers-Ramanujan continued fraction" | https://arxiv.org/abs/2405.06678 | Recent eta-quotient / generalized-eta modularity constraints for products built from the RR continued fraction; published in `Ramanujan J. 68 (2025)` | No exact hero-pattern source identity extracted; widens the active RR-product landscape and supports the eta-recognition trunk |
| Recent RR order-10 modular-function lane | Aricheta, Guadalupe, "Ramanujan's continued fractions of order 10 as modular functions" | https://arxiv.org/abs/2404.05756 | Recent modular-function and modular-equation identities around the RR orbit; published in `J. Number Theory 278 (2026)` | No exact hit recovered in the scanned formulas; confirms nearby RR structure remains active |
| Recent RR modular-equation lane | Guadalupe, "A remark on modular equations involving Rogers-Ramanujan continued fraction via 5-dissections" | https://arxiv.org/abs/2410.14149 | Modern modular-equation identities for the RR continued fraction | No direct coefficient-level match extracted in this pass |
| Recent RR multi-power modular-identity lane | Deka Baruah, Talukdar, "Identities for the Rogers-Ramanujan Continued Fraction" | https://arxiv.org/abs/2410.17110 | New modular identities relating `R(q)`, `R(q^2)`, `R(q^3)`, `R(q^4)`, `R(q^5)`, `R(q^6)`, `R(q^12)`, and `R(q^20)` | Supports keeping explicit RR modular-equation coordinates active; no exact hero reduction has been extracted yet |
| Recent Ramanujan periodic-point lane | Akkarapakam, Morton, "Periodic points of algebraic functions related to a continued fraction of Ramanujan" | https://nyjm.albany.edu/j/2024/30-36.html | A Ramanujan continued fraction whose special values become periodic points of a fixed algebraic function; more importantly for this project, it gives a Weber-Schlafli coordinate change `2*p(8τ) = 1/v(τ) - v(τ)` and low-degree exact relations in the deeper coordinates `p` and `b` | The current Morton-inspired sample box, the first Weber-Schlafli `P` lane, and the direct companion `B` lane all still give `0` hits on the tail-family ladder; the later Weber class-invariant return bridge `Q_XK_ws`, `L_XK_ws` also stays outside the named `GG` modular-equation boxes, so this paper still anchors the right orbit while the current honest reading remains "derived return bridge," not a positive source identification |
| Recent machine-discovery RR lane | Yamamoto, "Proof and Generalization of Conjectures of Ramanujan Machine" | https://arxiv.org/abs/2403.09729 | Modern discovery/proof workflow on RR-adjacent conjectural identities | Relevant for novelty burden and recognition posture, but no exact hero source identity |
| Recent RR-adjacent sign / product behavior | Baruah, Sarma, "Sign Patterns and Congruences of certain infinite products involving the Rogers-Ramanujan continued fraction" | https://arxiv.org/abs/2503.08517 | Modern RR-adjacent infinite-product behavior around the RR continued fraction | Not an exact source match, but another reminder that RR-adjacent work is current |
| Recent RR-adjacent conjecture follow-up | Ghoshal, Jana, "Discussion on some conjectures regarding the periodicity of sign patterns of certain infinite products involving the Rogers-Ramanujan Continued Fractions" | https://arxiv.org/abs/2503.17950 | Recent follow-up on RR continued-fraction conjectures and sign-pattern periodicity | Not a source match; used only as evidence that the area remains active |
| Survey of RR continued fraction properties, values, identities, generalizations | Berndt, Rebaka, "The Rogers-Ramanujan continued fraction" | https://arxiv.org/abs/2512.19952 | Recent survey-level snapshot of the RR continued-fraction literature | No exact hit seen, but this widens the novelty burden because the literature is still active |

## Coverage Read

What is reasonably covered now:

- the immediate page-43 / RR / cubic / Bauer-Muir neighborhood
- several older `1 + q^n` denominator examples
- a modern structural Hirschhorn paper and a direct GG power-relation paper
- a base GG modular-equation paper, an odd-prime modular-equation / modular-unit paper, and a later GG property paper
- a Ramanujan-Weber / classical-Weber dictionary paper that upgrades the
  current `g12_ws / p12_ws / G_f2_ws` shell from anonymous notation to a
  source-backed classical Weber trio
- a classical Weber modular-function singular-value paper that supports the modular-curve side of the same orbit
- a recent RR multi-power modular-identity paper and a recent Ramanujan periodic-point / algebraic-function paper
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
- a tail-family-first `GG/Weber` modular-equation note over
  `U(t^2)`, `U(t^3)`, `U(t^4)` and their gap-normalized residuals through depth
  `3`
- an exact tail-family-first quotient-coordinate pass on
  `Q_3 = GG(t^3)/GG(t)` and `Q_4 = GG(t^4)/GG(t)` with obstruction witnesses
- a coded Weber-Schlafli coordinate lane together with its direct companion
  `P/B` obstruction package and a later Ramanujan-Weber class-invariant
  compression / return-bridge pass through `g12_ws`, `p12_ws`, `G_X_ws`,
  `Q_XR_ws`, `K_XR_ws`, `Q_XK_ws`, and `L_XK_ws`
- a literature-closure check showing that the later return-bridge objects
  `Q_XK_ws` and `L_XK_ws` still do not re-enter the named `GG`
  modular-equation coordinates, now with exact Chan--Huang first-failure
  witness quartets rather than only hit-count summaries
- a notebook-level survey source plus a verified post-2020 RR bibliography slice covering modular products, modular functions, modular equations, and conjecture/proof papers

What is **not** covered well enough yet:

- a broad systematic bibliography across all modern Ramanujan-style `q`-continued fractions
- deeper nontrivial substitution / contraction / transform chains outside the current nearby families
- reciprocation / quotient / transform templates in the Gordon/Hirschhorn orbit beyond the current powered family ladders, quotient ladders, quotient-ladder two-layer boxes, mixed quotient-basis boxes, mixed quotient-basis two-layer boxes, and two-layer within-family corrections
- a theorem-grade quotient-coordinate elimination layer for the `GG` orbit that upgrades the current exact `Q_3` / `Q_4` residual witnesses into reusable modular-curve style obstructions
- a deeper named Weber / modular-function coordinate that turns the current
  exact obstruction quartets for `Q_XK_ws`, `L_XK_ws` into either a positive
  identification or a more structural uniqueness theorem
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
