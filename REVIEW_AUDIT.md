# Review Candidate Audit

This note records a conservative first-pass audit of the current `review` candidates produced by `Ramanujan Discovery Lab`.

Audit date: `2026-03-13`

## Method

- Compare each `review` template to the nearest built-in benchmark at high precision.
- Deduplicate exploratory outputs at the structural-family level before writing public artifacts.
- Cross-check against broad primary-source literature on generalized Rogers-Ramanujan and Ramanujan-style `q`-continued fractions.
- Use public-claim discipline:
  - no exact hit in a quick search is **not** enough to call something new
  - a candidate that sits very close to a well-studied family is treated as family-adjacent until proven otherwise

## Broad literature context

- Rogers-Ramanujan and related continued fractions already admit broad generalizations and arithmetic frameworks:
  - Griffin, Ono, Warnaar, 2014: https://arxiv.org/abs/1401.7718
  - Ciolan, Neiss, 2015: https://arxiv.org/abs/1504.06482
  - Berndt, Rebak, 2025 survey: https://arxiv.org/abs/2512.19952
- Ramanujan-style continued fractions of higher order remain active:
  - Rajkhowa, Saikia, 2023: https://arxiv.org/abs/2305.14988
  - Aricheta, Guadalupe, 2024: https://arxiv.org/abs/2404.05756
- Machine-generated continued fraction conjectures are already being proved and generalized:
  - Yamamoto, 2024: https://arxiv.org/abs/2403.09729

## Candidate classification

| Candidate | Closest benchmark | Local assessment | Public claim status |
| --- | --- | --- | --- |
| `bef31ddceea8` | `rogers_ramanujan_q4_normalized` (7 digits) | Strongest member of an off-diagonal single-numerator ladder near the `q^4` Rogers-Ramanujan family. Most likely a family-adjacent generalized variant. | Do not claim as new |
| `cb60fd71d1d7` | `rogers_ramanujan_q3_normalized` (7 digits) | Structured denominator perturbation of the `q^3` Rogers-Ramanujan family. Interesting enough for a second pass, but still too close to a known family for a novelty claim. | Keep under review |
| `e42a0d5f2679` | `rogers_ramanujan_q4_normalized` (6 digits) | Second member of the same off-diagonal single-numerator ladder near the `q^4` Rogers-Ramanujan family. Best treated as family-adjacent rather than new. | Do not claim as new |
| `1125ffe48b3b` | `shifted_rr_fixture` (6 digits) | Family-adjacent perturbation mixing a `q^2` Rogers-Ramanujan numerator with a higher-order extra term and denominator shift. Worth a second pass, but not a publishable claim. | Keep under review |
| `9dbafd59364c` | `rogers_ramanujan_q4_normalized` (5 digits) | Weaker third member of the off-diagonal single-numerator ladder near the `q^4` Rogers-Ramanujan family. Stable, but not a novelty signal. | Do not claim as new |
| `e2cc74240b6f` | `rogers_ramanujan_q3_normalized` (5 digits) | Analogous off-diagonal single-numerator perturbation near the `q^3` Rogers-Ramanujan family. Family-adjacent and not ready for a public claim. | Do not claim as new |

## Result

Current conclusion:

- after benchmark-relative equivalence dedupe, the public release keeps six `review` candidates
- four of the six `review` candidates form a short off-diagonal single-numerator ladder near the Rogers-Ramanujan `q^3` and `q^4` families
- none of the six `review` candidates is ready for a public “new formula” claim
- `cb60fd71d1d7` and `1125ffe48b3b` remain the only two that currently look worth a second pass
- the next serious step is still not publicity; it is a stronger equivalence reduction plus targeted literature search around the nearest classical families

## Notes

- On `2026-03-13`, quick targeted web searches for the exact exponent-pattern strings behind the strongest representatives of these families did not produce a clear exact match.
- That absence is only weak evidence, because the surrounding literature on generalized Ramanujan-style continued fractions is already broad.
