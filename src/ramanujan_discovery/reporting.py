from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ramanujan_discovery.storage import read_candidates


def _generated_at() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def build_report(input_path: str, output_path: str) -> None:
    records = read_candidates(input_path)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    generated_at = _generated_at()

    lines = [
        "# Ramanujan Discovery Report",
        "",
        f"Generated: {generated_at}",
        "",
        "## Summary",
        "",
        f"- Verified candidates: {len(records)}",
        "- Status values: `known`, `known_variant`, `fixture`, `review`",
        "",
        "| id | target | kind | digits | stability | status |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ]

    for record in records:
        lines.append(
            f"| {record.id} | {record.matched_target} | {record.benchmark_kind} | "
            f"{record.digits_agree} | {record.stability_score} | {record.novelty_status} |"
        )

    if not records:
        lines.extend(["", "No verified candidates were found."])

    for record in records:
        lines.extend(
            [
                "",
                f"## Candidate `{record.id}`",
                "",
                f"- Target: `{record.matched_target}`",
                f"- Benchmark kind: `{record.benchmark_kind}`",
                f"- Digits agree: `{record.digits_agree}`",
                f"- Stability score: `{record.stability_score}`",
                f"- Novelty status: `{record.novelty_status}`",
                f"- Closest benchmark: `{record.closest_benchmark}` ({record.closest_benchmark_digits} digits)",
                f"- Family bucket: `{record.family_bucket}`",
                f"- Equivalence key: `{record.equivalence_key}`",
                f"- Notes: {record.notes}",
                f"- Template signature: `{record.template.signature()}`",
                f"- Template LaTeX: `{record.template.latex()}`",
                f"- Sample q values: `{record.q_values}`",
                f"- Sample estimates: `{record.value_estimates}`",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_site(input_path: str, output_dir: str, title: str = "Ramanujan Discovery Lab") -> None:
    records = read_candidates(input_path)
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    generated_at = _generated_at()

    public_records = []
    for record in records:
        payload = record.to_dict()
        payload["template_signature"] = record.template.signature()
        payload["template_latex"] = record.template.latex()
        payload["complexity_score"] = record.template.complexity_score()
        public_records.append(payload)

    payload = {
        "title": title,
        "generated_at": generated_at,
        "records": public_records,
        "status_counts": {
            "known": sum(record.novelty_status == "known" for record in records),
            "known_variant": sum(record.novelty_status == "known_variant" for record in records),
            "fixture": sum(record.novelty_status == "fixture" for record in records),
            "review": sum(record.novelty_status == "review" for record in records),
        },
        "benchmark_families": [
            {
                "name": "Rogers-Ramanujan",
                "count": 4,
                "description": "Product-backed baselines at q, q^2, q^3, and q^4.",
            },
            {
                "name": "Ramanujan cubic",
                "count": 3,
                "description": "Independent cubic benchmarks at q, q^2, and q^3.",
            },
            {
                "name": "Internal fixtures",
                "count": 2,
                "description": "Regression-only shifted and denominator-perturbed control templates.",
            },
        ],
        "pipeline": [
            {
                "step": "discover",
                "detail": "Enumerate a bounded q-continued-fraction grid and keep numerically stable candidates.",
            },
            {
                "step": "verify",
                "detail": "Re-score survivors at higher precision against the benchmark catalog.",
            },
            {
                "step": "report",
                "detail": "Publish a markdown artifact with signatures, buckets, and benchmark alignment.",
            },
            {
                "step": "site",
                "detail": "Ship the same snapshot to GitHub Pages without claiming unexplained templates are new formulas.",
            },
        ],
    }

    review_records = [record for record in public_records if record["novelty_status"] == "review"]
    top_leads = [record for record in review_records if record["id"] in {"cb60fd71d1d7", "1125ffe48b3b"}]
    if len(top_leads) < 2:
        seen_ids = {record["id"] for record in top_leads}
        for record in review_records:
            if record["id"] in seen_ids:
                continue
            top_leads.append(record)
            seen_ids.add(record["id"])
            if len(top_leads) >= 2:
                break
    payload["top_leads"] = top_leads

    (directory / ".nojekyll").write_text("", encoding="utf-8")
    (directory / "results.json").write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    (directory / "styles.css").write_text(_site_css(), encoding="utf-8")
    (directory / "index.html").write_text(_site_html(payload), encoding="utf-8")


def _site_html(payload: dict[str, object]) -> str:
    data_json = json.dumps(payload, ensure_ascii=True)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{payload["title"]}</title>
  <link rel="stylesheet" href="styles.css">
</head>
<body>
  <main class="page">
    <section class="hero">
      <p class="eyebrow">Ramanujan-style q-continued fractions</p>
      <h1>{payload["title"]}</h1>
      <p class="lede">A static snapshot of benchmark rediscoveries and exploratory candidates generated by the local-first CLI pipeline.</p>
      <div class="meta">
        <span>Generated: {payload["generated_at"]}</span>
        <span>Records: {len(payload["records"])}</span>
        <a href="review-audit.html">Review audit</a>
      </div>
    </section>
    <section class="stats" id="stats"></section>
    <section class="panel spotlight">
      <div class="panel-header">
        <h2>Top Leads</h2>
        <p>The current release only highlights two weak leads for deeper manual audit. They are not public novelty claims.</p>
      </div>
      <div class="cards" id="top-leads"></div>
    </section>
    <section class="info-grid">
      <section class="panel">
        <div class="panel-header">
          <h2>Benchmark Families</h2>
          <p>The catalog stays intentionally small and independently checkable.</p>
        </div>
        <div class="stack" id="benchmark-families"></div>
      </section>
      <section class="panel">
        <div class="panel-header">
          <h2>Pipeline Discipline</h2>
          <p>This release is designed to surface candidates, not to auto-announce discoveries.</p>
        </div>
        <div class="stack" id="pipeline-list"></div>
      </section>
    </section>
    <section class="panel">
      <div class="panel-header">
        <h2>Candidate Gallery</h2>
        <p>Known benchmarks and review-worthy unexplained templates are rendered side by side.</p>
      </div>
      <div class="cards" id="cards"></div>
    </section>
  </main>
  <script>
    const payload = {data_json};
    const stats = document.getElementById("stats");
    const topLeads = document.getElementById("top-leads");
    const benchmarkFamilies = document.getElementById("benchmark-families");
    const pipelineList = document.getElementById("pipeline-list");
    const cards = document.getElementById("cards");

    const renderCard = (record, className = "") => {{
      const card = document.createElement("article");
      card.className = `card status-${{record.novelty_status}} ${{className}}`.trim();
      card.innerHTML = `
        <div class="card-top">
          <span class="status-pill">${{record.novelty_status}}</span>
          <span class="digits">${{record.digits_agree}} digits</span>
        </div>
        <h3>${{record.id}}</h3>
        <p class="target">${{record.matched_target}} <span>${{record.benchmark_kind}}</span></p>
        <p class="notes">${{record.notes}}</p>
        <dl>
          <div><dt>Stability</dt><dd>${{record.stability_score}}</dd></div>
          <div><dt>Closest benchmark</dt><dd>${{record.closest_benchmark}} (${{record.closest_benchmark_digits}} digits)</dd></div>
          <div><dt>Family bucket</dt><dd><code>${{record.family_bucket}}</code></dd></div>
          <div><dt>Complexity</dt><dd>${{record.complexity_score}}</dd></div>
        </dl>
        <details>
          <summary>Show details</summary>
          <p><strong>Signature</strong><br><code>${{record.template_signature}}</code></p>
          <p><strong>q values</strong><br><code>${{record.q_values.join(", ")}}</code></p>
          <p><strong>Estimates</strong><br><code>${{record.value_estimates.join("<br>")}}</code></p>
          <p><strong>LaTeX</strong><br><code>${{record.template_latex}}</code></p>
        </details>
      `;
      return card;
    }};

    Object.entries(payload.status_counts).forEach(([key, value]) => {{
      const item = document.createElement("article");
      item.className = "stat";
      item.innerHTML = `<h3>${{key}}</h3><p>${{value}}</p>`;
      stats.appendChild(item);
    }});

    payload.top_leads.forEach((record) => {{
      topLeads.appendChild(renderCard(record, "lead-card"));
    }});

    payload.benchmark_families.forEach((family) => {{
      const item = document.createElement("article");
      item.className = "info-card";
      item.innerHTML = `
        <h3>${{family.name}}</h3>
        <p class="digits">${{family.count}} catalog entries</p>
        <p class="notes">${{family.description}}</p>
      `;
      benchmarkFamilies.appendChild(item);
    }});

    payload.pipeline.forEach((entry) => {{
      const item = document.createElement("article");
      item.className = "info-card";
      item.innerHTML = `
        <h3>${{entry.step}}</h3>
        <p class="notes">${{entry.detail}}</p>
      `;
      pipelineList.appendChild(item);
    }});

    payload.records.forEach((record) => {{
      cards.appendChild(renderCard(record));
    }});
  </script>
</body>
</html>
"""


def _site_css() -> str:
    return """\
:root {
  --bg: #f6f1e8;
  --panel: #fffaf2;
  --ink: #1f1b16;
  --muted: #6e6456;
  --accent: #bb4d00;
  --accent-soft: #f3c9a2;
  --review: #0e6f7a;
  --known: #256029;
  --shadow: rgba(38, 28, 16, 0.12);
  --border: rgba(70, 50, 30, 0.12);
  --font-sans: "IBM Plex Sans", "Segoe UI", sans-serif;
  --font-serif: "IBM Plex Serif", Georgia, serif;
}

* {
  box-sizing: border-box;
}

body {
  margin: 0;
  background:
    radial-gradient(circle at top right, rgba(187, 77, 0, 0.12), transparent 28%),
    linear-gradient(180deg, #f8f2e7 0%, #f2eadc 100%);
  color: var(--ink);
  font-family: var(--font-sans);
}

.page {
  width: min(1120px, calc(100% - 32px));
  margin: 0 auto;
  padding: 32px 0 64px;
}

.hero {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 28px;
  padding: 32px;
  box-shadow: 0 18px 40px var(--shadow);
}

.eyebrow {
  margin: 0;
  text-transform: uppercase;
  letter-spacing: 0.18em;
  color: var(--accent);
  font-size: 0.8rem;
}

.hero h1,
.panel h2 {
  font-family: var(--font-serif);
}

.hero h1 {
  margin: 12px 0;
  font-size: clamp(2.2rem, 6vw, 4.2rem);
}

.lede {
  max-width: 700px;
  color: var(--muted);
  font-size: 1.05rem;
  line-height: 1.6;
}

.meta {
  display: flex;
  flex-wrap: wrap;
  gap: 12px;
  margin-top: 20px;
  color: var(--muted);
}

.stats {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 14px;
  margin: 24px 0;
}

.info-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 18px;
  margin-bottom: 24px;
}

.stat,
.card,
.panel {
  background: var(--panel);
  border: 1px solid var(--border);
  box-shadow: 0 10px 24px var(--shadow);
}

.stat {
  border-radius: 20px;
  padding: 18px 20px;
}

.stat h3 {
  margin: 0;
  font-size: 0.95rem;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--muted);
}

.stat p {
  margin: 8px 0 0;
  font-size: 2rem;
  font-family: var(--font-serif);
}

.panel {
  border-radius: 28px;
  padding: 24px;
}

.spotlight {
  margin-bottom: 24px;
}

.panel-header p {
  color: var(--muted);
}

.cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 18px;
}

.stack {
  display: grid;
  gap: 14px;
}

.card {
  border-radius: 22px;
  padding: 18px;
  display: flex;
  flex-direction: column;
  gap: 12px;
}

.card-top {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 8px;
}

.status-pill {
  border-radius: 999px;
  padding: 6px 12px;
  font-size: 0.82rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  background: var(--accent-soft);
  color: var(--accent);
}

.status-known .status-pill,
.status-known_variant .status-pill,
.status-fixture .status-pill {
  background: rgba(37, 96, 41, 0.12);
  color: var(--known);
}

.status-review .status-pill {
  background: rgba(14, 111, 122, 0.12);
  color: var(--review);
}

.lead-card {
  border-color: rgba(14, 111, 122, 0.25);
  background:
    linear-gradient(180deg, rgba(14, 111, 122, 0.06), transparent 28%),
    var(--panel);
}

.card h3 {
  margin: 0;
  font-family: var(--font-serif);
  font-size: 1.4rem;
}

.target {
  margin: 0;
  color: var(--muted);
}

.target span {
  display: inline-block;
  margin-left: 6px;
  padding: 2px 8px;
  border-radius: 999px;
  background: rgba(0, 0, 0, 0.04);
  font-size: 0.78rem;
}

.notes {
  margin: 0;
  line-height: 1.5;
}

dl {
  margin: 0;
}

dl div {
  display: grid;
  grid-template-columns: 92px 1fr;
  gap: 8px;
  padding: 6px 0;
  border-bottom: 1px dashed rgba(0, 0, 0, 0.08);
}

dt {
  color: var(--muted);
}

dd {
  margin: 0;
}

.info-card {
  border-radius: 20px;
  padding: 16px;
  border: 1px solid var(--border);
  background: rgba(255, 255, 255, 0.55);
}

.info-card h3 {
  margin: 0 0 8px;
  font-family: var(--font-serif);
}

details {
  margin-top: auto;
}

summary {
  cursor: pointer;
  color: var(--accent);
  font-weight: 700;
}

code {
  font-family: "Cascadia Code", "Fira Code", monospace;
  font-size: 0.92rem;
  white-space: pre-wrap;
  word-break: break-word;
}

@media (max-width: 640px) {
  .page {
    width: min(100% - 20px, 1120px);
  }

  .hero,
  .panel {
    padding: 20px;
    border-radius: 20px;
  }
}
"""
