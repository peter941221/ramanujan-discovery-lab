from pathlib import Path
import json

from ramanujan_discovery.cli import main


def test_cli_end_to_end(tmp_path: Path):
    candidates = tmp_path / "candidates.jsonl"
    verified = tmp_path / "verified.jsonl"
    report = tmp_path / "report.md"
    analysis = tmp_path / "analysis.md"
    site_dir = tmp_path / "site"

    assert main(
        [
            "discover",
            "--depth",
            "36",
            "--precision",
            "80",
            "--budget-hours",
            "0.1",
            "--out",
            str(candidates),
        ]
    ) == 0
    assert candidates.exists()

    assert main(
        [
            "verify",
            "--in",
            str(candidates),
            "--precision",
            "160",
            "--out",
            str(verified),
        ]
    ) == 0
    assert verified.exists()

    assert main(
        [
            "report",
            "--in",
            str(verified),
            "--out",
            str(report),
        ]
    ) == 0
    assert report.exists()
    assert "Ramanujan Discovery Report" in report.read_text(encoding="utf-8")
    assert "Equivalence key" in report.read_text(encoding="utf-8")

    assert main(
        [
            "analyze",
            "--in",
            str(verified),
            "--candidate-id",
            "cb60fd71d1d7",
            "--out",
            str(analysis),
        ]
    ) == 0
    assert analysis.exists()
    assert "Hero Case Analysis" in analysis.read_text(encoding="utf-8")

    assert main(
        [
            "site",
            "--in",
            str(verified),
            "--out-dir",
            str(site_dir),
        ]
    ) == 0
    assert (site_dir / "index.html").exists()
    assert (site_dir / "results.json").exists()
    payload = json.loads((site_dir / "results.json").read_text(encoding="utf-8"))
    assert all("equivalence_key" in record for record in payload["records"])
    assert any(note["title"] == "RR(q^3) Neighborhood" for note in payload["audit_notes"])
    index_text = (site_dir / "index.html").read_text(encoding="utf-8")
    assert "Audit Notes" in index_text
    assert "E2CC Plain-Step Audit" in index_text
