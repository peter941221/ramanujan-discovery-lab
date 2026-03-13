from pathlib import Path

from ramanujan_discovery.cli import main


def test_cli_end_to_end(tmp_path: Path):
    candidates = tmp_path / "candidates.jsonl"
    verified = tmp_path / "verified.jsonl"
    report = tmp_path / "report.md"
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
