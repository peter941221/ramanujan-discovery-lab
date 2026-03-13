from __future__ import annotations

import json
from pathlib import Path

from ramanujan_discovery.models import CandidateRecord


def write_candidates(path: str | Path, records: list[CandidateRecord]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record.to_dict(), ensure_ascii=True))
            handle.write("\n")


def read_candidates(path: str | Path) -> list[CandidateRecord]:
    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    records: list[CandidateRecord] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            records.append(CandidateRecord.from_dict(json.loads(stripped)))
    return records
