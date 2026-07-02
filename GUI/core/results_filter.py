"""Filter SRARec results files by sample (SRA / BAM id)."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def extract_sra_id(line: str) -> str | None:
    line = line.strip()
    if not line or not line.startswith("Sra:"):
        return None
    return line.split("Sra:")[1].split("\t")[0].strip() or None


def list_sras(results_path: Path) -> list[str]:
    """Return unique SRA/sample ids in file order."""
    path = Path(results_path)
    if not path.is_file():
        return []
    seen: set[str] = set()
    ordered: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        sra = extract_sra_id(line)
        if sra and sra not in seen:
            seen.add(sra)
            ordered.append(sra)
    return ordered


def filter_lines(results_path: Path, sra_ids: Iterable[str] | None = None) -> list[str]:
    """Return result lines, optionally restricted to *sra_ids*."""
    path = Path(results_path)
    if not path.is_file():
        raise FileNotFoundError(f"Results file not found: {path}")

    allowed = set(sra_ids) if sra_ids is not None else None
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text:
            continue
        sra = extract_sra_id(text)
        if sra is None:
            continue
        if allowed is None or sra in allowed:
            lines.append(text)
    return lines


def write_filtered_results(
    results_path: Path,
    output_path: Path,
    sra_ids: Iterable[str] | None = None,
) -> Path:
    lines = filter_lines(results_path, sra_ids)
    if not lines:
        raise ValueError("No result lines match the selected samples.")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path
