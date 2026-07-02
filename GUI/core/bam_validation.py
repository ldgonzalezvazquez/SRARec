"""Validate BAM files against the SRARec Desktop input contract."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import pysam


@dataclass
class BamValidationResult:
    ok: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    index_created: bool = False

    @property
    def message(self) -> str:
        parts = []
        if self.errors:
            parts.append("Errors:\n" + "\n".join(f"  - {e}" for e in self.errors))
        if self.warnings:
            parts.append("Warnings:\n" + "\n".join(f"  - {w}" for w in self.warnings))
        return "\n\n".join(parts)


def _header_sort_order(bam: pysam.AlignmentFile) -> Optional[str]:
    header = bam.header.to_dict()
    hd = header.get("HD") or {}
    if isinstance(hd, list) and hd:
        hd = hd[0]
    return hd.get("SO")


def validate_bam(
    bam_path: Path,
    *,
    create_index: bool = False,
) -> BamValidationResult:
    """
    Validate that *bam_path* meets the SRARec Desktop BAM contract.

    If *create_index* is True and the `.bai` is missing, create it with pysam.
    """
    result = BamValidationResult(ok=True)
    bam_path = Path(bam_path).resolve()

    if not bam_path.is_file():
        result.ok = False
        result.errors.append(f"BAM file not found: {bam_path}")
        return result

    bai_path = Path(str(bam_path) + ".bai")
    if not bai_path.is_file():
        if create_index:
            try:
                pysam.index(str(bam_path))
                result.index_created = True
                result.warnings.append(f"Created BAM index: {bai_path.name}")
            except Exception as exc:
                result.ok = False
                result.errors.append(f"Failed to create BAM index: {exc}")
                return result
        else:
            result.ok = False
            result.errors.append(
                "BAM index (.bai) not found. Create it or allow the app to index the BAM."
            )
            return result

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            sort_order = _header_sort_order(bam)
            if sort_order and sort_order != "coordinate":
                result.ok = False
                result.errors.append(
                    f"BAM must be coordinate-sorted (@HD SO:coordinate). Found SO:{sort_order}."
                )

            if not bam.references:
                result.ok = False
                result.errors.append("BAM has no reference sequences in the header.")
            else:
                mapped = 0
                dup_count = 0
                sample_limit = 50_000
                for i, read in enumerate(bam.fetch(until_eof=True)):
                    if read.is_unmapped:
                        continue
                    mapped += 1
                    if read.is_duplicate:
                        dup_count += 1
                    if i >= sample_limit:
                        break

                if mapped == 0:
                    result.warnings.append("No mapped reads found in the BAM (sampled).")
                elif dup_count > 0 and dup_count / max(mapped, 1) > 0.05:
                    result.warnings.append(
                        f"Many duplicate reads detected (~{100 * dup_count / mapped:.1f}% in sample). "
                        "Deduplication before analysis is recommended."
                    )

    except Exception as exc:
        result.ok = False
        result.errors.append(f"Cannot read BAM: {exc}")

    return result
