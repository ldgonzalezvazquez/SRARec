"""Generate pileup files using pysam (embedded HTSlib, no external samtools binary)."""

from __future__ import annotations

from pathlib import Path

import pysam


def write_mpileup(
    bam_path: Path,
    ref_fasta: Path,
    pileup_path: Path,
    *,
    min_base_qual: int = 20,
    include_anomalous: bool = True,
) -> None:
    """
    Write a tab-delimited pileup file compatible with SRARec's parser.

    Flags mirror pasos_ravi: -A -B -Q <min_base_qual> -f <ref>
    """
    bam_path = Path(bam_path).resolve()
    ref_fasta = Path(ref_fasta).resolve()
    pileup_path = Path(pileup_path)
    pileup_path.parent.mkdir(parents=True, exist_ok=True)

    args = [
        "-B",
        "-Q",
        str(min_base_qual),
        "-f",
        str(ref_fasta),
        str(bam_path),
    ]
    if include_anomalous:
        args = ["-A"] + args

    pysam.samtools.mpileup(*args, save_stdout=str(pileup_path))

    if not pileup_path.is_file() or pileup_path.stat().st_size == 0:
        raise RuntimeError(f"mpileup produced an empty or missing file: {pileup_path}")
