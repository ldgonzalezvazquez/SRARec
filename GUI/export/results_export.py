"""Export analysis results and figures."""

from __future__ import annotations

import shutil
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Flag, auto
from pathlib import Path
from typing import Iterator, List, Optional

import matplotlib.pyplot as plt

from GUI.core.results_filter import filter_lines, list_sras, write_filtered_results
from GUI.plotting.profile_plot import _load_graphs_module


class ExportFormat(Flag):
    GAMETES = auto()
    RECOMBINATION = auto()
    COVERAGE = auto()
    BEGFIN = auto()
    RAW = auto()
    ALL = GAMETES | RECOMBINATION | COVERAGE | BEGFIN | RAW


@dataclass
class ExportOutcome:
    paths: list[Path]
    message: str


def _parse_rec_events(line: str) -> List[dict]:
    if "Rec:" not in line:
        return []
    sra = line.split("Sra:")[1].split("\t")[0] if "Sra:" in line else "unknown"
    rec_field = line.split("Rec:")[1].split("No_Rec:")[0].strip()
    recs = [r.strip() for r in rec_field.split(";") if r.strip()]
    rows = []

    for idx, rec in enumerate(recs, start=1):
        if "/" in rec:
            sub_events = rec.split("/")
            parts = sub_events[0].split(",", 6)
            if len(parts) < 7:
                continue
            start, end, count = parts[0], parts[1], parts[2]
            h1, h2, mr, mf = parts[3], parts[4], parts[5], parts[6]
            rows.append(_row_from_parts(sra, idx, start, end, count, h1, h2, mr, mf))
            for sub in sub_events[1:]:
                sub_parts = sub.split(",", 3)
                if len(sub_parts) < 4:
                    continue
                h1, h2, mr, mf = sub_parts[0], sub_parts[1], sub_parts[2], sub_parts[3]
                rows.append(_row_from_parts(sra, idx, start, end, count, h1, h2, mr, mf))
        else:
            parts = rec.split(",", 6)
            if len(parts) < 7:
                continue
            start, end, count, h1, h2, mr, mf = parts
            rows.append(_row_from_parts(sra, idx, start, end, count, h1, h2, mr, mf))

    return rows


def _row_from_parts(sra, idx, start, end, count, h1, h2, mr, mf):
    mr_val = mr.split(":")[1] if ":" in mr else mr
    try:
        mf_val = float(mf.split(":")[1].split()[0]) * 100
    except (ValueError, IndexError):
        mf_val = 0.0
    a1, a2 = h1.split(".")
    b1, b2 = h2.split(".")
    gametes = " ".join([a1 + b1, a1 + b2, a2 + b1, a2 + b2])
    return {
        "SRA": sra,
        "Event": idx,
        "Start": start,
        "End": end,
        "Count": count,
        "Gametes": gametes,
        "Minor_Reads": mr_val,
        "Minor_Freq_pct": f"{mf_val:.2f}",
    }


@contextmanager
def _filtered_results_path(
    results_path: Path,
    sra_ids: list[str] | None,
) -> Iterator[Path]:
    if not sra_ids:
        yield Path(results_path)
        return
    all_sras = list_sras(results_path)
    if set(sra_ids) >= set(all_sras):
        yield Path(results_path)
        return
    tmp = results_path.parent / f".{results_path.stem}_filtered_tmp.txt"
    try:
        write_filtered_results(results_path, tmp, sra_ids)
        yield tmp
    finally:
        if tmp.is_file():
            tmp.unlink(missing_ok=True)


def export_formatted_tsv(
    results_path: Path,
    output_path: Path,
    *,
    sra_ids: list[str] | None = None,
) -> None:
    """Export recombination events (gamete summary) to TSV."""
    lines = filter_lines(results_path, sra_ids)
    rows: list[dict] = []
    for line in lines:
        rows.extend(_parse_rec_events(line))

    header = "\t".join(
        ["SRA", "Event", "Start", "End", "Count", "Gametes", "Minor_Reads", "Minor_Freq_pct"]
    )
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write(header + "\n")
        for row in rows:
            fout.write(
                "\t".join(
                    str(row[k])
                    for k in [
                        "SRA",
                        "Event",
                        "Start",
                        "End",
                        "Count",
                        "Gametes",
                        "Minor_Reads",
                        "Minor_Freq_pct",
                    ]
                )
                + "\n"
            )


def export_raw_lines(
    results_path: Path,
    output_path: Path,
    *,
    sra_ids: list[str] | None = None,
) -> None:
    lines = filter_lines(results_path, sra_ids)
    if not lines:
        raise ValueError("No result lines match the selected samples.")
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_raw_line(results_path: Path, output_path: Path, *, last_line_only: bool = True) -> None:
    """Copy raw results (last line or full file)."""
    text = Path(results_path).read_text(encoding="utf-8")
    if last_line_only:
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            raise ValueError("Results file is empty.")
        Path(output_path).write_text(lines[-1] + "\n", encoding="utf-8")
    else:
        shutil.copy2(results_path, output_path)


def _genome_len_for_results(results_path: Path, annotation_path: Path | None) -> int | None:
    graphs = _load_graphs_module()
    genome_len = None
    if annotation_path and Path(annotation_path).is_file():
        _, gb_len = graphs._load_genes(str(annotation_path), genome_len=None)
        if gb_len is not None:
            genome_len = gb_len
    if genome_len is None:
        try:
            events = graphs._load_events(str(results_path))
            if not events.empty:
                genome_len = int(max(events["pA"].max(), events["pB"].max())) + 1
        except Exception:
            genome_len = None
    return genome_len


def export_recombination_values(
    results_path: Path,
    output_path: Path,
    *,
    sra_ids: list[str] | None = None,
    annotation_path: Path | None = None,
    alpha: float = 1.0,
) -> None:
    graphs = _load_graphs_module()
    with _filtered_results_path(results_path, sra_ids) as filtered:
        genome_len = _genome_len_for_results(filtered, annotation_path)
        events = graphs._load_events(str(filtered))
        positions, rho = graphs.compute_recombination_profile(
            events,
            genome_len=genome_len,
            alpha=alpha,
            exclusive_positions=False,
            aggregate_by_patient=False,
        )
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("Position\tRecombination_Ratio\n")
        for pos, value in zip(positions, rho):
            fout.write(f"{int(pos)}\t{value:.6f}\n")


def export_coverage_values(
    results_path: Path,
    output_path: Path,
    *,
    sra_ids: list[str] | None = None,
    annotation_path: Path | None = None,
) -> None:
    graphs = _load_graphs_module()
    with _filtered_results_path(results_path, sra_ids) as filtered:
        genome_len = _genome_len_for_results(filtered, annotation_path)
        if genome_len is None:
            raise ValueError("Could not determine genome length for coverage export.")
        coverage = graphs._compute_coverage_from_begfin(str(filtered), genome_len=genome_len)
        events = graphs._load_events(str(filtered))
        positions, _ = graphs.compute_recombination_profile(
            events,
            genome_len=genome_len,
            exclusive_positions=False,
            aggregate_by_patient=False,
        )
    if coverage is None or not coverage.size:
        raise ValueError("No BegFin_Reads coverage data available.")
    cov_len = min(len(coverage), len(positions))
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("Position\tCoverage\n")
        for pos, cov in zip(positions[:cov_len], coverage[:cov_len]):
            fout.write(f"{int(pos)}\t{cov:.6f}\n")


def export_begfin_concentration(
    results_path: Path,
    output_path: Path,
    *,
    sra_ids: list[str] | None = None,
    annotation_path: Path | None = None,
) -> None:
    graphs = _load_graphs_module()
    with _filtered_results_path(results_path, sra_ids) as filtered:
        genome_len = _genome_len_for_results(filtered, annotation_path)
        if genome_len is None:
            raise ValueError("Could not determine genome length for BegFin export.")
        concentration = graphs._compute_begfin_concentration(str(filtered), genome_len=genome_len)
        events = graphs._load_events(str(filtered))
        positions, _ = graphs.compute_recombination_profile(
            events,
            genome_len=genome_len,
            exclusive_positions=False,
            aggregate_by_patient=False,
        )
    if concentration is None or not concentration.size:
        raise ValueError("No BegFin_Reads concentration data available.")
    conc_len = min(len(concentration), len(positions))
    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("Position\tRead_Start_End_Count\n")
        for pos, value in zip(positions[:conc_len], concentration[:conc_len]):
            fout.write(f"{int(pos)}\t{value:.6f}\n")


def _output_stem(prefix: str, sra_ids: list[str] | None) -> str:
    if not sra_ids:
        return prefix
    all_count = len(sra_ids)
    if all_count <= 3:
        tag = "_".join(sra_ids)
    else:
        tag = f"{all_count}_samples"
    return f"{prefix}_{tag}"


def export_results_bundle(
    results_path: Path,
    output_dir: Path,
    *,
    formats: ExportFormat,
    sra_ids: list[str] | None = None,
    file_prefix: str | None = None,
    annotation_path: Path | None = None,
) -> ExportOutcome:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = file_prefix or Path(results_path).stem
    stem = _output_stem(prefix, sra_ids)
    saved: list[Path] = []

    selected = formats
    if not selected:
        raise ValueError("Select at least one export format.")

    if ExportFormat.GAMETES in selected:
        path = output_dir / f"{stem}_gametes.tsv"
        export_formatted_tsv(results_path, path, sra_ids=sra_ids)
        saved.append(path)

    if ExportFormat.RAW in selected:
        path = output_dir / f"{stem}_raw.txt"
        export_raw_lines(results_path, path, sra_ids=sra_ids)
        saved.append(path)

    if ExportFormat.RECOMBINATION in selected:
        path = output_dir / f"{stem}_recombination_values.txt"
        export_recombination_values(
            results_path, path, sra_ids=sra_ids, annotation_path=annotation_path
        )
        saved.append(path)

    if ExportFormat.COVERAGE in selected:
        path = output_dir / f"{stem}_coverage_values.txt"
        export_coverage_values(
            results_path, path, sra_ids=sra_ids, annotation_path=annotation_path
        )
        saved.append(path)

    if ExportFormat.BEGFIN in selected:
        path = output_dir / f"{stem}_begfin_concentration.txt"
        export_begfin_concentration(
            results_path, path, sra_ids=sra_ids, annotation_path=annotation_path
        )
        saved.append(path)

    if not saved:
        raise ValueError("Select at least one export format.")

    return ExportOutcome(
        paths=saved,
        message=f"Saved {len(saved)} file(s) to {output_dir.resolve()}",
    )


def save_figure(fig: plt.Figure, output_path: Path, *, dpi: int = 300) -> None:
    """Save matplotlib figure to PNG, PDF, or SVG based on extension."""
    ext = output_path.suffix.lower().lstrip(".")
    fmt = ext if ext in {"png", "pdf", "svg", "jpg", "jpeg"} else "png"
    fig.savefig(output_path, format=fmt, dpi=dpi, bbox_inches="tight")


def read_last_results_line(results_path: Path) -> Optional[str]:
    lines = [ln.strip() for ln in Path(results_path).read_text(encoding="utf-8").splitlines() if ln.strip()]
    return lines[-1] if lines else None
