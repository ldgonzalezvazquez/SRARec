"""BAM-only recombination analysis for SRARec Desktop."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Callable, Optional

import pysam
from tqdm import tqdm

from GUI.config import MIN_BASE_QUAL
from GUI.core.bam_validation import validate_bam
from GUI.core.pileup_pysam import write_mpileup
from GUI.core.srarec_loader import load_srarec_core
from GUI.paths import get_reference_fasta


LogCallback = Optional[Callable[[str], None]]


def _to_one_based(pos: int) -> int:
    """Convert internal 0-based coordinates to 1-based for results output."""
    return int(pos) + 1


class AnalysisCancelled(Exception):
    """Raised when the user cancels a running analysis."""


def _log(cb: LogCallback, message: str) -> None:
    if cb:
        cb(message)


def run_bam_analysis(
    bam_path: Path,
    output_dir: Path,
    reference_key: str,
    platform: str,
    min_reads: int,
    min_gamete_reads: int,
    min_gamete_freq: float,
    *,
    no_indels: bool = False,
    verbose: bool = False,
    log_callback: LogCallback = None,
    cancel_event=None,
    clear_results: bool = True,
) -> Path:
    """
    Run read-level recombination analysis from a preprocessed BAM file.

    Returns the path to ``results_ravi.txt`` in *output_dir*.
    """
    bam_path = Path(bam_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    validation = validate_bam(bam_path, create_index=True)
    if not validation.ok:
        raise ValueError(validation.message or "BAM validation failed.")
    for warning in validation.warnings:
        _log(log_callback, f"[WARN] {warning}")

    refseq = get_reference_fasta(reference_key)
    sample_id = bam_path.stem
    work_dir = output_dir / sample_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # Use the BAM in the work directory (symlink to avoid copying large files)
    work_bam = work_dir / bam_path.name
    if not work_bam.exists():
        try:
            os.symlink(bam_path, work_bam)
        except OSError:
            shutil.copy2(bam_path, work_bam)

    bai = Path(str(work_bam) + ".bai")
    if not bai.is_file():
        pysam.index(str(work_bam))

    core = load_srarec_core()
    core.VERBOSE = verbose
    core.MIN_READS_LIST = [int(min_reads)]
    core.READS_GAMETOS_MINORITARIOS_LIST = [int(min_gamete_reads)]
    core.THRESHOLD_GAMETO_FREQ_LIST = [float(min_gamete_freq)]

    original_cwd = os.getcwd()
    result_filename = output_dir / "results_ravi.txt"
    if clear_results and result_filename.is_file():
        result_filename.unlink()

    modo = "ravi"
    pileup_file = work_dir / f"{sample_id}_ravi.pileup"

    try:
        os.chdir(work_dir)
        _log(log_callback, f"[INFO] Sample: {sample_id}")
        _log(log_callback, f"[INFO] Reference: {refseq}")
        _log(log_callback, f"[INFO] BAM: {work_bam}")

        if not pileup_file.is_file():
            _log(log_callback, "[STEP] Generating pileup (pysam)...")
            write_mpileup(
                work_bam,
                refseq,
                pileup_file,
                min_base_qual=MIN_BASE_QUAL,
                include_anomalous=True,
            )
        else:
            _log(log_callback, "[INFO] Reusing existing pileup file.")

        if cancel_event and cancel_event.is_set():
            raise AnalysisCancelled()

        param_combinations = [
            (int(min_reads), int(min_gamete_reads), float(min_gamete_freq))
        ]

        for mr, rgm, freq_val in tqdm(param_combinations, desc="Thresholds"):
            if cancel_event and cancel_event.is_set():
                raise AnalysisCancelled()

            _log(
                log_callback,
                f"[ANALYSIS] MIN_READS={mr}, gamete_reads={rgm}, gamete_freq={freq_val}",
            )

            core.MIN_READS = mr
            core.reads_gametos_minoritarios = rgm
            core.THRESHOLD_GAMETO_FREQ = freq_val

            posiciones_pileup = core.parse_pileup_and_extract_polimorfismos(
                pileup_file, mr
            )
            _log(log_callback, f"[INFO] Polymorphic positions from pileup: {len(posiciones_pileup)}")

            if len(posiciones_pileup) < 2:
                with open(result_filename, "a", encoding="utf-8") as rf:
                    rf.write(
                        f"Sra:{sample_id}\tCoinf_detected:False\tMin_read_num:{mr}\t"
                        f"Min_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\t"
                        f"Platform:{platform}\tRec:\tNo_Rec:\tBegFin_Reads:\n"
                    )
                continue

            _log(log_callback, "[STEP] Scanning polymorphisms in BAM...")
            with pysam.AlignmentFile(str(work_bam), "rb") as bam:
                ref_names = bam.references
            ref_fa = pysam.FastaFile(str(refseq))
            dic_pol_pos, max_len_reads, reads_in_fin = core.scan_region_for_polymorphisms(
                str(work_bam),
                ref_names[0],
                posiciones_pileup,
                mr,
                ref_fa,
                core.min_base_qual,
                no_indels,
            )
            ref_fa.close()

            _log(log_callback, "[STEP] Four-gamete test...")
            _results, comparaciones_hechas_global, dist_por_pos = (
                core.pipeline_agrupacion_por_distribuciones(
                    dic_pol_pos,
                    max_len=max_len_reads,
                    threshold=0.1,
                    min_size=2,
                )
            )
            core.filtrar_rangos_mas_grandes(comparaciones_hechas_global)

            coinfection_detected = False
            for _pos, freq_dict in dist_por_pos.items():
                count_over_01 = sum(1 for f in freq_dict.values() if f > 0.1)
                if count_over_01 >= 2:
                    coinfection_detected = True
                    break

            true_list = []
            false_list = []
            for (pA, pB), info in comparaciones_hechas_global.items():
                if info["anotar"]:
                    if info["res"] is True:
                        bigger_count = core.count_true_superranges_in(
                            pA, pB, comparaciones_hechas_global
                        )
                        true_list.append((pA, pB, info["detail"], bigger_count))
                    elif info["res"] is False:
                        false_list.append((pA, pB, info["detail"]))

            true_pairs = ";".join(
                f"{_to_one_based(pA)},{_to_one_based(pB)},{bigger_count},{detail}"
                for (pA, pB, detail, bigger_count) in true_list
            )
            false_pairs = ";".join(
                f"{_to_one_based(a)},{_to_one_based(b)},{detail}" for (a, b, detail) in false_list
            )
            reads_in_fin_str = ";".join(
                f"{_to_one_based(st)},{_to_one_based(en)}" for (st, en) in reads_in_fin
            )

            with open(result_filename, "a", encoding="utf-8") as rf:
                rf.write(
                    f"Sra:{sample_id}\tCoinf_detected:{coinfection_detected}\t"
                    f"Min_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\t"
                    f"Platform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\t"
                    f"BegFin_Reads:{reads_in_fin_str}\n"
                )

        if not result_filename.is_file():
            raise RuntimeError("Analysis finished but no results file was written.")

        # Save last line for quick replotting
        lines = result_filename.read_text(encoding="utf-8").strip().splitlines()
        if lines:
            last_file = output_dir / "results_ravi_last.txt"
            last_file.write_text(lines[-1] + "\n", encoding="utf-8")

        _log(log_callback, f"[OK] Results written to {result_filename}")
        return result_filename

    finally:
        os.chdir(original_cwd)


def format_results_summary(results_line: str) -> str:
    """Human-readable summary of recombination events from one results line."""
    if "Rec:" not in results_line:
        return "No results to display."

    try:
        rec_part = results_line.split("Rec:")[1].split("No_Rec:")[0].strip()
    except IndexError:
        return "Could not parse results line."

    events = [e.strip() for e in rec_part.split(";") if e.strip()]
    if not events:
        return "No recombination-compatible events found."

    lines = ["RESULTS:", ""]
    for i, rec in enumerate(events, start=1):
        parts = rec.split(",")
        if len(parts) < 7:
            continue
        start, end, count = parts[0], parts[1], parts[2]
        hap1, hap2 = parts[3], parts[4]
        mireads = parts[5].split(":")[1] if ":" in parts[5] else parts[5]
        try:
            minfreq = float(parts[6].split(":")[1].split()[0]) * 100
        except (ValueError, IndexError):
            minfreq = 0.0

        a1, a2 = hap1.split(".")
        b1, b2 = hap2.split(".")
        gametes = " ".join([a1 + b1, a1 + b2, a2 + b1, a2 + b2])
        lines.append(
            f"- Event {i}: between {start} and {end}, found {count} times "
            f"forming gametes {gametes}, with minimum support of {mireads} reads "
            f"and {minfreq:.2f}% minor gamete frequency"
        )
        lines.append("")

    return "\n".join(lines)
