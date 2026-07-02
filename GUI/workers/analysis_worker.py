"""Background analysis worker (QThread)."""

from __future__ import annotations

import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal

from GUI.core.bam_analysis import AnalysisCancelled, format_results_summary, run_bam_analysis


class AnalysisWorker(QThread):
    log_message = Signal(str)
    sample_finished = Signal(object, str, int, int, str)  # path, summary, index, total, sample_id
    finished_ok = Signal(object, str)  # results_path, final message
    finished_error = Signal(str)

    def __init__(
        self,
        bam_paths: list[Path],
        output_dir: Path,
        reference_key: str,
        platform: str,
        min_reads: int,
        min_gamete_reads: int,
        min_gamete_freq: float,
        no_indels: bool,
        verbose: bool,
        parent=None,
    ):
        super().__init__(parent)
        self.bam_paths = [Path(p) for p in bam_paths]
        self.output_dir = output_dir
        self.reference_key = reference_key
        self.platform = platform
        self.min_reads = min_reads
        self.min_gamete_reads = min_gamete_reads
        self.min_gamete_freq = min_gamete_freq
        self.no_indels = no_indels
        self.verbose = verbose
        self.cancel_event = threading.Event()

    def cancel(self) -> None:
        self.cancel_event.set()

    def run(self) -> None:
        total = len(self.bam_paths)
        results_path: Path | None = None
        try:
            for index, bam_path in enumerate(self.bam_paths, start=1):
                if self.cancel_event.is_set():
                    raise AnalysisCancelled()

                sample_id = bam_path.stem
                self.log_message.emit(
                    f"[BATCH] Sample {index}/{total}: {sample_id} ({bam_path})"
                )

                results_path = run_bam_analysis(
                    bam_path,
                    self.output_dir,
                    self.reference_key,
                    self.platform,
                    self.min_reads,
                    self.min_gamete_reads,
                    self.min_gamete_freq,
                    no_indels=self.no_indels,
                    verbose=self.verbose,
                    log_callback=self.log_message.emit,
                    cancel_event=self.cancel_event,
                    clear_results=(index == 1),
                )

                lines = results_path.read_text(encoding="utf-8").strip().splitlines()
                summary = format_results_summary(lines[-1]) if lines else "No results."
                self.sample_finished.emit(
                    results_path, summary, index, total, sample_id
                )

            if results_path is None:
                raise RuntimeError("No BAM files were processed.")

            self.finished_ok.emit(
                results_path,
                f"Batch complete: {total} sample(s). Results: {results_path}",
            )
        except AnalysisCancelled:
            self.finished_error.emit("Analysis cancelled by user.")
        except Exception as exc:
            self.finished_error.emit(str(exc))
