"""Main application window."""

from __future__ import annotations

import subprocess
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
    QApplication,
)

from GUI.config import APP_NAME
from GUI.core.bam_validation import validate_bam
from GUI.dialogs.about_dialog import AboutDialog
from GUI.dialogs.bam_requirements_dialog import BamRequirementsDialog
from GUI.dialogs.custom_ref_help import CustomRefHelpDialog
from GUI.dialogs.graph_export_dialog import GraphExportDialog
from GUI.dialogs.save_results_dialog import SaveResultsDialog
from GUI.export.results_export import (
    export_results_bundle,
    export_raw_line,
    read_last_results_line,
    save_figure,
)
from GUI.paths import get_annotation_gb
from GUI.plotting.profile_plot import build_profile_figure, generate_profile_html
from GUI.widgets.console_panel import ConsolePanel
from GUI.widgets.graph_panel import GraphPanel
from GUI.widgets.input_panel import InputPanel
from GUI.widgets.params_panel import ParamsPanel
from GUI.workers.analysis_worker import AnalysisWorker


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1400, 900)
        self._worker: AnalysisWorker | None = None
        self._results_path: Path | None = None

        self._build_ui()
        self._build_menus()
        self.statusBar().showMessage("Ready")

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left.setMaximumWidth(420)

        self.input_panel = InputPanel()
        self.params_panel = ParamsPanel()
        left_layout.addWidget(self.input_panel)
        left_layout.addWidget(self.params_panel)

        self.run_btn = QPushButton("Run analysis")
        self.run_btn.clicked.connect(self._run_analysis)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._cancel_analysis)
        left_layout.addWidget(self.run_btn)
        left_layout.addWidget(self.cancel_btn)

        export_row = QHBoxLayout()
        self.save_tsv_btn = QPushButton("Save results...")
        self.save_tsv_btn.clicked.connect(self._save_tsv)
        self.export_raw_btn = QPushButton("Export raw")
        self.export_raw_btn.setEnabled(False)
        self.export_raw_btn.clicked.connect(self._export_raw)
        self.save_graph_btn = QPushButton("Save graph")
        self.save_graph_btn.clicked.connect(
            lambda: self._graph_export_dialog("save")
        )
        export_row.addWidget(self.save_tsv_btn)
        export_row.addWidget(self.export_raw_btn)
        export_row.addWidget(self.save_graph_btn)
        left_layout.addLayout(export_row)
        left_layout.addStretch()

        self.input_panel.custom_ref_btn.clicked.connect(self._show_custom_ref_help)

        right_split = QSplitter(Qt.Vertical)
        self.console_panel = ConsolePanel()
        self.graph_panel = GraphPanel()
        self.graph_panel.set_html_callback(self._open_html_graph)
        self.graph_panel.set_regenerate_callback(
            lambda: self._graph_export_dialog("regenerate")
        )
        right_split.addWidget(self.console_panel)
        right_split.addWidget(self.graph_panel)
        right_split.setSizes([280, 520])

        root.addWidget(left)
        root.addWidget(right_split, stretch=1)

    def _build_menus(self) -> None:
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        act_save_tsv = QAction("Save results...", self)
        act_save_tsv.triggered.connect(self._save_tsv)
        file_menu.addAction(act_save_tsv)
        act_export_raw = QAction("Export raw...", self)
        act_export_raw.triggered.connect(self._export_raw)
        file_menu.addAction(act_export_raw)
        act_save_graph = QAction("Save graph...", self)
        act_save_graph.triggered.connect(lambda: self._graph_export_dialog("save"))
        file_menu.addAction(act_save_graph)
        file_menu.addSeparator()
        act_quit = QAction("Quit", self)
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        run_menu = menubar.addMenu("&Run")
        act_run = QAction("Run analysis", self)
        act_run.triggered.connect(self._run_analysis)
        run_menu.addAction(act_run)
        act_cancel = QAction("Cancel", self)
        act_cancel.triggered.connect(self._cancel_analysis)
        run_menu.addAction(act_cancel)

        help_menu = menubar.addMenu("&Help")
        act_bam = QAction("BAM requirements", self)
        act_bam.triggered.connect(self._show_bam_requirements)
        help_menu.addAction(act_bam)
        act_custom = QAction("Custom reference (CLI)", self)
        act_custom.triggered.connect(self._show_custom_ref_help)
        help_menu.addAction(act_custom)
        help_menu.addSeparator()
        act_about = QAction("About", self)
        act_about.triggered.connect(self._show_about)
        help_menu.addAction(act_about)

    def _set_running(self, running: bool) -> None:
        self.run_btn.setEnabled(not running)
        self.cancel_btn.setEnabled(running)

    def _run_analysis(self) -> None:
        bam_paths = self.input_panel.bam_paths()
        if not bam_paths:
            QMessageBox.warning(
                self,
                "Missing BAM",
                "Please add at least one BAM file.",
            )
            return

        if not self._validate_bam_batch(bam_paths):
            return

        output_dir = self.input_panel.output_dir()
        results_file = self.input_panel.results_file_path()

        self.console_panel.clear()
        self.graph_panel.reset()
        self._results_path = None
        self._set_raw_export_enabled(False)
        self._set_running(True)

        self.console_panel.append(f"[INFO] Output folder: {output_dir.resolve()}")
        self.console_panel.append(f"[INFO] Combined results: {results_file.resolve()}")
        self.console_panel.append(
            f"[INFO] Processing {len(bam_paths)} BAM file(s) in series..."
        )
        self.statusBar().showMessage(
            f"Running batch (0/{len(bam_paths)})..."
        )

        self._worker = AnalysisWorker(
            bam_paths=bam_paths,
            output_dir=output_dir,
            reference_key=self.input_panel.reference_key(),
            platform=self.input_panel.platform(),
            min_reads=self.params_panel.min_reads.value(),
            min_gamete_reads=self.params_panel.min_gamete_reads.value(),
            min_gamete_freq=self.params_panel.min_gamete_freq_fraction(),
            no_indels=self.params_panel.no_indels.isChecked(),
            verbose=self.params_panel.verbose.isChecked(),
        )
        self._worker.log_message.connect(self.console_panel.append)
        self._worker.sample_finished.connect(self._on_sample_finished)
        self._worker.finished_ok.connect(self._on_analysis_ok)
        self._worker.finished_error.connect(self._on_analysis_error)
        self._worker.start()

    def _validate_bam_batch(self, bam_paths: list[Path]) -> bool:
        missing_files = [p for p in bam_paths if not p.is_file()]
        if missing_files:
            QMessageBox.warning(
                self,
                "Missing BAM",
                "The following BAM files were not found:\n"
                + "\n".join(str(p) for p in missing_files),
            )
            return False

        need_index: list[Path] = []
        for bam in bam_paths:
            validation = validate_bam(bam, create_index=False)
            if not validation.ok:
                need_index.append(bam)
            for w in validation.warnings:
                self.console_panel.append(f"[WARN] {bam.name}: {w}")

        if need_index:
            names = "\n".join(p.name for p in need_index)
            reply = QMessageBox.question(
                self,
                "BAM index missing",
                f"Index (.bai) missing for:\n{names}\n\nCreate index now?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return False
            for bam in need_index:
                validation = validate_bam(bam, create_index=True)
                if not validation.ok:
                    QMessageBox.critical(
                        self,
                        "BAM validation failed",
                        f"{bam.name}: {validation.message}",
                    )
                    return False
                for w in validation.warnings:
                    self.console_panel.append(f"[WARN] {bam.name}: {w}")

        return True

    def _refresh_graph(self) -> None:
        if not self._results_path:
            return
        try:
            self.graph_panel.load_from_results(
                self._results_path,
                self.input_panel.reference_key(),
            )
        except Exception as exc:
            self.console_panel.append(f"[WARN] Could not refresh graph: {exc}")

    def _on_sample_finished(
        self,
        results_path,
        summary: str,
        index: int,
        total: int,
        sample_id: str,
    ) -> None:
        self._results_path = Path(results_path)
        self.console_panel.append(f"\n--- {sample_id} ({index}/{total}) ---")
        self.console_panel.append(summary)
        self.statusBar().showMessage(
            f"Sample {index}/{total} finished: {sample_id}"
        )
        self._refresh_graph()
        self._set_raw_export_enabled(True)

    def _on_analysis_ok(self, results_path, message: str) -> None:
        self._set_running(False)
        self._results_path = Path(results_path)
        self.console_panel.append(f"[OK] {message}")
        self._refresh_graph()
        self._set_raw_export_enabled(True)
        self.statusBar().showMessage("Analysis finished successfully.")

    def _cancel_analysis(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self.console_panel.append("[INFO] Cancellation requested...")

    def _on_analysis_error(self, message: str) -> None:
        self._set_running(False)
        self.console_panel.append(f"[ERROR] {message}")
        QMessageBox.critical(self, "Analysis failed", message)
        self.statusBar().showMessage("Analysis failed.")

    def _set_raw_export_enabled(self, enabled: bool) -> None:
        self.export_raw_btn.setEnabled(enabled)

    def _graph_export_dialog(self, mode: str) -> None:
        dlg = GraphExportDialog(
            self,
            mode=mode,
            default_results=self._results_path,
            default_output_dir=self.input_panel.output_dir(),
        )
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        opts = dlg.options()
        ann = get_annotation_gb(self.input_panel.reference_key())
        saved: list[str] = []
        self._results_path = opts.results_path

        try:
            if opts.update_display or (opts.save_png and not opts.save_html):
                fig = build_profile_figure(
                    opts.results_path,
                    ann,
                    sra_ids=opts.sra_ids,
                )
                if opts.update_display:
                    self.graph_panel.show_figure(fig)
                    self.graph_panel._results_path = opts.results_path
                    self.graph_panel._reference_key = self.input_panel.reference_key()
                if opts.save_png and opts.png_path and not opts.save_html:
                    save_figure(fig, opts.png_path)
                    saved.append(str(opts.png_path))

            if opts.save_html and opts.html_path:
                html_path = generate_profile_html(
                    opts.results_path,
                    ann,
                    opts.html_path,
                    sra_ids=opts.sra_ids,
                    output_png=opts.png_path if opts.save_png else None,
                )
                saved.append(str(html_path))
                if opts.save_png and opts.png_path and opts.html_path:
                    saved.append(str(opts.png_path))

            if saved:
                unique_saved = list(dict.fromkeys(saved))
                self.console_panel.append("[OK] Graph export:")
                for path in unique_saved:
                    self.console_panel.append(f"  {path}")
                QMessageBox.information(
                    self,
                    "Graph saved",
                    "Files written:\n" + "\n".join(unique_saved),
                )
            elif opts.update_display:
                self.console_panel.append("[OK] Graph regenerated in panel.")
        except Exception as exc:
            QMessageBox.critical(self, "Graph export failed", str(exc))

    def _save_tsv(self) -> None:
        dlg = SaveResultsDialog(
            self,
            default_results=self._results_path,
            default_output_dir=self.input_panel.output_dir(),
        )
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return

        opts = dlg.options()
        ann = get_annotation_gb(self.input_panel.reference_key())
        try:
            outcome = export_results_bundle(
                opts.results_path,
                opts.output_dir,
                formats=opts.formats,
                sra_ids=opts.sra_ids,
                file_prefix=opts.file_prefix,
                annotation_path=ann,
            )
            self.console_panel.append(f"[OK] {outcome.message}")
            for path in outcome.paths:
                self.console_panel.append(f"  {path}")
            QMessageBox.information(
                self,
                "Results exported",
                outcome.message + "\n\n" + "\n".join(str(p) for p in outcome.paths),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _export_raw(self) -> None:
        if not self._results_path:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export raw results", "", "Text files (*.txt)"
        )
        if path:
            export_raw_line(self._results_path, Path(path), last_line_only=True)
            QMessageBox.information(self, "Saved", f"Raw line saved to:\n{path}")

    def _open_html_graph(self) -> None:
        if not self._results_path:
            return
        try:
            ann = get_annotation_gb(self.input_panel.reference_key())
            out_dir = self.input_panel.output_dir()
            out_dir.mkdir(parents=True, exist_ok=True)
            html_path = out_dir / f"{self._results_path.stem}_profile.html"
            html_path = generate_profile_html(self._results_path, ann, html_path)
            if not html_path.is_file():
                raise FileNotFoundError(f"HTML file was not created: {html_path}")
            subprocess.run(["open", str(html_path)], check=True)
            self.console_panel.append(f"[OK] HTML graph: {html_path.resolve()}")
        except Exception as exc:
            QMessageBox.critical(self, "HTML graph failed", str(exc))

    def _show_about(self) -> None:
        AboutDialog(self).exec()

    def _show_custom_ref_help(self) -> None:
        CustomRefHelpDialog(self).exec()

    def _show_bam_requirements(self) -> None:
        BamRequirementsDialog(self).exec()

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.cancel()
            self._worker.wait(3000)
        event.accept()
        app = QApplication.instance()
        if app is not None:
            app.quit()
