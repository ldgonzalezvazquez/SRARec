"""Dialog for exporting SRARec results in multiple formats."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from GUI.core.results_filter import list_sras
from GUI.export.results_export import ExportFormat


@dataclass
class SaveResultsOptions:
    results_path: Path
    output_dir: Path
    file_prefix: str
    sra_ids: list[str]
    formats: ExportFormat


class SaveResultsDialog(QDialog):
    def __init__(
        self,
        parent=None,
        *,
        default_results: Path | None = None,
        default_output_dir: Path | None = None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Save results")
        self.setMinimumWidth(520)
        self._default_output_dir = default_output_dir or Path.home()

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Choose the results file, samples to include, output folder, "
            "and which tables to export."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        results_row = QHBoxLayout()
        self.results_edit = QLineEdit()
        browse_results = QPushButton("Browse...")
        browse_results.clicked.connect(self._browse_results)
        results_row.addWidget(self.results_edit, stretch=1)
        results_row.addWidget(browse_results)
        layout.addWidget(QLabel("Results file:"))
        layout.addLayout(results_row)

        layout.addWidget(QLabel("Samples to include:"))
        self.sample_list = QListWidget()
        self.sample_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.sample_list.setMaximumHeight(140)
        layout.addWidget(self.sample_list)

        sample_btns = QHBoxLayout()
        select_all = QPushButton("Select all")
        select_all.clicked.connect(self.sample_list.selectAll)
        clear_sel = QPushButton("Clear selection")
        clear_sel.clicked.connect(self.sample_list.clearSelection)
        sample_btns.addWidget(select_all)
        sample_btns.addWidget(clear_sel)
        sample_btns.addStretch()
        layout.addLayout(sample_btns)

        out_row = QHBoxLayout()
        self.output_edit = QLineEdit(str(self._default_output_dir))
        browse_out = QPushButton("Browse...")
        browse_out.clicked.connect(self._browse_output)
        out_row.addWidget(self.output_edit, stretch=1)
        out_row.addWidget(browse_out)
        layout.addWidget(QLabel("Output folder:"))
        layout.addLayout(out_row)

        self.prefix_edit = QLineEdit()
        layout.addWidget(QLabel("File name prefix:"))
        layout.addWidget(self.prefix_edit)

        formats_box = QGroupBox("Export formats")
        formats_layout = QVBoxLayout(formats_box)
        self.gametes_check = QCheckBox("Gamete events summary (TSV)")
        self.recomb_check = QCheckBox("Recombination ratio per position")
        self.coverage_check = QCheckBox("Coverage per position (BegFin_Reads)")
        self.begfin_check = QCheckBox("Read start/end concentration per position")
        self.raw_check = QCheckBox("Raw SRARec result lines")
        self.all_check = QCheckBox("All of the above")
        self.gametes_check.setChecked(True)
        self.all_check.toggled.connect(self._on_all_toggled)
        for box in (
            self.gametes_check,
            self.recomb_check,
            self.coverage_check,
            self.begfin_check,
            self.raw_check,
        ):
            formats_layout.addWidget(box)
        formats_layout.addWidget(self.all_check)
        layout.addWidget(formats_box)

        if default_results:
            self.results_edit.setText(str(default_results))
            self.prefix_edit.setText(default_results.stem)
            self._reload_samples()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Export")
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.results_edit.editingFinished.connect(self._on_results_changed)

    def _on_results_changed(self) -> None:
        self._reload_samples()
        path = Path(self.results_edit.text().strip())
        if path.is_file():
            self.prefix_edit.setText(path.stem)

    def _browse_results(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select results file",
            self.results_edit.text() or str(self._default_output_dir),
            "Results (*.txt *.tsv);;All files (*)",
        )
        if path:
            self.results_edit.setText(path)
            self._reload_samples()
            self.prefix_edit.setText(Path(path).stem)

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Select output folder",
            self.output_edit.text() or str(self._default_output_dir),
        )
        if path:
            self.output_edit.setText(path)

    def _reload_samples(self) -> None:
        self.sample_list.clear()
        path = Path(self.results_edit.text().strip())
        if not path.is_file():
            return
        for sra in list_sras(path):
            self.sample_list.addItem(sra)
        self.sample_list.selectAll()

    def _on_all_toggled(self, checked: bool) -> None:
        for box in (
            self.gametes_check,
            self.recomb_check,
            self.coverage_check,
            self.begfin_check,
            self.raw_check,
        ):
            box.setChecked(checked)

    def _selected_formats(self) -> ExportFormat:
        if self.all_check.isChecked():
            return ExportFormat.ALL
        fmt = ExportFormat(0)
        if self.gametes_check.isChecked():
            fmt |= ExportFormat.GAMETES
        if self.recomb_check.isChecked():
            fmt |= ExportFormat.RECOMBINATION
        if self.coverage_check.isChecked():
            fmt |= ExportFormat.COVERAGE
        if self.begfin_check.isChecked():
            fmt |= ExportFormat.BEGFIN
        if self.raw_check.isChecked():
            fmt |= ExportFormat.RAW
        return fmt

    def _validate_and_accept(self) -> None:
        path = Path(self.results_edit.text().strip())
        if not path.is_file():
            QMessageBox.warning(self, "Missing file", "Select a valid results file.")
            return

        selected = [item.text() for item in self.sample_list.selectedItems()]
        if not selected:
            QMessageBox.warning(self, "No samples", "Select at least one sample.")
            return

        out_dir = Path(self.output_edit.text().strip())
        if not out_dir:
            QMessageBox.warning(self, "Output folder", "Select an output folder.")
            return

        prefix = self.prefix_edit.text().strip() or path.stem
        formats = self._selected_formats()
        if not formats:
            QMessageBox.warning(self, "No format", "Select at least one export format.")
            return

        self._options = SaveResultsOptions(
            results_path=path,
            output_dir=out_dir,
            file_prefix=prefix,
            sra_ids=selected,
            formats=formats,
        )
        self.accept()

    def options(self) -> SaveResultsOptions:
        return self._options
