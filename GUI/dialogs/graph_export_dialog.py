"""Dialog for regenerating or saving recombination profile graphs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import Qt
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


@dataclass
class GraphExportOptions:
    results_path: Path
    sra_ids: list[str]
    save_png: bool
    save_html: bool
    png_path: Path | None
    html_path: Path | None
    update_display: bool


class GraphExportDialog(QDialog):
    def __init__(
        self,
        parent=None,
        *,
        mode: str = "save",
        default_results: Path | None = None,
        default_output_dir: Path | None = None,
    ):
        super().__init__(parent)
        self._mode = mode
        self._default_output_dir = default_output_dir or Path.home()

        title = "Regenerate graph" if mode == "regenerate" else "Save graph"
        self.setWindowTitle(title)
        self.setMinimumWidth(520)

        layout = QVBoxLayout(self)

        intro = QLabel(
            "Select the results file and samples (BAMs) to include, "
            "then choose which graph files to write."
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

        formats_box = QGroupBox("Output formats")
        formats_layout = QVBoxLayout(formats_box)
        self.png_check = QCheckBox("PNG image")
        self.html_check = QCheckBox("Interactive HTML")
        if mode == "regenerate":
            self.png_check.setChecked(True)
            self.html_check.setChecked(True)
        else:
            self.png_check.setChecked(True)
            self.html_check.setChecked(False)
        self.png_check.toggled.connect(self._sync_path_enabled)
        self.html_check.toggled.connect(self._sync_path_enabled)
        formats_layout.addWidget(self.png_check)
        formats_layout.addWidget(self.html_check)
        layout.addWidget(formats_box)

        png_row = QHBoxLayout()
        self.png_edit = QLineEdit()
        browse_png = QPushButton("Browse...")
        browse_png.clicked.connect(self._browse_png)
        png_row.addWidget(self.png_edit, stretch=1)
        png_row.addWidget(browse_png)
        self.png_label = QLabel("PNG path:")
        layout.addWidget(self.png_label)
        layout.addLayout(png_row)

        html_row = QHBoxLayout()
        self.html_edit = QLineEdit()
        browse_html = QPushButton("Browse...")
        browse_html.clicked.connect(self._browse_html)
        html_row.addWidget(self.html_edit, stretch=1)
        html_row.addWidget(browse_html)
        self.html_label = QLabel("HTML path:")
        layout.addWidget(self.html_label)
        layout.addLayout(html_row)

        if default_results:
            self.results_edit.setText(str(default_results))
            self._reload_samples()
            self._suggest_output_paths()

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        ok_label = "Regenerate" if mode == "regenerate" else "Save"
        buttons.button(QDialogButtonBox.Ok).setText(ok_label)
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.results_edit.editingFinished.connect(self._on_results_changed)
        self._sync_path_enabled()

    def _on_results_changed(self) -> None:
        self._reload_samples()
        self._suggest_output_paths()

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
            self._suggest_output_paths()

    def _reload_samples(self) -> None:
        self.sample_list.clear()
        path = Path(self.results_edit.text().strip())
        if not path.is_file():
            return
        for sra in list_sras(path):
            self.sample_list.addItem(sra)
        self.sample_list.selectAll()

    def _suggest_output_paths(self) -> None:
        path = Path(self.results_edit.text().strip())
        if not path.is_file():
            return
        stem = path.stem
        out_dir = self._default_output_dir
        self.png_edit.setText(str(out_dir / f"{stem}_profile.png"))
        self.html_edit.setText(str(out_dir / f"{stem}_profile.html"))

    def _browse_png(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save PNG graph",
            self.png_edit.text() or str(self._default_output_dir),
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)",
        )
        if path:
            self.png_edit.setText(path)

    def _browse_html(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save HTML graph",
            self.html_edit.text() or str(self._default_output_dir),
            "HTML (*.html)",
        )
        if path:
            self.html_edit.setText(path)

    def _sync_path_enabled(self) -> None:
        png_on = self.png_check.isChecked()
        html_on = self.html_check.isChecked()
        self.png_edit.setEnabled(png_on)
        self.png_label.setEnabled(png_on)
        self.html_edit.setEnabled(html_on)
        self.html_label.setEnabled(html_on)

    def _validate_and_accept(self) -> None:
        path = Path(self.results_edit.text().strip())
        if not path.is_file():
            QMessageBox.warning(self, "Missing file", "Select a valid results file.")
            return

        selected = [item.text() for item in self.sample_list.selectedItems()]
        if not selected:
            QMessageBox.warning(self, "No samples", "Select at least one sample.")
            return

        if not self.png_check.isChecked() and not self.html_check.isChecked():
            QMessageBox.warning(self, "No format", "Select PNG and/or HTML.")
            return

        png_path = None
        html_path = None
        if self.png_check.isChecked():
            if not self.png_edit.text().strip():
                QMessageBox.warning(self, "PNG path", "Enter a path for the PNG file.")
                return
            png_path = Path(self.png_edit.text().strip())
        if self.html_check.isChecked():
            if not self.html_edit.text().strip():
                QMessageBox.warning(self, "HTML path", "Enter a path for the HTML file.")
                return
            html_path = Path(self.html_edit.text().strip())

        self._options = GraphExportOptions(
            results_path=path,
            sra_ids=selected,
            save_png=self.png_check.isChecked(),
            save_html=self.html_check.isChecked(),
            png_path=png_path,
            html_path=html_path,
            update_display=(self._mode == "regenerate"),
        )
        self.accept()

    def options(self) -> GraphExportOptions:
        return self._options
