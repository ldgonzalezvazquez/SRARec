"""Input selection panel (BAM, platform, reference, output)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QVBoxLayout,
)

from GUI.config import PLATFORMS, REFERENCES
from GUI.paths import default_output_dir


class InputPanel(QGroupBox):
    def __init__(self, parent=None):
        super().__init__("Input", parent)

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("BAM files (one or more, processed in series):"))
        self.bam_list = QListWidget()
        self.bam_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.bam_list.setMaximumHeight(120)
        layout.addWidget(self.bam_list)

        bam_btn_row = QHBoxLayout()
        add_bams = QPushButton("Add BAMs...")
        add_bams.clicked.connect(self._add_bams)
        remove_bams = QPushButton("Remove selected")
        remove_bams.clicked.connect(self._remove_selected_bams)
        clear_bams = QPushButton("Clear")
        clear_bams.clicked.connect(self._clear_bams)
        bam_btn_row.addWidget(add_bams)
        bam_btn_row.addWidget(remove_bams)
        bam_btn_row.addWidget(clear_bams)
        layout.addLayout(bam_btn_row)

        self.bam_count_label = QLabel("0 BAM file(s) selected")
        self.bam_count_label.setStyleSheet("color: #666;")
        layout.addWidget(self.bam_count_label)

        # Platform
        layout.addWidget(QLabel("Sequencing platform:"))
        self.platform_combo = QComboBox()
        for label, value in PLATFORMS:
            self.platform_combo.addItem(label, value)
        layout.addWidget(self.platform_combo)

        # Reference
        ref_row = QHBoxLayout()
        layout.addWidget(QLabel("Reference:"))
        self.ref_combo = QComboBox()
        for label, value in REFERENCES:
            self.ref_combo.addItem(label, value)
        ref_row.addWidget(self.ref_combo, stretch=1)
        self.custom_ref_btn = QPushButton("Custom reference?")
        ref_row.addWidget(self.custom_ref_btn)
        layout.addLayout(ref_row)

        # Output
        out_row = QHBoxLayout()
        self.output_edit = QLineEdit(str(default_output_dir()))
        browse_out = QPushButton("Browse...")
        browse_out.clicked.connect(self._browse_output)
        out_row.addWidget(self.output_edit, stretch=1)
        out_row.addWidget(browse_out)
        layout.addWidget(QLabel("Output folder:"))
        layout.addLayout(out_row)

        self.results_hint = QLabel(
            f"Combined results: {default_output_dir() / 'results_ravi.txt'}"
        )
        self.results_hint.setWordWrap(True)
        self.results_hint.setStyleSheet("color: #666; font-size: 11px;")
        layout.addWidget(self.results_hint)
        self.output_edit.textChanged.connect(self._update_results_hint)

    def _update_results_hint(self) -> None:
        self.results_hint.setText(
            f"Combined results: {self.output_dir() / 'results_ravi.txt'}"
        )

    def _add_bams(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select BAM file(s)",
            "",
            "BAM files (*.bam);;All files (*)",
        )
        if not paths:
            return
        existing = {self.bam_list.item(i).text() for i in range(self.bam_list.count())}
        for path in paths:
            if path not in existing:
                self.bam_list.addItem(path)
        self._update_bam_count()

    def _remove_selected_bams(self) -> None:
        for item in self.bam_list.selectedItems():
            row = self.bam_list.row(item)
            self.bam_list.takeItem(row)
        self._update_bam_count()

    def _clear_bams(self) -> None:
        self.bam_list.clear()
        self._update_bam_count()

    def _update_bam_count(self) -> None:
        n = self.bam_list.count()
        self.bam_count_label.setText(
            f"{n} BAM file{'s' if n != 1 else ''} selected"
        )

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Select output folder",
            self.output_edit.text() or str(default_output_dir()),
        )
        if path:
            self.output_edit.setText(path)

    def bam_paths(self) -> list[Path]:
        paths: list[Path] = []
        for i in range(self.bam_list.count()):
            text = self.bam_list.item(i).text().strip()
            if text:
                paths.append(Path(text))
        return paths

    def bam_path(self) -> Path | None:
        """First selected BAM, for backward compatibility."""
        paths = self.bam_paths()
        return paths[0] if paths else None

    def output_dir(self) -> Path:
        text = self.output_edit.text().strip()
        return Path(text) if text else default_output_dir()

    def platform(self) -> str:
        return self.platform_combo.currentData()

    def reference_key(self) -> str:
        return self.ref_combo.currentData()

    def results_file_path(self) -> Path:
        return self.output_dir() / "results_ravi.txt"
