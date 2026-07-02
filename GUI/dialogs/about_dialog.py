"""About dialog for SRARec Desktop."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QVBoxLayout,
)

from GUI.config import APP_NAME, APP_VERSION, GITHUB_URL, PAPER_STATUS


class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"About {APP_NAME}")
        self.setMinimumWidth(420)

        layout = QVBoxLayout(self)

        title = QLabel(f"<h2>{APP_NAME}</h2>")
        title.setTextFormat(Qt.RichText)
        layout.addWidget(title)

        version = QLabel(f"Version {APP_VERSION}")
        layout.addWidget(version)

        desc = QLabel(
            "Read-level recombination analysis from preprocessed BAM files. "
            "Part of the SRARec project."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        github = QLabel(
            f'<a href="{GITHUB_URL}">GitHub repository</a>'
        )
        github.setTextFormat(Qt.RichText)
        github.setOpenExternalLinks(True)
        layout.addWidget(github)

        paper = QLabel(f"<b>Publication:</b> {PAPER_STATUS}")
        paper.setWordWrap(True)
        paper.setTextFormat(Qt.RichText)
        layout.addWidget(paper)

        license_lbl = QLabel("Licensed under the MIT License.")
        layout.addWidget(license_lbl)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
