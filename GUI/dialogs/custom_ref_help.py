"""Help dialog for custom reference indexing via CLI."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QTextEdit,
    QVBoxLayout,
)


class CustomRefHelpDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Custom reference (CLI)")
        self.setMinimumWidth(520)

        layout = QVBoxLayout(self)

        intro = QLabel(
            "SRARec Desktop only includes bundled SARS-CoV-2 and HIV-1 references. "
            "To use a custom reference, index it with the full CLI from the repository root:"
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        cmd = QTextEdit()
        cmd.setReadOnly(True)
        cmd.setPlainText(
            "python3 SRARec1.3.py --index_only --custom_ref /path/to/reference.fasta\n\n"
            "# Then run analysis with the CLI (not this desktop app):\n"
            "python3 SRARec1.3.py SAMPLE_FOLDER --bam --platform ILLUMINA --virus my_virus"
        )
        cmd.setMaximumHeight(120)
        layout.addWidget(cmd)

        note = QLabel(
            "The index is stored under custom_reference_index/ in the working directory. "
            "Custom-reference analyses must be run with SRARec1.3.py on the command line."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
