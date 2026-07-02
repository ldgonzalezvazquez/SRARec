"""BAM input requirements help dialog."""

from __future__ import annotations

from PySide6.QtWidgets import QDialog, QDialogButtonBox, QTextEdit, QVBoxLayout


BAM_REQUIREMENTS_TEXT = """BAM input requirements for SRARec Desktop
========================================

Your BAM file must be ready for read-level analysis:

1. Coordinate-sorted
   Header must contain: @HD SO:coordinate

2. Indexed
   A .bai index file must exist next to the BAM (the app can create it if missing).

3. Deduplicated (recommended)
   PCR/optical duplicates should be removed before analysis.

4. Paired-end data
   Mate information should be correct (fixmate applied during preprocessing).

5. Mapped reads
   The file should contain reads aligned to the same reference you select in the GUI.

Preprocessing (trimming, mapping, deduplication) must be done before using this app.
Use the full SRARec1.3.py CLI or your preferred pipeline to produce the BAM.
"""


class BamRequirementsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("BAM requirements")
        self.setMinimumSize(480, 360)

        layout = QVBoxLayout(self)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(BAM_REQUIREMENTS_TEXT)
        layout.addWidget(text)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
