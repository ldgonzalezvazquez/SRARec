"""Console log panel."""

from __future__ import annotations

from PySide6.QtGui import QColor, QTextCharFormat, QTextCursor
from PySide6.QtWidgets import QGroupBox, QTextEdit, QVBoxLayout


class ConsolePanel(QGroupBox):
    def __init__(self, parent=None):
        super().__init__("Console", parent)
        layout = QVBoxLayout(self)
        self.text = QTextEdit()
        self.text.setReadOnly(True)
        self.text.setStyleSheet(
            "QTextEdit { background-color: #1e1e1e; color: #f0f0f0; font-family: Menlo, monospace; }"
        )
        layout.addWidget(self.text)

    def append(self, message: str) -> None:
        cursor = self.text.textCursor()
        cursor.movePosition(QTextCursor.End)

        fmt = QTextCharFormat()
        if message.startswith("[ERROR]") or "Error" in message[:20]:
            fmt.setForeground(QColor("#ff6b6b"))
        elif message.startswith("[WARN]"):
            fmt.setForeground(QColor("#ffd166"))
        elif message.startswith("[OK]") or message.startswith("[STEP]"):
            fmt.setForeground(QColor("#7bd88f"))
        elif message.startswith("[INFO]") or message.startswith("[ANALYSIS]"):
            fmt.setForeground(QColor("#7ec8e3"))
        else:
            fmt.setForeground(QColor("#f0f0f0"))

        cursor.insertText(message.rstrip() + "\n", fmt)
        self.text.setTextCursor(cursor)
        self.text.ensureCursorVisible()

    def clear(self) -> None:
        self.text.clear()
