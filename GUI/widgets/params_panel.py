"""Analysis parameters panel."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QSpinBox,
    QVBoxLayout,
)

from GUI.config import (
    DEFAULT_MIN_GAMETE_FREQ_PCT,
    DEFAULT_MIN_GAMETE_READS,
    DEFAULT_MIN_READS,
)


class ParamsPanel(QGroupBox):
    def __init__(self, parent=None):
        super().__init__("Parameters", parent)
        form = QFormLayout(self)

        self.min_reads = QSpinBox()
        self.min_reads.setRange(1, 10000)
        self.min_reads.setValue(DEFAULT_MIN_READS)
        form.addRow("Min reads per polymorphism:", self.min_reads)

        self.min_gamete_reads = QSpinBox()
        self.min_gamete_reads.setRange(1, 10000)
        self.min_gamete_reads.setValue(DEFAULT_MIN_GAMETE_READS)
        form.addRow("Min reads per gamete:", self.min_gamete_reads)

        self.min_gamete_freq = QDoubleSpinBox()
        self.min_gamete_freq.setRange(0.0, 100.0)
        self.min_gamete_freq.setDecimals(2)
        self.min_gamete_freq.setSuffix(" %")
        self.min_gamete_freq.setValue(DEFAULT_MIN_GAMETE_FREQ_PCT)
        form.addRow("Min gamete frequency:", self.min_gamete_freq)

        self.no_indels = QCheckBox("Ignore indels when scanning polymorphisms")
        form.addRow(self.no_indels)

        self.verbose = QCheckBox("Verbose output")
        form.addRow(self.verbose)

    def min_gamete_freq_fraction(self) -> float:
        return self.min_gamete_freq.value() / 100.0


class AdvancedOptionsDialog(QDialog):
    """Placeholder for future options; verbose is on main panel."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced options")
        layout = QVBoxLayout(self)
        layout.addWidget(
            __import__("PySide6.QtWidgets", fromlist=["QLabel"]).QLabel(
                "Additional options will appear here in future releases."
            )
        )
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
