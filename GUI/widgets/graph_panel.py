"""Embedded matplotlib graph panel."""

from __future__ import annotations

from typing import Callable, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PySide6.QtWidgets import QGroupBox, QHBoxLayout, QPushButton, QVBoxLayout

from GUI.plotting.profile_plot import build_profile_figure
from GUI.paths import get_annotation_gb

plt.ioff()


class GraphPanel(QGroupBox):
    def __init__(self, parent=None):
        super().__init__("Recombination profile", parent)
        self._figure: Optional[plt.Figure] = None
        self._canvas: Optional[FigureCanvasQTAgg] = None
        self._toolbar: Optional[NavigationToolbar2QT] = None
        self._results_path = None
        self._reference_key = None
        self._on_open_html: Callable | None = None
        self._on_regenerate: Callable | None = None

        outer = QVBoxLayout(self)
        btn_row = QHBoxLayout()
        self.regen_btn = QPushButton("Regenerate graph")
        self.html_btn = QPushButton("Open HTML graph")
        self.html_btn.setEnabled(False)
        btn_row.addWidget(self.regen_btn)
        btn_row.addWidget(self.html_btn)
        btn_row.addStretch()
        outer.addLayout(btn_row)

        self.plot_container = QVBoxLayout()
        outer.addLayout(self.plot_container, stretch=1)

        self.regen_btn.clicked.connect(self._regenerate_clicked)
        self.html_btn.clicked.connect(self._open_html_clicked)

    def set_html_callback(self, callback) -> None:
        self._on_open_html = callback

    def set_regenerate_callback(self, callback) -> None:
        self._on_regenerate = callback

    def _clear_canvas(self) -> None:
        if self._toolbar is not None:
            self._toolbar.deleteLater()
            self._toolbar = None
        if self._canvas is not None:
            self._canvas.deleteLater()
            self._canvas = None
        if self._figure is not None:
            plt.close(self._figure)
            self._figure = None

    def show_figure(self, figure: plt.Figure) -> None:
        self._clear_canvas()
        self._figure = figure
        if figure.canvas is None:
            FigureCanvasAgg(figure)
        self._canvas = FigureCanvasQTAgg(figure)
        self._toolbar = NavigationToolbar2QT(self._canvas, self)
        self.plot_container.addWidget(self._toolbar)
        self.plot_container.addWidget(self._canvas)
        self._canvas.draw()
        self.html_btn.setEnabled(True)

    def load_from_results(
        self,
        results_path,
        reference_key: str,
        *,
        sra_ids: list[str] | None = None,
    ) -> None:
        self._results_path = results_path
        self._reference_key = reference_key
        ann = get_annotation_gb(reference_key)
        fig = build_profile_figure(results_path, ann, sra_ids=sra_ids)
        self.show_figure(fig)

    def _regenerate_clicked(self) -> None:
        if self._on_regenerate:
            self._on_regenerate()
            return
        if self._results_path and self._reference_key:
            self.load_from_results(self._results_path, self._reference_key)

    def _open_html_clicked(self) -> None:
        if self._on_open_html:
            self._on_open_html()

    def current_figure(self) -> Optional[plt.Figure]:
        return self._figure

    def reset(self) -> None:
        self._results_path = None
        self._reference_key = None
        self._clear_canvas()
        self.html_btn.setEnabled(False)
