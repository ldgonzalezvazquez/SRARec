"""Startup splash screen for SRARec Desktop."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import QApplication, QLabel, QVBoxLayout, QWidget

from GUI.paths import app_icon_path


class SplashScreen(QWidget):
    """Frameless splash shown only during application startup."""

    def __init__(self) -> None:
        super().__init__(
            None,
            Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint | Qt.Tool,
        )
        self.setAttribute(Qt.WA_DeleteOnClose, True)

        icon_path = app_icon_path()
        pixmap = QPixmap(str(icon_path)) if icon_path.is_file() else QPixmap(256, 256)
        if pixmap.isNull():
            pixmap = QPixmap(256, 256)
            pixmap.fill(Qt.white)

        scaled = pixmap.scaled(220, 220, Qt.KeepAspectRatio, Qt.SmoothTransformation)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 20)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignCenter)

        icon_label = QLabel()
        icon_label.setPixmap(scaled)
        icon_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(icon_label)

        title = QLabel("SRARec")
        title_font = QFont()
        title_font.setPointSize(22)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        credits = QLabel("Luis Daniel González Vázquez & Darren P. Martin")
        credits.setAlignment(Qt.AlignCenter)
        credits.setWordWrap(True)
        layout.addWidget(credits)

        year = QLabel("2026")
        year.setAlignment(Qt.AlignCenter)
        layout.addWidget(year)

        self.setFixedSize(360, 380)
        self.setStyleSheet("background-color: #e8f4fc; color: #1a1a1a;")

        screen = QApplication.primaryScreen()
        if screen is not None:
            geo = screen.availableGeometry()
            self.move(
                geo.x() + (geo.width() - self.width()) // 2,
                geo.y() + (geo.height() - self.height()) // 2,
            )

    def finish(self, window) -> None:
        self.hide()
        self.close()
        self.deleteLater()
        if window is not None:
            window.raise_()
            window.activateWindow()


def create_splash_screen() -> SplashScreen:
    return SplashScreen()
