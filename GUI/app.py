#!/usr/bin/env python3
"""SRARec Desktop — macOS GUI entry point."""

from __future__ import annotations

import multiprocessing

multiprocessing.freeze_support()

import os
import sys
from pathlib import Path

# Ensure repository root is on sys.path (dev and PyInstaller bundle).
if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
    _ROOT = Path(sys._MEIPASS)
else:
    _ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# PySide6 6.11 on macOS (especially x86_64 Python under Rosetta) may not locate
# platform plugins unless QT_PLUGIN_PATH is set before importing QtWidgets.
import PySide6

_qt_plugins = Path(PySide6.__file__).resolve().parent / "Qt" / "plugins"
if _qt_plugins.is_dir():
    os.environ["QT_PLUGIN_PATH"] = str(_qt_plugins)

# Agg backend: headless figures (no extra plot windows) and no Qt conflict with PySide6.
import matplotlib

matplotlib.use("Agg")

from PySide6.QtGui import QIcon
from PySide6.QtNetwork import QLocalServer, QLocalSocket
from PySide6.QtWidgets import QApplication

from GUI.config import APP_NAME
from GUI.paths import app_icon_path
from GUI.widgets.splash_screen import create_splash_screen

_INSTANCE_KEY = "org.srarec.desktop.single_instance"


def _another_instance_running() -> bool:
    sock = QLocalSocket()
    sock.connectToServer(_INSTANCE_KEY)
    if sock.waitForConnected(300):
        sock.close()
        return True
    return False


def _claim_single_instance() -> QLocalServer | None:
    server = QLocalServer()
    if server.listen(_INSTANCE_KEY):
        return server
    server.removeServer(_INSTANCE_KEY)
    if server.listen(_INSTANCE_KEY):
        return server
    return None


def main() -> int:
    if _another_instance_running():
        return 0

    app = QApplication(sys.argv)
    instance_server = _claim_single_instance()

    def _release_instance_server() -> None:
        nonlocal instance_server
        if instance_server is not None:
            instance_server.close()
            instance_server = None

    app.aboutToQuit.connect(_release_instance_server)
    app.setQuitOnLastWindowClosed(True)

    app.setApplicationName(APP_NAME)
    app.setOrganizationName("SRARec")

    icon_path = app_icon_path()
    if icon_path.is_file():
        app.setWindowIcon(QIcon(str(icon_path)))

    splash = create_splash_screen()
    splash.show()
    app.processEvents()

    from GUI.widgets.main_window import MainWindow

    window = MainWindow()
    if icon_path.is_file():
        window.setWindowIcon(QIcon(str(icon_path)))

    splash.finish(window)
    window.show()
    app.processEvents()

    if instance_server is not None:
        app._instance_server = instance_server  # type: ignore[attr-defined]

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
