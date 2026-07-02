# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for SRARec Desktop (macOS)
#
# Build from repo root:
#   GUI/.venv/bin/pyinstaller GUI/packaging/SRARecDesktop.spec

import sys
from pathlib import Path

import PyInstaller

from PyInstaller.utils.hooks import collect_all, collect_submodules

REPO = Path(SPECPATH).resolve().parent.parent
GUI = REPO / "GUI"
BUILD_ROOT = Path("/Users/cmelab/Desktop/SRARec-desktop-build")
_MP_RTHOOK = Path(PyInstaller.__path__[0]) / "hooks" / "rthooks" / "pyi_rth_multiprocessing.py"

block_cipher = None

datas = [
    (str(GUI / "resources"), "GUI/resources"),
    (str(GUI / "icon" / "icon.jpeg"), "GUI/icon"),
    (str(REPO / "SRARec1.3.py"), "."),
    (str(REPO / "generate_graphs.py"), "."),
]

binaries = []
hiddenimports = collect_submodules("GUI") + [
    "pysam",
    "pysam.libchtslib",
    "Bio",
    "Bio.SeqIO",
    "plotly",
    "plotly.graph_objects",
    "plotly.subplots",
    "kaleido",
    "tqdm",
]

for pkg in ("PySide6", "matplotlib", "numpy", "pandas"):
    try:
        tmp = collect_all(pkg)
        datas += tmp[0]
        binaries += tmp[1]
        hiddenimports += tmp[2]
    except Exception:
        pass

a = Analysis(
    [str(GUI / "app.py")],
    pathex=[str(REPO)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[str(_MP_RTHOOK)],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="SRARec Desktop",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="SRARec Desktop",
)

app = BUNDLE(
    coll,
    name="SRARec Desktop.app",
    icon=str(GUI / "icon" / "icon.icns"),
    bundle_identifier="org.srarec.desktop",
)
