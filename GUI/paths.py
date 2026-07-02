"""Resolve bundled and repository paths (dev vs frozen app)."""

from __future__ import annotations

import sys
from pathlib import Path

def _frozen_base() -> Path | None:
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS)
    return None


_FROZEN_BASE = _frozen_base()
_GUI_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _FROZEN_BASE if _FROZEN_BASE is not None else _GUI_DIR.parent


def gui_dir() -> Path:
    if _FROZEN_BASE is not None:
        candidate = _FROZEN_BASE / "GUI"
        if candidate.is_dir():
            return candidate
    return _GUI_DIR


def repo_root() -> Path:
    return _REPO_ROOT


def resource_base() -> Path:
    """Root for bundled resources (PyInstaller _MEIPASS or GUI/resources)."""
    if hasattr(sys, "_MEIPASS"):
        candidate = Path(sys._MEIPASS) / "GUI" / "resources"
        if candidate.is_dir():
            return candidate
        candidate = Path(sys._MEIPASS) / "resources"
        if candidate.is_dir():
            return candidate
    bundled = _GUI_DIR / "resources"
    if bundled.is_dir():
        return bundled
    return _GUI_DIR / "resources"


def srarec_script_path() -> Path:
  path = repo_root() / "SRARec1.3.py"
  if path.is_file():
      return path
  if hasattr(sys, "_MEIPASS"):
      alt = Path(sys._MEIPASS) / "SRARec1.3.py"
      if alt.is_file():
          return alt
  return path


def generate_graphs_path() -> Path:
    path = repo_root() / "generate_graphs.py"
    if path.is_file():
        return path
    if hasattr(sys, "_MEIPASS"):
        alt = Path(sys._MEIPASS) / "generate_graphs.py"
        if alt.is_file():
            return alt
    return path


def get_reference_fasta(reference_key: str) -> Path:
    key = reference_key.lower().replace("_", "-")
    resources = resource_base() / "reference_indexes"

    mapping = {
        "sars-cov-2": resources / "sars_cov_2" / "NC_045512.2.fna",
        "hiv-1": resources / "hiv_1" / "NC_001802.1.fna",
    }
    if key not in mapping:
        raise ValueError(f"Unknown reference key: {reference_key}")

    path = mapping[key]
    if path.is_file():
        return path

    # Fallback to repo root reference_indexes during development
    fallback = {
        "sars-cov-2": _REPO_ROOT / "reference_indexes" / "sars_cov_2" / "NC_045512.2.fna",
        "hiv-1": _REPO_ROOT / "reference_indexes" / "hiv_1" / "NC_001802.1.fna",
    }
    fb = fallback[key]
    if fb.is_file():
        return fb
    raise FileNotFoundError(f"Reference FASTA not found: {path}")


def get_annotation_gb(reference_key: str) -> Path:
    key = reference_key.lower().replace("_", "-")
    resources = resource_base() / "annotations"
    mapping = {
        "sars-cov-2": resources / "SC2.gb",
        "hiv-1": resources / "HIV-1.gb",
    }
    if key not in mapping:
        raise ValueError(f"Unknown reference key: {reference_key}")

    path = mapping[key]
    if path.is_file():
        return path

    fallback = {
        "sars-cov-2": _REPO_ROOT / "SC2.gb",
        "hiv-1": _REPO_ROOT / "HIV-1.gb",
    }
    fb = fallback[key]
    if fb.is_file():
        return fb
    raise FileNotFoundError(f"Annotation GenBank not found: {path}")


def default_output_dir() -> Path:
    return Path.home() / "SRARec_runs"


def app_icon_path() -> Path:
    """Path to the application icon (JPEG) for the window icon."""
    candidates = [
        gui_dir() / "icon" / "icon.jpeg",
        _GUI_DIR / "icon" / "icon.jpeg",
    ]
    if _FROZEN_BASE is not None:
        candidates.insert(0, _FROZEN_BASE / "GUI" / "icon" / "icon.jpeg")
    for path in candidates:
        if path.is_file():
            return path
    return _GUI_DIR / "icon" / "icon.jpeg"
