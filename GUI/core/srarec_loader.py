"""Load SRARec1.3 as an importable module (cached)."""

from __future__ import annotations

import importlib.util
import sys
from functools import lru_cache
from types import ModuleType

from GUI.paths import srarec_script_path


@lru_cache(maxsize=1)
def load_srarec_core() -> ModuleType:
    script = srarec_script_path()
    if not script.is_file():
        raise FileNotFoundError(f"SRARec core script not found: {script}")

    module_name = "srarec_core"
    spec = importlib.util.spec_from_file_location(module_name, str(script))
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module spec from {script}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module
