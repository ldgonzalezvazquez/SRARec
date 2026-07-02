"""Wrappers around generate_graphs.py for the desktop GUI."""

from __future__ import annotations

import importlib.util
import os
import sys
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path
from types import ModuleType

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from GUI.paths import generate_graphs_path

plt.ioff()


@lru_cache(maxsize=1)
def _load_graphs_module() -> ModuleType:
    path = generate_graphs_path()
    if not path.is_file():
        raise FileNotFoundError(f"generate_graphs.py not found: {path}")

    name = "generate_graphs"
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@contextmanager
def _use_output_workdir(workdir: Path):
    """Run graph generation with cwd set to a writable output directory."""
    workdir.mkdir(parents=True, exist_ok=True)
    prev = Path.cwd()
    os.chdir(workdir)
    try:
        yield workdir
    finally:
        os.chdir(prev)


@contextmanager
def _patch_plt_for_embedding():
    """Block native plot windows and build figures in memory only."""
    _orig_close = plt.close
    _orig_show = plt.show
    _orig_subplots = plt.subplots
    created: list[Figure] = []

    def _no_close(_fig=None):
        return None

    def _no_show(*_args, **_kwargs):
        return None

    def _agg_subplots(*_args, figsize=None, **kwargs):
        for key in (
            "nrows", "ncols", "num", "sharex", "sharey",
            "squeeze", "subplot_kw", "gridspec_kw",
        ):
            kwargs.pop(key, None)
        size = figsize if figsize is not None else (12, 6)
        fig = Figure(figsize=size)
        ax = fig.add_subplot(111)
        FigureCanvasAgg(fig)
        created.append(fig)
        return fig, ax

    plt.close = _no_close
    plt.show = _no_show
    plt.subplots = _agg_subplots
    try:
        yield created
    finally:
        plt.close = _orig_close
        plt.show = _orig_show
        plt.subplots = _orig_subplots


def build_profile_figure(
    results_path: Path,
    annotation_path: Path | None = None,
    *,
    alpha: float = 1.0,
    sra_ids: list[str] | None = None,
) -> Figure:
    """Build a matplotlib Figure for embedding in the Qt graph panel."""
    graphs = _load_graphs_module()
    results_path = Path(results_path)

    from GUI.core.results_filter import filter_lines
    import tempfile

    if sra_ids:
        lines = filter_lines(results_path, sra_ids)
        if not lines:
            raise ValueError("No result lines match the selected samples.")
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8"
        )
        tmp.write("\n".join(lines) + "\n")
        tmp.close()
        work_path = Path(tmp.name)
    else:
        work_path = results_path

    try:
        genome_len = None
        genes_df = None
        if annotation_path and Path(annotation_path).is_file():
            genes_df, gb_len = graphs._load_genes(str(annotation_path), genome_len=None)
            if gb_len is not None:
                genome_len = gb_len

        events = graphs._load_events(str(work_path))
        positions, rho = graphs.compute_recombination_profile(
            events,
            genome_len=genome_len,
            alpha=alpha,
            exclusive_positions=False,
            aggregate_by_patient=False,
        )

        coverage_for_plot = None
        begfin_for_plot = None
        try:
            cov = graphs._compute_coverage_from_begfin(str(work_path), genome_len=genome_len)
            if cov is not None and cov.size:
                coverage_for_plot = cov / 1e6
        except Exception:
            pass
        try:
            conc = graphs._compute_begfin_concentration(str(work_path), genome_len=genome_len)
            if conc is not None and conc.size:
                begfin_for_plot = conc / 1e6
        except Exception:
            pass

        num_sras = graphs._count_unique_sras(str(work_path))
        title = f"Genome-wide recombination ratio ({num_sras} sample(s))"

        with _patch_plt_for_embedding() as created:
            graphs._generate_matplotlib_figure(
                positions,
                rho,
                coverage=coverage_for_plot,
                begfin_concentration=begfin_for_plot,
                genes_df=genes_df if genes_df is not None and not genes_df.empty else None,
                output_png=None,
                title=title,
            )
            if not created:
                raise RuntimeError("Matplotlib figure was not created.")
            return created[-1]
    finally:
        if sra_ids and work_path != results_path and work_path.is_file():
            work_path.unlink(missing_ok=True)


def generate_profile_html(
    results_path: Path,
    annotation_path: Path | None,
    output_html: Path,
    *,
    alpha: float = 1.0,
    sra_ids: list[str] | None = None,
    output_png: Path | None = None,
) -> Path:
    """Generate interactive Plotly HTML and return its path."""
    from GUI.core.results_filter import filter_lines
    import shutil
    import tempfile

    graphs = _load_graphs_module()
    results_path = Path(results_path)
    ann = str(annotation_path) if annotation_path and Path(annotation_path).is_file() else None
    output_html = Path(output_html).resolve()
    workdir = output_html.parent
    workdir.mkdir(parents=True, exist_ok=True)

    if sra_ids:
        lines = filter_lines(results_path, sra_ids)
        if not lines:
            raise ValueError("No result lines match the selected samples.")
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8"
        )
        tmp.write("\n".join(lines) + "\n")
        tmp.close()
        input_path = Path(tmp.name)
    else:
        input_path = results_path.resolve()

    try:
        with _patch_plt_for_embedding():
            with _use_output_workdir(workdir):
                out = graphs.generate_recombination_html(
                    results_path=str(input_path),
                    output_html=output_html.name,
                    alpha=alpha,
                    annotation_path=ann,
                )

        out_path = Path(out)
        if not out_path.is_absolute():
            out_path = workdir / out_path
        out_path = out_path.resolve()

        if output_png is not None:
            auto_png = out_path.with_suffix(".png")
            if auto_png.is_file():
                output_png = Path(output_png)
                output_png.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(auto_png, output_png)

        return out_path
    finally:
        if sra_ids and input_path != results_path.resolve() and input_path.is_file():
            input_path.unlink(missing_ok=True)
