#!/usr/bin/env python3
"""
4_plot_phase_maps_from_table.py

Reviewed by Luis Daniel González Vázquez on 2026-04-29. High-level review only, but logic risk in this script is relatively low.
Figures are not publication-grade; if publication quality is needed, this should be redesigned or finalized manually. This script remains as a reproducible baseline.

Summary: Generates per-patient phase maps from phase_table_0.75.tsv (or only3indels_phase_table_0.75.tsv).

Recommended execution order in this pipeline:
0 -> 1 -> 2 -> 3 -> 4
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SPIKE_AA_LEN = 1273
SPIKE_FUNCTIONAL_REGIONS_AA = [
    ("SP", 1, 13),
    ("NTD", 14, 305),
    ("RBD", 319, 541),
    ("S1/S2", 676, 695),
    ("FP", 788, 806),
    ("HR1", 912, 984),
    ("HR2", 1163, 1213),
    ("TM", 1213, 1237),
    ("CT", 1237, 1273),
]
REGION_COLORS = [
    "#8dd3c7",
    "#ffffb3",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
]


def maybe_prefix_name(name: str, enabled: bool, prefix: str = "only3indels_"):
    if not enabled:
        return name
    if name.startswith(prefix):
        return name
    return prefix + name


def parse_side_tokens(side_text):
    tokens = []
    for tok in (side_text or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^(\d+)(.*)$", tok)
        if not m:
            continue
        pos = int(m.group(1))
        aa = (m.group(2) or "").strip()
        if not aa:
            aa = "?"
        tokens.append((pos, aa))
    # deduplicate while preserving order
    seen = set()
    out = []
    for t in tokens:
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def parse_cut_range(cut_range_text):
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", cut_range_text or "")
    if not m:
        return None
    a = int(m.group(1))
    b = int(m.group(2))
    return (min(a, b), max(a, b))


def load_rows(phase_table_path):
    by_patient = defaultdict(list)
    with phase_table_path.open() as fh:  # Read phase_table.tsv
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:  # Iterate phase_table.tsv rows
            patients = [p.strip() for p in (r.get("patients") or "").split(",") if p.strip()]
            if not patients:
                continue
            row_obj = {
                "timepoints": r.get("timepoints", ""),
                "cut_range": r.get("cut_range", ""),
                "functional_region": r.get("functional_region", "-"),
                "left_tokens": parse_side_tokens(r.get("5p_side", "")),
                "right_tokens": parse_side_tokens(r.get("3p_side", "")),
            }
            for p in patients:
                by_patient[p].append(row_obj)
    return by_patient


def compute_zoom_window(rows, pad=20, min_span=120):
    xs = []
    for row in rows:
        xs.extend(pos for pos, _ in row["left_tokens"])
        xs.extend(pos for pos, _ in row["right_tokens"])
        cut = parse_cut_range(row["cut_range"])
        if cut:
            xs.extend([cut[0], cut[1]])

    if not xs:
        return (1, SPIKE_AA_LEN)

    x_min = max(1, min(xs) - pad)
    x_max = min(SPIKE_AA_LEN, max(xs) + pad)
    span = x_max - x_min
    if span < min_span:
        center = (x_min + x_max) / 2.0
        half = min_span / 2.0
        x_min = max(1, int(center - half))
        x_max = min(SPIKE_AA_LEN, int(center + half))
    return (x_min, x_max)


def draw_patient_map(patient, rows, out_path):
    fig, ax = plt.subplots(figsize=(14, 4.2))
    # More aggressive zoom to separate positions and reduce lateral gaps
    x_min, x_max = compute_zoom_window(rows, pad=1, min_span=18)
    ax.set_xlim(x_min, x_max)
    y_rows_top = max(2.0, len(rows) * 0.8 + 0.2)
    y_plot_top = y_rows_top + 1.55
    ax.set_ylim(-1.55, y_plot_top)

    ax.set_xlabel("Spike amino-acid position")
    ax.set_title(f"{patient} - Phase map from phase_table")
    ax.set_yticks([])
    ax.tick_params(axis="y", left=False, labelleft=False)

    # Functional region band at the bottom, above the X axis
    region_y0 = -1.25
    region_y1 = -0.70
    for i, (name, a0, a1) in enumerate(SPIKE_FUNCTIONAL_REGIONS_AA):
        rs = max(a0, x_min)
        re = min(a1, x_max)
        if rs > re:
            continue
        ax.add_patch(
            Rectangle(
                (rs, region_y0),
                re - rs,
                region_y1 - region_y0,
                facecolor=REGION_COLORS[i % len(REGION_COLORS)],
                edgecolor="none",
                alpha=0.85,
                zorder=0,
            )
        )
        ax.text((rs + re) / 2.0, (region_y0 + region_y1) / 2.0, name, ha="center", va="center", fontsize=7, color="#222")

    for i, row in enumerate(rows, start=1):
        y = i * 0.8

        lt = row["left_tokens"]
        rt = row["right_tokens"]
        # Replace points with per-row labels (e.g., 454E/K)
        for pos, aa in lt:
            ax.text(
                pos,
                y,
                f"{pos}{aa}",
                color="#1f77b4",
                fontsize=6.5,
                ha="center",
                va="center",
                zorder=5,
            )
        for pos, aa in rt:
            ax.text(
                pos,
                y,
                f"{pos}{aa}",
                color="#2ca02c",
                fontsize=6.5,
                ha="center",
                va="center",
                zorder=5,
            )

    ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=250)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot per-patient phase maps from phase_table TSV.")
    ap.add_argument("--phase-table", default="phase_table_0.75.tsv")
    ap.add_argument("--out-dir", default="phase_maps_from_table")
    ap.add_argument(
        "--only-3n-indels",
        action="store_true",
        help="Use only3indels_ prefixed input/output names (expects filtered phase table).",
    )
    args = ap.parse_args()

    phase_table_arg = Path(args.phase_table)
    out_dir_arg = Path(args.out_dir)
    if args.only_3n_indels:
        pref_phase_name = maybe_prefix_name(phase_table_arg.name, enabled=True)
        candidate = phase_table_arg.with_name(pref_phase_name)
        phase_table = candidate if candidate.exists() else phase_table_arg
        out_dir = out_dir_arg.with_name(maybe_prefix_name(out_dir_arg.name, enabled=True))
    else:
        phase_table = phase_table_arg
        out_dir = out_dir_arg
    out_dir.mkdir(parents=True, exist_ok=True)

    by_patient = load_rows(phase_table)
    for patient, rows in sorted(by_patient.items()):
        out = out_dir / f"{patient}_phase_map.png"
        draw_patient_map(patient, rows, out)
        print(f"[OK] {out}")


if __name__ == "__main__":
    main()
