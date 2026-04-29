#!/usr/bin/env python3
"""
2_plot_s_recomb_rate_with_regions.py

Code structure was not fully reviewed, but the resulting figure was validated against a parallel reviewed version on 2026-04-29 by Luis Daniel González Vázquez.
Code adapted from the original version with AI assistance for output formatting.

Summary: Draws the recombination profile over Spike with functional region annotations from the results HTML.

Recommended execution order in this pipeline:
0 -> 1 -> 2 -> 3 -> 4
"""

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt

SPIKE_START_1B = 21563
SPIKE_END_1B = 25384
REGION_PALETTE = [
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
    "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
    "#ccebc5",
]

# Regions in Spike AA coordinates (Wuhan-Hu-1 / NC_045512.2). https://www.nature.com/articles/s41401-020-0485-4
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


def build_spike_aa_profile_from_nt(pos_ratio, indexing_mode: str):
    """
    Build Spike AA-level profile from nucleotide-level pos_ratio.
    Uses max(codon_nt_values) for each amino-acid position.
    """
    spike_nt_positions = list(range(SPIKE_START_1B, SPIKE_END_1B + 1))
    spike_nt_values = [value_at_position_1based(pos_ratio, p, indexing_mode) for p in spike_nt_positions]
    aa_values = []
    aa_positions = []
    aa_idx = 1
    for i in range(0, len(spike_nt_values), 3):
        codon_vals = spike_nt_values[i : i + 3]
        if not codon_vals:
            continue
        aa_positions.append(aa_idx)
        aa_values.append(max(codon_vals))
        aa_idx += 1
    return aa_positions, aa_values


def load_pos_ratio(source_path: Path):
    text = source_path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"const DATA=(\{.*?\});</script>", text, re.S)
    if not m:
        raise RuntimeError("Could not find DATA block in source file")
    data = json.loads(m.group(1))
    arr = data.get("pos_ratio")
    if not isinstance(arr, list):
        raise RuntimeError("DATA.pos_ratio is not a list")
    if len(arr) < SPIKE_END_1B:
        raise RuntimeError(
            f"DATA.pos_ratio too short ({len(arr)}), expected at least {SPIKE_END_1B}"
        )
    return data, arr


def detect_indexing_mode(data: dict, pos_ratio: list):
    """
    Detect whether HTML stores position values as 0-based or 1-based indexing.

    Strategy:
    1) If an explicit positions array exists, inspect its first value (0 or 1).
    2) Otherwise infer from array lengths:
       - len(pos_ratio) == genome_len     -> values indexed by (position_1based - 1), i.e. "0-based array".
       - len(pos_ratio) == genome_len + 1 -> values indexed by position_1based directly, i.e. "1-based array".
    """
    for key in ("positions", "position", "genome_pos", "pos", "x", "nt_pos"):
        vals = data.get(key)
        if isinstance(vals, list) and vals:
            first = vals[0]
            if first == 0:
                return "0-based"
            if first == 1:
                return "1-based"

    genome_len = data.get("genome_len")
    if isinstance(genome_len, int):
        if len(pos_ratio) == genome_len:
            return "0-based"
        if len(pos_ratio) == genome_len + 1:
            return "1-based"

    raise RuntimeError(
        "Could not infer indexing mode from HTML (neither explicit positions nor consistent lengths)."
    )


def value_at_position_1based(pos_ratio, pos1, indexing_mode: str):
    if indexing_mode == "0-based":
        idx = pos1 - 1
    elif indexing_mode == "1-based":
        idx = pos1
    else:
        raise ValueError(f"Unknown indexing mode: {indexing_mode}")
    if idx < 0 or idx >= len(pos_ratio):
        raise IndexError(f"Position out of range: pos1={pos1}, idx={idx}, len={len(pos_ratio)}")
    return float(pos_ratio[idx])


def write_profile_tsv(out_path: Path, positions_1based, values):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["position_1based", "recomb_rate"])
        for p, v in zip(positions_1based, values):
            w.writerow([p, f"{float(v):.12g}"])


def main():
    ap = argparse.ArgumentParser(description="Plot recombination rate across Spike with functional regions")
    ap.add_argument("--source", default="Recombination_SC2_9feb.html")
    ap.add_argument("--out", default="s_recomb_rate_with_regions.png")
    ap.add_argument(
        "--out-genome-tsv",
        default="recomb_rate_profile_whole_genome_1based.tsv",
        help="Output TSV with per-position values for whole genome (1-based positions).",
    )
    ap.add_argument(
        "--out-spike-tsv",
        default="recomb_rate_profile_spike_1based.tsv",
        help="Output TSV with per-position values for Spike region (1-based positions).",
    )
    args = ap.parse_args()

    source_path = Path(args.source)
    out_figure = Path(args.out)
    data, pos_ratio = load_pos_ratio(source_path)
    indexing_mode = detect_indexing_mode(data, pos_ratio)

    # Profiles exported explicitly in 1-based coordinates.
    if indexing_mode == "0-based":
        genome_positions = list(range(1, len(pos_ratio) + 1))
    else:
        # 1-based array (index 0 unused or dummy)
        genome_positions = list(range(1, len(pos_ratio)))
    genome_values = [value_at_position_1based(pos_ratio, p, indexing_mode) for p in genome_positions]
    spike_positions = list(range(SPIKE_START_1B, SPIKE_END_1B + 1))
    spike_values = [value_at_position_1based(pos_ratio, p, indexing_mode) for p in spike_positions]

    write_profile_tsv(Path(args.out_genome_tsv), genome_positions, genome_values)
    write_profile_tsv(Path(args.out_spike_tsv), spike_positions, spike_values)

    xs, ys = build_spike_aa_profile_from_nt(pos_ratio, indexing_mode)
    ymax = max(ys) if ys else 1.0

    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.fill_between(xs, ys, color="#b34700", alpha=0.28, linewidth=0)
    ax.plot(xs, ys, color="#7a2f00", linewidth=1.0)
    if xs:
        ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(0, ymax * 1.15 if ymax > 0 else 1.0)
    ax.set_xlabel("Spike amino-acid position")
    ax.set_ylabel("Recombination rate")
    ax.grid(alpha=0.18, axis="y")

    y0 = ymax * 1.02 if ymax > 0 else 0.02
    y1 = ymax * 1.10 if ymax > 0 else 0.10
    ylim_top = ax.get_ylim()[1] or 1.0
    ymin_frac = y0 / ylim_top
    ymax_frac = y1 / ylim_top
    for i, (name, aa_s, aa_e) in enumerate(SPIKE_FUNCTIONAL_REGIONS_AA):
        ns = aa_s
        ne = aa_e
        ax.axvspan(
            ns,
            ne,
            ymin=ymin_frac,
            ymax=ymax_frac,
            color=REGION_PALETTE[i % len(REGION_PALETTE)],
            alpha=0.85,
            lw=0,
        )
        ax.text((ns + ne) / 2, y1, name, ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_figure, dpi=260)
    plt.close(fig)
    print(f"[OK] Indexing detected in HTML DATA.pos_ratio: {indexing_mode}")
    print(f"[OK] Figure generated: {out_figure}")
    print(f"[OK] Whole-genome TSV (1-based): {args.out_genome_tsv}")
    print(f"[OK] Spike TSV (1-based): {args.out_spike_tsv}")


if __name__ == "__main__":
    main()
