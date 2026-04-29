#!/usr/bin/env python3
"""
3_plot_global_circos_annotated_zoom.py

Reviewed by Luis Daniel González Vázquez on 2026-04-29.
The coordinate drawing details were not reviewed exhaustively, but output behavior was validated, with focus on ensuring arcs match real gamete and reconstruction evidence.

Summary: Generates global and per-patient circos plots using gamete results and reconstructed data to color arcs (gray/green/blue/orange).

Recommended execution order in this pipeline:
0 -> 1 -> 2 -> 3 -> 4
"""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


SPIKE_START_1B = 21563
SPIKE_END_1B = 25384
SPIKE_AA_LEN = 1273
COLOR_BG_UNRESOLVED = "#8a9096"
COLOR_RESOLVED_NONSYN = "#1f77b4"
COLOR_RESOLVED_NONSYN_AMB = "#7a1f2a"
COLOR_RESOLVED_OTHER = "#2ca02c"
COLOR_RESOLVED_INDEL = "#e7a14a"

SPIKE_FUNCTIONAL_REGIONS = [
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


def nt_to_spike_aa_index(pos_nt_0b: int):
    """
    Maps a 0-based nt position (results-slow) to a 1-based AA index in Spike.
    If the polymorphism range supporting a recombination signal is outside Spike,
    that evidence is discarded.
    """
    pos_nt_1b = pos_nt_0b + 1
    if pos_nt_1b < SPIKE_START_1B or pos_nt_1b > SPIKE_END_1B:
        return None
    return ((pos_nt_1b - SPIKE_START_1B) // 3) + 1


def theta_from_idx(idx):
    frac = (idx - 1) / SPIKE_AA_LEN
    return math.pi / 2 - 2 * math.pi * frac


def pol2cart(r, theta):
    return r * math.cos(theta), r * math.sin(theta)


def draw_chord(ax, th1, th2, color, lw=1.0, alpha=0.6, ctrl_r=0.20):
    x1, y1 = pol2cart(1.0, th1)
    x2, y2 = pol2cart(1.0, th2)
    cx1, cy1 = pol2cart(ctrl_r, th1)
    cx2, cy2 = pol2cart(ctrl_r, th2)
    verts = [(x1, y1), (cx1, cy1), (cx2, cy2), (x2, y2)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    patch = PathPatch(MplPath(verts, codes), facecolor="none", edgecolor=color, lw=lw, alpha=alpha)
    ax.add_patch(patch)


def draw_region_track(ax):
    palette = [
        "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
        "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
        "#ccebc5", "#ffed6f",
    ]
    r_inner = 1.22
    r_outer = 1.34
    for i, (name, start, end) in enumerate(SPIKE_FUNCTIONAL_REGIONS):
        n = max(12, int((end - start + 1) / 3))
        thetas = [theta_from_idx(start + (end - start) * j / (n - 1)) for j in range(n)]
        outer = [pol2cart(r_outer, th) for th in thetas]
        inner = [pol2cart(r_inner, th) for th in reversed(thetas)]
        poly = outer + inner
        ax.fill([p[0] for p in poly], [p[1] for p in poly], color=palette[i % len(palette)], alpha=0.9, linewidth=0)
        th_m = theta_from_idx((start + end) / 2.0)
        xt, yt = pol2cart(1.42, th_m)
        ax.text(xt, yt, name, fontsize=12, ha="center", va="center")


def load_spike_ratio_profile_from_html(html_path: Path):
    text = html_path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"const DATA=(\{.*?\});</script>", text, re.S)
    if not m:
        return None
    data = json.loads(m.group(1))
    pos_ratio = data.get("pos_ratio")
    if not isinstance(pos_ratio, list) or len(pos_ratio) < SPIKE_END_1B:
        return None
    # Index detection:
    # - if positions array starts at 0 -> 0-based
    # - if it starts at 1 -> 1-based
    # - if it does not start at 0/1 (or does not exist) -> default to 1-based
    indexing_mode = "1-based"
    for key in ("positions", "position", "genome_pos", "pos", "x", "nt_pos"):
        vals = data.get(key)
        if isinstance(vals, list) and vals:
            first = vals[0]
            if first == 0:
                indexing_mode = "0-based"
            elif first == 1:
                indexing_mode = "1-based"
            else:
                indexing_mode = "1-based"
            break
    aa_profile = []
    for aa_idx in range(1, SPIKE_AA_LEN + 1):
        codon_start_1b = SPIKE_START_1B + (aa_idx - 1) * 3
        vals = []
        for nt in (codon_start_1b, codon_start_1b + 1, codon_start_1b + 2):
            idx = nt - 1 if indexing_mode == "0-based" else nt  # nt is in biological coordinates; in 0-based arrays subtract 1.
            if 0 <= idx < len(pos_ratio):
                vals.append(float(pos_ratio[idx]))
        aa_profile.append(max(vals) if vals else 0.0)
    return aa_profile


def draw_recomb_rate_track(ax, aa_profile):
    if not aa_profile:
        return
    vmax = max(aa_profile)
    if vmax <= 0:
        return
    r_base = 1.01
    amp = 0.18
    thetas = [theta_from_idx(i) for i in range(1, SPIKE_AA_LEN + 1)]
    outer = [pol2cart(r_base + amp * (v / vmax), th) for th, v in zip(thetas, aa_profile)]
    inner = [pol2cart(r_base, th) for th in reversed(thetas)]
    poly = outer + inner
    ax.fill([p[0] for p in poly], [p[1] for p in poly], color="#b34700", alpha=0.45, linewidth=0)
    ax.plot([p[0] for p in outer], [p[1] for p in outer], color="#7a2f00", linewidth=1.0, alpha=0.95)

def sample_belongs_to_patient(sample: str, patient_filter: str):
    if not patient_filter:
        return True
    return sample.startswith(f"{patient_filter}-")

def parse_results_line(line: str):
    return dict(part.split(":", 1) for part in line.split("\t") if ":" in part)


def detail_indel_status(detail: str):
    """Return (has_indel, has_non3n_indel) from results Detail field."""
    txt = detail or ""
    matches = list(re.finditer(r"([+-])(\d+)", txt))
    if matches:
        lens = [int(m.group(2)) for m in matches]
        has_non3 = any((n % 3) != 0 for n in lens)
        return True, has_non3
    # Fallback for formats without explicit length
    if "+" in txt or "-" in txt:
        return True, True
    return False, False


def maybe_prefix_output(path: Path, enabled: bool, prefix: str = "only3indels_"):
    if not enabled:
        return path
    if path.name.startswith(prefix):
        return path
    return path.with_name(prefix + path.name)


def pair_key(i1: int, i2: int):
    return (i1, i2) if i1 <= i2 else (i2, i1)

def aa_pair_key_from_fields(fields):
    a = fields.get("PosA")
    b = fields.get("PosB")
    if not a or not b:
        return None
    try:
        p1 = int(a)
        p2 = int(b)
    except ValueError:
        return None
    aa1 = nt_to_spike_aa_index(p1)
    aa2 = nt_to_spike_aa_index(p2)
    if aa1 is None or aa2 is None:
        return None
    return tuple(sorted((aa1, aa2)))


def load_resolved_event_keys(path: Path):
    if not path.exists():
        return set()
    keys = set()
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            s, a, b = r.get("sample", ""), r.get("posA", ""), r.get("posB", "")
            if s and a and b:
                keys.add((s, a, b))
    return keys


def load_pair_both_nonsyn_flags(path: Path, patient_filter: str = None):
    if not path.exists(): 
        return {}, {}
    aa_sets_by_patient = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            patient = (r.get("patient") or "").strip()
            if patient_filter and patient != patient_filter:
                continue
            for gk in ("gamete1_aa", "gamete2_aa", "gamete3_aa", "gamete4_aa"):
                sig = r.get(gk, "")
                parts = sig.split("|")
                if len(parts) < 2:
                    continue
                try:
                    i1s, a1 = parts[0].split(":", 1)
                    i2s, a2 = parts[1].split(":", 1)
                    i1 = int(i1s)
                    i2 = int(i2s)
                except Exception:
                    continue
                key = pair_key(i1, i2)
                if a1 != "*":
                    aa_sets_by_patient[patient][key][i1].add(a1)
                if a2 != "*":
                    aa_sets_by_patient[patient][key][i2].add(a2)
    
    #AQUIUIIII
    
    flags = {}
    amb_flags = {}
    has_clean_resolution = defaultdict(bool)
    has_amb_resolution = defaultdict(bool)
    for _patient, pairs_map in aa_sets_by_patient.items():
        for key, pos_map in pairs_map.items():
            p1, p2 = key
            s1 = pos_map.get(p1, set())
            s2 = pos_map.get(p2, set())
            pair_both_nonsyn = len(s1) >= 2 and len(s2) >= 2
            if not pair_both_nonsyn:
                continue
            flags[key] = True
            if ("amb" in s1) or ("amb" in s2):
                has_amb_resolution[key] = True
            else:
                has_clean_resolution[key] = True
    for key in flags:
        # amb only dominates when no better (non-amb) resolution exists in any patient.
        amb_flags[key] = has_amb_resolution.get(key, False) and not has_clean_resolution.get(key, False)
    return flags, amb_flags


def iter_result_pairs(
    results_path: Path,
    patient_filter: str = None,
    only_3n_indels: bool = False,
):
    """Yield (sample, posA, posB, aa_pair_key, has_indel) from results file."""
    with results_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = parse_results_line(line)
            sample = fields.get("Sra")
            if not sample:
                continue
            if patient_filter and not sample_belongs_to_patient(sample, patient_filter):
                continue
            has_indel, has_non3 = detail_indel_status(fields.get("Detail", ""))
            if only_3n_indels and has_non3:
                continue
            key = aa_pair_key_from_fields(fields)
            if key is None:
                continue
            a = fields.get("PosA")
            b = fields.get("PosB")
            if not a or not b:
                continue
            yield sample, a, b, key, has_indel


def collect_aggregated_pairs(
    results_path: Path,
    resolved_keys,
    patient_filter: str = None,
    only_3n_indels: bool = False,
):
    agg = defaultdict(lambda: {"resolved": 0, "unresolved": 0})
    for sample, a, b, key, _has_indel in iter_result_pairs(
        results_path,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    ):
        is_resolved = (sample, a, b) in resolved_keys
        if is_resolved:
            agg[key]["resolved"] += 1
        else:
            agg[key]["unresolved"] += 1
    return agg

def collect_indel_pair_flags(
    results_path: Path,
    patient_filter: str = None,
    only_3n_indels: bool = False,
):
    """Mark AA pairs where event detail includes indels."""
    flags = defaultdict(bool)
    for _sample, _a, _b, key, has_indel in iter_result_pairs(
        results_path,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    ):
        if has_indel:
            flags[key] = True
    return flags

def collect_all_resolved_are_indel_flags(
    results_path: Path,
    resolved_keys,
    patient_filter: str = None,
    only_3n_indels: bool = False,
):
    stats = defaultdict(lambda: {"resolved": 0, "resolved_indel": 0})
    for sample, a, b, key, has_indel in iter_result_pairs(
        results_path,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    ):
        if (sample, a, b) not in resolved_keys:
            continue
        stats[key]["resolved"] += 1
        if has_indel:
            stats[key]["resolved_indel"] += 1
    out = {}
    for key, d in stats.items():
        out[key] = d["resolved"] > 0 and d["resolved"] == d["resolved_indel"]
    return out

def make_global_annotated(
    results_path,
    reconstructed_path,
    html_path,
    out_png,
    patient_filter=None,
    only_3n_indels: bool = False,
):
    resolved_keys = load_resolved_event_keys(reconstructed_path)
    agg = collect_aggregated_pairs(
        results_path,
        resolved_keys,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    )
    indel_flags = collect_indel_pair_flags(
        results_path,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    )
    all_resolved_are_indel = collect_all_resolved_are_indel_flags(
        results_path,
        resolved_keys,
        patient_filter=patient_filter,
        only_3n_indels=only_3n_indels,
    )
    both_nonsyn, nonsyn_has_amb = load_pair_both_nonsyn_flags(reconstructed_path, patient_filter=patient_filter)
    aa_profile = load_spike_ratio_profile_from_html(html_path)

    fig, ax = plt.subplots(figsize=(10.8, 10.8))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, lw=1.2, color="#666"))
    draw_recomb_rate_track(ax, aa_profile)
    draw_region_track(ax)

    # ticks
    tick_positions = [1] + list(range(100, SPIKE_AA_LEN + 1, 100))
    for idx in tick_positions:
        th = theta_from_idx(idx)
        x0, y0 = pol2cart(1.0, th)
        x1, y1 = pol2cart(1.04, th)
        ax.plot([x0, x1], [y0, y1], color="#999", lw=0.6)
        xt, yt = pol2cart(1.10, th)
        ax.text(xt, yt, str(idx), fontsize=10.5, ha="center", va="center", color="#666")

    items = []
    for (a, b), d in agg.items():
        total = d["resolved"] + d["unresolved"]
        items.append((a, b, d["resolved"], d["unresolved"], total))
    items.sort(key=lambda x: x[4], reverse=True)
    chord_lw = 1.8
    chord_alpha = 0.78
    chord_ctrl_r = 0.22

    # draw all pairs
    for a, b, nres, _nunres, _total in items:
        th1 = theta_from_idx(a)
        th2 = theta_from_idx(b)
        is_indel_pair = indel_flags.get((a, b), False)
        if nres > 0 and both_nonsyn.get((a, b), False):
            if all_resolved_are_indel.get((a, b), False):
                color = COLOR_RESOLVED_INDEL
            else:
                color = COLOR_RESOLVED_NONSYN_AMB if nonsyn_has_amb.get((a, b), False) else COLOR_RESOLVED_NONSYN
        elif is_indel_pair:
            color = COLOR_RESOLVED_INDEL
        elif nres > 0:
            color = COLOR_RESOLVED_OTHER
        else:
            color = COLOR_BG_UNRESOLVED
        draw_chord(ax, th1, th2, color=color, lw=chord_lw, alpha=chord_alpha, ctrl_r=chord_ctrl_r)

    # Clean bounds for upper circos panel
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    fig.tight_layout(pad=0.15)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Global Spike circos figure.")
    ap.add_argument("--results", default="results_gametes.txt")
    ap.add_argument("--reconstructed", default="reconstructed_gamete_aas_0.75.tsv")
    ap.add_argument("--html", default="Recombination_SC2_9feb.html")
    ap.add_argument("--out-global", default="global_circos_all_patients_annotated_nodes.png")
    ap.add_argument("--patient", default="", help="Optional patient ID filter, e.g., 024-A")
    ap.add_argument(
        "--only-3n-indels",
        action="store_true",
        help="Ignore events containing indels whose length is not multiple of 3; prefix output with only3indels_.",
    )
    args = ap.parse_args()

    out_png = maybe_prefix_output(Path(args.out_global), enabled=args.only_3n_indels)

    make_global_annotated(
        Path(args.results),
        Path(args.reconstructed),
        Path(args.html),
        out_png,
        patient_filter=(args.patient.strip() or None),
        only_3n_indels=args.only_3n_indels,
    )
    print(f"[OK] Annotated global plot: {out_png}")


if __name__ == "__main__":
    main()
