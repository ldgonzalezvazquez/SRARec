#!/usr/bin/env python3
"""
1_export_phase_table_tsv.py

Reviewed by Luis Daniel González Vázquez on 2026-04-29.

Summary: Builds phasing tables (phase_table_*.tsv and only3indels variants) from reconstructed_gamete_aas_*.tsv.

Recommended execution order in this pipeline:
0 -> 1 -> 2 -> 3 -> 4
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

SPIKE_START_1B = 21563
SPIKE_END_1B = 25384
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
    Returns None when a reported polymorphic breakpoint position is outside Spike.
    """
    pos_nt_1b = pos_nt_0b + 1
    if pos_nt_1b < SPIKE_START_1B or pos_nt_1b > SPIKE_END_1B:
        return None
    return ((pos_nt_1b - SPIKE_START_1B) // 3) + 1


def parse_results_line(line: str):
    return dict(part.split(":", 1) for part in line.split("\t") if ":" in part)


def detail_indel_status(detail: str):
    txt = detail or ""
    matches = list(re.finditer(r"([+-])(\d+)", txt))
    if matches:
        lens = [int(m.group(2)) for m in matches]
        return True, any((n % 3) != 0 for n in lens)  # (has_indel, has_non3n_indel)
    if "+" in txt or "-" in txt:
        return True, True
    return False, False


def maybe_prefix_output(path: Path, enabled: bool, prefix: str = "only3indels_"):
    if not enabled:
        return path
    if path.name.startswith(prefix):
        return path
    return path.with_name(prefix + path.name)


def load_allowed_event_keys_from_results(results_path: Path, only_3n_indels: bool):
    """Build (sample, posA, posB) keys allowed by indel-length filter."""
    if not only_3n_indels:
        return None
    allowed = set()
    with results_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = parse_results_line(line)
            sample = (fields.get("Sra") or "").strip()
            a = (fields.get("PosA") or "").strip()  # 0-based positions from results_gametes.txt
            b = (fields.get("PosB") or "").strip()  # 0-based positions from results_gametes.txt
            if not sample or not a or not b:
                continue
            _has_indel, has_non3 = detail_indel_status(fields.get("Detail", ""))
            if has_non3:
                continue
            allowed.add((sample, a, b))
    return allowed


def region_for_aa(aa_idx):
    for name, start, end in SPIKE_FUNCTIONAL_REGIONS:
        if start <= aa_idx <= end:
            return name
    return "-"


def parse_gamete_sig(sig):
    parts = (sig or "").split("|")
    if len(parts) < 2:
        return None
    try:
        i1s, a1 = parts[0].split(":", 1)
        i2s, a2 = parts[1].split(":", 1)
        return int(i1s), a1, int(i2s), a2
    except Exception:
        return None


def infer_cut_ranges(intervals):
    """Split into the minimum number of non-empty intersection ranges."""
    if not intervals:
        return tuple()
    ints = sorted((min(a, b), max(a, b)) for a, b in intervals)
    ints.sort(key=lambda x: x[1])
    ranges = []
    cur_l, cur_r = ints[0]
    for a, b in ints[1:]:
        nl = max(cur_l, a)
        nr = min(cur_r, b)
        if nl <= nr:
            cur_l, cur_r = nl, nr
        else:
            ranges.append((cur_l, cur_r))
            cur_l, cur_r = a, b
    ranges.append((cur_l, cur_r))
    return tuple(ranges)

def overlaps(interval, cut_range):
    a, b = interval
    l, r = cut_range
    return not (b < l or a > r)


def parse_side_token(tok):
    i = 0
    while i < len(tok) and tok[i].isdigit():
        i += 1
    if i == 0 or i == len(tok):
        return None
    try:
        pos = int(tok[:i])
    except Exception:
        return None
    aas = [x for x in tok[i:].split("/") if x]
    return pos, set(aas)

def normalize_side_text(txt):
    """Canonicalize a side-text without merging with others."""
    pos_map = defaultdict(set)
    for tok in (txt or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        parsed = parse_side_token(tok)
        if not parsed:
            continue
        p, aas = parsed
        pos_map[p].update(aas)
    out = []
    for p in sorted(pos_map):
        opts = sorted(pos_map[p])
        if opts:
            out.append(f"{p}{'/'.join(opts)}")
    return ",".join(out)


def aa_side_token(pos, aas):
    opts = sorted({aa for aa in aas if aa and aa != "*"})
    if not opts:
        return None
    return f"{pos}{'/'.join(opts)}"

def load_reconstructed_rows(path, allowed_event_keys=None):
    rows = []
    with path.open() as fh:  # reconstructed_gamete_aas_*.tsv
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            patient = (r.get("patient") or "").strip()
            sample = (r.get("sample") or "").strip()
            if not patient or not sample:
                continue
            try:
                a = int(r.get("posA", ""))  # 0-based positions from reconstructed_gamete_aas_*.tsv
                b = int(r.get("posB", ""))
            except Exception:
                continue
            if allowed_event_keys is not None:
                if (sample, str(a), str(b)) not in allowed_event_keys:
                    continue
            aa1 = nt_to_spike_aa_index(a)
            aa2 = nt_to_spike_aa_index(b)
            if aa1 is None or aa2 is None:
                continue
            parsed = []
            for gk in ("gamete1_aa", "gamete2_aa", "gamete3_aa", "gamete4_aa"):
                p = parse_gamete_sig(r.get(gk, ""))
                if p:
                    parsed.append(p)
            if not parsed:
                continue
            rows.append(
                {
                    "patient": patient,
                    "sample": sample,
                    "timepoint": (r.get("timepoint") or "").strip(),
                    "aa_pair": tuple(sorted((aa1, aa2))),
                    "gametes": parsed,
                }
            )
    return rows


def patient_cut_groups(rows_for_patient):
    by_sample = defaultdict(list)
    for r in rows_for_patient:
        by_sample[r["sample"]].append(r)

    groups = defaultdict(list)  # cut_signature -> [sample]
    for sample, srows in by_sample.items():
        intervals = [row["aa_pair"] for row in srows]
        cut_sig = infer_cut_ranges(intervals)
        groups[cut_sig].append(sample)
    return by_sample, groups


def phase_rows_for_patient(rows_for_patient):
    by_sample, groups = patient_cut_groups(rows_for_patient)
    out_rows = []

    for cut_sig, samples in sorted(groups.items(), key=lambda x: (len(x[0]), x[0])):
        if not cut_sig:
            continue
        all_rows = [r for s in samples for r in by_sample[s]]
        timepoints = sorted({(r.get("timepoint") or "") for r in all_rows if r.get("timepoint")})

        for cut_range in cut_sig:
            left = defaultdict(set)   # pos -> {aa}
            right = defaultdict(set)  # pos -> {aa}
            for row in all_rows:
                a, b = row["aa_pair"]
                if not overlaps((a, b), cut_range):
                    continue
                for i1, aa1, i2, aa2 in row["gametes"]:
                    if i1 <= i2:
                        left[i1].add(aa1)
                        right[i2].add(aa2)
                    else:
                        left[i2].add(aa2)
                        right[i1].add(aa1)

            left_tokens = [aa_side_token(p, left[p]) for p in sorted(left)]
            right_tokens = [aa_side_token(p, right[p]) for p in sorted(right)]
            left_tokens = [t for t in left_tokens if t]
            right_tokens = [t for t in right_tokens if t]
            if not left_tokens or not right_tokens:
                continue

            left_txt = ",".join(left_tokens)
            right_txt = ",".join(right_tokens)
            cut_l, cut_r = cut_range
            cut_mid = (cut_l + cut_r) / 2.0
            out_rows.append(
                {
                    "patients": rows_for_patient[0]["patient"],
                    "timepoints": ",".join(timepoints),
                    "cut_range": f"{cut_l}-{cut_r}",
                    "functional_region": region_for_aa(int(round(cut_mid))),  # Region assigned from cut_range midpoint
                    "5p_side": left_txt,
                    "3p_side": right_txt,
                    "source": "cut_group",
                }
            )
    return out_rows


def export_one(
    reconstructed_path: Path,
    out_tsv: Path,
    patient_filter: str = None,
    allowed_event_keys=None,
):
    rows = load_reconstructed_rows(reconstructed_path, allowed_event_keys=allowed_event_keys)
    by_patient = defaultdict(list)
    for r in rows:
        by_patient[r["patient"]].append(r)

    if patient_filter:
        if patient_filter not in by_patient:
            raise ValueError(
                f"Patient filter '{patient_filter}' not found in reconstructed input."
            )
        patients = [patient_filter]
    else:
        patients = sorted(by_patient.keys())

    patient_phase_rows = []
    for p in patients:
        patient_phase_rows.extend(phase_rows_for_patient(by_patient[p]))

    # Patient mode: direct output by cut/group
    if patient_filter:
        out_rows = []
        for r in patient_phase_rows:
            out_rows.append(
                {
                    "patients": r["patients"],
                    "timepoints": r["timepoints"],
                    "cut_range": r["cut_range"],
                    "functional_region": r["functional_region"],
                    "5p_side": r["5p_side"],
                    "3p_side": r["3p_side"],
                    "source": r["source"],
                }
            )
    else:
        # Global mode: group by exact phasing (range + canonical sides)
        grouped = {}
        for r in patient_phase_rows:
            left_norm = normalize_side_text(r["5p_side"])
            right_norm = normalize_side_text(r["3p_side"])
            key = (r["cut_range"], left_norm, right_norm)
            if key not in grouped:
                grouped[key] = {
                    "patients": set(),
                    "timepoints": defaultdict(set),
                    "cut_range": r["cut_range"],
                    "functional_region": r["functional_region"],
                    "left": left_norm,
                    "right": right_norm,
                    "source": r["source"],
                }
            g = grouped[key]
            g["patients"].add(r["patients"])
            if r["timepoints"]:
                for tp in r["timepoints"].split(","):
                    if tp:
                        g["timepoints"][r["patients"]].add(tp)

        out_rows = []
        for g in grouped.values():
            patient_list = sorted(g["patients"])
            tp_txt_parts = []
            for p in patient_list:
                tps = sorted(g["timepoints"].get(p, set()))
                if tps:
                    tp_txt_parts.append(f"{p}:{'/'.join(tps)}")
            out_rows.append(
                {
                    "patients": ",".join(patient_list),
                    "timepoints": ";".join(tp_txt_parts),
                    "cut_range": g["cut_range"],
                    "functional_region": g["functional_region"],
                    "5p_side": g["left"],
                    "3p_side": g["right"],
                    "source": g["source"],
                }
            )

    out_rows.sort(key=lambda r: (r["patients"], r["cut_range"], r["5p_side"], r["3p_side"]))
    with out_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "patients",
                "timepoints",
                "cut_range",
                "functional_region",
                "5p_side",
                "3p_side",
                "source",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(out_rows)
    return len(out_rows)


def main():
    ap = argparse.ArgumentParser(description="Export phase-table rows to TSV.")
    ap.add_argument("--reconstructed", default="reconstructed_gamete_aas_0.75.tsv")
    ap.add_argument("--out", default="phase_table_0.75.tsv")
    ap.add_argument("--results", default="results_gametes.txt")
    ap.add_argument("--patient", default="", help="Optional patient ID filter, e.g., 024-A")
    ap.add_argument(
        "--only-3n-indels",
        action="store_true",
        help="Ignore events containing indels whose length is not multiple of 3; prefix output with only3indels_.",
    )
    args = ap.parse_args()

    out_tsv = maybe_prefix_output(Path(args.out), enabled=args.only_3n_indels)
    reconstructed_path = Path(args.reconstructed)
    if not reconstructed_path.exists():
        raise FileNotFoundError(
            f"Missing reconstructed input: {reconstructed_path}. "
            "Run script 0 first (0_build_reconstructed_gamete_aas.py)."
        )
    allowed = load_allowed_event_keys_from_results(Path(args.results), only_3n_indels=args.only_3n_indels)
    n = export_one(
        reconstructed_path,
        out_tsv,
        patient_filter=(args.patient.strip() or None),
        allowed_event_keys=allowed,
    )
    print(f"[OK] Rows exported: {n}")
    print(f"[OK] TSV: {out_tsv}")


if __name__ == "__main__":
    main()
