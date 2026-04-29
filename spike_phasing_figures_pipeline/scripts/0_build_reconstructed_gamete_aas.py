#!/usr/bin/env python3
"""
0_build_reconstructed_gamete_aas.py

Reviewed on 2026-04-28 by Luis Daniel González Vázquez.

Summary: Generates reconstructed_gamete_aas_*.tsv from results_gametes.txt (0-based) joined with *_genes.tsv to reconstruct nucleotide and amino-acid gametes.

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
DEFAULT_THRESHOLDS = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

GENE_CELL_RE = re.compile(
    r"^(?P<alt>\d+)\/(?P<cov>\d+)\/A(?P<A>\d+)C(?P<C>\d+)T(?P<T>\d+)G(?P<G>\d+)N(?P<N>\d+)\/\+(?P<ins_raw>.+)-(?P<del_raw>.+)$"
)


def parse_results_line(line: str):
    return dict(part.split(":", 1) for part in line.split("\t") if ":" in part)


def parse_sample(sample: str):
    m = re.match(r"^(.*)-([^\\-]+)_pTrim$", sample or "")
    if m:
        return m.group(1), m.group(2)
    if "-" in (sample or ""):
        left, right = sample.rsplit("-", 1)
        return left, right
    return sample or "", ""


def nt0_to_nt1(pos_nt0: int):
    return pos_nt0 + 1


def nt1_to_aa_idx_in_spike(pos_nt1: int):
    """
    Returns None when a reported polymorphic breakpoint position is outside Spike.
    """
    if pos_nt1 < SPIKE_START_1B or pos_nt1 > SPIKE_END_1B:
        return None
    return ((pos_nt1 - SPIKE_START_1B) // 3) + 1


def codon_start_1b(pos_nt1: int):
    return pos_nt1 - ((pos_nt1 - SPIKE_START_1B) % 3)


def translate_codon(codon: str):
    codon = codon.upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "amb"
    return CODON_TABLE.get(codon, "amb")


def parse_gametes(detail: str):
    out = []
    if not detail:
        return out
    m = re.search(r"gametes:([^\\n]+)$", detail)
    if not m:
        return out
    for tok in m.group(1).split(";"):
        parts = tok.split(">")
        if len(parts) < 4:
            continue
        try:
            out.append((parts[0].strip(), parts[1].strip(), int(parts[2]), float(parts[3])))
        except Exception:
            continue
    return out


def detail_indel_has_non3n(detail: str):
    matches = list(re.finditer(r"([+-])(\d+)", detail or ""))
    if not matches:
        return False
    return any((int(m.group(2)) % 3) != 0 for m in matches)


def parse_indel_allele(allele: str):
    m = re.match(r"^([+-])(\d+)([A-Za-z]+)?$", allele or "")
    if not m:
        return None
    sign = m.group(1)
    n = int(m.group(2))
    seq = (m.group(3) or "").upper()
    return sign, n, seq


def aa_from_indel_allele(allele: str, pos1: int | None = None):
    p = parse_indel_allele(allele)
    if not p:
        return None
    sign, n, seq = p
    if n % 3 != 0:
        return f"{sign}@"
    # 3n indel but not at codon start: translation depends on local context
    # (neighboring bases), so translating "seq" in isolation is not reliable.
    if pos1 is not None:
        phase = (pos1 - codon_start_1b(pos1)) % 3
        if phase != 0:
            return "amb"
    if not seq:
        # without sequence we cannot reconstruct the affected peptide
        return "amb"
    # traducimos la secuencia del indel en marcos de 3
    aas = []
    for i in range(0, len(seq) - (len(seq) % 3), 3):
        aa = translate_codon(seq[i:i + 3])
        aas.append(aa if aa != "amb" else "X")
    if not aas:
        return "amb"
    return f"{sign}{''.join(aas)}"


def parse_gene_cell(cell: str):
    if not cell or cell == "NA":
        return None
    m = GENE_CELL_RE.match(cell.strip())
    if not m:
        return None
    d = m.groupdict()
    # ins_raw can be "0" or "TC500"; extract trailing numeric count
    im = re.search(r"(\d+)$", d["ins_raw"])
    ins_count = int(im.group(1)) if im else 0
    # del_raw can be "0" or "TC2"; extract trailing numeric count
    dm = re.search(r"(\d+)$", d["del_raw"])
    del_count = int(dm.group(1)) if dm else 0
    return {
        "cov": int(d["cov"]),
        "A": int(d["A"]),
        "C": int(d["C"]),
        "T": int(d["T"]),
        "G": int(d["G"]),
        "N": int(d["N"]),
        "ins": ins_count,
        "del": del_count,
    }


def consensus_base(counts, threshold: float):
    if counts is None:
        return None, 0.0
    total = counts.get("cov", 0)
    if total <= 0:
        # conservative fallback if cov is missing or invalid
        total = counts["A"] + counts["C"] + counts["T"] + counts["G"] + counts["N"]
    if total <= 0:
        return None, 0.0
    ordered = sorted(
        [("A", counts["A"]), ("C", counts["C"]), ("T", counts["T"]), ("G", counts["G"])],
        key=lambda x: x[1],
        reverse=True,
    )
    b, c = ordered[0]
    freq = c / total
    if freq >= threshold:
        return b, freq
    return None, freq


def load_results_events(results_path: Path, only_3n_indels: bool):
    events = []
    by_patient_positions = defaultdict(set)
    by_patient_tps = defaultdict(set)
    for line in results_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip():
            continue
        fields = parse_results_line(line)
        sample = (fields.get("Sra") or "").strip()
        if not sample:
            continue
        try:
            posA0 = int(fields.get("PosA", "").strip())
            posB0 = int(fields.get("PosB", "").strip())
        except Exception:
            continue
        posA1 = nt0_to_nt1(posA0)
        posB1 = nt0_to_nt1(posB0)
        aaA = nt1_to_aa_idx_in_spike(posA1)
        aaB = nt1_to_aa_idx_in_spike(posB1)
        if aaA is None or aaB is None:
            continue
        detail = fields.get("Detail", "")
        if only_3n_indels and detail_indel_has_non3n(detail):
            continue
        gam = parse_gametes(detail)
        if len(gam) < 4:
            continue
        patient, tp = parse_sample(sample)
        event = {
            "patient": patient,
            "sample": sample,
            "timepoint": tp,
            "posA0": posA0,
            "posB0": posB0,
            "posA1": posA1,
            "posB1": posB1,
            "aaA": aaA,
            "aaB": aaB,
            "gametes": gam[:4],
        }
        events.append(event)
        by_patient_tps[patient].add(tp)
        for p1 in (posA1, posB1):
            cs = codon_start_1b(p1)
            by_patient_positions[patient].update([cs, cs + 1, cs + 2])
    return events, by_patient_positions, by_patient_tps


def load_patient_gene_data(genes_dir: Path, needed_positions, needed_timepoints):
    """Return dict[pos1][tp] = parsed gene cell counts."""
    out = defaultdict(dict)
    if not genes_dir.exists():
        return out
    with genes_dir.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if not rd.fieldnames:
            return out
        tp_cols = [c for c in rd.fieldnames if c not in ("position", "gene", "codon_pos", "aa_idx")]
        tp_cols = [c for c in tp_cols if c in needed_timepoints]
        for r in rd:
            try:
                pos1 = int(r.get("position", ""))
            except Exception:
                continue
            if pos1 not in needed_positions:
                continue
            for tp in tp_cols:
                out[pos1][tp] = parse_gene_cell(r.get(tp, ""))
    return out


def codon_aa_from_patient_data(pos1, allele, tp, threshold, patient_data):
    indel_aa = aa_from_indel_allele(allele, pos1=pos1)
    if indel_aa is not None:
        return indel_aa
    if len(allele) != 1 or allele.upper() not in "ACGTN":
        return "amb"

    cs = codon_start_1b(pos1)
    codon_positions = [cs, cs + 1, cs + 2]
    bases = []
    for p in codon_positions:
        if p == pos1:
            b = allele.upper()
            if b == "N":
                return "amb"
            bases.append(b)
        else:
            b, _ = consensus_base(patient_data.get(p, {}).get(tp), threshold)
            if not b:
                return "amb"
            bases.append(b)
    return translate_codon("".join(bases))


def build_rows(results_path: Path, genes_dir: Path, thresholds, only_3n_indels: bool):
    events, needed_positions, needed_tps = load_results_events(results_path, only_3n_indels)
    patient_gene = {}
    for patient in sorted(needed_positions):
        fp = genes_dir / f"{patient}_genes.tsv"
        if not fp.exists():
            patient_gene[patient] = defaultdict(dict)
            continue
        patient_gene[patient] = load_patient_gene_data(fp, needed_positions[patient], needed_tps[patient])

    rows_by_thr = {t: [] for t in thresholds}
    for ev in events:
        patient = ev["patient"]
        tp = ev["timepoint"]
        pdata = patient_gene.get(patient, defaultdict(dict))
        for thr in thresholds:
            out = {
                "patient": patient,
                "sample": ev["sample"],
                "timepoint": tp,
                "posA": str(ev["posA0"]),  # Keep 0-based to match results_gametes.txt.
                "posB": str(ev["posB0"]),  # Keep 0-based to match results_gametes.txt.
                "posA1b": str(ev["posA1"]),
                "posB1b": str(ev["posB1"]),
            }
            for i, (a1, a2, cnt, frq) in enumerate(ev["gametes"], start=1):
                aa1 = codon_aa_from_patient_data(ev["posA1"], a1, tp, thr, pdata)
                aa2 = codon_aa_from_patient_data(ev["posB1"], a2, tp, thr, pdata)
                out[f"gamete{i}_nt"] = f"{ev['aaA']}:{a1}|{ev['aaB']}:{a2}"
                out[f"gamete{i}_aa"] = f"{ev['aaA']}:{aa1}|{ev['aaB']}:{aa2}"
                out[f"gamete{i}_reads"] = str(cnt)
                out[f"gamete{i}_freq"] = f"{frq:.6f}"
            rows_by_thr[thr].append(out)
    return rows_by_thr


def write_rows(rows, out_path: Path):
    fields = [
        "patient", "sample", "timepoint", "posA", "posB",
        "posA1b", "posB1b",
        "gamete1_nt", "gamete2_nt", "gamete3_nt", "gamete4_nt",
        "gamete1_aa", "gamete2_aa", "gamete3_aa", "gamete4_aa",
        "gamete1_reads", "gamete2_reads", "gamete3_reads", "gamete4_reads",
        "gamete1_freq", "gamete2_freq", "gamete3_freq", "gamete4_freq",
    ]
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(
        description="Build reconstructed_gamete_aas_*.tsv from results_gametes.txt + *_genes.tsv (results 0-based)."
    )
    ap.add_argument("--results", default="results_gametes.txt")
    ap.add_argument("--genes-dir", default="LINES_CHICUADRADO/patient_genes_tsv")
    ap.add_argument("--thresholds", default="0.70,0.75,0.80,0.85,0.90,0.95")
    ap.add_argument("--out-prefix", default="reconstructed_gamete_aas_")
    ap.add_argument("--only-3n-indels", action="store_true")
    args = ap.parse_args()

    thresholds = []
    for t in args.thresholds.split(","):
        t = t.strip()
        if not t:
            continue
        thresholds.append(float(t))
    if not thresholds:
        thresholds = list(DEFAULT_THRESHOLDS)

    rows_by_thr = build_rows(
        Path(args.results),
        Path(args.genes_dir),
        thresholds=thresholds,
        only_3n_indels=args.only_3n_indels,
    )
    for t in thresholds:
        name = f"{args.out_prefix}{t:.2f}.tsv"
        out_path = Path(name)
        if args.only_3n_indels and not out_path.name.startswith("only3indels_"):
            out_path = out_path.with_name("only3indels_" + out_path.name)
        write_rows(rows_by_thr[t], out_path)
        print(f"[OK] {out_path} rows={len(rows_by_thr[t])}")


if __name__ == "__main__":
    main()
