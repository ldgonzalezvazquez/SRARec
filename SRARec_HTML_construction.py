# -*- coding: utf-8 -*-
import os
# forzar locale básico y evitar nl_langinfo(CODESET) errors
os.environ.setdefault('LC_ALL', 'C')
os.environ.setdefault('LANG',   'C.UTF-8')

"""
rec_dashboard_ldhat_v5_fix_parallel.py
─────────────────────────────────────
• separate Gene-selector  +  nsp-selector
• full heat-map grid (light-grey = below thresholds)
• modular JS  (plotly_white)
• NO f-string in HTML  →  no "{{ … }}" escaping headaches
• NEW METHODOLOGY - genome profiles per sample and mean ρ calculation
• RESUME FUNCTIONALITY - proper cache handling
• PARALLEL PROCESSING - patient genome profiles calculated in parallel using multiprocessing
"""

import argparse, gzip, json, pickle, re, warnings, sys
from pathlib import Path
from functools import lru_cache
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm.auto import tqdm

warnings.filterwarnings("ignore", category=UserWarning)
BAR = dict(bar_format=" {l_bar}{bar:25}{r_bar}")

# ═════════ helpers ═════════════════════════════════════════

def process_patient_batch(args_tuple):
    """
    Procesa un lote de pacientes en paralelo
    args_tuple: (patient_batch, events_df, genome_len, alpha)
    """
    patient_batch, events_df, genome_len, alpha = args_tuple
    results = {}
    
    for sample in patient_batch:
        pos_t = np.zeros(genome_len, dtype=np.float32)
        pos_f = np.zeros(genome_len, dtype=np.float32)
        
        # Filtrar eventos para esta muestra específica
        sample_events = events_df[events_df["sample"] == sample]
        
        for _, ev in sample_events.iterrows():
            a, b = sorted((int(ev.pA), int(ev.pB)))
            if b <= a: continue
            w = 1.0/(b - a)
            if ev.true == 1:
                pos_t[a:b] += w
            else:
                pos_f[a:b] += w
        
        tot = pos_t + pos_f
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_arr = np.where(tot>0, alpha * pos_t/tot, 0.0)
        
        results[sample] = ratio_arr
    
    return results

def calculate_gene_rho_batch(args_tuple):
    """
    Calcula ρ por genes para un lote de pacientes en paralelo
    args_tuple: (patient_batch, sample_genome_profiles, genes, cds)
    """
    patient_batch, sample_genome_profiles, genes, cds = args_tuple
    results = []
    for s in patient_batch:
        sample_profile = sample_genome_profiles[s]
        for g in genes:
            # Encontrar las posiciones que flanquea este gen
            gene_positions = cds[cds["lbl"] == g]
            if len(gene_positions) > 0:
                start_pos = gene_positions.iloc[0]["start"]
                end_pos = gene_positions.iloc[0]["end"]
                # Calcular media de ρ en las posiciones del gen (incluyendo ceros)
                gene_rho_values = sample_profile[start_pos:end_pos]
                rho_val = np.mean(gene_rho_values)
                results.append((s, g, rho_val))
    return results

def calculate_nsp_rho_batch(args_tuple):
    """
    Calcula ρ por NSPs para un lote de pacientes en paralelo
    args_tuple: (patient_batch, sample_genome_profiles, nsps, pep)
    """
    patient_batch, sample_genome_profiles, nsps, pep = args_tuple
    results = []
    for s in patient_batch:
        sample_profile = sample_genome_profiles[s]
        for r in nsps:
            # Encontrar las posiciones que flanquea este NSP
            nsp_positions = pep[pep["lbl"] == r]
            if len(nsp_positions) > 0:
                start_pos = nsp_positions.iloc[0]["start"]
                end_pos = nsp_positions.iloc[0]["end"]
                # Calcular media de ρ en las posiciones del NSP (incluyendo ceros)
                nsp_rho_values = sample_profile[start_pos:end_pos]
                rho_val = np.mean(nsp_rho_values)
                results.append((s, r, rho_val))
    return results

def calculate_protein_rho_batch(args_tuple):
    """
    Calcula ρ por proteínas para un lote de pacientes en paralelo
    args_tuple: (patient_batch, sample_genome_profiles, all_proteins_list, proteins, nsps, all_proteins, pep)
    """
    patient_batch, sample_genome_profiles, all_proteins_list, proteins, nsps, all_proteins, pep = args_tuple
    results = []
    for s in patient_batch:
        sample_profile = sample_genome_profiles[s]
        for p in all_proteins_list:
            # Buscar en genes primero, luego en NSPs
            protein_positions = None
            if p in proteins:
                protein_positions = all_proteins[all_proteins["lbl"] == p]
            elif p in nsps:
                protein_positions = pep[pep["lbl"] == p]
            
            if protein_positions is not None and len(protein_positions) > 0:
                start_pos = protein_positions.iloc[0]["start"]
                end_pos = protein_positions.iloc[0]["end"]
                # Calcular media de ρ en las posiciones de la proteína (incluyendo ceros)
                protein_rho_values = sample_profile[start_pos:end_pos]
                rho_val = np.mean(protein_rho_values)
                results.append((s, p, rho_val))
    return results

def split_pid_time(s: str):
    p = s.split("-")
    return ("-".join(p[:2]) if len(p) >= 3 else s,
            "-".join(p[2:])  if len(p) >= 3 else "1")

_mire = re.compile(r"mireads:(\d+)")
_mfre = re.compile(r"minfreq:([0-9.]+)")

def parse_pairs(field: str, expect_big=True):
    if not field or field.strip() in ("NA", ""): return []
    out = []
    for tok in field.split(";"):
        if tok.strip() == "": continue
        p = tok.split(",")
        try:
            pA, pB = int(p[0]), int(p[1])
        except ValueError:
            continue
        big = int(p[2]) if expect_big and len(p) >= 3 else 1
        mireads = int(_mire.search(tok).group(1)) if _mire.search(tok) else 0
        minfreq = float(_mfre.search(tok).group(1)) if _mfre.search(tok) else 0.0
        out.append((pA, pB, big, mireads, minfreq))
    return out

def load_ann(gb: Path):
    rec = SeqIO.read(str(gb), "genbank")
    cds, pep, all_proteins, seen = [], [], [], {}
    
    for f in rec.features:
        if f.type not in ("CDS", "mat_peptide"):
            continue
            
        prod = f.qualifiers.get("product", [""])[0] or ""
        gene_qual = f.qualifiers.get("gene", [""])[0] or ""
        gene_name = gene_qual or prod
        
        # Para genes: usar gene_name (incluye ORF1ab)
        if f.type == "CDS":
            lbl = gene_name
        else:
            # Para mat_peptide: usar prod y limpiar
            lbl = prod.replace("non-structural protein ", "nsp").strip() or "unknown"
        
        st, en = int(f.location.start), int(f.location.end)
        prev = seen.get(lbl)
        if prev and (en - st) <= (prev[2] - prev[1]):
            continue
        seen[lbl] = (lbl, st, en, f.type)
    
    # Para proteínas: excluir poliproteínas como en interactive_heatmap_upset.py
    for lbl, st, en, ftype in seen.values():
        if ftype == "mat_peptide":
            pep.append((lbl, st, en))
        else:
            cds.append((lbl, st, en))
            # Para all_proteins: excluir SOLO ORF1ab y ORF1a (poliproteínas), incluir ORF10
            if "orf1ab" not in lbl.lower() and "orf1a" not in lbl.lower():
                all_proteins.append((lbl, st, en))
    
    # Agregar NSPs a all_proteins
    for lbl, st, en in pep:
        all_proteins.append((lbl, st, en))
    
    cds_df = (pd.DataFrame(cds, columns=["lbl", "start", "end"])
                .sort_values("start").reset_index(drop=True))
    pep_df = (pd.DataFrame(pep, columns=["lbl", "start", "end"])
                .sort_values("start").reset_index(drop=True))
    all_proteins_df = (pd.DataFrame(all_proteins, columns=["lbl", "start", "end"])
                        .sort_values("start").reset_index(drop=True))
    return len(rec.seq), cds_df, pep_df, all_proteins_df

def anno_at(pos: int, table_df):
    """Optimized version that works with DataFrame directly"""
    if table_df is None or table_df.empty: return None
    hit = table_df.query("start <= @pos and end > @pos")
    return hit.lbl.iloc[0] if not hit.empty else None

def ensure_cols(ev: pd.DataFrame, cds, pep, all_proteins=None):
    """
    Optimized version: NO calcula anotaciones que se recalculan después
    Solo calcula sample y mid que son necesarios para el procesamiento
    """
    ev = ev.copy()
    if {"patient", "time"}.issubset(ev.columns) is False:
        pt = ev["sra"].apply(split_pid_time).apply(pd.Series)
        pt.columns = ["patient", "time"]
        ev = pd.concat([ev, pt], axis=1)
    
    ev["sample"] = ev["patient"] + " | " + ev["time"]
    
    # NO calcular mid aquí - no se usa en ningún lugar del HTML
    # Esto elimina ~6M operaciones innecesarias adicionales
    
    return ev

# ═════════ SIMPLIFIED PROCESSING ═══════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results",  help="Results file (if not provided, will use existing cache)")
    ap.add_argument("--cache",    help="Cache file (.cache.pkl.gz) - if not provided, will be auto-generated")
    ap.add_argument("--genbank",  required=True)
    ap.add_argument("--out",      default="rec_dashboard.html")
    ap.add_argument("--alpha",    type=float, default=1.0)
    ap.add_argument("--min_reads", type=int, default=None, help="Minimum reads threshold to filter events")
    ap.add_argument("--min_freq", type=float, default=None, help="Minimum frequency threshold to filter events")
    ap.add_argument("--cores", type=int, default=None, help="Number of cores to use for parallel processing (default: 75% of available cores)")
    args = ap.parse_args()
    
    # Modificar nombre de salida si hay filtros aplicados
    out_name = args.out
    if args.min_reads is not None or args.min_freq is not None:
        base_name = Path(args.out).stem
        ext = Path(args.out).suffix
        filters = []
        if args.min_reads is not None:
            filters.append(f"minR{args.min_reads}")
        if args.min_freq is not None:
            filters.append(f"minF{args.min_freq}")
        out_name = f"{base_name}_{'_'.join(filters)}{ext}"
    
    out = Path(out_name)
    
    # Determine cache file path
    if args.cache:
        cache = Path(args.cache)
    else:
        cache = out.with_suffix(".cache.pkl.gz")

    # ── cache load / build ──────────────────────────────────
    cache_data = None
    if cache.exists():
        print("→ Loading cache from:", cache.name)
        try:
            cache_data = pickle.load(gzip.open(cache, "rb"))
            events = pd.DataFrame(cache_data["events"])
            cds = pd.DataFrame(cache_data["cds"])
            pep = pd.DataFrame(cache_data["pep"])
            genome_len = cache_data["len"]
            print("✔ Cache loaded successfully")
        except Exception as e:
            print(f"⚠️ Error loading cache: {e}")
            print("→ Rebuilding from scratch...")
            cache_data = None
    
    # Always load annotations from GenBank file
    genome_len, cds, pep, all_proteins = load_ann(Path(args.genbank))
    
    if cache_data is None:
        if not args.results:
            print("❌ Error: No results file provided and no valid cache found")
            print("Please provide --results argument or ensure cache file exists")
            sys.exit(1)
        
        print("→ Parsing results file...")
        rows = []
        total_events = 0
        with open(args.results) as fh:
            total_lines = sum(1 for _ in fh); fh.seek(0)
            progress_bar = tqdm(fh, total=total_lines, desc="Reading file lines", unit="line")
            for ln in progress_bar:
                if "Sra:" not in ln: continue
                kv = {x.split(":", 1)[0]: x.split(":", 1)[1]
                      for x in ln.strip().split("\t") if ":" in x}
                sra = kv["Sra"].strip()
                patient, time = split_pid_time(sra)
                for lst, flag in ((parse_pairs(kv.get("Rec", "")), 1),
                                  (parse_pairs(kv.get("No_Rec", ""), False), 0)):
                    for pA, pB, bg, mr, mf in lst:
                       rows.append((patient, sra, time,
                                    pA, pB, abs(pB - pA),
                                    flag, bg, mr, mf))
                       total_events += 1
                       # Update progress bar with event count
                       if total_events % 10000 == 0:
                           progress_bar.set_postfix(events=f"{total_events:,}")
        events = pd.DataFrame(rows, columns=
            ["patient","sra","time","pA","pB","dist",
             "true","bigger","mireads","minfreq"])
        
        print(f"→ Parsed {total_events:,} events from {total_lines:,} lines")
        
        # Aplicar filtros si se especificaron
        if args.min_reads is not None or args.min_freq is not None:
            original_count = len(events)
            # Los eventos FALSE nunca se filtran por reads/freq
            true_events = events[events['true'] == 1]
            false_events = events[events['true'] == 0]
            
            if args.min_reads is not None:
                true_events = true_events[true_events['mireads'] >= args.min_reads]
            if args.min_freq is not None:
                true_events = true_events[true_events['minfreq'] >= args.min_freq]
            
            # Reconstruir events con true filtrados + false sin filtrar
            events = pd.concat([true_events, false_events], ignore_index=True)
            filtered_count = len(events)
            print(f"→ Filtered TRUE events: {len(events[events['true']==1]):,} (FALSE events always kept: {len(false_events):,})")
            if args.min_reads is not None:
                print(f"   min_reads >= {args.min_reads} (only for TRUE events)")
            if args.min_freq is not None:
                print(f"   min_freq >= {args.min_freq} (only for TRUE events)")

        events = ensure_cols(events, cds, pep, all_proteins)

        # ── SAVE TO CACHE ─────────────────────────────────────────────────────────
        cache_data = {
            "events": events,
            "cds": cds,
            "pep": pep,
            "len": genome_len
        }
        pickle.dump(cache_data, gzip.open(cache, "wb"))
        print("✔ Data saved to cache:", cache.name)

        # lists  (now that 'sample' column exists)
        genes   = cds.lbl.tolist()
        nsps    = pep.lbl.tolist()
        proteins = all_proteins.lbl.tolist()  # All proteins (genes + NSPs)
        # Compatible sort without key parameter (for older pandas compatibility)
        samples_series = events["sample"].drop_duplicates()
        samples_df = samples_series.to_frame()
        samples_df['_sort_key'] = (samples_series.str.split(" \\| ").str[0] + 
                                    samples_series.str.split(" \\| ").str[1])
        samples = samples_df.sort_values('_sort_key')['sample'].tolist()
        patients = sorted(set(events.patient))

        print(f"→ Processing {len(events)} events with {len(samples)} samples")
        print(f"→ Genes: {len(genes)}, NSPs: {len(nsps)}")

    # Definir lista de todas las proteínas (genes + NSPs, sin duplicados)
    all_proteins_list = list(dict.fromkeys(proteins + nsps))



    # ── PREPARE FOR ρ CALCULATION ──────────────────────────────────────────────
    print("→ Preparing for ρ calculation...")


    # ── GENOME PROFILE PER SAMPLE (PARALLEL) ──────────────────────────────────
    print("→ Calculating genome profiles per sample (parallel)...")
    # Crear perfiles genómicos individuales por paciente usando paralelización
    
    # Determinar número de cores
    if args.cores is not None:
        n_cores = max(1, args.cores)
    else:
        # Usar 75% de los disponibles para no sobrecargar
        n_cores = max(1, int(cpu_count() * 0.75))
    print(f"→ Using {n_cores} cores for parallel processing (available: {cpu_count()})")
    
    # Dividir pacientes en lotes para procesamiento paralelo
    batch_size = max(1, len(samples) // n_cores)
    patient_batches = [samples[i:i+batch_size] for i in range(0, len(samples), batch_size)]
    
    # Preparar argumentos para cada proceso
    batch_args = [(batch, events, genome_len, args.alpha) for batch in patient_batches]
    
    # Procesar lotes en paralelo
    sample_genome_profiles = {}
    with Pool(processes=n_cores) as pool:
        # Usar tqdm para mostrar progreso de los lotes
        batch_results = list(tqdm(
            pool.imap(process_patient_batch, batch_args),
            total=len(batch_args),
            desc="Processing patient batches",
            unit="batch"
        ))
    
    # Combinar resultados de todos los lotes
    for batch_result in batch_results:
        sample_genome_profiles.update(batch_result)
    
    print(f"→ Processed {len(sample_genome_profiles)} samples successfully")
    
    # Perfil genómico global (suma de todos los pacientes)
    print("→ Calculating global genome profile...")
    pos_t_global = np.zeros(genome_len, dtype=np.float32)
    pos_f_global = np.zeros(genome_len, dtype=np.float32)
    
    # IMPORTANTE: Para el perfil genómico usar TODOS los eventos SIN FILTRAR
    # (igual que en el backup) para mantener consistencia con versiones anteriores
    # Usar eventos filtrados (ya incluye solo TRUE filtrados + FALSE sin filtrar)
    print("  → Using filtered events (TRUE filtered + FALSE unfiltered) for genome profile")
    events_for_genome = events
    
    for _, ev in events_for_genome.iterrows():  # events YA FILTRADOS
        a, b = sorted((int(ev.pA), int(ev.pB)))
        if b <= a: continue
        w = 1.0/(b - a)
        if ev.true == 1:
            pos_t_global[a:b] += w
        else:
            pos_f_global[a:b] += w
    tot_global = pos_t_global + pos_f_global
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_arr_global = np.where(tot_global>0, args.alpha * pos_t_global/tot_global, 0.0)
    pos_ratio = ratio_arr_global.tolist()
    
    # ── CALCULATE ρ MATRICES FROM GENOME PROFILES ───────────────────────────────
    print("→ Calculating ρ matrices from genome profiles...")
    
    # Genes - calcular ρ como media de sitios que flanquean cada gen (PARALELO)
    print("→ Calculating gene ρ matrices (parallel)...")
    batch_size = max(1, len(samples) // n_cores)
    patient_batches = [samples[i:i+batch_size] for i in range(0, len(samples), batch_size)]
    batch_args = [(batch, sample_genome_profiles, genes, cds) for batch in patient_batches]
    
    rows = []
    with Pool(processes=n_cores) as pool:
        batch_results = list(tqdm(
            pool.imap(calculate_gene_rho_batch, batch_args),
            total=len(batch_args),
            desc="Processing gene ρ batches",
            unit="batch"
        ))
    
    for batch_result in batch_results:
        rows.extend([(g, s, rho_val) for s, g, rho_val in batch_result])
    rho_df = pd.DataFrame(rows, columns=["gene", "sample", "rho"])
    mat = (rho_df.pivot(index="gene", columns="sample", values="rho")
             .reindex(genes).reindex(columns=samples).fillna(0))
    
    # Calcular ρ promedio para "Whole CDS" por muestra (media del perfil genómico)
    genome_wide_rho = []
    
    for j, s in enumerate(samples):
        # Calcular media del perfil genómico para esta muestra (incluyendo ceros)
        sample_profile = sample_genome_profiles[s]
        rho_genome = np.mean(sample_profile)
        genome_wide_rho.append(float(rho_genome))
    
    # Añadir fila "Whole CDS" al principio
    genes_with_genome = ["Whole CDS"] + genes
    z = [genome_wide_rho] + mat.values.tolist()
    
    # NSPs - calcular ρ como media de sitios que flanquean cada NSP (PARALELO)
    print("→ Calculating NSP ρ matrices (parallel)...")
    batch_args_r = [(batch, sample_genome_profiles, nsps, pep) for batch in patient_batches]
    
    rows_r = []
    with Pool(processes=n_cores) as pool:
        batch_results_r = list(tqdm(
            pool.imap(calculate_nsp_rho_batch, batch_args_r),
            total=len(batch_args_r),
            desc="Processing NSP ρ batches",
            unit="batch"
        ))
    
    for batch_result in batch_results_r:
        rows_r.extend([(r, s, rho_val) for s, r, rho_val in batch_result])
    rho_r_df = pd.DataFrame(rows_r, columns=["region", "sample", "rho"])
    mat_r = (rho_r_df.pivot(index="region", columns="sample", values="rho")
                 .reindex(nsps).reindex(columns=samples).fillna(0))
    zR = mat_r.values.tolist()

    # Proteins - calcular ρ como media de sitios que flanquean cada proteína (PARALELO)
    print("→ Calculating protein ρ matrices (parallel)...")
    batch_args_p = [(batch, sample_genome_profiles, all_proteins_list, proteins, nsps, all_proteins, pep) for batch in patient_batches]
    
    rows_p = []
    with Pool(processes=n_cores) as pool:
        batch_results_p = list(tqdm(
            pool.imap(calculate_protein_rho_batch, batch_args_p),
            total=len(batch_args_p),
            desc="Processing protein ρ batches",
            unit="batch"
        ))
    
    for batch_result in batch_results_p:
        rows_p.extend([(p, s, rho_val) for s, p, rho_val in batch_result])
    rho_p_df = pd.DataFrame(rows_p, columns=["protein", "sample", "rho"])
    mat_p = (rho_p_df.pivot(index="protein", columns="sample", values="rho")
                 .reindex(all_proteins_list).reindex(columns=samples).fillna(0))
    
    # Calcular ρ promedio para "Whole CDS" en protein view
    protein_wide_rho = []
    for j, s in enumerate(samples):
        sample_profile = sample_genome_profiles[s]
        rho_prot = np.mean(sample_profile)
        protein_wide_rho.append(float(rho_prot))
    
    proteins_with_whole = ["Whole CDS"] + all_proteins_list
    zP_with_whole = [protein_wide_rho] + mat_p.values.tolist()
    zP = zP_with_whole
    
    # Precalculate genome statistics for hover (OPTIMIZED - vectorized)
    print("→ Precalculating genome statistics (vectorized)...")
    
    # Solo eventos TRUE para las estadísticas
    true_events = events[events["true"] == 1].copy()
    
    # Crear un DataFrame con todas las posiciones que cada evento cubre
    print("  → Expanding events to positions...")
    position_data = []
    for _, event in tqdm(true_events.iterrows(), total=len(true_events), desc="Expanding to positions", unit="event"):
        pA, pB = int(event["pA"]), int(event["pB"])
        start_pos, end_pos = sorted([pA, pB])
        for pos in range(start_pos, end_pos):
            position_data.append({
                "pos": pos,
                "event_id": f"{pA}_{pB}",
                "patient": event["patient"],
                "sample": event["sample"],
                "time": event["time"]
            })
    
    pos_df = pd.DataFrame(position_data)
    
    # Calcular estadísticas por posición usando groupby (MUCHO MÁS RÁPIDO)
    print("  → Calculating statistics per position...")
    pos_stats = pos_df.groupby("pos").agg({
        "event_id": "nunique",  # Número de eventos únicos
        "patient": "nunique",   # Número de pacientes únicos
        "sample": "nunique"     # Número de muestras únicas (timepoints)
    }).rename(columns={
        "event_id": "event_count",
        "patient": "patient_count", 
        "sample": "timepoint_count"
    })
    
    # Obtener contribuidores por posición
    print("  → Calculating contributors...")
    contributors_dict = {}
    for pos, group in tqdm(pos_df.groupby("pos"), desc="Contributors", unit="pos"):
        contrib = []
        for patient, patient_group in group.groupby("patient"):
            times = sorted(patient_group["time"].unique())
            contrib.append([patient, times])
        contributors_dict[pos] = contrib
    
    # Obtener anotaciones de genes/NSPs por posición
    print("  → Adding gene/NSP annotations...")
    gene_annotations = {}
    nsp_annotations = {}
    
    for _, gene_row in cds.iterrows():
        for pos in range(gene_row["start"], gene_row["end"]):
            gene_annotations[pos] = gene_row["lbl"]
    
    for _, nsp_row in pep.iterrows():
        for pos in range(nsp_row["start"], nsp_row["end"]):
            nsp_annotations[pos] = nsp_row["lbl"]
    
    # Construir genome_stats final
    print("  → Building final genome stats...")
    genome_stats = []
    for pos in tqdm(range(genome_len), desc="Final stats", unit="pos"):
        if pos in pos_stats.index:
            row = pos_stats.loc[pos]
            event_count = int(row["event_count"])
            patient_count = int(row["patient_count"])
            timepoint_count = int(row["timepoint_count"])
            contrib = contributors_dict.get(pos, [])
        else:
            event_count = 0
            patient_count = 0
            timepoint_count = 0
            contrib = []
        
        gene_name = gene_annotations.get(pos, "")
        nsp_name = nsp_annotations.get(pos, "")
        
        genome_stats.append({
            "pos": pos,
            "ratio": float(pos_ratio[pos]),
            "event_count": event_count,
            "patient_count": patient_count,
            "timepoint_count": timepoint_count,
            "contributors": contrib,
            "gene": gene_name,
            "nsp": nsp_name
        })

    # ── BUILD PAYLOAD ────────────────────────────────────────────────────────
    print("→ Building payload...")
    # Usar eventos ya filtrados (los filtros se aplicaron según args.min_reads y args.min_freq)
    
    # Expandir eventos para incluir TODOS los genes/proteínas que flanquean cada evento
    # Solo procesar eventos TRUE (los FALSE no aportan información útil)
    print("→ Expanding TRUE events to include all flanking genes/proteins...")
    true_events_for_expansion = events[events["true"] == 1]
    print(f"→ Processing {len(true_events_for_expansion):,} TRUE events (out of {len(events):,} total)")
    expanded_events = []
    
    for _, ev in tqdm(true_events_for_expansion.iterrows(), total=len(true_events_for_expansion), desc="Expanding TRUE events", unit="event"):
        pA, pB = int(ev["pA"]), int(ev["pB"])
        start_pos, end_pos = sorted([pA, pB])
        
        # Encontrar todos los genes que flanquea este evento
        flanking_genes = []
        for _, gene_row in cds.iterrows():
            gene_start, gene_end = gene_row["start"], gene_row["end"]
            # Evento flanquea gen si se superpone
            if not (end_pos <= gene_start or start_pos >= gene_end):
                flanking_genes.append(gene_row["lbl"])
        
        # Encontrar todos los NSPs que flanquea este evento
        flanking_nsps = []
        for _, nsp_row in pep.iterrows():
            nsp_start, nsp_end = nsp_row["start"], nsp_row["end"]
            # Evento flanquea NSP si se superpone
            if not (end_pos <= nsp_start or start_pos >= nsp_end):
                flanking_nsps.append(nsp_row["lbl"])
        
        # Encontrar todas las proteínas que flanquea este evento
        flanking_proteins = []
        for _, prot_row in all_proteins.iterrows():
            prot_start, prot_end = prot_row["start"], prot_row["end"]
            # Evento flanquea proteína si se superpone
            if not (end_pos <= prot_start or start_pos >= prot_end):
                flanking_proteins.append(prot_row["lbl"])
        
        # Crear una entrada para cada gen que flanquea
        for gene in flanking_genes:
            expanded_events.append({
                "pA": ev["pA"], "pB": ev["pB"],
                "true": ev["true"], "gene": gene, "region": None, "protein": None,
                "sample": ev["sample"], "mireads": ev["mireads"], 
                "minfreq": ev["minfreq"], "bigger": ev["bigger"]
            })
        
        # Crear una entrada para cada NSP que flanquea
        for nsp in flanking_nsps:
            expanded_events.append({
                "pA": ev["pA"], "pB": ev["pB"],
                "true": ev["true"], "gene": None, "region": nsp, "protein": nsp,
                "sample": ev["sample"], "mireads": ev["mireads"], 
                "minfreq": ev["minfreq"], "bigger": ev["bigger"]
            })
        
        # Crear una entrada para cada proteína que flanquea (que no sea NSP)
        for protein in flanking_proteins:
            if protein not in flanking_nsps:  # Evitar duplicados con NSPs
                expanded_events.append({
                    "pA": ev["pA"], "pB": ev["pB"],
                    "true": ev["true"], "gene": None, "region": None, "protein": protein,
                    "sample": ev["sample"], "mireads": ev["mireads"], 
                    "minfreq": ev["minfreq"], "bigger": ev["bigger"]
                })
    
    # NO eliminar duplicados - cada evento puede aparecer múltiples veces
    # (una vez por cada gen/NSP/proteína que flanquea)
    # Esto es necesario para que los popups funcionen correctamente
    genome_events = expanded_events.copy()
    
    # Todos los eventos en genome_events son TRUE (ya filtrados arriba)
    events_list = genome_events
    print(f"→ Original events: {len(events):,}")
    print(f"→ Expanded events: {len(expanded_events):,}")
    print(f"→ Genome events (deduplicated): {len(genome_events):,}")
    print(f"→ Events for HTML (TRUE only): {len(events_list):,}")
    print(f"→ Expansion factor: {len(expanded_events)/len(events):.1f}x")

    payload = dict(
        genes       = genes_with_genome,  # Incluir "Whole genome" como primera fila
        nsps        = nsps,
        proteins    = proteins_with_whole,  # All proteins for protein view
        samples     = samples,
        z           = z,  # Ya incluye fila "Whole genome"
        nsps_z      = zR,
        proteins_z  = zP,  # Protein ρ data
        genome_len  = genome_len,
        events      = events_list,  # Incluir muestra representativa de eventos
        genes_start = [0] + cds.start.tolist(),  # Añadir start=0 para "Whole genome"
        genes_end   = [genome_len] + cds.end.tolist(),  # Añadir end=genome_len para "Whole genome"
        nsps_start  = pep.start.tolist(),
        nsps_end    = pep.end.tolist(),
        proteins_start = [0] + all_proteins.start.tolist(),  # Añadir start=0 para "Whole CDS"
        proteins_end   = [genome_len] + all_proteins.end.tolist(),  # Añadir end=genome_len para "Whole CDS"
        patients    = patients,
        pos_ratio   = pos_ratio,
        genome_stats = genome_stats,  # Estadísticas precalculadas del genoma
        alpha       = args.alpha,
        contact     = "luisdaniel.gonzalez@uvigo.es"
    )

    # ── HTML TEMPLATE ────────────────────────────────────────────────────────
    print("→ Generating HTML...")
    html_tpl = r"""
<!DOCTYPE html>
    <!-- Generated at 2025-09-18 -->
    <html><head><meta charset="utf-8">
    <title>SARS-CoV-2 Recombination Analysis</title>
    <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
     * {
       margin: 0;
       padding: 0;
       box-sizing: border-box;
     }
     
     body {
       font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
       background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
       min-height: 100vh;
       color: #333;
       line-height: 1.6;
     }
     
     .container {
       max-width: 1400px;
       margin: 0 auto;
       padding: 20px;
     }
     
     .header {
       background: rgba(255, 255, 255, 0.95);
       backdrop-filter: blur(10px);
       border-radius: 15px;
       padding: 25px;
       margin-bottom: 30px;
       box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
       border: 1px solid rgba(255, 255, 255, 0.2);
       position: sticky;
       top: 20px;
       z-index: 1000;
     }
     
     .header h1 {
       color: #2c3e50;
       font-size: 2.2em;
       margin-bottom: 15px;
       text-align: center;
       font-weight: 300;
     }
     
     .controls {
       display: flex;
       flex-wrap: wrap;
       gap: 20px;
       align-items: center;
       justify-content: center;
     }
     
     .control-group {
       display: flex;
       flex-direction: column;
       align-items: center;
       gap: 8px;
     }
     
     .control-group.hidden {
       display: none;
     }
     
     .control-group label {
       font-weight: 600;
       color: #34495e;
       font-size: 0.9em;
       text-transform: uppercase;
       letter-spacing: 0.5px;
     }
     
     select {
       padding: 10px 15px;
       border: 2px solid #e0e6ed;
       border-radius: 8px;
       background: white;
       font-size: 14px;
       min-width: 150px;
       transition: all 0.3s ease;
       box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
     }
     
     select:focus {
       outline: none;
       border-color: #667eea;
       box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
     }
     
     input[type="number"] {
       padding: 8px 12px;
       border: 2px solid #e0e6ed;
       border-radius: 8px;
       background: white;
       font-size: 14px;
       width: 80px;
       transition: all 0.3s ease;
       box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
     }
     
     input[type="number"]:focus {
       outline: none;
       border-color: #667eea;
       box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
     }
     
     button {
       padding: 12px 24px;
       border: none;
       border-radius: 8px;
       font-weight: 600;
       cursor: pointer;
       transition: all 0.3s ease;
       font-size: 14px;
       text-transform: uppercase;
       letter-spacing: 0.5px;
     }
     
     #toggle {
       background: linear-gradient(135deg, #667eea, #764ba2);
       color: white;
       box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
     }
     
     #toggle:hover {
       transform: translateY(-2px);
       box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
     }
     
     #go {
       background: linear-gradient(135deg, #27ae60, #2ecc71);
       color: white;
       box-shadow: 0 4px 15px rgba(39, 174, 96, 0.3);
     }
     
     #go:hover {
       transform: translateY(-2px);
       box-shadow: 0 6px 20px rgba(39, 174, 96, 0.4);
     }
     
     #reset {
       background: linear-gradient(135deg, #e67e22, #f39c12);
       color: white;
       box-shadow: 0 4px 15px rgba(230, 126, 34, 0.3);
     }
     
     #reset:hover {
       transform: translateY(-2px);
       box-shadow: 0 6px 20px rgba(230, 126, 34, 0.4);
     }
     
     .chart-section {
       background: rgba(255, 255, 255, 0.95);
       backdrop-filter: blur(10px);
       border-radius: 15px;
       padding: 25px;
       margin-bottom: 30px;
       box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
       border: 1px solid rgba(255, 255, 255, 0.2);
     }
     
     .chart-section h2 {
       color: #2c3e50;
       font-size: 1.5em;
       margin-bottom: 20px;
       font-weight: 400;
       border-bottom: 2px solid #e0e6ed;
       padding-bottom: 10px;
     }
     
     .chart-container {
       margin-bottom: 30px;
       position: relative;
       overflow: hidden;
       min-height: auto;
     }
     
     .chart-container:last-child {
       margin-bottom: 0;
     }
     
     .chart-container > div {
       width: 100% !important;
       height: auto !important;
       min-height: auto;
     }
     
    /* === Control único del espacio tras los pies de figura ===
       Ajusta este valor para ampliar/reducir el hueco final. */
    :root { --footnote-gap: 2px; }  /* p.ej.: 0px, 4px, 8px */

    /* Reset de contador para numerar figuras según aparición */
    body { counter-reset: figure; }

    /* Pie de figura */
    .footnote {
      font-size: 14px;
      margin-top: 6px;
      margin-bottom: 0;            /* ← clave para evitar hueco debajo */
      color: #7f8c8d;
      font-style: italic;
      text-align: center;
      padding: 2px 8px;            /* ← padding vertical reducido */
      background: rgba(236, 240, 241, 0.5);
      border-radius: 6px;
      line-height: 1.2;            /* ← line-height reducido */
    }
    .footnote:empty { display: none; } /* evita reservar espacio si viene vacía */
    /* Si el pie de figura es el último elemento del contenedor, elimina margen extra */
    .footnote:last-child { margin-bottom: 0; }

    .footnote::before {
      counter-increment: figure;
      content: "Figure " counter(figure) ". ";
      font-weight: 600;
      color: #2c3e50;
    }

    /* ===== Colapsar huecos del CONTENEDOR cuando termina en footnote =====
       Aplica a wrappers típicos de tus figuras. Si usas otros, añádelos aquí. */
    .figure,
    .figure-wrapper,
    .plot-card,
    .chart-card,
    .card {
      /* Por si algún estilo global metía mucho padding al final */
      padding-bottom: var(--footnote-gap);
      /* Espacio externo bajo la tarjeta/figura (tócalo si quieres menos hueco entre figuras) */
      margin-bottom: var(--footnote-gap);
    }
    /* Si el ÚLTIMO hijo es el footnote, no dejes margen extra dentro */
    .figure > .footnote:last-child,
    .figure-wrapper > .footnote:last-child,
    .plot-card > .footnote:last-child,
    .chart-card > .footnote:last-child,
    .card > .footnote:last-child {
      margin-bottom: 0 !important;
    }

    /* Versión con :has() (moderna). Si tu navegador lo soporta, asegura aún más el colapso.
       Puedes dejarlo sin problema; Safari/Chromium modernos lo soportan. */
    .figure:has(> .footnote:last-child),
    .figure-wrapper:has(> .footnote:last-child),
    .plot-card:has(> .footnote:last-child),
    .chart-card:has(> .footnote:last-child),
    .card:has(> .footnote:last-child) {
      padding-bottom: var(--footnote-gap) !important;
      margin-bottom: var(--footnote-gap) !important;
    }

    /* Blindaje extra por si algún contenedor reservaba altura fija */
    .figure-footer,
    .footnote-container {
      height: auto !important;
      min-height: 0 !important;
      padding-bottom: 0 !important;
    }

    /* Si usas Plotly y su contenedor añade hueco inferior, elimínalo */
    .js-plotly-plot { margin-bottom: 0 !important; }

     
     .popup-overlay {
       position: fixed;
       top: 0;
       left: 0;
       width: 100%;
       height: 100%;
       background: rgba(0, 0, 0, 0.5);
       display: flex;
       align-items: center;
       justify-content: center;
       z-index: 10000;
     }
     
     .popup-content {
       background: white;
       padding: 20px;
       border-radius: 10px;
       box-shadow: 0 10px 30px rgba(0,0,0,0.3);
       max-width: 90vw;
       max-height: 80vh;
       overflow-y: auto;
       position: relative;
     }
     
     .popup-header {
       display: flex;
       justify-content: space-between;
       align-items: center;
       margin-bottom: 20px;
       padding-bottom: 10px;
       border-bottom: 1px solid #ecf0f1;
     }
     
     .popup-header h3 {
       margin: 0;
       color: #2c3e50;
     }
     
     .close-btn {
       position: absolute;
       top: 10px;
       right: 15px;
       background: #e74c3c;
       color: white;
       border: none;
       border-radius: 50%;
       width: 30px;
       height: 30px;
       cursor: pointer;
       font-size: 16px;
       display: flex;
       align-items: center;
       justify-content: center;
       transition: all 0.3s ease;
     }
     
     .close-btn:hover {
       background: #c0392b;
       transform: scale(1.1);
     }
     
     .info-btn {
       margin-left: 6px;
       font-size: 11px;
       padding: 0 6px;
       height: 20px;
       line-height: 18px;
       background: #e9ecef;
       color: #333;
       border: 1px solid #ccc;
       border-radius: 10px;
       cursor: pointer;
     }
     
     .view-toggle {
       text-align: center;
       margin: 20px 0;
     }
     
     .view-toggle button {
       background: #6c757d;
       margin: 0 10px;
     }
     
     .view-toggle button.active {
       background: linear-gradient(135deg, #667eea, #764ba2);
     }

     /* Contenedor de figura + pie: quita alturas fijas/reservas */
     .figure-card,
     .figure-container,
     .chart-card,
     [data-figure],
     .section-figure {
       height: auto !important;
       min-height: auto !important;   /* ← quita "reservas" de altura */
       padding-bottom: 0 !important;  /* ← evita colchón inferior */
     }

     /* Si las figuras están en una GRID y ves filas "altas" de más */
     .figure-grid {
       grid-auto-rows: minmax(min-content, auto);
       align-content: start;
     }

     /* Grupo de filtro en columna: título encima del control */
     .filter-group {
       display: flex;
       flex-direction: column;
       align-items: center;
     }

     .filter-group > label {
       font-weight: 600;
       margin-bottom: 6px;
       display: block;
       text-align: center;
     }
</style></head><body>
    <div class="container">
      <div class="header">
        <h1>🧬 SARS-CoV-2 Recombination Analysis</h1>
        <div class="controls">
          <div class="control-group">
            <button id="code-btn" style="background: #2c3e50; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; margin-right: 10px;">💻 Code</button>
            <button id="methodology-btn" style="background: #3498db; color: white; border: none; border-radius: 6px; padding: 8px 16px; cursor: pointer; font-size: 14px; margin-right: 10px;">📖 Methodology</button>
</div>
<div class="control-group">
            <label for="selP">👥 Patients</label>
<select id="selP" multiple size="5"></select>
</div>
          <div class="filter-group" id="gene-filter">
            <label for="selG">🧬 Genes</label>
<select id="selG" multiple size="5"></select>
</div>
          <div class="filter-group" id="protein-filter" style="display:none">
            <label for="selPr">🦠 Proteins</label>
<select id="selPr" multiple size="5"></select>
</div>
          <div class="control-group" style="display: flex; flex-direction: column; gap: 8px;">
</div>
<div class="control-group">
            <button id="toggle">🔄 Protein View</button>
</div>
          <div class="control-group" style="display: flex; gap: 10px;">
            <button id="go">✅ Apply Filters</button>
            <button id="reset">🔄 Reset</button>
</div>
</div>
</div>

<div id="gene-plots">
        <div class="chart-section">
          <h2>📊 Recombination Analysis - Genes</h2>
          <div class="chart-container">
          <h3>ρ (genes) - <span id="gene-title">All genes</span></h3>
            <div id="hm"></div>
            <div class="footnote">Heatmap displaying log<sub>10</sub> recombination rate (ρ) for each gene across all samples. Color intensity represents recombination activity: darker purple indicates lower rates, while bright green/yellow shows higher recombination. Light grey cells indicate values below filtering thresholds. ρ is the population recombination rate, calculated as the ratio of recombination events to total events scaled by α. Use the gene filter above to focus on specific genes.</div>
          </div>
</div>

        <div class="chart-section">
          <h2>📈 Genome Profile - Genes</h2>
          <div class="chart-container">
<h3>Genome recombination ratio (by position)</h3>
            <div id="prof_ratio"></div>
            <div class="footnote">Genome-wide recombination profile showing recombination ratio = α·(true events)/(total events) at each genomic position. The plot reveals hotspots and coldspots of recombination activity across the genome. Hover over data points to see detailed gene and NSP annotations with patient information. Gene annotations are displayed as colored bars below the plot. Click on points to explore specific genomic regions.</div>
          </div>
</div>

        <div class="chart-section">
          <h2>⏰ Time Series - Genes</h2>
          <div class="chart-container">
<h3>Time series – genes - <span id="ts-gene-title">All genes</span></h3>
            <div id="ts_genes"></div>
            <div class="footnote">Temporal evolution of recombination activity showing mean log<sub>10</sub> ρ over selected genes plotted against sampling time. Each line represents a different sample, allowing comparison of recombination patterns over time. Use the gene filter to focus on specific genes and observe how recombination rates change throughout the sampling period.</div>
          </div>
</div>

</div>

<div id="protein-plots" style="display:none">
        <div class="chart-section">
          <h2>📊 Recombination Analysis - Proteins</h2>
          <div class="chart-container">
            <h3>ρ (proteins) - <span id="protein-title">All proteins</span></h3>
            <div id="cds"></div>
            <div class="footnote">Heatmap displaying log<sub>10</sub> recombination rate (ρ) for each protein across all samples. Color intensity represents recombination activity: darker purple indicates lower rates, while bright green/yellow shows higher recombination. Light grey cells indicate values below filtering thresholds. ρ is the population recombination rate, calculated as the ratio of recombination events to total events scaled by α. Use the protein filter above to focus on specific proteins.</div>
          </div>
</div>

        <div class="chart-section">
          <h2>📈 Genome Profile - Proteins</h2>
          <div class="chart-container">
            <h3>Genome recombination ratio (by position)</h3>
            <div id="prof_ratio_nsps"></div>
            <div class="footnote">Genome-wide recombination profile showing recombination ratio = α·(true events)/(total events) at each genomic position. The plot reveals hotspots and coldspots of recombination activity across the genome. Hover over data points to see detailed gene and protein annotations with patient information. Protein annotations are displayed as colored bars below the plot. Click on points to explore specific genomic regions.</div>
          </div>
</div>

        <div class="chart-section">
          <h2>⏰ Time Series - Proteins</h2>
          <div class="chart-container">
            <h3>Time series – proteins - <span id="ts-protein-title">All proteins</span></h3>
            <div id="ts_nsps"></div>
            <div class="footnote">Temporal evolution of recombination activity showing mean log<sub>10</sub> ρ over selected proteins plotted against sampling time. Each line represents a different sample, allowing comparison of recombination patterns over time. Use the protein filter to focus on specific proteins and observe how recombination rates change throughout the sampling period.</div>
          </div>
</div>

</div>
    </div>
    
    
<!-- Modal for contributors -->
<div id="contrib-modal" style="display:none; position:fixed; top:0; left:0; right:0; bottom:0; background:rgba(0,0,0,0.4); align-items:center; justify-content:center; z-index:2000;">
  <div style="background:white; padding:20px; border-radius:8px; max-width:700px; max-height:80vh; overflow:auto; margin: auto; position:relative; padding-top:40px;">
    <button id="contrib-close" style="position:absolute; top:5px; right:10px; background:none; border:none; font-size:28px; cursor:pointer; color:#666; line-height:1;">&times;</button>
    <h3 style="margin-top:0;">Contributors at position <span id="contrib-pos"></span></h3>
    <div id="contrib-content"></div>
  </div>
  </div>

<!-- Modal for events -->
<div id="events-modal" style="display:none; position:fixed; top:0; left:0; right:0; bottom:0; background:rgba(0,0,0,0.4); align-items:center; justify-content:center; z-index:2000;">
  <div style="background:white; padding:20px; border-radius:8px; max-width:900px; max-height:80vh; overflow:auto; margin: auto; position:relative;">
    <button id="events-close" style="position:absolute; top:10px; right:10px; background:none; border:none; font-size:24px; cursor:pointer; color:#666;">&times;</button>
    <h3>Event Details</h3>
    <div id="events-content"></div>
  </div>
</div>
    
<!-- Modal for filter info -->
<div id="info-modal" style="display:none; position:fixed; top:0; left:0; right:0; bottom:0; background:rgba(0,0,0,0.4); align-items:center; justify-content:center; z-index:2000;">
  <div style="background:white; padding:20px; border-radius:8px; max-width:700px; max-height:80vh; overflow:auto;">
    <h3 id="info-title">Filter info</h3>
    <div id="info-content"></div>
    <div style="text-align:right; margin-top:15px;"><button id="info-close">Close</button></div>
  </div>
  </div>
    
    <!-- Modal for methodology -->
    <div id="methodology-modal" style="display:none; position:fixed; top:0; left:0; right:0; bottom:0; background:rgba(0,0,0,0.4); align-items:center; justify-content:center; z-index:2000;">
      <div style="background:white; padding:20px; border-radius:8px; max-width:800px; max-height:80vh; overflow-y:auto;">
        <h3>Methodology</h3>
        <div id="methodology-content" style="text-align: left; line-height: 1.6; font-size: 14px;">
          <p><strong>Methodology (from BAM to recombination figures)</strong></p>
          <p>Analyses began from coordinate-sorted, indexed BAM files (one per sample), processed with standard tooling for alignment data handling (sorting/indexing and pileup generation with SAMtools/BCFtools). All computations operated directly on mapped reads; no local re-assembly or consensus calling was used to infer phase.</p>
          <p>Candidate polymorphic sites were identified by scanning alignment pileups for multi-allelic evidence. For each such position, the pipeline assembled—at read level—the sets of reads supporting each observed allele (single-nucleotide changes and short indels), applying base-quality screening and excluding reads that appeared to support more than one allele at the same position to avoid ambiguous evidence. Deletions were left-anchored so that extended gaps were not double-counted across adjacent coordinates.</p>
          <p>To determine which site pairs could be phased on the same molecules, positions were organized into windows whose span did not exceed the maximum mapped read length observed in a sample. Within each window, positions were grouped when their allele-frequency profiles were similar (frequencies within a narrow tolerance of a running group average for at least one component). This yielded compact sets of co-observed positions suitable for read-level recombination testing.</p>
          <p>Recombination was assessed with a read-level four-gamete test: for each co-covered site pair, the four possible haplotypes (AB, Ab, aB, ab) were counted from the intersection of read sets supporting each allele. A pair was classified as four-gamete-positive when all four haplotypes were observed directly on reads (or read pairs) with sufficient support; pairs lacking any haplotype—or supported exclusively by indel-only intersections—were treated as non-positive. Testing proceeded (i) within groups—starting from the most distant positions and refining toward the median—and (ii) between adjacent groups within the same window, while tracking containment to retain only minimal (non-contained) positive ranges.</p>
          <p>For each sample, the pipeline produced a compact record encoding all tested site pairs. For positive pairs, the record includes the genomic coordinates, whether a pair is contained within a broader positive range, and summary support metrics (minimum per-haplotype read count and minimum per-haplotype frequency among co-covering reads). Non-positive pairs were also retained. These per-sample lines constitute the sole inputs to downstream aggregation and visualization.</p>
          <p>Pair records were then annotated to coding features (genes, non-structural proteins, individual proteins) using the reference GenBank. Each pair contributes a uniform weight across the genomic interval spanned by its two sites; specifically, a pair spanning (b−a) bases contributes 1/(b−a) to every covered position. For any sample × region (gene/NSP/protein), contributions are summed separately for four-gamete-positive (T) and non-positive (F) evidence, and a region-level recombination score is computed as ρ = T / (T + F). The same T/(T+F) ratio is computed genome-wide at nucleotide resolution to produce a positional profile. The self-contained interactive HTML renders heatmaps of log10(ρ) by genes and by proteins/NSPs, a genome-wide positional track of ρ, and optional multivariate summaries, built directly from these per-sample records.</p>
          <p><strong>References (for this section only)</strong><br>
          Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.<br>
          Hudson RR, Kaplan NL. Statistical properties of the number of recombination events in the history of a sample of DNA sequences. Genetics. 1985 Aug;111(1):147–64. doi: 10.1093/genetics/111.1.147. PMID: 4029609; PMCID: PMC1202594.</p>
        </div>
        <div style="text-align:right; margin-top:15px;"><button id="methodology-close">Close</button></div>
      </div>
    </div>
    
<script>const DATA=__DATA__;</script>
<script>
const selP=document.getElementById('selP'),
      selG=document.getElementById('selG'),
      selPr=document.getElementById('selPr');
// fill selectors BEFORE redraw
DATA.patients.forEach(p=>selP.add(new Option(p,p)));

/* numeric filter inputs - REMOVED */
const iSam=document.getElementById('minSam');

// Event listeners para recalcular heatmap automáticamente - FILTROS NUMÉRICOS REMOVIDOS

/* log warning if empty */
if(DATA.nsps.length===0){
  console.warn("No nsp labels shipped – selector will appear empty.");
}
document.addEventListener('keydown',e=>{
 if(e.key==='a'){[...selP.options].forEach(o=>o.selected=true);}
 if(e.key==='n'){[...selP.options].forEach(o=>o.selected=false);}
});
function idx(a,x){return a.indexOf(x);}
function resetFilters() {
  // Reset all selectors
  [...selP.options].forEach(o => o.selected = false);
  [...selG.options].forEach(o => o.selected = false);
  [...selPr.options].forEach(o => o.selected = false);
  
  // Reset numeric inputs - FILTROS NUMÉRICOS REMOVIDOS
  
  // Redraw with reset values
  setTimeout(redraw, 0);
}

function updateTitles(pickG, pickPr) {
  // Update gene titles
  const geneTitle = pickG.length > 0 ? pickG.join(', ') : 'All genes';
  const tsGeneTitle = pickG.length > 0 ? pickG.join(', ') : 'All genes';
  
  // Update protein titles
  const proteinTitle = pickPr.length > 0 ? pickPr.join(', ') : 'All proteins';
  const tsProteinTitle = pickPr.length > 0 ? pickPr.join(', ') : 'All proteins';
  
  // Update DOM elements
  const geneTitleEl = document.getElementById('gene-title');
  const tsGeneTitleEl = document.getElementById('ts-gene-title');
  const proteinTitleEl = document.getElementById('protein-title');
  const tsProteinTitleEl = document.getElementById('ts-protein-title');
  
  if (geneTitleEl) geneTitleEl.textContent = geneTitle;
  if (tsGeneTitleEl) tsGeneTitleEl.textContent = tsGeneTitle;
  if (proteinTitleEl) proteinTitleEl.textContent = proteinTitle;
  if (tsProteinTitleEl) tsProteinTitleEl.textContent = tsProteinTitle;
}

// --- Helpers: soporte por evento y mapas pos->gene/protein (seguros) ---
function getEventSupport(ev) {
  // Usa la primera clave disponible que refleje el "soporte" del evento (como en las tablas)
  return ev.support ?? ev.events ?? ev.n ?? ev.n_events ?? ev.support_count ?? 0;
}
function eventPasses(ev, minMr, minMf, minEvt) {
  if (minMr && ev.mireads < minMr) return false;
  if (minMf && ev.minfreq < minMf) return false;
  if (minEvt && getEventSupport(ev) < minEvt) return false;
  return true;
}
function geneAt(pos){ const a = DATA.pos_gene || DATA.gene_by_pos; return a ? (a[pos] ?? '') : ''; }
function nspAt(pos){  const a = DATA.pos_nsp  || DATA.nsp_by_pos;  return a ? (a[pos] ?? '') : ''; }

// Recalcular perfil del genoma (alpha*T/(T+F)) usando TODOS los eventos (true y false)
function recalculateGenomeProfile(minEvt, minMr, minMf, selectedGenes, selectedSamples) {
  const a = DATA.alpha || 1.0;
  const N = DATA.genome_len;
  const selSamples = (selectedSamples && selectedSamples.length) ? new Set(selectedSamples) : null;
  const pos_t = new Float32Array(N), pos_f = new Float32Array(N);
  for (const ev of DATA.events) {
    if (minMr && ev.mireads < minMr) continue;
    if (minMf && ev.minfreq < minMf) continue;
    if (selSamples && !selSamples.has(ev.sample)) continue;
    const A = Math.min(ev.pA, ev.pB)|0, B = Math.max(ev.pA, ev.pB)|0;
    if (B<=A) continue;
    const w = 1.0/(B-A);
    const arr = (ev.true===1) ? pos_t : pos_f;
    for (let p=A; p<B && p<N; p++) if (p>=0) arr[p]+=w;
  }
  const out = Array(N);
  for (let p=0; p<N; p++){
    const T=pos_t[p], F=pos_f[p], tot=T+F;
    out[p] = (tot>0) ? a*(T/tot) : 0;
  }
  
  // DEBUG: Imprimir valores del perfil del genoma
  console.log('=== DEBUG: Perfil del genoma recalculado ===');
  const nonZeroValues = out.filter(v => v > 0);
  console.log(`Valores no cero: ${nonZeroValues.length} de ${out.length}`);
  if (nonZeroValues.length > 0) {
    console.log(`Rango: ${Math.min(...nonZeroValues).toFixed(6)} a ${Math.max(...nonZeroValues).toFixed(6)}`);
    console.log(`Primeros 10 valores no cero:`, nonZeroValues.slice(0, 10).map(v => v.toFixed(6)));
  }
  
  return out;
}

// Recalcula el heatmap de genes usando la MEDIA del ratio por posición dentro del gen.
// Crucial: el filtrado por Min events se aplica POR EVENTO (per-event), no como recuento por gen.
function recalculateHeatmapWithFilters(data, minEvt, minMr, minMf, selectedGenes, selectedSamples) {
  const N      = DATA.genome_len;
  const alpha  = DATA.alpha || 1.0;
  const genes  = (selectedGenes && selectedGenes.length) ? selectedGenes : DATA.genes;
  const samples= (selectedSamples && selectedSamples.length) ? selectedSamples : DATA.samples;

  // Indexar eventos por muestra aplicando UMBRALES INDIVIDUALES (per-event)
  const evBySample = {};
  for (const s of samples) evBySample[s] = [];
  for (const ev of DATA.events) {
    if (!evBySample.hasOwnProperty(ev.sample)) continue;
    if (!eventPasses(ev, minMr, minMf, minEvt)) continue; // <-- per-event
    evBySample[ev.sample].push(ev);
  }

  // Acumular T/F por posición para cada muestra (ponderación 1/(b-a))
  const TF = {};
  for (const s of samples) {
    TF[s] = { t: new Float32Array(N), f: new Float32Array(N) };
    for (const ev of evBySample[s]) {
      const a = Math.min(ev.pA, ev.pB) | 0;
      const b = Math.max(ev.pA, ev.pB) | 0;
      if (b <= a) continue;
      const w = 1.0 / (b - a);
      const arr = (ev.true === 1) ? TF[s].t : TF[s].f;
      for (let p = a; p < b && p < N; p++) if (p >= 0) arr[p] += w;
    }
  }
  
           // DEBUG: Imprimir valores de T/F para algunas muestras
           console.log('=== DEBUG: Valores T/F por muestra ===');
           for (const s of samples.slice(0, 3)) {
             const t_sum = Array.from(TF[s].t).reduce((a,b) => a+b, 0);
             const f_sum = Array.from(TF[s].f).reduce((a,b) => a+b, 0);
             console.log(`Muestra ${s}: T_sum=${t_sum.toFixed(6)}, F_sum=${f_sum.toFixed(6)}`);
           }
           
           // DEBUG: Mostrar en el HTML también
           const debugDiv = document.createElement('div');
           debugDiv.innerHTML = '<h3>DEBUG: Valores T/F</h3>';
           for (const s of samples.slice(0, 3)) {
             const t_sum = Array.from(TF[s].t).reduce((a,b) => a+b, 0);
             const f_sum = Array.from(TF[s].f).reduce((a,b) => a+b, 0);
             debugDiv.innerHTML += `<p>Muestra ${s}: T_sum=${t_sum.toFixed(6)}, F_sum=${f_sum.toFixed(6)}</p>`;
           }
           document.body.appendChild(debugDiv);
           
  const zFlat = [];
  const allVals = [];
  for (const gene of genes) {
    const gi = DATA.genes.indexOf(gene);
    if (gi === -1) continue;
    const g0 = DATA.genes_start[gi], g1 = DATA.genes_end[gi];
    for (const s of samples) {
      // Media del ratio por posición dentro del gen (solo posiciones con soporte)
      const t = TF[s].t, f = TF[s].f;
      let sum = 0, cnt = 0;
      for (let p = g0; p < g1 && p < N; p++) {
        const T = t[p], F = f[p], tot = T + F;
        if (tot > 0) { sum += alpha * (T / tot); cnt++; }
      }
      const val = (cnt > 0) ? (sum / cnt) : -9;
      
      zFlat.push(val);
      if (val !== -9 && !isNaN(val)) allVals.push(val);
    }
  }
  const valid = allVals.length ? allVals : [-9];
  const zmin = (valid[0] === -9) ? -9 : Math.min(...valid);
  const zmax = (valid[0] === -9) ? -9 : Math.max(...valid);
  
  // DEBUG: Imprimir valores del heatmap calculados
  console.log('=== DEBUG: Valores del heatmap ===');
  const heatmapValues = zFlat.filter(v => v !== -9);
  console.log(`Valores válidos del heatmap: ${heatmapValues.length}`);
  if (heatmapValues.length > 0) {
    console.log(`Rango: ${Math.min(...heatmapValues).toFixed(6)} a ${Math.max(...heatmapValues).toFixed(6)}`);
    console.log(`Primeros 10 valores:`, heatmapValues.slice(0, 10).map(v => v.toFixed(6)));
    console.log(`Últimos 10 valores:`, heatmapValues.slice(-10).map(v => v.toFixed(6)));
  }
  
  return { z: zFlat, zmin, zmax };
}

// Función para mostrar popup con eventos detallados
function showEventsPopup(gene, sample, rhoValue, events, minEvt, minMr, minMf) {
  const modal = document.getElementById('events-modal');
  const content = document.getElementById('events-content');
  
  // rhoValue ya viene como el valor correcto de ρ (no es log10)
  let html = `
    <div style="margin-bottom: 15px;">
      <strong>Gene/Protein:</strong> ${gene}<br>
      <strong>Sample:</strong> ${sample}<br>
      <strong>ρ (rho):</strong> ${rhoValue.toFixed(6)}
    </div>
  `;
  
  if (events.length === 0) {
    html += '<p>No events found for this combination.</p>';
  } else {
    html += `<h4>Events (${events.length} total)</h4>`;
    html += '<table id="events-table" style="width: 100%; border-collapse: collapse; margin-bottom: 20px;">';
    html += '<tr style="background: #f8f9fa;"><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTable(0)">Event ID ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTable(1)">Times detected ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTable(2)">Min Reads ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTable(3)">Min Freq (%) ↕</th></tr>';
    
    events.forEach(eventData => {
      const event = eventData.event;
      const eventId = `${event.pA},${event.pB}`;
      html += `<tr>
        <td style="border: 1px solid #ddd; padding: 8px;">${eventId}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${event.bigger || 1}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${event.mireads || 0}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${((event.minfreq || 0) * 100).toFixed(2)}%</td>
      </tr>`;
    });
    html += '</table>';
  }
  
  content.innerHTML = html;
  modal.style.display = 'flex';
}

// Popup para eventos de múltiples pacientes (cuando se hace clic en "Average")
function showEventsPopupWithPatients(gene, sample, rhoValue, events) {
  const modal = document.getElementById('events-modal');
  const content = document.getElementById('events-content');
  
  let html = `
    <div style="margin-bottom: 15px;">
      <strong>Gene/Protein:</strong> ${gene}<br>
      <strong>Sample:</strong> ${sample}<br>
      <strong>ρ (rho):</strong> ${rhoValue.toFixed(6)}
    </div>
  `;
  
  if (events.length === 0) {
    html += '<p>No events found for this combination.</p>';
  } else {
    html += `<h4>Events (${events.length} total)</h4>`;
    html += '<table id="events-table-with-patients" style="width: 100%; border-collapse: collapse; margin-bottom: 20px;">';
    html += '<tr style="background: #f8f9fa;"><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTableWithPatients(0)">Patient ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTableWithPatients(1)">Event ID ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTableWithPatients(2)">Times detected ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTableWithPatients(3)">Min Reads ↕</th><th style="border: 1px solid #ddd; padding: 8px; cursor: pointer;" onclick="sortTableWithPatients(4)">Min Freq (%) ↕</th></tr>';
    
    events.forEach(eventData => {
      const event = eventData.event;
      const eventId = `${event.pA},${event.pB}`;
      const patient = event.sample || 'N/A';
      html += `<tr>
        <td style="border: 1px solid #ddd; padding: 8px;">${patient}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${eventId}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${event.bigger || 1}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${event.mireads || 0}</td>
        <td style="border: 1px solid #ddd; padding: 8px;">${((event.minfreq || 0) * 100).toFixed(2)}%</td>
      </tr>`;
    });
    html += '</table>';
  }
  
  content.innerHTML = html;
  modal.style.display = 'flex';
}

function redraw(){
  const pickP=[...selP.selectedOptions].map(o=>o.value);
  const pickG=[...selG.selectedOptions].map(o=>o.value);
  const pickPr=[...selPr.selectedOptions].map(o=>o.value);
  
  // Update titles
  updateTitles(pickG, pickPr);
  const minEvt=0,  // Filtros numéricos removidos
        minMr =0,
        minMf =0;
  // filas ordenadas por paciente y luego tiempo numérico
  const cols = DATA.samples
    .filter(s => !pickP.length || pickP.includes(s.split(' | ')[0]))
    .sort((a,b) => {
      const [pa,ta] = a.split(' | ');  const [pb,tb] = b.split(' | ');
      if(pa!==pb) return pa.localeCompare(pb);
      return (parseFloat(ta.replace(/[^0-9.]/g,''))||0)
           - (parseFloat(tb.replace(/[^0-9.]/g,''))||0);
    });
  
  // Calcular heatmap: VISUALIZACIÓN en log10, pero mantener valores originales para hover/popup
  const allRhoValues = [...DATA.z.flat(), ...DATA.nsps_z.flat(), ...DATA.proteins_z.flat()].filter(v => v > 0);
  const allRhoLog10 = allRhoValues.map(v => Math.log10(v + 1e-9));
  const zmax_log = allRhoLog10.length > 0 ? Math.max(...allRhoLog10) : 0;
  const zmin_log = allRhoLog10.length > 0 ? Math.min(...allRhoLog10) : -9;
  const idxC = cols.map(c => DATA.samples.indexOf(c));
  
  // genes: swap axes (x=genes, y=samples) and add first row with mean across samples (EXACTAMENTE COMO Rec_A.py)
  const idxG = [];
  DATA.genes.forEach((g,i)=>{
    if(pickG.length && !pickG.includes(g)) return;
    idxG.push(i);
  });
  const xGenes = idxG.map(i=>DATA.genes[i]);
  
  // Mean row (includes zeros; uses current patient selection via idxC)
  const meanRowGenes = xGenes.map((_,k)=>{
    const i = idxG[k];
    const vals = idxC.map(j=>DATA.z[i][j]).filter(v=>v!=null);  // Incluir valores 0
    if(vals.length===0) return zmin_log;  // log10(0) -> valor mínimo
    const s = vals.reduce((a,b)=>a+b,0);
    const avgRho = s/vals.length;
    return Math.log10(avgRho + 1e-9);  // Convertir a log10 para visualización
  });
  
  // Sample rows (convertir a log10 para visualización)
  const sampleRowsGenes = cols.map(s=>{
    const j = idx(DATA.samples, s);
    return xGenes.map((_,k)=>{
      const i = idxG[k];
      const rhoVal = DATA.z[i][j];
      return rhoVal > 0 ? Math.log10(rhoVal + 1e-9) : zmin_log;
    });
  });
  
  const zG = [meanRowGenes, ...sampleRowsGenes];
  const yRowsGenes = ['Average', ...cols];
  const geneHeight = Math.max(450, yRowsGenes.length * 20);
  
  Plotly.newPlot('hm',[{
    z:      zG,
    x:      xGenes,
    y:      yRowsGenes,
    type:   'heatmap',
    colorscale:'Viridis',
    zmin:   zmin_log,
    zmax:   zmax_log,
    colorbar: { len: 0.4, thickness: 18, y: 1.0, yanchor: 'top', title: 'log₁₀(ρ)' }
  }],{
    template:'plotly_white',height: geneHeight,
    xaxis:{tickangle:45,automargin:true, side:'top'},
    yaxis:{automargin:true, type:'category', categoryorder:'array', categoryarray:yRowsGenes, autorange:'reversed'}
  });
  
  // Event listener para clic en el heatmap
  document.getElementById('hm').on('plotly_click', function(data) {
    const point = data.points[0];
    const gene = point.x;
    const sample = point.y;
    
    // Si se hace clic en la fila "Average", mostrar todos los eventos de ese gen para TODOS los pacientes
    if (sample === 'Average') {
      const geneIdx = DATA.genes.indexOf(gene);
      if (geneIdx === -1) return;
      
      // Calcular ρ promedio para este gen (incluyendo valores 0)
      const geneRhoValues = DATA.z[geneIdx].filter(v => v != null);
      const avgRho = geneRhoValues.length > 0 ? geneRhoValues.reduce((a,b) => a+b, 0) / geneRhoValues.length : 0;
      
      // Obtener TODOS los eventos de este gen (de todas las muestras)
      const eventsForPopup = [];
      for (let i = 0; i < DATA.events.length; i++) {
        const event = DATA.events[i];
        const matchesGene = (gene === "Whole CDS") || (event.gene === gene);
        if (matchesGene && event.true === 1) {
          eventsForPopup.push({ event, meetsFilters: true });
        }
      }
      
      // Mostrar popup con columna de paciente
      showEventsPopupWithPatients(gene, "All samples", avgRho, eventsForPopup);
      return;
    }
    
    // Encontrar el índice del gen y muestra
    const geneIdx = DATA.genes.indexOf(gene);
    const sampleIdx = DATA.samples.indexOf(sample);
    
    if (geneIdx === -1 || sampleIdx === -1) return;
    
    // Obtener el valor de ρ directamente del heatmap (ya está calculado correctamente)
    const rhoValue = DATA.z[geneIdx][sampleIdx] || 0;
    
    // Obtener eventos de recombinación para mostrar en el popup
    const eventsForPopup = [];
    for (let i = 0; i < DATA.events.length; i++) {
      const event = DATA.events[i];
      // Si es "Whole CDS", incluir TODOS los eventos del paciente
      const matchesGene = (gene === "Whole CDS") || (event.gene === gene);
      if (matchesGene && event.sample === sample && event.true === 1) {
        eventsForPopup.push({ event, meetsFilters: true });
      }
    }
    
    // Crear popup con los eventos de recombinación
    showEventsPopup(gene, sample, rhoValue, eventsForPopup, minEvt, minMr, minMf);
  });
  // protein ρ heat-map (EXACTAMENTE COMO Rec_A.py - misma estructura que gene view)
  const idxRows = [];
  if (currentView === 'protein') {
    // Proteins heatmap with swapped axes and mean row on top (IGUAL QUE gene view)
    const labels     = DATA.proteins;
    const allZ       = DATA.proteins_z;

    const labelToIdx = Object.fromEntries(labels.map((l,i)=>[l,i]));
    const xProt = (pickPr.length ? pickPr : labels).filter(l => l in labelToIdx);
    const protIdxs = xProt.map(l => labelToIdx[l]);

    // Mean row ignoring thresholds (convertir a log10 para visualización)
    const meanRowProt = xProt.map((_,k)=>{
      const i = protIdxs[k];
      const vals = idxC.map(j=>allZ[i][j]).filter(v=>v!=null && v>0);  // Solo valores > 0
      if(vals.length===0) return zmin_log;  // log10(0) -> valor mínimo
      const s = vals.reduce((a,b)=>a+b,0);
      const avgRho = s/vals.length;
      return Math.log10(avgRho + 1e-9);  // Convertir a log10 para visualización
    });
    
    // Sample rows (convertir a log10 para visualización)
    const sampleRowsProt = cols.map(s=>{
      const j = idx(DATA.samples, s);
      return xProt.map((_,k)=>{
        const i = protIdxs[k];
        const rhoVal = allZ[i][j];
        return rhoVal > 0 ? Math.log10(rhoVal + 1e-9) : zmin_log;
      });
    });
    
    const zNR = [meanRowProt, ...sampleRowsProt];
    const yRowsProt = ['Average', ...cols];
    const proteinHeight = Math.max(450, yRowsProt.length * 20);
    Plotly.newPlot('cds',[{
      z:         zNR,
      x:         xProt,
      y:         yRowsProt,
      type:      'heatmap',
      colorscale:'Viridis',
      zmin:      zmin_log,
      zmax:      zmax_log,
      colorbar: { len: 0.4, thickness: 18, y: 1.0, yanchor: 'top', title: 'log₁₀(ρ)' }
    }],{
      template:'plotly_white',
      height: proteinHeight,
      xaxis:{tickangle:45,automargin:true, side:'top'},
      yaxis:{automargin:true, type:'category', categoryorder:'array', categoryarray:yRowsProt, autorange:'reversed'}
    });
    
    // Event listener para clic en el heatmap de proteínas
    document.getElementById('cds').on('plotly_click', function(data) {
      const point = data.points[0];
      const protein = point.x;
      const sample = point.y;
      
      // Si se hace clic en la fila "Average", mostrar todos los eventos de esa proteína para TODOS los pacientes
      if (sample === 'Average') {
        const proteinIdx = DATA.proteins.indexOf(protein);
        if (proteinIdx === -1) return;
        
        // Calcular ρ promedio para esta proteína (incluyendo valores 0)
        const proteinRhoValues = DATA.proteins_z[proteinIdx].filter(v => v != null);
        const avgRho = proteinRhoValues.length > 0 ? proteinRhoValues.reduce((a,b) => a+b, 0) / proteinRhoValues.length : 0;
        
        // Obtener TODOS los eventos de esta proteína (de todas las muestras)
        const eventsForPopup = [];
        for (let i = 0; i < DATA.events.length; i++) {
          const event = DATA.events[i];
          let matchesProtein = false;
          if (protein === "Whole CDS") {
            matchesProtein = true; // Todos los eventos
          } else {
            // Buscar en protein, region (NSPs) y gene (genes)
            matchesProtein = (event.protein === protein || event.region === protein || event.gene === protein);
          }
          if (matchesProtein && event.true === 1) {
            eventsForPopup.push({ event, meetsFilters: true });
          }
        }
        
        // Mostrar popup con columna de paciente
        showEventsPopupWithPatients(protein, "All samples", avgRho, eventsForPopup);
        return;
      }
      
      // Encontrar el índice de la proteína y muestra
      const proteinIdx = DATA.proteins.indexOf(protein);
      const sampleIdx = DATA.samples.indexOf(sample);
      
      if (proteinIdx === -1 || sampleIdx === -1) return;
      
      // Obtener el valor de ρ directamente del heatmap (ya está calculado correctamente)
      const rhoValue = DATA.proteins_z[proteinIdx][sampleIdx] || 0;
      
      // Obtener eventos de recombinación para mostrar en el popup
      const eventsForPopup = [];
      for (let i = 0; i < DATA.events.length; i++) {
        const event = DATA.events[i];
        // Verificar primero si es "Whole CDS" (todos los eventos del paciente)
        // Si no, filtrar por proteína específica
        let matchesProtein = false;
        if (protein === "Whole CDS") {
          matchesProtein = true; // Todos los eventos del paciente
        } else {
          // Buscar en protein, region (NSPs) y gene (genes)
          matchesProtein = (event.protein === protein || event.region === protein || event.gene === protein);
        }
        if (matchesProtein && event.sample === sample && event.true === 1) {
          eventsForPopup.push({ event, meetsFilters: true });
        }
      }
      
      // Crear popup con los eventos de recombinación
      showEventsPopup(protein, sample, rhoValue, eventsForPopup, 0, 0, 0);
    });
  } else {
    // For gene view, show all NSPs by default, filter only if specific ones selected (EXACTAMENTE COMO Rec_A.py)
    DATA.nsps.forEach((r,i)=>{
      if(pickG.length && !pickG.includes(r)) return;
      idxRows.push(i);
    });
    
    // Mean row ignoring thresholds (convertir a log10 para visualización)
    const meanRowNsps = idxRows.map(i=>{
      const vals = idxC.map(j=>DATA.nsps_z[i][j]).filter(v=>v!=null && v>0);  // Solo valores > 0
      if(vals.length===0) return zmin_log;  // log10(0) -> valor mínimo
      const avgRho = vals.reduce((a,b)=>a+b,0)/vals.length;
      return Math.log10(avgRho + 1e-9);  // Convertir a log10 para visualización
    });
    
    // Sample rows (convertir a log10 para visualización)
    const sampleRowsNsps = cols.map(s=>{
      const j = idx(DATA.samples, s);
      return idxRows.map(i => {
        const rhoVal = DATA.nsps_z[i][j];
        return rhoVal > 0 ? Math.log10(rhoVal + 1e-9) : zmin_log;
      });
    });
    
    const zNR = [meanRowNsps, ...sampleRowsNsps];
    const yRowsNsps = ['Average', ...cols];
    const nspHeight = Math.max(450, yRowsNsps.length * 20);
    
    Plotly.newPlot('cds',[{
      z:         zNR,
      x:         idxRows.map(i=>DATA.nsps[i]),
      y:         yRowsNsps,
      type:      'heatmap',
      colorscale:'Viridis',
      zmin:      zmin_log,
      zmax:      zmax_log,
      colorbar: { len: 0.4, thickness: 18, y: 1.0, yanchor: 'top', title: 'log₁₀(ρ)' }
    }],{
      template:'plotly_white',
      height: nspHeight,
      xaxis:{tickangle:45,automargin:true, side:'top'},
      yaxis:{automargin:true, type:'category', categoryorder:'array', categoryarray:yRowsNsps, autorange:'reversed'}
    });
    
    // Event listener para clic en el heatmap de NSPs (gene view)
    document.getElementById('cds').on('plotly_click', function(data) {
      const point = data.points[0];
      const nsp = point.x;
      const sample = point.y;
      
      // Skip if clicking on Average row
      if (sample === 'Average') return;
      
      // Encontrar el índice del NSP y muestra
      const nspIdx = DATA.nsps.indexOf(nsp);
      const sampleIdx = DATA.samples.indexOf(sample);
      
      if (nspIdx === -1 || sampleIdx === -1) return;
      
      // Obtener el valor de ρ directamente del heatmap (ya está calculado correctamente)
      const rhoValue = DATA.nsps_z[nspIdx][sampleIdx] || 0;
      
      // Obtener eventos de recombinación para mostrar en el popup
      const eventsForPopup = [];
      for (let i = 0; i < DATA.events.length; i++) {
        const event = DATA.events[i];
        if (event.region === nsp && event.sample === sample && event.true === 1) {
          eventsForPopup.push({ event, meetsFilters: true });
        }
      }
      
      // Crear popup con los eventos de recombinación
      showEventsPopup(nsp, sample, rhoValue, eventsForPopup, 0, 0, 0);
    });
  }
    // ── Genome recombination ratio by position (PRECALCULATED) ── COPIADO EXACTAMENTE DE Rec_A.py
    const N          = DATA.genome_len;
    const alpha      = DATA.alpha;
    const samplesSel = cols;
    
    // Use precalculated data instead of calculating on the fly
    const ratioArr = DATA.pos_ratio;  // Already precalculated in Python
    const genomeStats = DATA.genome_stats;  // Precalculated statistics
    
    // Create arrays for compatibility with existing code
    const evtCount = genomeStats.map(stat => stat.event_count);
    const patCount = genomeStats.map(stat => stat.patient_count);

    // Helpers gene/protein by position
    const gene_at = p => {
      for(let i=0;i<DATA.genes.length;i++){
        if(p>=DATA.genes_start[i] && p<DATA.genes_end[i]) return DATA.genes[i];
      }
      return '';
    };
    const nsp_at = p => {
      for(let i=0;i<DATA.nsps.length;i++){
        if(p>=DATA.nsps_start[i] && p<DATA.nsps_end[i]) return DATA.nsps[i];
      }
      return '';
    };

    // Prepare palette; shapes will be computed after yRange is defined
    const palette = ['#636efa','#EF553B','#00cc96','#ab63fa','#FFA15A',
                     '#19d3f3','#FF6692','#B6E880','#FF97FF','#FECB52'];
    let shapes = [];

    // Hover text simplified using precalculated data
    const hoverInfo = Array(N).fill('').map((_, p) => {
      const stat = genomeStats[p];
      return `Position: ${p}<br>Gene: ${stat.gene}<br>Protein: ${stat.nsp}<br>Patients: ${stat.patient_count}<br>Timepoints: ${stat.timepoint_count}<br>Events: ${stat.event_count}<br>Ratio: ${stat.ratio.toFixed(3)}`;
    });
    // Traces: Ratio without markers, with simplified hover
    const traceRatio = {
      name: 'Ratio', showlegend: true, mode: 'lines',
      x: [...Array(N).keys()], y: ratioArr,
      text: hoverInfo,
      hovertemplate: '%{text}<extra></extra>'
    };
    const traceGenes = DATA.genes.map((g,i)=>({
      name:g, showlegend:true,
      x:[null], y:[null], mode:'markers',
      marker:{symbol:'square', size:10, color:palette[i%palette.length]}
    }));
    // Calculate adaptive range for genome profile
    let xRange = [0, N-1];
    if (currentView === 'gene' && pickG.length > 0) {
      const selectedGenes = pickG;
      let minPos = N, maxPos = 0;
      
      selectedGenes.forEach(gene => {
        const geneIndex = DATA.genes.indexOf(gene);
        if (geneIndex !== -1) {
          const start = DATA.genes_start[geneIndex];
          const end = DATA.genes_end[geneIndex];
          if (start !== undefined && end !== undefined) {
            minPos = Math.min(minPos, start);
            maxPos = Math.max(maxPos, end);
          }
        }
      });
      
      if (minPos < N && maxPos > 0) {
        xRange = [Math.max(0, minPos - 1000), Math.min(N-1, maxPos + 1000)];
      }
    }
    
    // Calculate Y axis range to include bars at bottom (adaptive and thin bars)
    const maxRatio = Math.max(...ratioArr);
    const lower = maxRatio > 0 ? -0.1 * maxRatio : -0.05;
    const upper = maxRatio > 0 ? maxRatio * 1.05 : 0.05;
    const yRange = [lower, upper];
    const ySpan  = yRange[1] - yRange[0];
    const barH   = Math.max(upper * 0.01, ySpan * 0.02);
    const y0bar  = yRange[0] + ySpan * 0.015;
    // Now compute gene bars
    shapes  = DATA.genes_start.map((st,i)=>({
      type:'rect', x0:st, x1:DATA.genes_end[i],
      y0:y0bar, y1:y0bar+barH,
      fillcolor:palette[i%palette.length], line:{width:0}
    }));
    
    // Draw ratio with weights
    Plotly.newPlot('prof_ratio', [traceRatio, ...traceGenes], {
      template:'plotly_white', height:450,
      xaxis:{title:'Position', range: xRange},
      yaxis:{title:'Ratio',    range: yRange},
      shapes
    });

    // Click handler to show contributors modal
    const profDiv = document.getElementById('prof_ratio');
    profDiv.on('plotly_click', function(ev){
      const pt = ev.points && ev.points[0];
      if(!pt) return;
      const pos = Math.round(pt.x);
      const stat = genomeStats[pos];
      if(!stat) return;
      const modal = document.getElementById('contrib-modal');
      const posEl = document.getElementById('contrib-pos');
      const content = document.getElementById('contrib-content');
      posEl.textContent = pos.toString();
      if (stat.contributors && stat.contributors.length) {
        const html = stat.contributors.map(([pat, times])=>{
          // Ordenar tiempos numéricamente (extraer números y ordenar)
          const sortedTimes = times.slice().sort((a, b) => {
            const numA = parseFloat(a.toString().replace(/[^0-9.]/g, '')) || 0;
            const numB = parseFloat(b.toString().replace(/[^0-9.]/g, '')) || 0;
            return numA - numB;
          });
          return `<div><b>${pat}</b>: ${sortedTimes.join(', ')}</div>`;
        }).join('');
        content.innerHTML = html;
      } else {
        content.innerHTML = '<i>No positive contributors at this position.</i>';
      }
      modal.style.display = 'flex';
    });

    // ── GENOME PROFILE FOR PROTEINS ─────────────────────────────── COPIADO EXACTAMENTE DE Rec_A.py
    // Same ratio data but with protein annotations and bars
    const traceRatioNsps = {
      name: 'Ratio', showlegend: true, mode: 'lines',
      x: [...Array(N).keys()], y: ratioArr,
      text: hoverInfo,
      hovertemplate: '%{text}<extra></extra>'
    };
    const traceProteins = DATA.proteins.map((p,i)=>({
      name:p, showlegend:true,
      x:[null], y:[null], mode:'markers',
      marker:{symbol:'square', size:10, color:palette[i%palette.length]}
    }));
    
    // Create protein-specific shapes (bars) for protein view (use same yRange)
    const shapesProteins = [];
    DATA.proteins.forEach((protein, i) => {
      const start = DATA.proteins_start[i];
      const end = DATA.proteins_end[i];
      if (start !== undefined && end !== undefined) {
        shapesProteins.push({
          type: 'rect',
          x0: start, x1: end,
          y0: y0bar, y1: y0bar + barH,
          fillcolor: palette[i % palette.length],
          line: { width: 0 },
          layer: 'below'
        });
      }
    });
    
    // Calculate adaptive range for protein genome profile (use all proteins list)
    let xRangeProteins = [0, N-1];
    if (currentView === 'protein' && pickPr.length > 0) {
      const selectedProteins = pickPr;
      let minPos = N, maxPos = 0;
      
      selectedProteins.forEach(protein => {
        const proteinIndex = DATA.proteins.indexOf(protein);
        if (proteinIndex !== -1) {
          const start = DATA.proteins_start[proteinIndex];
          const end = DATA.proteins_end[proteinIndex];
          if (start !== undefined && end !== undefined) {
            minPos = Math.min(minPos, start);
            maxPos = Math.max(maxPos, end);
          }
        }
      });
      
      if (minPos < N && maxPos > 0) {
        xRangeProteins = [Math.max(0, minPos - 1000), Math.min(N-1, maxPos + 1000)];
      }
    }
    
    // Draw ratio with protein weights (use same Y range)
    Plotly.newPlot('prof_ratio_nsps', [traceRatioNsps, ...traceProteins], {
      template:'plotly_white', height:450,
      xaxis:{title:'Position', range: xRangeProteins},
      yaxis:{title:'Ratio',    range: yRange},
      shapes: shapesProteins
    });

    // Click handler to show contributors modal for protein view
    const profDivNsps = document.getElementById('prof_ratio_nsps');
    profDivNsps.on('plotly_click', function(ev){
      const pt = ev.points && ev.points[0];
      if(!pt) return;
      const pos = Math.round(pt.x);
      const stat = genomeStats[pos];
      if(!stat) return;
      const modal = document.getElementById('contrib-modal');
      const posEl = document.getElementById('contrib-pos');
      const content = document.getElementById('contrib-content');
      posEl.textContent = pos.toString();
      if (stat.contributors && stat.contributors.length) {
        const html = stat.contributors.map(([pat, times])=>{
          // Ordenar tiempos numéricamente (extraer números y ordenar)
          const sortedTimes = times.slice().sort((a, b) => {
            const numA = parseFloat(a.toString().replace(/[^0-9.]/g, '')) || 0;
            const numB = parseFloat(b.toString().replace(/[^0-9.]/g, '')) || 0;
            return numA - numB;
          });
          return `<div><b>${pat}</b>: ${sortedTimes.join(', ')}</div>`;
        }).join('');
        content.innerHTML = html;
      } else {
        content.innerHTML = '<i>No positive contributors at this position.</i>';
      }
      modal.style.display = 'flex';
    });


  // ── Prepara agrupación TIME SERIES – GENES y NSPS ─────────────────────
  const tsPats = pickP.length ? pickP : DATA.patients;
  // Siempre usar genes para el primer gráfico (ts_genes) independientemente de la vista
  const regsG  = pickG.length ? pickG : DATA.genes;
  // Para el segundo gráfico, usar proteínas solo cuando hay selección de proteínas
  // Si no hay selección, usar genes también para consistencia
  const regsR  = currentView === 'protein' ? (pickPr.length ? pickPr : DATA.genes) : (pickG.length ? pickG : DATA.nsps);

  const tracesG = tsPats.map(p => {
    const byT = {};
    DATA.samples
      .filter(s => s.startsWith(p + " |"))
      .forEach(s => {
        const num = parseFloat(s.split(" | ")[1].replace(/[^0-9.]/g, '')) || 0;
        (byT[num] ??= []).push(s);
      });
    const xs = Object.keys(byT).map(Number).sort((a, b) => a - b);
    const ys = xs.map(t => {
      const all = byT[t].flatMap(s => {
        const col = idx(DATA.samples, s);
        // AHORA usa regsG y DATA.z, convertir a log10 para visualización
        return regsG.map(g => {
          const rho = DATA.z[DATA.genes.indexOf(g)][col];
          return rho != null && rho > 0 ? Math.log10(rho + 1e-9) : null;
        }).filter(v => v != null);
      });
      return all.length ? all.reduce((a, b) => a + b, 0) / all.length : null;
    });
    return { x: xs, y: ys, name: p, mode: 'lines+markers' };
  });
  Plotly.newPlot('ts_genes', tracesG, {
    template: 'plotly_white', height: 300,
    xaxis: { title: 'Time' },
    yaxis: { title: 'mean log₁₀ρ – genes', zeroline: false }
  });

  // ── Time series – NSPS/PROTEINS ───────────────────────────────
  const tracesR = tsPats.map(p => {
    const byT = {};
    DATA.samples.filter(s => s.startsWith(p + " |")).forEach(s => {
      const num = parseFloat(s.split(" | ")[1].replace(/[^0-9.]/g, '')) || 0;
      (byT[num] ??= []).push(s);
    });
    const xs = Object.keys(byT).map(Number).sort((a, b) => a - b);
    const ys = xs.map(t => {
      const all = byT[t].flatMap(s => {
        const col = idx(DATA.samples, s);
        if (currentView === 'protein' && pickPr.length > 0) {
          // Solo usar proteins_z cuando hay proteínas seleccionadas
          return regsR.map(r => {
            const rho = DATA.proteins_z[DATA.proteins.indexOf(r)][col];
            return rho != null && rho > 0 ? Math.log10(rho + 1e-9) : null;
          }).filter(v => v != null);
        } else if (currentView === 'protein' && pickPr.length === 0) {
          // Si estamos en protein view pero no hay selección, usar genes (DATA.z)
          return regsR.map(r => {
            const rho = DATA.z[DATA.genes.indexOf(r)][col];
            return rho != null && rho > 0 ? Math.log10(rho + 1e-9) : null;
          }).filter(v => v != null);
        } else {
          // Gene view: usar nsps_z
          return regsR.map(r => {
            const rho = DATA.nsps_z[DATA.nsps.indexOf(r)][col];
            return rho != null && rho > 0 ? Math.log10(rho + 1e-9) : null;
          }).filter(v => v != null);
        }
      });
      return all.length ? all.reduce((a, b) => a + b, 0) / all.length : null;
    });
    return { x: xs, y: ys, name: p, mode: 'lines+markers' };
  });
  Plotly.newPlot('ts_nsps', tracesR, {
    template: 'plotly_white', height: 300,
    xaxis: { title: 'Time' },
    yaxis: { title: currentView === 'protein' ? 'mean log₁₀ρ – proteins' : 'mean log₁₀ρ – nsps', zeroline: false }
  });

 // ── BIPLOT PCA GENES ──────────────────────────────────────────
 // ── PCA REMOVED ───────────────────────────────────────────────────────────
 // PCA functionality has been completely removed as requested

  // ── ACTUALIZAR GRÁFICO DE VARIABILIDAD DEL GENOMA ─────────────────────────
  // Usar el perfil PRECALCULADO (como en Rec_A) para evitar perfiles binarios (0/α)
  // y mantener la suavidad de la curva. No recalculamos desde DATA.events aquí.
  const filteredRatioArr = ratioArr.slice();

  // Actualizar gráficos (+ ajustar eje Y para el nuevo rango)
  Plotly.restyle('prof_ratio', { y: [filteredRatioArr] }, [0]);
  const newMax = Math.max(...filteredRatioArr.filter(v => isFinite(v)));
  if (isFinite(newMax) && newMax > 0) {
    Plotly.relayout('prof_ratio', { 'yaxis.range': [-newMax*0.15, newMax*1.05] });
    Plotly.relayout('prof_ratio_nsps', { 'yaxis.range': [-newMax*0.15, newMax*1.05] });
  }
  Plotly.restyle('prof_ratio_nsps', { y: [filteredRatioArr] }, [0]);
}

// ── VIEW TOGGLE FUNCTIONALITY ───────────────────────────────────────────
let currentView = 'gene';

function switchView(view) {
  currentView = view;
  
  // Update button text
  const toggleBtn = document.getElementById('toggle');
  toggleBtn.textContent = view === 'gene' ? '🔄 Protein View' : '🔄 Gene View';
  
  // Show/hide plot sections
  document.getElementById('gene-plots').style.display = view === 'gene' ? 'block' : 'none';
  document.getElementById('protein-plots').style.display = view === 'protein' ? 'block' : 'none';
  
  // Show/hide filters based on view
  document.getElementById('gene-filter').style.display = view === 'gene' ? 'block' : 'none';
  document.getElementById('protein-filter').style.display = view === 'protein' ? 'block' : 'none';
  
  // Populate filters based on current view
  if (view === 'gene') {
    // Populate gene filter
    const selG = document.getElementById('selG');
    selG.innerHTML = '';
    DATA.genes.forEach(g => {
      const option = document.createElement('option');
      option.value = g;
      option.textContent = g;
      selG.appendChild(option);
    });
  } else {
    // Populate protein filter
    const selPr = document.getElementById('selPr');
    selPr.innerHTML = '';
    DATA.proteins.forEach(p => {
      const option = document.createElement('option');
      option.value = p;
      option.textContent = p;
      selPr.appendChild(option);
    });
  }
  
  // Redraw plots for the current view
  redraw();
}

// Event listeners for view toggle buttons
document.addEventListener('DOMContentLoaded', function() {
  // Toggle button
  document.getElementById('toggle').onclick = () => {
    // Toggle based on the global currentView to avoid DOM-dependent logic
    switchView(currentView === 'gene' ? 'protein' : 'gene');
  };

  // Code button
  document.getElementById('code-btn').onclick = () => {
    window.open('https://github.com/ldgonzalezvazquez/SRARec', '_blank');
  };

  // Methodology button
  document.getElementById('methodology-btn').onclick = () => {
    document.getElementById('methodology-modal').style.display = 'flex';
  };

  // Methodology close button
  document.getElementById('methodology-close').onclick = () => {
    document.getElementById('methodology-modal').style.display = 'none';
  };

  // Other buttons
  document.getElementById('go').onclick = redraw;
  document.getElementById('reset').onclick = resetFilters;
  document.getElementById('contrib-close').onclick = () => {
    document.getElementById('contrib-modal').style.display = 'none';
  };
  document.getElementById('events-close').onclick = () => {
    document.getElementById('events-modal').style.display = 'none';
  };
  document.getElementById('info-close').onclick = () => {
    document.getElementById('info-modal').style.display = 'none';
  };

  // Initialize gene filter
  const selG = document.getElementById('selG');
  DATA.genes.forEach(g => {
    const option = document.createElement('option');
    option.value = g;
    option.textContent = g;
    selG.appendChild(option);
  });
  
  // Initialize protein filter (hidden initially)
  const selPr = document.getElementById('selPr');
  DATA.proteins.forEach(p => {
    const option = document.createElement('option');
    option.value = p;
    option.textContent = p;
    selPr.appendChild(option);
  });
  
  // Initialize patient filter
  const selP = document.getElementById('selP');
  DATA.samples.forEach(s => {
    const option = document.createElement('option');
    option.value = s;
    option.textContent = s;
    selP.appendChild(option);
  });
  
  // Set initial view to gene view
  switchView('gene');
  
  // Initialize plots
  redraw();
});

// Info modal logic
const INFO_TEXT = {
  'min-events': '<b>Min events</b>: Minimum times the algorithm detects a recombination event at the same position for the same patient and time. Cells below this are masked (TO CHANGE).',
  'min-reads': '<b>Min reads</b>: Minimum reads of the less-supported gamete required to consider an event valid. The cell uses the maximum value across its events (TO CHANGE).',
  'min-freq': '<b>Min freq %</b>: Minimum frequency (%) of the less-supported gamete (TO CHANGE).'
};

document.querySelectorAll('.info-btn').forEach(btn=>{
  btn.addEventListener('click', ()=>{
    const key = btn.getAttribute('data-info');
    const modal = document.getElementById('info-modal');
    const title = document.getElementById('info-title');
    const content = document.getElementById('info-content');
    const titles = { 'min-events':'Min events', 'min-reads':'Min reads', 'min-freq':'Min freq %' };
    title.textContent = titles[key] || 'Filter info';
    content.innerHTML = INFO_TEXT[key] || '';
    modal.style.display = 'flex';
  });
});

// Funciones de ordenamiento para tablas
function sortTable(columnIndex) {
  const table = document.getElementById('events-table');
  if (!table) return;
  
  const tbody = table.querySelector('tbody') || table;
  const rows = Array.from(tbody.querySelectorAll('tr')).slice(1); // Excluir header
  
  const isAscending = table.getAttribute('data-sort-direction') !== 'asc';
  table.setAttribute('data-sort-direction', isAscending ? 'asc' : 'desc');
  
  rows.sort((a, b) => {
    const aVal = a.cells[columnIndex].textContent.trim();
    const bVal = b.cells[columnIndex].textContent.trim();
    
    // Para columnas numéricas (Times detected, Min Reads, Min Freq)
    if (columnIndex > 0) {
      const aNum = parseFloat(aVal.replace('%', ''));
      const bNum = parseFloat(bVal.replace('%', ''));
      return isAscending ? aNum - bNum : bNum - aNum;
    } else {
      // Para Event ID (string)
      return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
    }
  });
  
  // Reordenar filas
  rows.forEach(row => tbody.appendChild(row));
  
  // Actualizar indicadores de ordenamiento
  const headers = table.querySelectorAll('th');
  headers.forEach((header, index) => {
    if (index === columnIndex) {
      header.textContent = header.textContent.replace(/[↑↓↕]/g, '') + (isAscending ? ' ↑' : ' ↓');
    } else {
      header.textContent = header.textContent.replace(/[↑↓↕]/g, '') + ' ↕';
    }
  });
}

function sortTableWithPatients(columnIndex) {
  const table = document.getElementById('events-table-with-patients');
  if (!table) return;
  
  const tbody = table.querySelector('tbody') || table;
  const rows = Array.from(tbody.querySelectorAll('tr')).slice(1); // Excluir header
  
  const isAscending = table.getAttribute('data-sort-direction') !== 'asc';
  table.setAttribute('data-sort-direction', isAscending ? 'asc' : 'desc');
  
  rows.sort((a, b) => {
    const aVal = a.cells[columnIndex].textContent.trim();
    const bVal = b.cells[columnIndex].textContent.trim();
    
    // Para columnas numéricas (Times detected, Min Reads, Min Freq)
    if (columnIndex > 1) {
      const aNum = parseFloat(aVal.replace('%', ''));
      const bNum = parseFloat(bVal.replace('%', ''));
      return isAscending ? aNum - bNum : bNum - aNum;
    } else {
      // Para Patient y Event ID (strings)
      return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
    }
  });
  
  // Reordenar filas
  rows.forEach(row => tbody.appendChild(row));
  
  // Actualizar indicadores de ordenamiento
  const headers = table.querySelectorAll('th');
  headers.forEach((header, index) => {
    if (index === columnIndex) {
      header.textContent = header.textContent.replace(/[↑↓↕]/g, '') + (isAscending ? ' ↑' : ' ↓');
    } else {
      header.textContent = header.textContent.replace(/[↑↓↕]/g, '') + ' ↕';
    }
  });
}
</script>
</body></html>"""

    # ── WRITE HTML ───────────────────────────────────────────────────────────
    html = html_tpl.replace("__MAIL__", payload["contact"]) \
                   .replace("__DATA__", json.dumps(payload))
    
    with open(args.out, "w") as fh:
        fh.write(html)
    
    print(f"✔ HTML written: {args.out}")
    print(f"✔ Processing completed successfully")

if __name__ == "__main__":
    main()