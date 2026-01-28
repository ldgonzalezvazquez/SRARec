# -*- coding: utf-8 -*-
"""
SRARec_RecombinationHTML
------------------------

Módulo ligero para generar un HTML interactivo con la tasa de recombinación
a lo largo del genoma (perfil global) a partir de un results_txt de SRARec.

Pensado para ser llamado desde la GUI, sin dependencias raras:
- Usa únicamente: pandas, numpy, plotly.
- No define callbacks ni popups personalizados; sólo zoom/pan/hover nativos.

Ejemplo de uso desde la GUI:

    from SRARec_RecombinationHTML import generate_recombination_html
    import webbrowser, os

    html_path = generate_recombination_html("results_ravi.txt")
    webbrowser.open(f"file://{os.path.abspath(html_path)}")

Soporta dos formatos de entrada:
1) Tabla de eventos plana con columnas: pA, pB, true
2) Tabla "SRARec" con columnas Rec / No_Rec (strings codificando pares)
"""

import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional, List, Tuple
from tqdm import tqdm

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError as e:
    raise ImportError(
        "Missing dependency 'plotly'. Install it with:\n"
        "    pip install plotly\n"
        "or añádelo a requirements.txt del entorno de la GUI."
    ) from e

# Matplotlib es opcional, para generar figuras PNG
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None
    mpatches = None

# Biopython es opcional, sólo para GenBank
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

# ---------------------------------------------------------------------------
# Parsing helpers (basados en tu script original)
# ---------------------------------------------------------------------------

# Ejemplo de campo:
#   "100,200,1[mireads=50][minfreq=0.02];250,400,1[mireads=30][minfreq=0.01]"
_mire = re.compile(r"mireads=(\d+)")
_mfre = re.compile(r"minfreq=([0-9]*\.?[0-9]+)")

def _parse_pairs(field: str, expect_big: bool = True):
    """
    Parsea un campo tipo 'Rec' o 'No_Rec' en una lista de:
        (pA, pB, bigger, mireads, minfreq)

    Si el formato no cuadra, simplemente ignora entradas problemáticas.
    """
    if field is None:
        return []
    field = str(field).strip()
    if not field or field in ("NA", "nan"):
        return []

    out = []
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok:
            continue

        parts = tok.split(",")
        try:
            pA = int(parts[0])
            pB = int(parts[1])
        except (ValueError, IndexError):
            continue

        bigger = 1
        if expect_big and len(parts) >= 3:
            try:
                bigger = int(parts[2])
            except ValueError:
                bigger = 1

        mireads = 0
        m = _mire.search(tok)
        if m:
            try:
                mireads = int(m.group(1))
            except ValueError:
                mireads = 0

        minfreq = 0.0
        m = _mfre.search(tok)
        if m:
            try:
                minfreq = float(m.group(1))
            except ValueError:
                minfreq = 0.0

        out.append((pA, pB, bigger, mireads, minfreq))

    return out


def _load_events(results_path: str) -> pd.DataFrame:
    """
    Carga eventos de recombinación desde un results_txt.

    Soporta tres formatos:
    1) Tabla con columnas: pA, pB, true
    2) Tabla con columnas: Rec, No_Rec
    3) Formato SRARec clásico: línea con campos tipo:
         Sra:...  Rec:...  No_Rec:...
       separados por tabuladores.

    Devuelve DataFrame con columnas: pA (int), pB (int), true (0/1), sra (str, opcional).
    """
    if not os.path.exists(results_path):
        raise FileNotFoundError(f"Results file not found: {results_path}")

    # --- Pequeño "sniff" para evitar leer con pandas los resultados SRARec gigantes ---
    # Si detectamos un formato clásico de SRARec (Sra:/Rec:/No_Rec:), NO usamos read_csv,
    # y pasamos directamente al parser línea a línea (mucho más ligero en RAM).
    is_srarec_style = False
    try:
        with open(results_path, "r", encoding="utf-8", errors="ignore") as fh:
            for _line in fh:
                _line = _line.strip()
                if not _line or _line.startswith("#"):
                    continue
                if _line.startswith("Sra:") or "Rec:" in _line or "No_Rec:" in _line:
                    is_srarec_style = True
                break
    except OSError:
        # Si falla el sniff, simplemente seguimos como antes
        pass

    df = None
    if not is_srarec_style:
        # --- Intento 1: CSV/TSV estructurado ---
        try:
            df = pd.read_csv(results_path, sep=None, engine="python")
        except Exception:
            df = None

    if df is not None and not df.empty:
        # Caso 1: columnas pA, pB, true
        if {"pA", "pB", "true"}.issubset(df.columns):
            events = df[["pA", "pB", "true"]].copy()
            events["pA"] = pd.to_numeric(events["pA"], errors="coerce")
            events["pB"] = pd.to_numeric(events["pB"], errors="coerce")
            events["true"] = pd.to_numeric(events["true"], errors="coerce").fillna(0).astype(int)
            events = events.dropna(subset=["pA", "pB"])
            return events

        # Caso 2: columnas Rec / No_Rec
        if {"Rec", "No_Rec"}.issubset(df.columns):
            rows = []
            for _, row in df.iterrows():
                rec_field = row.get("Rec", "")
                norec_field = row.get("No_Rec", "")

                for (pA, pB, bigger, mireads, minfreq) in _parse_pairs(rec_field, expect_big=True):
                    rows.append((pA, pB, 1))

                for (pA, pB, bigger, mireads, minfreq) in _parse_pairs(norec_field, expect_big=False):
                    rows.append((pA, pB, 0))

            if rows:
                return pd.DataFrame(rows, columns=["pA", "pB", "true"])

    # --- Intento 2: formato SRARec tipo "clave:valor" por línea ---
    rows = []
    with open(results_path, "r", encoding="utf-8", errors="ignore") as f:
        # Contar líneas primero para la barra de progreso
        total_lines = sum(1 for _ in f)
        f.seek(0)  # Volver al inicio
        
        for line in tqdm(f, total=total_lines, desc="Loading events"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Campos separados por tabuladores: Sra:..., Coinf_detected:..., Rec:..., No_Rec:...
            parts = line.split("\t")
            rec_field = None
            norec_field = None
            sra_id = None
            for p in parts:
                p = p.strip()
                if p.startswith("Sra:"):
                    sra_id = p[len("Sra:"):].strip()
                if p.startswith("Rec:"):
                    rec_field = p[len("Rec:"):].strip()
                elif p.startswith("No_Rec:"):
                    norec_field = p[len("No_Rec:"):].strip()

            if rec_field:
                for (pA, pB, bigger, mireads, minfreq) in _parse_pairs(rec_field, expect_big=True):
                    rows.append((pA, pB, 1, sra_id))

            if norec_field:
                for (pA, pB, bigger, mireads, minfreq) in _parse_pairs(norec_field, expect_big=False):
                    rows.append((pA, pB, 0, sra_id))

    if rows:
        return pd.DataFrame(rows, columns=["pA", "pB", "true", "sra"])

    # Si nada ha funcionado:
    raise ValueError(
        "Unsupported results format.\n"
        "Expected one of:\n"
        " - columns: pA, pB, true\n"
        " - columns: Rec, No_Rec\n"
        " - or SRARec-style lines with 'Rec:' and 'No_Rec:' fields."
    )


def _count_unique_sras(results_path: str) -> int:
    """Cuenta el número de SRAs únicos en el archivo de resultados.
    
    Parameters
    ----------
    results_path : str
        Fichero de resultados de SRARec.
    
    Returns
    -------
    int
        Número de SRAs únicos encontrados.
    """
    sras = set()
    try:
        with open(results_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                for p in parts:
                    p = p.strip()
                    if p.startswith("Sra:"):
                        sra_id = p[len("Sra:"):].strip()
                        if sra_id:
                            sras.add(sra_id)
                        break
    except OSError:
        pass
    return len(sras)


def _count_unique_sras(results_path: str) -> int:
    """Cuenta el número de SRAs únicos en el archivo de resultados.
    
    Parameters
    ----------
    results_path : str
        Fichero de resultados de SRARec.
    
    Returns
    -------
    int
        Número de SRAs únicos encontrados.
    """
    sras = set()
    try:
        with open(results_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                for p in parts:
                    p = p.strip()
                    if p.startswith("Sra:"):
                        sra_id = p[len("Sra:"):].strip()
                        if sra_id:
                            sras.add(sra_id)
                        break
    except OSError:
        pass
    return len(sras)


def _compute_coverage_from_begfin(results_path: str, genome_len: int) -> np.ndarray:
    """Construye un perfil de profundidad a partir del campo BegFin_Reads.

    Parameters
    ----------
    results_path : str
        Fichero de resultados de SRARec.
    genome_len : int
        Longitud del genoma utilizada para el perfil de recombinación.

    Returns
    -------
    cov : np.ndarray
        Array de longitud genome_len con el número de lecturas que cubren cada posición.
        Si no se encuentra ningún campo BegFin_Reads, devuelve un array de ceros.
    """
    if genome_len is None or genome_len <= 0:
        return np.zeros(0, dtype=float)

    cov = np.zeros(int(genome_len) + 2, dtype=np.int64)

    try:
        with open(results_path, "r", encoding="utf-8", errors="ignore") as f:
            # Contar líneas para la barra de progreso
            total_lines = sum(1 for _ in f)
            f.seek(0)
            
            for line in tqdm(f, total=total_lines, desc="Computing coverage"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                begfin_field = None
                for p in parts:
                    p = p.strip()
                    if p.startswith("BegFin_Reads:"):
                        begfin_field = p[len("BegFin_Reads:"):].strip()
                        break

                if not begfin_field:
                    continue

                for tok in begfin_field.split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    try:
                        start_str, end_str = tok.split(",", 1)
                        start = int(start_str)
                        end = int(end_str)
                    except ValueError:
                        continue

                    # Asegurar límites dentro del genoma
                    if end < 0 or start >= genome_len:
                        continue
                    start = max(start, 0)
                    end = min(end, genome_len - 1)
                    if start > end:
                        continue

                    cov[start] += 1
                    cov[end + 1] -= 1
    except OSError as e:
        print(f"[WARN] Could not read BegFin_Reads from {results_path}: {e}")
        return np.zeros(0, dtype=float)

    cov = np.cumsum(cov[:-1])
    return cov.astype(float)


def _compute_begfin_concentration(results_path: str, genome_len: int) -> np.ndarray:
    """Construye un perfil de concentración de inicios y fin de reads a partir del campo BegFin_Reads.

    Cuenta cuántos reads empiezan o terminan en cada posición, en lugar de cuántos reads
    cubren esa posición (como hace _compute_coverage_from_begfin).

    Parameters
    ----------
    results_path : str
        Fichero de resultados de SRARec.
    genome_len : int
        Longitud del genoma utilizada para el perfil de recombinación.

    Returns
    -------
    concentration : np.ndarray
        Array de longitud genome_len con el número de reads que empiezan o terminan en cada posición.
        Si no se encuentra ningún campo BegFin_Reads, devuelve un array de ceros.
    """
    if genome_len is None or genome_len <= 0:
        return np.zeros(0, dtype=float)

    concentration = np.zeros(int(genome_len), dtype=np.int64)

    try:
        with open(results_path, "r", encoding="utf-8", errors="ignore") as f:
            # Contar líneas para la barra de progreso
            total_lines = sum(1 for _ in f)
            f.seek(0)
            
            for line in tqdm(f, total=total_lines, desc="Computing begfin concentration"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                begfin_field = None
                for p in parts:
                    p = p.strip()
                    if p.startswith("BegFin_Reads:"):
                        begfin_field = p[len("BegFin_Reads:"):].strip()
                        break

                if not begfin_field:
                    continue

                for tok in begfin_field.split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    try:
                        start_str, end_str = tok.split(",", 1)
                        start = int(start_str)
                        end = int(end_str)
                    except ValueError:
                        continue

                    # Asegurar límites dentro del genoma
                    if end < 0 or start >= genome_len:
                        continue
                    start = max(start, 0)
                    end = min(end, genome_len - 1)
                    if start > end:
                        continue

                    # Contar inicios y fin de reads
                    concentration[start] += 1  # Inicio del read
                    concentration[end] += 1     # Fin del read
    except OSError as e:
        print(f"[WARN] Could not read BegFin_Reads from {results_path}: {e}")
        return np.zeros(0, dtype=float)

    return concentration.astype(float)


# ---------------------------------------------------------------------------
# Carga de genes desde anotación NCBI (GenBank o GFF3)
# ---------------------------------------------------------------------------

def _load_genes(annotation_path: str, genome_len: Optional[int] = None):
    """
    Carga anotaciones génicas desde un archivo GenBank.

    Usa features de tipo 'CDS' (Coding Sequence) para obtener solo las regiones codificantes.
    Si un CDS tiene múltiples exones (join()), cada exón se carga como una región separada.
    Devuelve:
      - genes_df: DataFrame con columnas: start, end, name, track
      - genome_length: largo del genoma del GenBank

    El campo 'track' indica en qué fila se debe dibujar la región para evitar solapamientos.
    """
    if not annotation_path:
        return pd.DataFrame(), None
    if not os.path.exists(annotation_path):
        raise FileNotFoundError(f"Annotation file not found: {annotation_path}")

    apath = os.path.abspath(annotation_path)
    lower = apath.lower()

    genes_list = []
    genome_length = None

    # ---- GenBank ----
    if lower.endswith((".gb", ".gbk", ".gbff")):
        if SeqIO is None:
            raise ImportError(
                "Biopython is required to parse GenBank files. "
                "Instala con: pip install biopython"
            )
        try:
            record = SeqIO.read(apath, 'genbank')
            genome_length = len(record.seq)  # Largo real del genoma
            print(f"[DEBUG] Parsing GenBank record: {record.id}, length={genome_length}")
            
            for feat in record.features:
                if feat.type == 'CDS':  # Usar CDS en lugar de 'gene' para obtener solo regiones codificantes
                    try:
                        # Obtener el nombre del gen
                        name = feat.qualifiers.get('gene', [''])[0]
                        if not name or not name.strip():
                            # Intentar con locus_tag como alternativa
                            name = feat.qualifiers.get('locus_tag', [''])[0]
                            if not name or not name.strip():
                                name = f"CDS_{int(feat.location.start) + 1}"
                        
                        # Parsear las regiones codificantes (puede haber múltiples en un join())
                        # feat.location puede ser una SimpleLocation o CompoundLocation (join)
                        if hasattr(feat.location, 'parts'):
                            # Es un join() con múltiples partes
                            for part in feat.location.parts:
                                start = int(part.start) + 1  # 1-based como en Graphs.py
                                end = int(part.end)
                                
                                if end <= start:
                                    continue
                                
                                genes_list.append((start, end, name))
                                print(f"[DEBUG] Added CDS region: {name} ({start}-{end})")
                        else:
                            # Es una región simple
                            start = int(feat.location.start) + 1  # 1-based como en Graphs.py
                            end = int(feat.location.end)
                            
                            if end <= start:
                                print(f"[DEBUG] Skipping CDS with end <= start: {start}, {end}")
                                continue
                            
                            genes_list.append((start, end, name))
                            print(f"[DEBUG] Added CDS: {name} ({start}-{end})")
                    except Exception as e:
                        print(f"[DEBUG] Skipping CDS due to location error: {e}")
                        continue
                
        except Exception as e:
            print(f"[ERROR] Failed to parse GenBank file: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame(), None

    else:
        # Formato no soportado
        raise ValueError(
            f"Unsupported annotation format for genes: {annotation_path}\n"
            "Use GenBank (.gb/.gbk/.gbff)."
        )

    if not genes_list:
        print(f"[WARN] No CDS regions parsed from annotation file: {annotation_path}")
        return pd.DataFrame(columns=["start", "end", "name", "track"]), genome_length

    # Agrupar regiones por nombre de gen
    genes_dict = {}
    for start, end, name in genes_list:
        if name not in genes_dict:
            genes_dict[name] = []
        genes_dict[name].append((start, end))
    
    # Ordenar las regiones de cada gen por posición de inicio
    for name in genes_dict:
        genes_dict[name].sort(key=lambda x: x[0])
    
    # Crear lista de genes con todas sus regiones para asignar tracks
    genes_with_regions = []
    for name, regions in genes_dict.items():
        # Encontrar el inicio mínimo y fin máximo de todas las regiones del gen
        min_start = min(r[0] for r in regions)
        max_end = max(r[1] for r in regions)
        genes_with_regions.append((min_start, max_end, name, regions))
    
    # Ordenar genes por posición de inicio mínima
    genes_with_regions.sort(key=lambda x: x[0])
    
    # Asignar pistas considerando cada gen como una unidad (todas sus regiones en el mismo track)
    genes_with_track = _asignar_pistas_por_gen(genes_with_regions)
    
    # Expandir las regiones de vuelta a la lista plana con el track asignado
    expanded_regions = []
    for name, track, regions in genes_with_track:
        for start, end in regions:
            expanded_regions.append((start, end, name, track))
    
    gdf = pd.DataFrame(expanded_regions, columns=["start", "end", "name", "track"])

    print(f"[INFO] Loaded {len(gdf)} CDS regions from {len(genes_dict)} genes in {annotation_path}")

    return gdf, genome_length


def _asignar_pistas_por_gen(genes_with_regions):
    """
    Asigna pistas (tracks) a genes considerando cada gen como una unidad.
    Todas las regiones de un mismo gen comparten el mismo track.
    
    Input: lista de tuplas (min_start, max_end, name, regions)
           donde regions es una lista de tuplas (start, end)
    Output: lista de tuplas (name, track, regions)
    """
    pistas_fin = []  # último "max_end" de cada pista
    genes_con_pista = []
    
    for min_start, max_end, name, regions in genes_with_regions:
        assigned = False
        for i, ultimo_end in enumerate(pistas_fin):
            if min_start > ultimo_end:  # el gen cabe en la pista i
                pistas_fin[i] = max_end
                genes_con_pista.append((name, i, regions))
                assigned = True
                break
        
        if not assigned:  # ninguna pista libre → crea una nueva
            pistas_fin.append(max_end)
            genes_con_pista.append((name, len(pistas_fin) - 1, regions))
    
    return genes_con_pista

# ---------------------------------------------------------------------------
# Cálculo del perfil de recombinación global (tipo gráfico 2)
# ---------------------------------------------------------------------------

def compute_recombination_profile(events: pd.DataFrame,
                                  genome_len: int = None,
                                  alpha: float = 1.0,
                                  exclusive_positions: bool = False,
                                  aggregate_by_patient: bool = False):
    """
    Calcula el perfil global de recombinación (ρ) a lo largo del genoma.

    Sigue la idea del script grande:
      - Para cada evento (pA, pB), con pA < pB:
          * contribuye con peso 1/(b-a) a todas las posiciones [a, b)
          * se separan eventos true (1) y false (0)
      - ρ(i) = alpha * true(i) / (true(i) + false(i)) cuando hay eventos.

    Parámetros
    ----------
    events : DataFrame con columnas pA, pB, true.
    genome_len : int, opcional
        Si no se da, se usa max(pA, pB) + 1.
    alpha : float
        Factor de escala (por defecto 1.0, usa el mismo que en tu script si quieres).
    exclusive_positions : bool
        Si True, para cada posición se toma el máximo valor entre rangos true/false
        y, si hay ambos, se conserva true y se descarta false.
    aggregate_by_patient : bool
        Si True, calcula valores por paciente (Sra) y luego suma T/F por posición.

    Returns
    -------
    positions : np.ndarray de ints
    rho       : np.ndarray de floats
    """
    if events.empty:
        raise ValueError("No events available to compute recombination profile.")

    # Limpiar y asegurar tipos
    ev = events.copy()
    ev["pA"] = pd.to_numeric(ev["pA"], errors="coerce")
    ev["pB"] = pd.to_numeric(ev["pB"], errors="coerce")
    ev["true"] = ev["true"].astype(int)
    ev = ev.dropna(subset=["pA", "pB"])
    ev = ev[ev["pB"] > ev["pA"]]

    if ev.empty:
        raise ValueError("No valid events after filtering pA < pB.")

    if genome_len is None:
        genome_len = int(max(ev["pA"].max(), ev["pB"].max())) + 1

    genome_len = max(genome_len, 1)

    pos_t = np.zeros(genome_len, dtype=float)
    pos_f = np.zeros(genome_len, dtype=float)

    if aggregate_by_patient and "sra" in ev.columns:
        # Calcular por paciente y sumar T/F por posición
        ev_clean = ev[ev["sra"].notna()]
        ev_sra_groups = ev_clean.groupby("sra")
        for sra_id, ev_sra in tqdm(ev_sra_groups, total=ev["sra"].nunique(),
                                   desc="Computing profile (per patient)"):
            s_pos_t = np.zeros(genome_len, dtype=float)
            s_pos_f = np.zeros(genome_len, dtype=float)
            for _, row in ev_sra.iterrows():
                a = int(row["pA"])
                b = int(row["pB"])
                if b <= a or a < 0 or b > genome_len:
                    continue
                w = 1.0 / float(b - a)
                if row["true"] == 1:
                    if exclusive_positions:
                        s_pos_t[a:b] = np.maximum(s_pos_t[a:b], w)
                    else:
                        s_pos_t[a:b] += w
                else:
                    if exclusive_positions:
                        s_pos_f[a:b] = np.maximum(s_pos_f[a:b], w)
                    else:
                        s_pos_f[a:b] += w
            if exclusive_positions:
                s_pos_f = np.where(s_pos_t > 0.0, 0.0, s_pos_f)
            pos_t += s_pos_t
            pos_f += s_pos_f
    else:
        # Sumar contribuciones globales
        for _, row in tqdm(ev.iterrows(), total=len(ev), desc="Computing profile"):
            a = int(row["pA"])
            b = int(row["pB"])
            if b <= a or a < 0 or b > genome_len:
                continue
            w = 1.0 / float(b - a)
            if row["true"] == 1:
                if exclusive_positions:
                    pos_t[a:b] = np.maximum(pos_t[a:b], w)
                else:
                    pos_t[a:b] += w
            else:
                if exclusive_positions:
                    pos_f[a:b] = np.maximum(pos_f[a:b], w)
                else:
                    pos_f[a:b] += w
        if exclusive_positions:
            # Si hay ambos estados en una posición, conservar true
            pos_f = np.where(pos_t > 0.0, 0.0, pos_f)

    tot = pos_t + pos_f
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(tot > 0.0, alpha * pos_t / tot, 0.0)

    positions = np.arange(genome_len, dtype=int)
    return positions, rho

# ---------------------------------------------------------------------------
# Generación de figura matplotlib
# ---------------------------------------------------------------------------

def _generate_matplotlib_figure(
    positions: np.ndarray,
    rho: np.ndarray,
    coverage: Optional[np.ndarray] = None,
    begfin_concentration: Optional[np.ndarray] = None,
    genes_df: Optional[pd.DataFrame] = None,
    output_png: str = None,
    title: str = "Genome-wide recombination ratio"
):
    """
    Genera una figura matplotlib similar al script proporcionado.
    
    Parameters
    ----------
    positions : np.ndarray
        Posiciones del genoma.
    rho : np.ndarray
        Valores de ρ (ratio de recombinación).
    coverage : np.ndarray, opcional
        Valores de cobertura (BegFin_Reads).
    genes_df : pd.DataFrame, opcional
        DataFrame con genes (columnas: start, end, name, track).
    output_png : str, opcional
        Ruta de salida para el PNG.
    title : str
        Título del gráfico.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("[WARN] Matplotlib not available, skipping PNG generation.")
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    # Asegurar que la rejilla quede por detrás de las líneas de datos
    ax.set_axisbelow(True)
    
    # Dibujar línea de recombinación
    # Usar positions[1:] para empezar desde posición 1 (como en el script original)
    x = positions[1:] if len(positions) > 1 else positions
    rho_plot = rho[1:] if len(rho) > 1 else rho
    
    max_rho = float(np.max(rho_plot)) if rho_plot.size else 0.0
    if max_rho <= 0:
        ymax = 0.05
    else:
        # Máximo dinámico sin sobrepasar 1.0, con +5% de margen visual
        ymax = min(max_rho, 1.0)
    y_axis_top = min(ymax * 1.05, 1.0)
    
    # Estilos consistentes con el HTML: azul encima, roja/gris visibles sin eclipsar
    PNG_BLUE_WIDTH = 1.1
    PNG_RED_WIDTH = 1.2
    PNG_GRAY_WIDTH = 1.0
    PNG_BLUE_ALPHA = 1.0
    PNG_RED_ALPHA = 0.45
    PNG_GRAY_ALPHA = 0.55

    # Línea azul de recombinación (por encima de las demás)
    ax.plot(x, rho_plot, color="#5C6BC0", linewidth=PNG_BLUE_WIDTH, alpha=PNG_BLUE_ALPHA, zorder=3)
    
    # Añadir línea de cobertura si está disponible (por encima de la línea azul)
    if coverage is not None and coverage.size > 0:
        # Escalar la cobertura para que se vea en el mismo rango
        coverage_plot = coverage[1:] if len(coverage) > 1 else coverage
        cov_len = min(len(coverage_plot), len(x))
        if cov_len > 0:
            coverage_scaled = coverage_plot[:cov_len]
            x_cov = x[:cov_len]
            max_cov = float(np.max(coverage_scaled)) if coverage_scaled.size > 0 else 0.0
            if max_cov > 0 and max_rho > 0:
                # Escalar para que el máximo de cobertura coincida con el máximo de rho
                scale_factor = max_rho / max_cov
                coverage_scaled = coverage_scaled * scale_factor
            ax.plot(x_cov, coverage_scaled, color="#555555", linewidth=PNG_GRAY_WIDTH, linestyle="--",
                   label="Analysis depth", alpha=PNG_GRAY_ALPHA, zorder=2)
    
    # Añadir línea de concentración de inicios/fin de reads si está disponible (por encima de la línea azul)
    if begfin_concentration is not None and begfin_concentration.size > 0:
        # Escalar la concentración para que se vea en el mismo rango
        concentration_plot = begfin_concentration[1:] if len(begfin_concentration) > 1 else begfin_concentration
        conc_len = min(len(concentration_plot), len(x))
        if conc_len > 0:
            concentration_scaled = concentration_plot[:conc_len]
            x_conc = x[:conc_len]
            max_conc = float(np.max(concentration_scaled)) if concentration_scaled.size > 0 else 0.0
            if max_conc > 0 and max_rho > 0:
                # Escalar para que el máximo de concentración coincida con el máximo de rho
                scale_factor = max_rho / max_conc
                concentration_scaled = concentration_scaled * scale_factor
            ax.plot(x_conc, concentration_scaled, color="#FF0000", linewidth=PNG_RED_WIDTH, linestyle="--",
                   label="Read start/end density", alpha=PNG_RED_ALPHA, zorder=2)
    
    # Dibujar genes si están disponibles
    # Usar altura fija para la pista de genes que no dependa del rango del eje Y
    # La altura se calcula como una fracción fija del rango del eje Y para mantener consistencia visual
    num_tracks = 0
    unique_genes = []
    gene_colors = {}
    # Altura fija de la pista de genes en unidades de datos (fracción fija del rango)
    gene_track_fraction = 0.05  # 5% del rango del eje Y para la altura total de genes (más gordas)
    gene_track_spacing_fraction = 0.02  # 2% adicional para el espaciado (menos separación)
    gene_bottom_offset_fraction = 0.05  # 5% adicional del rango del eje Y para separar los genes del 0
    
    if genes_df is not None and not genes_df.empty:
        # Obtener colores únicos para cada gen
        unique_genes = genes_df['name'].unique()
        try:
            # Usar la nueva API de matplotlib (3.7+)
            import matplotlib
            if hasattr(matplotlib, 'colormaps'):
                cmap_genes = matplotlib.colormaps['Set3'].resampled(len(unique_genes))
            else:
                cmap_genes = plt.cm.get_cmap("Set3", len(unique_genes))
        except (AttributeError, KeyError):
            # Fallback para versiones antiguas de matplotlib
            cmap_genes = plt.cm.get_cmap("Set3", len(unique_genes))
        gene_colors = {name: cmap_genes(i) for i, name in enumerate(unique_genes)}
        
        num_tracks = int(genes_df['track'].max()) + 1 if 'track' in genes_df.columns else 1
        # Calcular altura de cada track como fracción fija del rango del eje Y
        total_gene_height = ymax * gene_track_fraction
        track_height = total_gene_height / max(num_tracks, 1)
        track_spacing = ymax * gene_track_spacing_fraction
        
        # Agrupar regiones por gen y track para añadir líneas discontinuas
        genes_by_name_track = {}
        for idx, g in genes_df.iterrows():
            start = int(g["start"])
            end = int(g["end"])
            name = str(g.get("name", f"gene_{idx}")).strip() or f"gene_{idx}"
            track = int(g.get("track", 0))
            
            key = (name, track)
            if key not in genes_by_name_track:
                genes_by_name_track[key] = []
            genes_by_name_track[key].append((start, end))
        
        # Ordenar regiones por cada gen
        for key in genes_by_name_track:
            genes_by_name_track[key].sort(key=lambda x: x[0])
        
        # Dibujar genes agrupados
        for (name, track), regions in genes_by_name_track.items():
            color = gene_colors.get(name, "gray")
            bottom_offset = ymax * gene_bottom_offset_fraction
            y_bottom = -(bottom_offset + track_spacing + track * (track_height + track_spacing))
            y_center = y_bottom + (track_height / 2)
            
            # Dibujar cada región codificante
            for start, end in regions:
                width = end - start
                ax.broken_barh(
                    [(start, width)],
                    (y_bottom, track_height),
                    facecolors=color,
                    edgecolors="black",
                    linewidth=1.0,  # Más gordas las líneas
                    zorder=0  # Por debajo del cuadro
                )
            
            # Añadir líneas discontinuas entre regiones del mismo gen
            if len(regions) > 1:
                for i in range(len(regions) - 1):
                    # Línea desde el final de la región i hasta el inicio de la región i+1
                    end_prev = regions[i][1]
                    start_next = regions[i+1][0]
                    
                    # Línea discontinua
                    ax.plot([end_prev, start_next], [y_center, y_center], 
                           color=color, linestyle='--', linewidth=0.5, alpha=0.7, zorder=0)  # Por debajo del cuadro
    
    x_min = 1
    x_max = positions[-1] if len(positions) > 0 else 29903
    # Calcular y_min basado en la altura fija de los genes
    # Necesitamos calcular el espacio necesario para todos los tracks
    if num_tracks > 0:
        bottom_offset = ymax * gene_bottom_offset_fraction
        # Calcular el espacio necesario para el último track
        last_track_y_bottom = -(bottom_offset + track_spacing + (num_tracks - 1) * (track_height + track_spacing))
        # y_min debe incluir el último track completo más un pequeño margen
        y_min = last_track_y_bottom - track_height - (ymax * 0.01)  # Pequeño margen adicional
    else:
        y_min = -ymax * 0.02  # Pequeño margen si no hay genes
    
    ax.set_xlim(x_min, x_max)
    # IMPORTANTE: Establecer ylim ANTES de establecer los ticks para evitar que los ticks fuercen el rango
    ax.set_ylim(y_min, y_axis_top)
    
    # Ticks hasta el máximo visible (con tope 1.0)
    num_ticks = 6
    y_ticks = np.linspace(0.0, min(1.0, y_axis_top), num_ticks)
    ax.set_yticks(y_ticks)
    
    # Formatear los ticks para mostrar máximo 2 decimales
    from matplotlib.ticker import FuncFormatter
    def format_y_tick(value, pos):
        return f'{value:.2f}'
    ax.yaxis.set_major_formatter(FuncFormatter(format_y_tick))
    
    # Asegurar que el límite del eje Y se mantenga en ymax (no permitir que matplotlib lo ajuste)
    ax.set_ylim(y_min, y_axis_top)
    ax.set_xlabel("Position")
    ax.set_ylabel("FGT Ratio")
    ax.set_title(title)
    # Rejilla por detrás de las líneas (axisbelow=True)
    ax.grid(True, alpha=0.3)
    
    # Crear leyenda combinada: primero las líneas (recombinación y cobertura), luego los genes
    handles, labels = ax.get_legend_handles_labels()
    
    # Leyenda de genes si están disponibles
    if genes_df is not None and not genes_df.empty:
        gene_patches = []
        for gname in unique_genes:
            patch = mpatches.Patch(
                facecolor=gene_colors[gname],
                edgecolor="black",
                linewidth=0.3,
                label=gname
            )
            gene_patches.append(patch)
        
        # Combinar handles: primero las líneas, luego los genes
        all_handles = handles + gene_patches
        all_labels = labels + [gname for gname in unique_genes]
        
        # Leyenda debajo del eje X (position)
        gene_legend = ax.legend(
            handles=all_handles,
            labels=all_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),  # Subir un poco la leyenda
            ncol=len(all_handles),
            fontsize="small"
        )
        plt.subplots_adjust(top=0.90, bottom=0.25)  # Más espacio abajo para Position
    else:
        # Si no hay genes, solo mostrar leyenda de líneas debajo del eje X
        if handles:
            ax.legend(handles=handles, labels=labels, loc="upper center", 
                     bbox_to_anchor=(0.5, -0.10), fontsize="small", ncol=len(handles))
            plt.subplots_adjust(bottom=0.20)  # Más espacio abajo para Position
    
    if output_png:
        fig.savefig(output_png, dpi=300, bbox_inches="tight")
        print(f"[INFO] Matplotlib figure saved to: {output_png}")
    else:
        plt.show()
    
    plt.close(fig)


# ---------------------------------------------------------------------------
# Generación del HTML interactivo
# ---------------------------------------------------------------------------

def generate_recombination_html(
    results_path: str,
    output_html: str = None,
    alpha: float = 1.0,
    genome_len: int = None,
    title: str = "Genome-wide recombination ratio",
    annotation_path: Optional[str] = None,
    gene_track_label: str = "Genes",
    skip_coverage: bool = False,
    skip_begfin_concentration: bool = False,
    exclusive_positions: bool = False,
    aggregate_by_patient: bool = False
) -> str:
    """
    Genera un HTML con un gráfico interactivo (Plotly) de ρ a lo largo del genoma.

    - Soporta zoom, pan, selección de rango, exportar como PNG, etc.
    - No incluye callbacks de clic personalizados ni popups extra.

    Parámetros
    ----------
    results_path : str
        Ruta al results_txt de SRARec.
    output_html : str, opcional
        Ruta de salida. Si es None, usa "<basename>_recomb_profile.html".
    alpha : float
        Factor de escala para ρ (usa el mismo valor que en tu Graphs.py si quieres).
    genome_len : int, opcional
        Largo del genoma. Si None, se infiere del máximo pA/pB.
    title : str
        Título del gráfico.
    annotation_path : str, opcional
        Ruta a archivo NCBI de anotación (GenBank .gb/.gbff o GFF3 .gff/.gff3).
        Si se proporciona, se crea una pista de genes debajo del perfil de recombinación.
    gene_track_label : str
        Etiqueta mostrada para la pista de genes.
    exclusive_positions : bool
        Si True, no permite solapamientos true/false por posición y
        conserva el máximo valor por posición (true domina a false).
    aggregate_by_patient : bool
        Si True, agrega por paciente antes de sumar T/F por posición.

    Returns
    -------
    output_html : str
        Ruta al archivo HTML generado.
    """
    # Genes opcionales (si se proporciona anotación) - CARGAR PRIMERO para obtener genome_len
    genes_df = None
    gb_genome_len = None
    if annotation_path:
        try:
            genes_df, gb_genome_len = _load_genes(annotation_path, genome_len=genome_len)
            if genes_df is not None and genes_df.empty:
                print("[INFO] Annotation provided but no genes found; skipping gene track.")
                genes_df = None
        except Exception as e:
            # No abortar la figura si falla la anotación: sólo avisar por consola
            print(f"[WARN] Could not load genes from {annotation_path}: {e}")
            genes_df = None
    
    # Usar genome_len del GenBank si está disponible, sino el especificado
    if gb_genome_len is not None:
        genome_len = gb_genome_len
        print(f"[INFO] Using genome length from GenBank: {genome_len}")
    
    # Contar SRAs únicos para el título
    num_sras = _count_unique_sras(results_path)
    if title == "Genome-wide recombination ratio":
        title = f"Genome-wide recombination ratio ({num_sras} SRAs)"
    
    events = _load_events(results_path)
    if exclusive_positions and not aggregate_by_patient and "sra" in events.columns:
        if events["sra"].notna().any():
            aggregate_by_patient = True
            print("[INFO] exclusive_positions activo: agregando por paciente para evitar 0/1 global.")
    positions, rho = compute_recombination_profile(events,
                                                   genome_len=genome_len,
                                                   alpha=alpha,
                                                   exclusive_positions=exclusive_positions,
                                                   aggregate_by_patient=aggregate_by_patient)

    L = int(positions[-1]) + 1 if len(positions) else genome_len

    # Crear mapa de posición a genes para el hover (puede haber múltiples genes por posición)
    genes_at_pos = {}
    if genes_df is not None:
        for idx, g in genes_df.iterrows():
            start = int(g["start"])
            end = int(g["end"])
            name = str(g.get("name", f"gene_{idx}")).strip() or f"gene_{idx}"
            for pos in range(start, end):
                if pos < L:
                    if pos not in genes_at_pos:
                        genes_at_pos[pos] = []
                    genes_at_pos[pos].append(name)
    
    # Crear texto hover personalizado con anotación de genes (múltiples si aplica)
    hover_text = []
    for i, pos in enumerate(positions):
        gene_names = genes_at_pos.get(int(pos), [])
        if gene_names:
            genes_str = ", ".join(gene_names)
            hover_text.append(f"Position: {pos}<br>Gene: {genes_str}<br>ρ: {rho[i]:.4f}")
        else:
            hover_text.append(f"Position: {pos}<br>ρ: {rho[i]:.4f}")

    # Perfil de profundidad (coverage) a partir de BegFin_Reads
    cov_positions = None
    cov_millions = None
    cov_hover = None
    max_cov_millions = 0.0
    coverage_full = None  # Mantener el coverage completo para matplotlib y archivos

    if not skip_coverage:
        coverage_full = _compute_coverage_from_begfin(results_path, genome_len=genome_len)
    else:
        print("[INFO] Skipping coverage computation (--no-coverage flag enabled)")
    if coverage_full is not None and coverage_full.size:
        cov_len = min(len(coverage_full), len(positions))
        if cov_len > 0:
            coverage = coverage_full[:cov_len]
            cov_positions = positions[:cov_len]
            cov_millions = coverage / 1e6
            max_cov_millions = float(np.max(cov_millions)) if cov_millions.size else 0.0
            cov_hover = [
                f"Position: {int(pos)}<br>Coverage: {cov/1e6:.3f}M reads"
                for pos, cov in zip(cov_positions, coverage)
            ]

    # Perfil de concentración de inicios/fin de reads a partir de BegFin_Reads
    begfin_conc_positions = None
    begfin_conc_millions = None
    begfin_conc_hover = None
    max_begfin_conc_millions = 0.0
    begfin_concentration_full = None  # Mantener la concentración completa para matplotlib y archivos

    if not skip_begfin_concentration:
        begfin_concentration_full = _compute_begfin_concentration(results_path, genome_len=genome_len)
    else:
        print("[INFO] Skipping begfin concentration computation (--no-begfin-concentration flag enabled)")
    if begfin_concentration_full is not None and begfin_concentration_full.size:
        conc_len = min(len(begfin_concentration_full), len(positions))
        if conc_len > 0:
            begfin_concentration = begfin_concentration_full[:conc_len]
            begfin_conc_positions = positions[:conc_len]
            begfin_conc_millions = begfin_concentration / 1e6
            max_begfin_conc_millions = float(np.max(begfin_conc_millions)) if begfin_conc_millions.size else 0.0
            begfin_conc_hover = [
                f"Position: {int(pos)}<br>Read start/end density: {conc/1e6:.3f}M reads"
                for pos, conc in zip(begfin_conc_positions, begfin_concentration)
            ]

    # Figura principal (perfil de recombinación)
    ROW_HEIGHTS = [0.78, 0.22]
    VERTICAL_SPACING = 0.06

    num_tracks = 0
    use_gene_track = False
    if genes_df is not None and 'track' in genes_df.columns and not genes_df.empty:
        num_tracks = int(genes_df['track'].max()) + 1
        use_gene_track = num_tracks > 0

    if use_gene_track:
        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            row_heights=ROW_HEIGHTS,
            vertical_spacing=VERTICAL_SPACING
        )
    else:
        fig = make_subplots(rows=1, cols=1, shared_xaxes=True)

    # Ajustes de visualización para asegurar lectura clara en vista completa
    RED_OVERVIEW_THRESHOLD = 10000
    RED_OVERVIEW_TARGET_BINS = 2000
    GL_THRESHOLD = 12000
    BLUE_LINE_WIDTH = 2.0
    RED_LINE_WIDTH = 1.2
    RED_OPACITY = 0.45
    GRAY_LINE_WIDTH = 1.0
    GRAY_OPACITY = 0.55

    # Definir rango de Y (no negativos en el eje) y calcular posición de barras FUERA del gráfico
    max_rho = float(np.max(rho)) if rho.size else 0.0
    if max_rho <= 0:
        ymax = 0.05
        ymin = 0
    else:
        # El máximo del eje Y debe ser un poco más que el máximo de la línea de recombinación
        ymax = max_rho * 1.1
        ymin = 0

    # Añadir traza de profundidad si la hemos podido calcular
    if cov_positions is not None and cov_millions is not None and max_cov_millions > 0:
        # Escalar la cobertura para que ambos ejes muestren los mismos valores
        if max_rho > 0:
            scale_factor = max_rho / max_cov_millions
        else:
            scale_factor = 1.0
        cov_scaled = cov_millions * scale_factor
        cov_scatter = go.Scattergl if len(cov_positions) > GL_THRESHOLD else go.Scatter
        fig.add_trace(cov_scatter(
            x=cov_positions,
            y=cov_scaled,
            mode="lines",
            line=dict(width=GRAY_LINE_WIDTH, dash="dot", color="#555555"),
            opacity=GRAY_OPACITY,
            name="Analysis depth",
            legendrank=2,
            yaxis="y",  # Usar el mismo eje Y principal, sin segundo eje
            hovertemplate="%{text}<extra></extra>",
            text=cov_hover,
            showlegend=True
        ), row=1, col=1)
    
    # Añadir traza de concentración de inicios/fin de reads si la hemos podido calcular
    if begfin_conc_positions is not None and begfin_conc_millions is not None and max_begfin_conc_millions > 0:
        # Escalar la concentración para que ambos ejes muestren los mismos valores
        if max_rho > 0:
            scale_factor = max_rho / max_begfin_conc_millions
        else:
            scale_factor = 1.0
        begfin_conc_scaled = begfin_conc_millions * scale_factor
        red_scatter = go.Scattergl if len(begfin_conc_positions) > GL_THRESHOLD else go.Scatter

        fig.add_trace(red_scatter(
            x=begfin_conc_positions,
            y=begfin_conc_scaled,
            mode="lines",
            line=dict(width=RED_LINE_WIDTH, color="#FF0000"),
            opacity=RED_OPACITY,
            name="Read start/end density",
            legendrank=3,
            yaxis="y",  # Usar el mismo eje Y principal, sin segundo eje
            hovertemplate="%{text}<extra></extra>",
            text=begfin_conc_hover,
            showlegend=True
        ), row=1, col=1)
    
    # Barras de genes con múltiples pistas (tracks) para solapamientos
    shapes = []
    legend_traces = []
    palette = ['#636efa', '#EF553B', '#00cc96', '#ab63fa', '#FFA15A', '#19d3f3',
               '#FF6692', '#B6E880', '#FF97FF', '#FECB52', '#7A5195', '#1F77B4']

    gene_colors = {}  # Diccionario para asignar colores por nombre de gen
    genes_in_legend = set()  # Para evitar duplicados en la leyenda

    if use_gene_track and genes_df is not None and 'track' in genes_df.columns:
        # Obtener nombres únicos de genes y asignar colores
        unique_genes = genes_df['name'].unique()
        for i, gene_name in enumerate(unique_genes):
            gene_colors[gene_name] = palette[i % len(palette)]
        
        # Agrupar regiones por gen y track para añadir líneas discontinuas
        genes_by_name_track = {}
        for idx, g in genes_df.iterrows():
            start = int(g["start"])
            end = int(g["end"])
            name = str(g.get("name", f"gene_{idx}")).strip() or f"gene_{idx}"
            track = int(g.get("track", 0))
            
            key = (name, track)
            if key not in genes_by_name_track:
                genes_by_name_track[key] = []
            genes_by_name_track[key].append((start, end))
        
        # Ordenar regiones por cada gen
        for key in genes_by_name_track:
            genes_by_name_track[key].sort(key=lambda x: x[0])
        
        # Geometría fija en el eje de genes (y2), estable entre datasets
        GENE_BAND_TOP = 0.90
        GENE_BAND_BOTTOM = 0.10
        TRACK_SPACING_RATIO = 0.35
        total_height = GENE_BAND_TOP - GENE_BAND_BOTTOM
        track_height = total_height / max(num_tracks + (num_tracks - 1) * TRACK_SPACING_RATIO, 1)
        track_spacing = track_height * TRACK_SPACING_RATIO
        
        for (name, track), regions in genes_by_name_track.items():
            color = gene_colors.get(name, palette[0])
            # Bajar más las barras para que estén bien separadas de las líneas de datos
            # Layout: líneas arriba -> espacio -> genes -> espacio -> Position -> espacio -> leyenda
            y_top = GENE_BAND_TOP - (track * (track_height + track_spacing))
            y_bottom = y_top - track_height
            y_center = (y_top + y_bottom) / 2.0
            
            # Dibujar cada región codificante
            for start, end in regions:
                shapes.append(dict(
                    type="rect",
                    xref="x2",
                    yref="y2",
                    x0=start,
                    x1=end,
                    y0=y_bottom,
                    y1=y_bottom + track_height,
                    fillcolor=color,
                    line=dict(width=1.0, color="white"),  # Más gordas las líneas
                    layer="below"
                ))
            
            # Añadir líneas discontinuas entre regiones del mismo gen
            if len(regions) > 1:
                for i in range(len(regions) - 1):
                    # Línea desde el final de la región i hasta el inicio de la región i+1
                    end_prev = regions[i][1]
                    start_next = regions[i+1][0]
                    
                    # Línea discontinua simple
                    shapes.append(dict(
                        type="line",
                        xref="x2",
                        yref="y2",
                        x0=end_prev,
                        x1=start_next,
                        y0=y_center,
                        y1=y_center,
                        line=dict(width=0.5, color=color, dash="dash"),  # Línea discontinua más fina
                        layer="below"
                    ))

            # Leyenda: solo añadir una vez por gen (no por región)
            if name not in genes_in_legend:
                legend_traces.append(go.Scatter(
                    x=[None], y=[None],
                    mode="markers",
                    marker=dict(symbol="square", size=10, color=color),
                    name=name,
                    showlegend=True,
                    hoverinfo="skip"
                ))
                genes_in_legend.add(name)

    for trace in legend_traces:
        fig.add_trace(trace, row=1, col=1)

    # Azul siempre arriba: añadir al final para máxima visibilidad
    blue_scatter = go.Scattergl if len(positions) > GL_THRESHOLD else go.Scatter
    fig.add_trace(blue_scatter(
        x=positions,
        y=rho,
        mode="lines",
        line=dict(width=BLUE_LINE_WIDTH, color="#303F9F"),
        opacity=1.0,
        text=hover_text,
        hovertemplate="%{text}<extra></extra>",
        name="Recombination ratio",
        showlegend=True,
        legendrank=1
    ), row=1, col=1)

    # Layout estable y sin solapes
    X_TITLE_STANDOFF = 30
    Y_TITLE_STANDOFF = 15

    n_items = len(legend_traces) + 1  # azul + genes
    if cov_positions is not None and cov_millions is not None and max_cov_millions > 0:
        n_items += 1
    if begfin_conc_positions is not None and begfin_conc_millions is not None and max_begfin_conc_millions > 0:
        n_items += 1

    entrywidth_fraction = 0.11
    items_per_row = max(1, int(1.0 / entrywidth_fraction))
    legend_rows = max(1, int(np.ceil(n_items / items_per_row)))
    MARGIN = dict(l=80, r=30, t=70, b=140 + (legend_rows * 28))
    LEGEND_Y = -0.18

    # --- Compatibilidad Plotly ---
    # `ticklabelstandoff` no existe en algunas versiones (p.ej. 5.18.0),
    # así que lo intentamos y, si falla, reintentamos sin esa propiedad.
    xaxis_bottom_cfg = dict(
        title=dict(text="Position", standoff=X_TITLE_STANDOFF),
        automargin=True,
        ticklabelstandoff=6,
        showgrid=True,
        zeroline=False,
        rangemode="tozero"
    )
    xaxis_top_cfg = dict(
        title=None,
        showticklabels=False,
        showgrid=True,
        zeroline=False,
        rangemode="tozero"
    )
    yaxis_main_cfg = dict(
        title="FGT Ratio",
        title_standoff=Y_TITLE_STANDOFF,
        automargin=True,
        range=[ymin, ymax],
        showgrid=True,
        zeroline=True,
        rangemode="nonnegative"
    )
    yaxis_gene_cfg = dict(
        range=[0, 1],
        showgrid=False,
        zeroline=False,
        visible=False,
        fixedrange=True
    )

    layout_kwargs = dict(
        template="plotly_white",
        title=dict(text=title, x=0.5, font=dict(size=18)),
        hovermode="closest",
        margin=MARGIN,
        height=700,
        shapes=shapes,
        legend=dict(
            orientation="h",
            x=0,
            xanchor="left",
            y=LEGEND_Y,
            yanchor="top",
            entrywidth=entrywidth_fraction,
            entrywidthmode="fraction",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#999",
            borderwidth=1,
            font=dict(size=12)
        ),
        annotations=[],
    )

    try:
        fig.update_layout(**layout_kwargs)
        if use_gene_track:
            fig.update_yaxes(**yaxis_main_cfg, row=1, col=1)
            fig.update_yaxes(**yaxis_gene_cfg, row=2, col=1)
            fig.update_xaxes(**xaxis_top_cfg, row=1, col=1)
            fig.update_xaxes(**xaxis_bottom_cfg, row=2, col=1)
        else:
            fig.update_yaxes(**yaxis_main_cfg, row=1, col=1)
            fig.update_xaxes(**xaxis_bottom_cfg, row=1, col=1)
    except ValueError as e:
        if "ticklabelstandoff" in str(e) or "entrywidth" in str(e):
            xaxis_bottom_cfg.pop("ticklabelstandoff", None)
            xaxis_top_cfg.pop("ticklabelstandoff", None)
            layout_kwargs["legend"].pop("entrywidth", None)
            layout_kwargs["legend"].pop("entrywidthmode", None)
            layout_kwargs["legend"]["itemwidth"] = 100
            layout_kwargs["legend"]["itemsizing"] = "constant"
            fig.update_layout(**layout_kwargs)
            if use_gene_track:
                fig.update_yaxes(**yaxis_main_cfg, row=1, col=1)
                fig.update_yaxes(**yaxis_gene_cfg, row=2, col=1)
                fig.update_xaxes(**xaxis_top_cfg, row=1, col=1)
                fig.update_xaxes(**xaxis_bottom_cfg, row=2, col=1)
            else:
                fig.update_yaxes(**yaxis_main_cfg, row=1, col=1)
                fig.update_xaxes(**xaxis_bottom_cfg, row=1, col=1)
        else:
            raise

    # Determinar base para nombres de archivos
    base = os.path.splitext(os.path.basename(results_path))[0]
    
    # Añadir timestamp si no se especifica output_html (para identificar resultados más recientes)
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    use_timestamp = (output_html is None)
    
    if output_html is None:
        output_html = f"{base}_recomb_profile_{timestamp}.html"

    # full_html=True + include_plotlyjs='cdn' → archivo autocontenido (salvo CDN), perfecto para abrir desde la GUI
    fig.write_html(output_html, full_html=True, include_plotlyjs="cdn")
    
    # Generar figura matplotlib
    if use_timestamp:
        output_png = f"{base}_recomb_profile_{timestamp}.png"
    else:
        # Si se especificó output_html, usar el mismo patrón para PNG
        png_base = os.path.splitext(output_html)[0]
        output_png = f"{png_base}.png"
    # Preparar coverage para matplotlib (en millones, igual que en Plotly)
    coverage_for_matplotlib = None
    if coverage_full is not None and coverage_full.size > 0:
        coverage_for_matplotlib = coverage_full / 1e6  # Convertir a millones
    # Preparar concentración de begfin para matplotlib (en millones, igual que en Plotly)
    begfin_conc_for_matplotlib = None
    if begfin_concentration_full is not None and begfin_concentration_full.size > 0:
        begfin_conc_for_matplotlib = begfin_concentration_full / 1e6  # Convertir a millones
    _generate_matplotlib_figure(
        positions=positions,
        rho=rho,
        coverage=coverage_for_matplotlib,
        begfin_concentration=begfin_conc_for_matplotlib,
        genes_df=genes_df,
        output_png=output_png,
        title=title
    )
    
    # Guardar archivos con valores
    # 1) Valores de recombinación por posición
    if use_timestamp:
        recomb_output = f"{base}_recombination_values_{timestamp}.txt"
    else:
        recomb_output = f"{base}_recombination_values.txt"
    with open(recomb_output, "w") as f:
        f.write("Position\tRecombination_Ratio\n")
        for pos, r in zip(positions, rho):
            f.write(f"{pos}\t{r:.6f}\n")
    print(f"[INFO] Recombination values saved to: {recomb_output}")
    
    # 2) Valores de cobertura (BegFin_Reads) por posición
    if coverage_full is not None and coverage_full.size > 0:
        if use_timestamp:
            coverage_output = f"{base}_coverage_values_{timestamp}.txt"
        else:
            coverage_output = f"{base}_coverage_values.txt"
        # Asegurar que tenemos las posiciones correctas para el coverage
        cov_pos_for_file = positions[:len(coverage_full)] if len(coverage_full) <= len(positions) else positions
        coverage_for_file = coverage_full[:len(cov_pos_for_file)]
        with open(coverage_output, "w") as f:
            f.write("Position\tCoverage\n")
            for pos, cov in zip(cov_pos_for_file, coverage_for_file):
                f.write(f"{int(pos)}\t{cov:.6f}\n")
        print(f"[INFO] Coverage values saved to: {coverage_output}")
    
    return output_html


if __name__ == "__main__":
    # Uso rápido en CLI:
    import argparse
    script_dir = Path(__file__).resolve().parent
    default_results = script_dir / "results_slow.txt"
    default_annotation = script_dir / "NC_045512.2.gb"
    output_dir = script_dir / "output"

    p = argparse.ArgumentParser(
        description="Generate interactive genome-wide recombination profile HTML from SRARec results."
    )
    p.add_argument("results_txt", nargs="?", default=str(default_results),
                   help="Path to SRARec results txt (default: results_slow.txt in script dir)")
    p.add_argument("-o", "--output", help="Output HTML file (default: ./output/<results>_recomb_profile.html)")
    p.add_argument("--alpha", type=float, default=1.0,
                   help="Scaling factor for ρ (default: 1.0)")
    p.add_argument("--genome-len", type=int, default=None,
                   help="Genome length (if omitted, inferred from events)")
    p.add_argument("--annotation", type=str, default=str(default_annotation),
                   help="NCBI annotation file (GenBank .gb/.gbff or GFF3 .gff/.gff3) "
                        "to display genes below the recombination profile "
                        "(default: NC_045512.2.gb in script dir).")
    p.add_argument("--gene-track-label", type=str, default="Genes",
                   help="Label for the gene track (default: 'Genes').")
    p.add_argument("--no-coverage", action="store_true",
                   help="Skip coverage computation from BegFin_Reads (faster execution)")
    p.add_argument("--no-begfin-concentration", action="store_true",
                   help="Skip begfin concentration computation (concentration of read starts/ends)")
    p.add_argument("--exclusive-positions", action="store_true",
                   help="Use max-per-position and prefer true over false to avoid overlaps")
    p.add_argument("--aggregate-by-patient", action="store_true",
                   help="Aggregate per patient before summing T/F per position")
    args = p.parse_args()

    results_path = Path(args.results_txt)
    if not results_path.is_file():
        raise FileNotFoundError(f"[ERROR] Missing required file: {results_path}")

    annotation_path = Path(args.annotation) if args.annotation else None
    if annotation_path is None or not annotation_path.is_file():
        raise FileNotFoundError(f"[ERROR] Missing required file: {annotation_path}")

    if args.output:
        output_html = Path(args.output)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_html = output_dir / f"{results_path.stem}_recomb_profile.html"

    out = generate_recombination_html(
        results_path=str(results_path),
        output_html=str(output_html),
        alpha=args.alpha,
        genome_len=args.genome_len,
        annotation_path=str(annotation_path),
        gene_track_label=args.gene_track_label,
        skip_coverage=args.no_coverage,
        skip_begfin_concentration=args.no_begfin_concentration,
        exclusive_positions=args.exclusive_positions,
        aggregate_by_patient=args.aggregate_by_patient,
    )
    print(f"[INFO] Recombination profile HTML written to: {out}")