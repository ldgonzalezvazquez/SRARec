#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import subprocess
from pathlib import Path
import shutil
import pysam
from collections import defaultdict, Counter
import pprint
import tempfile
import time
import itertools
import random
from tqdm import tqdm
import argparse
import gzip
import re
from pathlib import Path

############################################
# Configuración Global
############################################
min_base_qual = 20
min_fraction_align_reads=0.9
############################
#Variables de ejecución
#########################
MAX_PIPELINE_ATTEMPTS = 3
MAX_DOWNLOAD_ATTEMPTS = 3
VERBOSE = False



# No guardar rangos inicio/fin por read (ahorra RAM/tiempo)
SAVE_READS_MINMAX = True
if sys.getrecursionlimit() < 20000:
    sys.setrecursionlimit(20000)


############################################
# Funciones Auxiliares
############################################

def print_header():
    # Códigos de colores ANSI
    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    MAGENTA = "\033[95m"
    CYAN = "\033[96m"
    BOLD = "\033[1m"
    RESET = "\033[0m"

    header = f"""
    {BOLD}{CYAN}███████╗██████═╗ █████╗ ███{GREEN}███═╗██████╗  ██████╗ 
    {BOLD}{CYAN}██╔════╝██╔══██║██╔══██╗██╔═{GREEN}═██║██╔══██╗██╔═══██╗   
    {BOLD}{CYAN}███████╗██████╔╝███████║████{GREEN}██╔╝███████║██║       
    {BOLD}{CYAN}╚════██║██╔═██╚╗██╔══██║██╔═{GREEN}██╚╗██╔════╗██║      
    {BOLD}{CYAN}███████║██║ ║██║██║  ██║██║ {GREEN}║██║║██████║╚██████╗ 
    {BOLD}{CYAN}╚══════╝╚═╝ ╚══╝╚═╝  ╚═╝╚═╝ {GREEN}╚══╝╚══════╝ ╚═════╝ 

    {BOLD}{YELLOW}Program:{RESET} {CYAN}SRA{GREEN}Rec{RESET}
    {BOLD}{YELLOW}Authors:{RESET} {MAGENTA}Luis Daniel González-Vázquez & Darren P. Martin{RESET}
    {BOLD}{YELLOW}Email:{RESET} {BLUE}luisdaniel.gonzalez@uvigo.es{RESET}
    {BOLD}{YELLOW}Year:{RESET} {RED}2025{RESET}
    {BOLD}{YELLOW}Version:{RESET} {CYAN}1.{GREEN}2{RESET}
    """

    print(header)

def print_step(message):
    """Imprime un paso principal en color blanco."""
    print(f"\033[1;37m{message}\033[0m")  # Blanco brillante

def print_substep(message):
    """Imprime un subpaso (comando ejecutado) en color azul."""
    print(f"\033[1;34m{message}\033[0m")  # Azul

def print_output(message):
    """Imprime la salida del programa ejecutado en color verde."""
    print(f"\033[1;32m{message}\033[0m")  # Verde

def print_error(message):
    """Imprime un error en color rojo."""
    print(f"\033[1;31m[ERROR] {message}\033[0m")

def print_verbose(message):
    """Imprime mensajes extras sólo si VERBOSE = True."""
    if VERBOSE:
        # Puedes elegir otro color, por ejemplo magenta (35)
        print(f"\033[1;35m[VERBOSE] {message}\033[0m")

def identify_platform_and_layout(srr):
    """
    Identifica la plataforma de secuenciación y el tipo de lectura (Single End o Paired End) para un SRR dado.
    Retorna un diccionario con el nombre de la plataforma y el layout.
    """
    try:
        cmd = f"esearch -db sra -query {srr} | efetch -format runinfo"
        ###print_verbose(f"Identificando plataforma con comando: {cmd}")
        runinfo = subprocess.check_output(cmd, shell=True, text=True).strip()

        lines = runinfo.split('\n')
        if len(lines) < 2:
            raise ValueError("No se encontró información en RunInfo.")

        headers = lines[0].split(',')
        data = lines[1].split(',')

        platform_index = headers.index("Platform") if "Platform" in headers else -1
        layout_index = headers.index("LibraryLayout") if "LibraryLayout" in headers else -1

        if platform_index == -1 or layout_index == -1:
            raise ValueError("No se encontraron las columnas 'Platform' o 'LibraryLayout' en RunInfo.")

        platform = data[platform_index].strip()
        layout = data[layout_index].strip()

        print(f"[INFO] Plataforma identificada: {platform}")
        print(f"[INFO] Tipo de lectura: {layout}")

        return {"platform": platform, "layout": layout}
    except Exception as e:
        print_error(f"[ERROR] Could not identify sequencing platform and layout: {e}")
        sys.exit(1)

def retry_command(max_attempts, cmd, shell=True):
    """
    Ejecuta un comando con reintentos en caso de fallo.
    """
    attempt = 1
    while attempt <= max_attempts:
        try:
            print_substep(f"Ejecutando: {cmd}")
            process = subprocess.Popen(
                cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = process.communicate() #espera a que el proceso termine y captura su salida (stdout) y errores (stderr)
            print_output(stdout)
            if process.returncode != 0: #process.returncode devuelve el código de salida del comando
                # Mostrar también el stderr para depuración
                print_output(f"[STDERR]: {stderr}") 
                raise subprocess.CalledProcessError(process.returncode, cmd) # Subporcess es el módulo que realiza el subproceso
            return True
        except subprocess.CalledProcessError:
            print_substep(f"Attempt {attempt} of {max_attempts} failed. Retrying...")
            attempt += 1
    print_error(f"Command failed after {max_attempts} attempts: {cmd}")
    return False

def filter_bam_by_alignment_fraction(input_bam, output_bam, min_fraction=0.7):
    """
    Filtra el BAM de entrada y escribe un nuevo BAM (output_bam) que solo contiene
    aquellas lecturas cuyo porcentaje de bases alineadas (query_alignment_length / query_length)
    sea al menos min_fraction (por defecto, 0.5 o 50%).

    Parámetros:
      input_bam (str): Ruta del BAM de entrada.
      output_bam (str): Ruta del BAM filtrado de salida.
      min_fraction (float): Fracción mínima de mapeo para conservar la lectura.
    """
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    for read in in_bam.fetch():
        # Si la lectura no está mapeada, se omite.
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        
        # read.query_alignment_length: longitud de la parte alineada (excluye inserciones y clips)
        # read.query_length: longitud total de la lectura
        if read.query_length > 0:
            fraction = float(read.query_alignment_length) / float(read.query_length)
        else:
            fraction = 0.0
        
        if fraction >= min_fraction:
            out_bam.write(read)
    
    in_bam.close()
    out_bam.close()


def run_bwa_mem(refseq, fastq1, fastq2, output_sam, output_bam, platform,min_fraction_align_reads):
    """
    Ejecuta bwa mem (o minimap2, según plataforma) para generar un archivo SAM
    y luego lo convierte a BAM. 
    """
    try:
        # Definir el read group correctamente
        read_group = f"@RG\\tID:{Path(fastq1).stem}\\tSM:{Path(fastq1).stem}\\tLB:Lib\\tPL:ILLUMINA" #Path(fastq1).stem es una expresión de Python que usa el módulo pathlib para extraer el nombre base del archivo sin extensión

        # Ajuste según la plataforma
        if platform == "ILLUMINA":
            if fastq2 is not None:
                map_cmd = f"bwa mem -M -R '{read_group}' {refseq} {fastq1} {fastq2} > {output_sam}"
            else:
                map_cmd = f"bwa mem -M -R '{read_group}' {refseq} {fastq1} > {output_sam}"
        elif platform in ["OXFORD NANOPORE", "ONT", "OXFORD_NANOPORE"]:
            if fastq2 is not None:  # Paired End
                map_cmd = f"minimap2 -ax map-ont -R '{read_group}' -t 8 {refseq} {fastq1} {fastq2} > {output_sam}"
            else:  # Single End
                map_cmd = f"minimap2 -ax map-ont -R '{read_group}' -t 8 {refseq} {fastq1} > {output_sam}"
        elif platform in ["PACBIO","PACBIO_SMRT"]:
            if fastq2 is not None:  # Paired End
                map_cmd = f"minimap2 -ax map-pb -R '{read_group}' -t 8 {refseq} {fastq1} {fastq2} > {output_sam}"
            else:  # Single End
                map_cmd = f"minimap2 -ax map-pb -R '{read_group}' -t 8 {refseq} {fastq1} > {output_sam}"
        else:
            raise ValueError(f"Plataforma no soportada: {platform}") #El comando raise en Python lanza una excepción manualmente. Se usa para indicar que ha ocurrido un error que impide que el programa continúe

        print_substep(f"Ejecutando: {map_cmd}")

        # Ejecutar el mapeo para generar el archivo .sam
        result = subprocess.run(map_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            # Mostramos también stderr si falla
            print_output(f"[STDERR]: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, map_cmd)

        print(f"[INFO] SAM file generated: {output_sam}")

        # Convertir el archivo .sam a .bam
        # Paso 1: Convertir SAM a BAM
        intermediate_bam = str(output_bam).replace(".bam", "_intermediate.bam")
        samtools_view_cmd = f"samtools view -q 30 -bS {output_sam} -o {intermediate_bam}"
        print_substep(f"Ejecutando: {samtools_view_cmd}")
        result_view = subprocess.run(samtools_view_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result_view.returncode != 0:
            print_output(f"[STDERR]: {result_view.stderr}")
            raise subprocess.CalledProcessError(result_view.returncode, samtools_view_cmd)

        #Paso 2: Filtrar read por calidad de alineamiento (cobertura del read con la referencia)
        intermediate_bam_00 = str(output_bam).replace(".bam", "_intermediate_00.bam")
        sort_cmd = f"samtools sort -o {intermediate_bam_00} {intermediate_bam} "
        retry_command(MAX_DOWNLOAD_ATTEMPTS, sort_cmd)
        subprocess.run(f"samtools index {intermediate_bam_00}", shell=True, check=True)
        intermediate_bam_0 = str(output_bam).replace(".bam", "_intermediate_0.bam")
        filter_bam_by_alignment_fraction(intermediate_bam_00,intermediate_bam_0,min_fraction_align_reads)


        # Paso 2: Ordenar el BAM por coordenada (importante para mpileup)
        # <-- ADAPTADO (quitamos "-n")
        if fastq2 is not None:
            intermediate_bam_1 = str(output_bam).replace(".bam", "_intermediate_1.bam")
            sort_paired_cmd = f"samtools sort -n -o {intermediate_bam_1} {intermediate_bam_0} "
            retry_command(MAX_DOWNLOAD_ATTEMPTS, sort_paired_cmd)
            rm_intermediate_bam_cmd=f"rm {intermediate_bam_0}"
            subprocess.run(rm_intermediate_bam_cmd, shell=True, stderr=subprocess.PIPE, text=True)
            fixmate_cmd = f"samtools fixmate -m {intermediate_bam_1} {intermediate_bam_0}"
            retry_command(MAX_DOWNLOAD_ATTEMPTS, fixmate_cmd)
            rm_intermediate_bam_cmd_2=f"rm {intermediate_bam_1}"
            subprocess.run(rm_intermediate_bam_cmd_2, shell=True, stderr=subprocess.PIPE, text=True)
        samtools_sort_cmd = f"samtools sort -o {output_bam} {intermediate_bam_0}"
        print_substep(f"Ejecutando: {samtools_sort_cmd}")
        result_sort = subprocess.run(samtools_sort_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result_sort.returncode != 0:
            print_output(f"[STDERR]: {result_sort.stderr}")
            raise subprocess.CalledProcessError(result_sort.returncode, samtools_sort_cmd)

        print_substep(f"[INFO] Sorted BAM file created: {output_bam}")

    except subprocess.CalledProcessError as e: #as e → Captura la excepción en la variable e para acceder a sus atributos
        print(f"[ERROR] Error executing: {e.cmd}\nExit code: {e.returncode}") #e.cmd → El comando que causó el error
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        sys.exit(1)

#BLOQUE FUNCIÓN SCAN_REGION_FOR_POLYMORPHISMS - HAY QUE OPTIMIZAR MUCHO ESTO. TENGO EJEMPLOS DE CÓMO HACERLO EN CHATGPT (SALTO DE REFERENCIA EN PILEUP):

def _get_base_or_indel_at_pos_indel(read, pos, ref_fasta=None, ref_name=None, noindel=False):
    """
    Devuelve (evento, calidad) para la posición 'pos' (0-based) de la referencia.
    
    - Si noindel=True, se IGNORAN inserciones y deleciones:
        * Solo se busca si existe una base en la lectura (qPos is not None).
        * Si la hay, se devuelve esa base y su calidad.
        * Si no la hay o hay deleción, se devuelve (None, None).
    
    - Si noindel=False, la lógica es la de siempre:
        * Primero buscar inserción "entre pos-1 y pos" (código igual a antes).
        * Luego buscar base/deleción en pos.
        * Combinar si ambas existen: p.ej. "A+3ACT", "-1A+3ACT".
    """

    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    # ----------------------------------------------------------------
    # MODO "noindel" => SOLO BASES
    # ----------------------------------------------------------------
    if noindel:
        for qPos, rPos in aligned_pairs:
            if rPos == pos:
                # Si en esta posición la lectura aporta base (qPos no es None),
                # devolvemos (base, calidad). Ignoramos inserciones/deleciones.
                if qPos is not None:
                    base_str  = read.query_sequence[qPos].upper()
                    base_qual = read.query_qualities[qPos]
                    return (base_str, base_qual)
                else:
                    # Hay deleción -> ignorarla
                    pass
        # Si llegamos aquí, no había base en la posición => devolvemos None
        return (None, None)

    # ----------------------------------------------------------------
    # MODO NORMAL => DETECTAR INSERCIONES + BASE/DELECIÓN
    # ----------------------------------------------------------------
    insertion_str  = None
    insertion_qual = None

    # (1) Buscar INSERCIÓN entre pos-1 y pos
    for i in range(len(aligned_pairs) - 1):
        q1, r1 = aligned_pairs[i]
        q2, r2 = aligned_pairs[i + 1]

        # "Entre pos-1 y pos": r1 == pos y r2 is None...
        if r1 is not None and r1 == pos and r2 is None and q2 is not None:
            # Detectar inserción
            ins_qpos = [q2]
            j = i + 1
            while (j + 1) < len(aligned_pairs):
                qq, rr = aligned_pairs[j + 1]
                if rr is None and qq == (ins_qpos[-1] + 1):
                    ins_qpos.append(qq)
                    j += 1
                else:
                    break

            ins_len = len(ins_qpos)
            ins_seq = ''.join(read.query_sequence[q].upper() for q in ins_qpos)
            insertion_str  = f"+{ins_len}{ins_seq}"

            # EN LUGAR DE USAR LA PRIMERA, TOMAMOS LA CALIDAD MÍNIMA DE TODAS LAS BASES
            ins_quals = [read.query_qualities[q] for q in ins_qpos]
            insertion_qual = min(ins_quals)

            break  # Tomamos solo la primera inserción que encontremos

    # (2) Buscar base o deleción en 'pos'
    base_str  = None
    base_qual = None

    for i, (qPos, rPos) in enumerate(aligned_pairs):
        if rPos == pos:
            if qPos is not None:
                # => hay base
                base_str  = read.query_sequence[qPos].upper()
                base_qual = read.query_qualities[qPos]
            else:
                # => deleción
                ref_positions = [pos]
                j = i
                while (j + 1) < len(aligned_pairs):
                    q2, r2 = aligned_pairs[j + 1]
                    if q2 is None and r2 == (ref_positions[-1] + 1):
                        ref_positions.append(r2)
                        j += 1
                    else:
                        break

                del_len = len(ref_positions)
                if ref_fasta and ref_name:
                    del_seq = ref_fasta.fetch(ref_name, ref_positions[0], ref_positions[-1] + 1).upper()
                    base_str = f"-{del_len}{del_seq}"
                else:
                    base_str = f"-{del_len}"
                base_qual = 0  # Calidad 0 (por convención en deleciones)

            break

    # (3) Combinar resultados según lo que se haya encontrado
    if not base_str and not insertion_str:
        return (None, None)

    if base_str and insertion_str:
        # Ej: "A+3ACT" o "-1A+3ACT"
        combined_str = base_str + insertion_str
        # Elegir calidad
        if base_qual == 0 and insertion_qual is not None:
            # Ej. si la base es una deleción => 0, 
            # entonces asignamos la calidad de la inserción
            final_qual = insertion_qual
        else:
            final_qual = base_qual
        return (combined_str, final_qual)

    if base_str and not insertion_str:
        return (base_str, base_qual)

    if insertion_str and not base_str:
        return (insertion_str, insertion_qual)

def parse_deletion(deletion_str):
    if '+' in deletion_str:
        return None, None

    if not deletion_str.startswith('-'):
        return None, None
    i = 1
    num_str = ''
    while i < len(deletion_str) and deletion_str[i].isdigit():
        num_str += deletion_str[i]
        i += 1
    if not num_str:
        return None, None
    num = int(num_str)
    seq = deletion_str[i:]
    return num, seq


def is_continuation(deletion_current, deletion_prev, delta):
    """
    Retorna True si 'deletion_current' es una continuación de 'deletion_prev'
    considerando que la diferencia de posiciones es delta.
    Es decir, si deletion_prev = "-5CTTCT" y delta = 1, se espera que deletion_current sea "-4TTCT".
    No hace falta gestionar las delecciones con inserciones después porque estas son de una sola base y, aúnque pueden estar repetidas con dos bases en la posición anterior la necesitamos para nombrar la inserción.
    """
    n_prev, seq_prev = parse_deletion(deletion_prev) #devuelve tamaño de la deleccion y secuencia
    n_cur, seq_cur = parse_deletion(deletion_current)
    if n_prev is None or n_cur is None:
        return False
    if n_prev != n_cur + delta: #El tamaño de la delección de ahora más delta tiene que ser igual al tamaño de la deleccion ya anotada
        return False
    if seq_prev[delta:] == seq_cur: #devuelve la secuencia a partir de delta, es decir lo que debería ser esa secuencia sin la diferencia (la parte no común)
        return True
    return False

def scan_region_for_polymorphisms(
    bam_file, 
    ref_name, 
    posiciones_a_analizar,
    min_reads=5, 
    ref_fasta=None, 
    min_base_qual=20,
    noindel=False
):
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = {}
    len_reads = []
    global_deletions = []
    read_minmax = {}

    for pos in tqdm(posiciones_a_analizar):
        base_reads = defaultdict(set)

        # 1) Recolectar lecturas
        for read in bam.fetch(ref_name, pos, pos + 1):
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            # Actualizar min/max solo si queremos guardarlo
            if SAVE_READS_MINMAX:
                r_start = read.reference_start
                r_end   = read.reference_end
                if read.query_name not in read_minmax:
                    read_minmax[read.query_name] = [r_start, r_end]
                else:
                    if r_start < read_minmax[read.query_name][0]:
                        read_minmax[read.query_name][0] = r_start
                    if r_end > read_minmax[read.query_name][1]:
                        read_minmax[read.query_name][1] = r_end

            base_str, base_qual = _get_base_or_indel_at_pos_indel(
                read, pos, ref_fasta, ref_name, noindel
            )
            if not base_str:
                continue

            # Filtrar por calidad si no es deleción
            if not base_str.startswith('-') and base_qual < min_base_qual:
                continue

            base_reads[base_str].add(read.query_name)
            len_reads.append(len(read.query_sequence))

        # 2) Eliminar lecturas conflictivas (las que estén en >1 variante)
        read_to_variants = defaultdict(list)
        for variant, rset in base_reads.items():
            for read_name in rset:
                read_to_variants[read_name].append(variant)

        for read_name, variants_list in read_to_variants.items():
            if len(variants_list) > 1:
                # Eliminar esa lectura de todas las variantes
                for v in variants_list:
                    base_reads[v].discard(read_name)

        # 3) Filtrar por min_reads
        passing_bases = {
            variant: rnames
            for variant, rnames in base_reads.items()
            if len(rnames) >= min_reads
        }

        # 4) Filtrar deleciones continuas
        keys_to_remove = []
        for variant in list(passing_bases.keys()):
            if variant.startswith('-'):
                remove_flag = False
                for (prev_pos, prev_del) in global_deletions:
                    delta = pos - prev_pos
                    if delta >= 1 and is_continuation(variant, prev_del, delta):
                        remove_flag = True
                        break
                if remove_flag:
                    keys_to_remove.append(variant)
                else:
                    global_deletions.append((pos, variant))

        for k in keys_to_remove:
            passing_bases.pop(k, None)

        # 5) Guardar solo si >= 2 variantes
        if len(passing_bases) >= 2:
            results[pos] = passing_bases

    bam.close()

    # Convertir read_minmax en lista (o vacío si no se guarda)
    if SAVE_READS_MINMAX:
        reads_in_fin = [(st, en) for (st, en) in read_minmax.values()]
    else:
        reads_in_fin = []

    # Devolver lo de siempre
    return results, max(len_reads), reads_in_fin


#FIN BLOQUE FUNCIÓN SCAN_REGION_FOR_POLYMORPHISMS
#INICIO BLOQUE FUNCIÓN COMPARAR_POSICIONES_CONSECUTIVAS:

def is_indel(variant):
    # Función para saber si un variant es un indel (inserción o deleción)
    return variant.startswith('-') or '+' in variant
    
    # Esta función determina cómo formatear una variante en el "detalle".
    #  - Si no es inserción, la devuelve igual.
    #  - Si es inserción, mira si hay un "X+Nseq" con X siendo 1 nucleótido.
    #    y X aparece como variante "sola" en variants_dict => quita X, dejando "+Nseq".

def format_variant(variants_dict, variant):
    # 1) Si no es inserción, lo devolvemos tal cual
    #    (deleción => se queda, base => se queda).
    if not '+' in variant:
        return variant
        
    # En "A+2CT"...
    #   - la "A" es la parte previa
    #   - comprobamos si "A" existe como variante "sola" en variants_dict
    #   - si sí, devolvemos "+2CT", si no, devolvemos "A+2CT"
    # LÓGICA:
    # Buscamos el primer '+' y separamos
    
    idx_plus = variant.find('+')   
    # base_previa es todo lo que esté antes del '+'
    base_previa = variant[:idx_plus]
    # Resto es la parte con '+' inclusive
    insercion = variant[idx_plus:]  # p.ej. "+2CT"

    # y existe como variante sola => devolvemos solo la parte de la inserción
    if base_previa in variants_dict:
        # hay una variante "A" sola => quitamos la base previa
        return insercion
    else:
        # mantenemos la base previa
        return variant

def check_comparacion_custom(posA, polyA, posB, polyB):
    """
    Realiza la prueba de 4 gametos entre posA y posB.
    
    Retorna (resultado, detalle), donde:
      - resultado = True/False/None
      - detalle:
         * Si True:  "vA1.vA2,vB1.vB2 / vA1'.vA2',vB1'.vB2' / ..." 
                    (para cada quadrupla que formó 4 gametos)
         * Si False: 
             - si TODAS las intersecciones son indel-based => "posAIndels,posBIndels"
             - de lo contrario => ""
         * Si None:  "" (sin intersecciones)
         
    EXTRA:
      - En inserciones tipo "A+2CT", si "A" existe sola en esa posición, 
        no se escribe ("+2CT" en su lugar). Si no existe sola, se deja "A+2CT".
    """
    # 1) Cobertura en cada posición (suma de todos los sets de variantes)
    coverageA = sum(len(rset) for rset in polyA.values())
    coverageB = sum(len(rset) for rset in polyB.values())
    coverage_min = min(coverageA, coverageB)    


    # Variables de control
    any_intersection = False       # Para saber si hubo al menos una intersección
    only_indels = True            # Se mantendrá True si en TODAS las intersecciones hay un indel en A o B
    posA_indels = set()           # Indels encontrados en posA
    posB_indels = set()           # Indels encontrados en posB

    # Almacena las quadruplas que SÍ forman los 4 gametos
    quadruplas_ok = []

    variantesA = list(polyA.keys())  # p.ej. ["A", "G", "A+2CT", "-3C", ...]
    variantesB = list(polyB.keys())

    n_gametos_true = 0
    for (vA1, vA2) in itertools.combinations(variantesA, 2):
        for (vB1, vB2) in itertools.combinations(variantesB, 2):
            min_freq=1
            min_inter=100000
            gametos = {
                (vA1, vB1): False,
                (vA1, vB2): False,
                (vA2, vB1): False,
                (vA2, vB2): False
            }

            # Revisar intersecciones
            for (gA, gB) in gametos:
                readsA = polyA[gA]
                readsB = polyB[gB]
                inter = readsA.intersection(readsB)
                
                inter_count=len(inter)
                if (inter_count != 0) and (inter_count < min_inter): #Siempre se va a poder hacer más restrictivo, pero no menos, porque si no los false también cambiarían
                    min_inter=inter_count

                # Filtro 1: umbral de lecturas mínimas
                enough_reads = (inter_count >= reads_gametos_minoritarios)
                
                # Filtro 2: umbral de frecuencia (respecto a coverage_min)
                # coverage_min puede ser 0; evitar división por cero
                if coverage_min > 0:
                    freq = inter_count / coverage_min
                else:
                    freq = 0.0

                if freq < min_freq:
                    min_freq = freq
                
                # Filtro por frecuencia menos restrictiva (ej. >= 0.1 => 10%)
                enough_freq = (freq >= THRESHOLD_GAMETO_FREQ)

                if enough_reads and enough_freq:
                    gametos[(gA, gB)] = True
                    any_intersection = True
                    # Revisar si hay un indel en gA o gB
                    if not (is_indel(gA) or is_indel(gB)):
                        only_indels = False
                    else:
                        # si gA es un indel => lo guardamos en posA_indels
                        if is_indel(gA):
                            posA_indels.add(gA)
                        # si gB es un indel => lo guardamos en posB_indels
                        if is_indel(gB):
                            posB_indels.add(gB)

            # Si se cumplieron los 4 gametos
            if all(gametos.values()):
                # Construir la salida "vA1.vA2,vB1.vB2"
                # pero formateando cada uno con format_variant
                fmtA1 = format_variant(polyA, vA1)
                fmtA2 = format_variant(polyA, vA2)
                fmtB1 = format_variant(polyB, vB1)
                fmtB2 = format_variant(polyB, vB2)

                combos_str = f"{fmtA1}.{fmtA2},{fmtB1}.{fmtB2},mireads:{min_inter},minfreq:{min_freq}"
                quadruplas_ok.append(combos_str)
            else:
                n_gametos_true_nuevo = sum(gametos.values())
                if n_gametos_true < n_gametos_true_nuevo:
                    n_gametos_true = n_gametos_true_nuevo
                    

    # Determinar el resultado
    if quadruplas_ok:
        # => tenemos al menos una quadrupla con 4 gametos
        detalle = "/".join(quadruplas_ok)
        return (True, detalle)
    else:
        # => no se formaron 4 gametos
        if any_intersection:
            # => hubo cruces
            if only_indels:
                # => todas las intersecciones son indel-based
                # formatear la lista de indels
                # y, en cada uno, aplicar format_variant
                if posA_indels:
                    posA_str = "/".join(sorted(format_variant(polyA, ind) for ind in posA_indels))
                else:
                    posA_str = ""

                if posB_indels:
                    posB_str = "/".join(sorted(format_variant(polyB, ind) for ind in posB_indels))
                else:
                    posB_str = ""

                return (False, f"{n_gametos_true},{posA_str},{posB_str}")
            else:
                # => hubo alguna intersección base-base => sin detalle
                return (False, f"{n_gametos_true}")
        else:
            # => ninguna intersección => None
            return (None, "")

def comparar_posiciones_consecutivas(dic_pol_pos):
    """
    Compara pares de posiciones consecutivas en dic_pol_pos usando 
    check_comparacion_custom, que devuelve (resultado, detalle).
    
    Parámetros:
      dic_pol_pos (dict): {pos -> {variante: set(read_names), ...}}
      
    Retorna:
      dict: {
          'True':  [(posA, posB, detalle), ...],
          'False': [(posA, posB, detalle), ...],
          'None':  [(posA, posB, detalle), ...]
      }
      donde:
        - detalle es la cadena devuelta por check_comparacion_custom (p.ej. "A.G,B.T/...").
    """

    comparaciones = {
        'True': [],
        'False': [],
        'None': []
    }

    # Obtener las posiciones ordenadas
    sorted_positions = sorted(dic_pol_pos.keys())

    # Iterar sobre pares de posiciones consecutivas
    for i in range(len(sorted_positions) - 1):
        posA = sorted_positions[i]
        posB = sorted_positions[i + 1]

        polyA = dic_pol_pos[posA]
        polyB = dic_pol_pos[posB]

        # Llamamos a check_comparacion_custom, que retorna (resultado, detalle)
        resultado, detalle = check_comparacion_custom(posA, polyA, posB, polyB)

        # Almacenar según el valor de 'resultado'
        if resultado is True:
            comparaciones['True'].append((posA, posB, detalle))
        elif resultado is False:
            comparaciones['False'].append((posA, posB, detalle))
        elif resultado is None:
            comparaciones['None'].append((posA, posB))

    return comparaciones

#FIN BLOQUE FUNCIÓN COMPARAR_POSICIONES_CONSECUTIVAS
#FUNCIÓN EXCLUSIVA IMPLEMENTACIÓN FAST:
def extract_polymorphic_fastq(bam_path, ref_name, posiciones_a_analizar, output_fastq, min_reads=5, ref_fasta_path=None,noindel=False):
    """
    Extrae los reads que mapean en posiciones polimórficas dentro de la región [start, end)
    del archivo BAM y genera un archivo FASTQ con dichos reads.
    
    Parámetros:
      bam_path: Ruta al archivo BAM.
      ref_name: Nombre de la referencia (ej: "NC_045512.2").
      start: Posición inicial (0-based).
      end: Posición final (0-based, no inclusiva).
      output_fastq: Ruta del archivo FASTQ resultante.
      min_reads: Mínimo número de reads necesarios para considerar una base/indel.
      ref_fasta_path: (Opcional) Ruta al archivo FASTA de la referencia, para extraer la secuencia real
                      en caso de deleciones.
    """
    # Abrir el FASTA de referencia si se proporciona la ruta.
    ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None

    # Obtener las posiciones polimórficas utilizando la función ya definida.
    polymorphic_positions,max_len_reads, reads_in_fin = scan_region_for_polymorphisms(bam_path, ref_name, posiciones_a_analizar, min_reads, ref_fasta, min_base_qual, noindel)
    
    # Cerrar el FASTA si se abrió.
    if ref_fasta:
        ref_fasta.close()
    
    # Reunir los nombres de todos los reads que aportan una base/indel en posiciones polimórficas.
    reads_of_interest = set()
    for pos, bases in polymorphic_positions.items():
        for rnames in bases.values():
            reads_of_interest.update(rnames)
    
    if not reads_of_interest:
        print("[INFO] No reads of interest found in the specified region.")
        return

    # Guardar la lista de reads en un archivo temporal (uno por línea)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp_list = tmp.name
        for rname in reads_of_interest:
            tmp.write(rname + "\n")
    
    # Construir el comando para extraer los reads de interés y generar un FASTQ.
    # El comando hace lo siguiente:
    #  1) samtools view -h -F 4: extrae las lecturas mapeadas (con cabecera)
    #  2) awk: filtra aquellas líneas cuyo primer campo (nombre del read) esté en el archivo temporal.
    #  3) samtools fastq -: convierte el SAM filtrado a FASTQ, escribiendo en el archivo de salida.
    cmd = (
        f"samtools view -h -F 4 {bam_path} | "
        f"awk 'BEGIN{{while((getline line < \"{tmp_list}\")>0) r[line]=1}} "
        f"/^@/{{print;next}} ($1 in r)' | "
        f"samtools fastq - > {output_fastq}"
    )
    
    print(f"Ejecutando: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        print(f"[INFO] FASTQ file with reads of interest generated: {output_fastq}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] No se pudo generar el archivo FASTQ de reads de interés: {e}")
        sys.exit(1)
    finally:
        # Eliminar el archivo temporal de la lista de reads.
        os.remove(tmp_list)

#FIN FUNCIÓN EXLUSIVA IMPLEMENTACIÓN FASTA

#INICIO FUNCIONES USO DE PILEUP

def quick_check_col5(bases_str,MIN_READS):

    base_counts = {}
    i = 0
    length = len(bases_str)

    while i < length:
        c = bases_str[i]

        if c in ('.', ','):
            b = "."
            base_counts[b] = base_counts.get(b, 0) + 1
            i += 1

        elif c == '^':
            # saltar '^' + mapQuality
            i += 2

        elif c == '$':
            i += 1

        elif c in ['+', '-']:
            sign = c
            i += 1
            num_str = ''
            # leer la longitud
            while i < length and bases_str[i].isdigit():
                num_str += bases_str[i]
                i += 1
            if num_str.isdigit():
                skip_len = int(num_str)
                if i+skip_len <= length:
                    indel_seq = bases_str[i:i+skip_len]
                else:
                    indel_seq = ''
                i += skip_len
                b = sign + indel_seq  # e.g. +AC
                base_counts[b] = base_counts.get(b, 0) + 1
            else:
                # no hay número, ignorar
                continue
        
        elif c.upper() in ('A', 'C', 'G', 'T'):
            b = c.upper()
            base_counts[b] = base_counts.get(b, 0) + 1
            i += 1

        else:
            # Caracter que no contemplamos -> lo ignoramos
            i += 1

    bases_cumplen = [b for b,cnt in base_counts.items() if cnt >= MIN_READS]
    return (len(bases_cumplen) >= 2)


def parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS): #HAY QUE VER SI ESTO ES 0-BASED O NO? POR EJEMPLO CON ESTE MÉTODO ESTAMOS SELECCIONANDO LA POSICIÓN 701 COMO LA PRIMERA, CON EL OTRO SELECCIONÁBAMOS LA POSICIÓN 700. -> AHORA LO AJUSTÉ A 1-BASED -> CONTÉ LAS POSICIONES SOBRE LA REFERENCIA Y LA 39 (1-BASED) ERA TAMBIÉN LA 39 EN EL PILEUP
    with open(pileup_file, 'r') as f:
        posiciones_polimorficas = []

        for line in f:
            fields = line.strip().split('\t')
            ####print_verbose(fields)
            if len(fields) < 6: #HAY QUE CAMBIAR ESTO SI DESPUÉS NO SACO LOS NOMBRES DE LAS LECTURAS
                continue
            pos = int(fields[1])
            #print(f"Posición: {pos}")
            bases = str(fields[4])
            if len(bases) < (MIN_READS*2):
                continue
            polimorfismo_si_o_no = quick_check_col5(bases,MIN_READS)
            if polimorfismo_si_o_no:
                ###print_verbose(f"Posición analizada: {pos}, bases analizadas: {bases}")
                posiciones_polimorficas.append(pos-1) #Ajusto a 0-based para leer con pysam

    return posiciones_polimorficas

#FIN FUNCIONES USO DE PILEUP
#INICIO FUNCIONES ALGORITMO DE BÚSQUEDA
def construir_distribuciones(dic_pol_pos):
    """
    Convierte dic_pol_pos (pos -> {base -> set(reads), ...})
    en un dict pos -> { base -> porcentaje }.

    Verifica que ningún read aparezca en varias variantes de la misma pos,
    calcula (#reads_variant / total_reads). 
    Imprime el resultado para depuración.
    """
    distributions_por_pos = {}

    print_substep("Construyendo distribuciones de polimorfismos ...")
    for pos, variant_dict in tqdm(dic_pol_pos.items()):
        all_reads = set()
        for variant, rset in variant_dict.items():
            for r in rset:
                if r in all_reads: #HAY QUE DESCOMENTAR ESTO CUANDO TERMINE DE ADAPTAR EL CÓDIGO Y VER POR QUÉ HAY READS REPETIDOS
                    raise ValueError(
                        f"Read '{r}' aparece en múltiples variantes en la posición {pos}."
                    )
                all_reads.add(r)

        total_reads = sum(len(rset) for rset in variant_dict.values())
        if total_reads == 0:
            continue

        dist = {}
        for variant, rset in variant_dict.items():
            dist[variant] = len(rset) / total_reads

        distributions_por_pos[pos] = dist

    ###print_verbose(f"\n[DEBUG] Distributions por posición: {distributions_por_pos}")
    return distributions_por_pos


def generar_todos_grupos_por_ventana(positions, max_read_length, min_size=2):
    """
    Genera TODAS las 'ventanas' (subconjuntos) donde 
    (pos_final - pos_inicial) <= max_read_length.
    
    Luego elimina aquellas ventanas que son un 
    subconjunto estricto de otra ventana mayor.
    """

    # == PARTE 1: Generar ventanas como antes ==
    sorted_positions = sorted(positions)
    todos_grupos = []
    n = len(sorted_positions)

    for i in range(n):
        grupo_actual = set()
        p_inicial = sorted_positions[i]
        grupo_actual.add(p_inicial)

        j = i + 1
        while j < n and (sorted_positions[j] - p_inicial) <= max_read_length:
            grupo_actual.add(sorted_positions[j])
            j += 1

        if len(grupo_actual) >= min_size:
            todos_grupos.append(grupo_actual)

    # == PARTE 2: Eliminar subconjuntos ==
    # 1) Ordenar descendentemente por tamaño
    todos_grupos_ordenados = sorted(todos_grupos, key=len, reverse=True)

    # 2) Recorremos en orden grande->pequeño
    ventanas_finales = []
    for v in todos_grupos_ordenados:
        # Si 'v' es subconjunto estricto de alguna en 'ventanas_finales', REVISAR Y ENTENDER BIEN ESTE CÓDIGO -> FUNCIONA, PERO HAY QUE ENTENDERLO
        # la descartamos
        # Sino, la añadimos
        if not any((v < x) for x in ventanas_finales):
            # 'v < x' significa v.issubset(x) y v != x
            ventanas_finales.append(v)

    # 3) Opcionalmente, podríamos reordenar 'ventanas_finales' 
    #    según tu preferencia (por tamaño asc/desc o 
    #    volver al orden en que se generaron).
    #    Aquí las dejamos en orden descendente por tamaño:
    return ventanas_finales

##############################################################################
# 1.A - Ignora la base: compara fracciones
##############################################################################

def agrupar_por_distribucion_individual(distributions_por_pos, threshold=0.05):
    """
    Agrupa posiciones basándose en las fracciones individuales de su distribución.
    
    Para cada posición (ej. 100, 101, etc.), para cada base (ej. 'A', 'G', etc.) y su porcentaje,
    se intenta agregar la posición a un grupo existente cuyo promedio (para esa sub-distribución)
    esté dentro del umbral (threshold) de la fracción de la posición.
    
    Cada grupo se representa como un diccionario con:
       - 'promedio': valor numérico (float) que es la media de las fracciones en el grupo.
       - 'posiciones': un set con las posiciones que han contribuido.
       
    Cuando se añade una nueva posición, se actualiza el promedio:
       nuevo_promedio = (promedio * N + porcentaje_actual) / (N+1),
       donde N es el número de posiciones ya en el grupo.
    
    Si la posición no encaja en ningún grupo para esa fracción, se crea un grupo nuevo.
    
    Al finalizar, se eliminan los grupos que sean estrictamente iguales o subconjuntos de otros.
    
    Retorna una lista de grupos.
    """
    grupos = []  # Cada grupo es un diccionario: {'promedio': float, 'posiciones': set(), 'base': base}
    
    ###print_verbose(f"\n------------------------NUEVA VENTANA: Distribuciones por pos:\n {distributions_por_pos}\n")

    # Procesamos las posiciones en orden creciente (aunque el orden puede no importar)
    for pos in sorted(distributions_por_pos.keys()):
        ###print_verbose(f"\nPOSICIÓN: {pos}")
        dist = list(set(distributions_por_pos[pos].values())) #Ejemplo [0.6, 0.4]
        ###print_verbose(f"DIST: {dist}")
        # Para cada base y su porcentaje en la posición:
        for porcentaje in dist:
            ####print_verbose(f"Porcentaje: {porcentaje}, grupos len: {len(grupos)}")
            if len(grupos) > 0 and any(abs(porcentaje - grupo['promedio']) <= threshold for grupo in grupos):
                for grupo in grupos:
                    if abs(porcentaje - grupo['promedio']) <= threshold:
                        grupo['posiciones'].add(pos)
                        ###print_verbose(f"Anterior media: {grupo['promedio']}, ajustando con el porcentaje: {porcentaje}")
                        N = len(grupo['posiciones'])
                        grupo['promedio'] = (grupo['promedio'] * (N - 1) + porcentaje) / N
                        ###print_verbose(f"Nueva media: {grupo['promedio']}")
            else:
                ###print_verbose("No existe el grupo: Creando ...")
                nuevo_grupo = {'promedio': porcentaje, 'posiciones': {pos}}
                ###print_verbose(f"Nuevo grupo: {nuevo_grupo}")
                grupos.append(nuevo_grupo)               
    
    ###print_verbose(f"Grupos no ordenados: {grupos}")

    grupos = sorted(grupos, key=lambda grupo: len(grupo['posiciones']), reverse=True)

    ###print_verbose(f"Grupos: {grupos}")

    grupos_finales = []
    for grupo in grupos:
        pos_set = grupo['posiciones']
        # Si existe en grupos_finales otro grupo cuyo conjunto de posiciones es mayor y contiene pos_set,
        # entonces descartamos este grupo.
        if not any(pos_set == otro['posiciones'] for otro in grupos_finales): #poner <= para eliminar también subgrupos
            grupos_finales.append(grupo)
    ###print_verbose(f"Grupos finales: {grupos_finales}") 

    return grupos_finales

##############################################################################
# LÓGICA REFINADA DENTRO DE GRUPO
##############################################################################

def existe_subtrue_mas_pequeno(p_first, p_last, comp_global):
    """
    Retorna True si dentro de comp_global existe un par (pA, pB)
    tal que:
      1) comp_global[(pA, pB)]["res"] == True
      2) pA >= p_first y pB <= p_last
      3) (pA, pB) != (p_first, p_last)
      4) anotar == True (opcional, si quieres exigir que además esté anotar=True)
    """
    for (pA, pB), info in comp_global.items():
        if info.get("res") is True and info.get("anotar", True) is True:
            if pA >= p_first and pB <= p_last and (pA, pB) != (p_first, p_last):
                return True
    return False

def desanotar_rangos_mas_grandes(p_first, p_last, comp_global):
    """
    Marca con anotar=False todos los rangos (pA, pB) en 'comp_global' que:
      - tengan res=True, anotar=True,
      - sean distintos de (p_first, p_last),
      - contengan a (p_first, p_last), es decir pA <= p_first y pB >= p_last.

    Esto sirve para 'descartar' rangos más grandes cuando preferimos quedarnos
    con el rango nuevo (más pequeño o más 'concreto').
    """
    for (pA, pB), info in comp_global.items():
        # Debe ser True y estar marcado para anotar
        if info.get("res") is True and info.get("anotar", True) is True:
            # Debe ser un rango distinto de (p_first, p_last)
            if (pA, pB) != (p_first, p_last):
                # Y contenerlo (inclusivo)
                if pA <= p_first and pB >= p_last:
                    ###print_verbose(f"   -> Marcando anotar=False para el rango grande {(pA, pB)}")
                    comp_global[(pA, pB)]["anotar"] = False

def refine_compare(pos_list, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global): #VER SI NO SE PUEDE PONER COMO UNA FUNCIÓN INDIVIDUAL
    
    resultados = []

    if len(pos_list) < 2:
        ###print_verbose("ONLY ONE POSITION")
        return

    p_first = pos_list[0]
    p_last  = pos_list[-1]
    ###print_verbose(f"Comparando: p_first = {p_first}, p_last = {p_last}")
    if p_first == p_last:
        return

    key = (p_first, p_last)


    if key in comparaciones_hechas_global:
        info = comparaciones_hechas_global[key]
        res = info["res"]
        detail = info["detail"]
        ###print_verbose(f"Ya se comparó {key}: {res}, {detail} (intersected={info['intersected']})")
        # Si ya fue marcado en grupos (intersected=True), no refinamos más
        if info["intersected"]:
            return
        else:
            # Si aparece por primera vez en grupos, actualizamos el flag a True
            comparaciones_hechas_global[key]["intersected"] = True
            ###print_verbose(f"Actualizando flag a True para {key} (ya se había hecho en intersección)")
    else:
        res, detail = check_comparacion_custom(p_first, dic_pol_pos[p_first], p_last, dic_pol_pos[p_last])
        # Aquí, en el contexto de grupos refinados, el flag se establece en True
        comparaciones_hechas_global[key] = {"res": res, "detail": detail, "intersected": True, "anotar": True}
        ###print_verbose(f"Result of {key} (group refinement): {res}, {detail}, intersected: True")

    ####print_verbose(f"Key={key} – Res={res}")


    if res is False:
        resultado_a_anotar=(p_first, p_last, False, detail)    
        if resultado_a_anotar not in resultados:
            resultados.append(resultado_a_anotar)
        return

    elif res is True:
        before_count = len(resultados)
        n = len(pos_list)
        if n == 2:
            resultado_a_anotar=(p_first, p_last, True, detail)
            if resultado_a_anotar not in resultados:
                resultados.append(resultado_a_anotar)
            desanotar_rangos_mas_grandes(p_first, p_last, comparaciones_hechas_global)
            return

        # Subdividir según paridad
        if n % 2 == 1: #impar
            mid_index = n // 2
            refine_compare(pos_list[:mid_index+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) #El +1 está bien, es porque se accede como a un range entoncss no se considera el último elementeo
            refine_compare(pos_list[mid_index:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) #Aquí el mid_index se refiere al primer elemento.
        else:
            half = n // 2
            refine_compare(pos_list[:half+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)
            refine_compare(pos_list[half-1:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

        # ---------------------------------------------
        # LÓGICA NUEVA DE ANOTACIÓN (SIN before_count):
        # ---------------------------------------------
        # Guardamos este resultado (res=True) en comparaciones_hechas_global,
        # pero no decidimos anotar hasta ver si ya existe un sub-True más pequeño:
        comparaciones_hechas_global[key] = {
            "res": True,
            "detail": detail,
            "intersected": True,
            # anotar se definirá según exista subTrue
            "anotar": True  # valor provisional
        }

        # Si ya existe un sub-True más pequeño, marcamos anotar=False
        if existe_subtrue_mas_pequeno(p_first, p_last, comparaciones_hechas_global):
            ###print_verbose(f"Ya existe un True más concreto dentro de {key}; se anota con anotar=False.")
            comparaciones_hechas_global[key]["anotar"] = False
        else:
            ###print_verbose(f"No existe sub-True más pequeño para {key}; se anota con anotar=True.")
            pass

        desanotar_rangos_mas_grandes(p_first, p_last, comparaciones_hechas_global)
        return

    else:  # res is None
        if len(pos_list) == 2:
            resultado_a_anotar=(p_first, p_last, None, detail)
            if resultado_a_anotar not in resultados:    
                resultados.append(resultado_a_anotar)
            return

        # Reducir p_last uno a uno (de derecha a izquierda)
        for idx in range(len(pos_list) - 2, 0, -1): #Esto para cada iteración devuelve el idx que es cada vez una posición menos desde la última posición y saltándose la primera, por ejemplo si fuese 3 sería 3,2,1 (len sería 4).
            list_candidate=pos_list[:idx+1]
            p_first = list_candidate[0]
            p_last  = list_candidate[-1]
            resultado_a_anotar=(p_first, p_last, None, detail)
            if resultado_a_anotar not in resultados:    
                resultados.append(resultado_a_anotar)
            refine_compare(list_candidate, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)
        return

##############################################################################
# compare_fn y comparar_grupos_interseccion
##############################################################################

def comparar_grupos_interseccion(grupos, check_comparacion_custom, dic_pol_pos, comparaciones_hechas_global):
    ###print_verbose("\nENTRAMOS EN LA FUNCIÓN comparar_grupos_interseccion\n")
    
    resultados = {'True': [], 'False': [], 'None': []}
    etiquetados = {}
    ###print_verbose(f"grupos: {grupos}")
    
    # Construir un diccionario que asocie cada posición a la(s) etiqueta(s) (los promedios de los grupos a los que pertenece)
    for grupo in grupos:
        for p in grupo["posiciones"]:
            if p in etiquetados:
                etiquetados[p] = f"{etiquetados[p]},{grupo['promedio']}" #Esto lo que hace es concatenar los valores promedios encontrados en esa posición, coge los promedios que ya había y le añade el nuevo separado con coma
            else:
                etiquetados[p] = str(grupo["promedio"])
    
    ###print_verbose(f"etiquetados: {etiquetados}")
    posiciones_ordenadas = sorted(etiquetados.keys())
    ###print_verbose(f"posiciones_ordenadas: {posiciones_ordenadas}")
    
    # Recorrer todas las parejas (pA, pB) con pA < pB para evitar repeticiones
    for i in range(len(posiciones_ordenadas)):
        for j in range(i + 1, len(posiciones_ordenadas)):
            pA = posiciones_ordenadas[i]
            pB = posiciones_ordenadas[j]
            # Convertimos las etiquetas a listas de floats para compararlas
            gA = sorted(map(float, etiquetados[pA].split(',')))
            gB = sorted(map(float, etiquetados[pB].split(',')))
            
            ###print_verbose(f"Comparing: pA: {pA}, gA: {gA} | pB: {pB}, gB: {gB}")
            
            # Solo comparamos si los grupos son diferentes
            if gA == gB:
                ###print_verbose(f"Los grupos para {pA} y {pB} son iguales; se omite la comparación.")
                continue
            
            # Usamos la tupla ordenada (mínimo, máximo) para evitar duplicados
            pair = (min(pA, pB), max(pA, pB))
            if pair not in comparaciones_hechas_global:
                res, detail = check_comparacion_custom(pA, dic_pol_pos[pA], pB, dic_pol_pos[pB])
                ####print_verbose(f"Key={pair} – Res={res}")
                comparaciones_hechas_global[pair] = {"res": res, "detail": detail, "intersected": False, "anotar": True}
                if existe_subtrue_mas_pequeno(pair[0], pair[1], comparaciones_hechas_global):
                    ###print_verbose(f"Ya existe un True más concreto dentro de {pair}; se anota con anotar=False.")
                    comparaciones_hechas_global[pair]["anotar"] = False
                else:
                    ###print_verbose(f"No existe sub-True más pequeño para {pair}; se anota con anotar=True.")
                    pass

                desanotar_rangos_mas_grandes(pair[0], pair[1], comparaciones_hechas_global)
                ###print_verbose(f"Comparaciones hechas global: {comparaciones_hechas_global}")
                polyA = dic_pol_pos[pA]
                polyB = dic_pol_pos[pB]
                ###print_verbose(f"Comparing position {pA} (reads: {polyA}) con {pB} (reads: {polyB})")
                res, detail = check_comparacion_custom(pA, polyA, pB, polyB)
                if res is True:
                    resultados['True'].append((pA, pB, detail))
                elif res is False:
                    resultados['False'].append((pA, pB, detail))
                else:
                    resultados['None'].append((pA, pB, detail))
            else:
                ###print_verbose(f"El par {pair} ya se comparó; se omite.")
                pass
    
    return resultados

##############################################################################
# 2. FUNCIÓN PRINCIPAL
##############################################################################

def pipeline_agrupacion_por_distribuciones(dic_pol_pos, max_len=5, threshold=0.05, min_size=2):
    """
    1) Construye distribuciones
    2) Genera ventanas
    3) Usa un set global 'comparaciones_hechas_global' para no repetir 
       comparaciones en TODAS las ventanas
    4) Para cada ventana:
       - Filtra las pos => sub_dist
       - Agrupa ignorando la base
       - Compara ENTRE grupos (adyacentes)
       - Compara DENTRO de cada grupo con la función refinada
    5) Devuelve un dict con la info
    """

    dist_por_pos = construir_distribuciones(dic_pol_pos)

    ###print_verbose(f"dist_por_pos: {dist_por_pos}")

    ventanas = generar_todos_grupos_por_ventana(dist_por_pos.keys(), max_len, min_size)

    ###print_verbose(f"ventanas: {ventanas}")

    results_all_windows = {}
    comparaciones_hechas_global = {} 

    print_substep("Analizando posiciones ...")
    for i, window in enumerate(ventanas):
        sub_dist = {pos: dist_por_pos[pos] for pos in window}
        grupos_subdist = agrupar_por_distribucion_individual(sub_dist, threshold)

        ###print_verbose(f"grupos_subdist: {grupos_subdist}")

        comp_entre = comparar_grupos_interseccion(grupos_subdist, check_comparacion_custom, dic_pol_pos, comparaciones_hechas_global)

        ###print_verbose(f"comp_entre: {comp_entre}")

        comp_en_grupos = {}
        for j, g in enumerate(grupos_subdist):
            ###print_verbose(f"G:{g}\n")
            # Usamos la versión refinada
            pos_sorted = sorted(g["posiciones"])
            comp_en_grupos[j] = refine_compare(pos_sorted, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

            ###print_verbose(f"comp_en_grupos[j]: {comp_en_grupos[j]}")

        results_all_windows[i] = {"window_positions": sorted(window),"grupos": grupos_subdist,"comparaciones_entre_grupos": comp_entre,"comparaciones_en_grupos": comp_en_grupos}

        ###print_verbose(f"results_all_windows[i]: {results_all_windows[i]}")

    return results_all_windows, comparaciones_hechas_global, dist_por_pos
# FIN FUNCIONES ALGORITMO DE BÚSQUEDA

def run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False):
    """
    Realiza el trimming de lecturas. Si 'filter_unpaired=True' y 'layout' es PAIRED,
    primero ejecuta repair.sh para descartar/ajustar pares huérfanos,
    luego el trimming normal con cutadapt.
    """

    print_step("Step 4: Trimming reads...")

    trimmed_reads = work_dir / "reads_trimmed.fastq"
    read1_trimmed = work_dir / "reads_trimmed_1.fastq"
    read2_trimmed = work_dir / "reads_trimmed_2.fastq"
    read1_trimmed_temp = work_dir / "reads_trimmed_1_temp.fastq"
    read2_trimmed_temp = work_dir / "reads_trimmed_2_temp.fastq"
    trimmed_reads_temp = work_dir / "reads_trimmed.fastq.temp"

    # Si ya existen los archivos resultantes, saltamos el proceso
    if (read1_trimmed.exists() and read2_trimmed.exists()) or trimmed_reads.exists():
        print_substep("[INFO] Trimmed read files already exist. Skipping trimming.")
        return True

    try:
        # ======================================
        # PAIRED-END
        # ======================================
        if layout == "PAIRED":
            # ------------------------------------------------------------------------
            # (1) Si filter_unpaired=True -> ejecutar repair.sh ANTES de trimming
            # ------------------------------------------------------------------------
            MAX_PIPELINE_ATTEMPTS=1
            if filter_unpaired:
                MAX_PIPELINE_ATTEMPTS=3
                # Rutas de salida del repair
                repaired_1 = work_dir / "reads_repaired_1.fastq"
                repaired_2 = work_dir / "reads_repaired_2.fastq"
                singletons = work_dir / "singletons.fastq"

                repair_cmd = (f"/mnt/netapp2/Store_uni/home/uvi/cm/lgv/DeNovo/software/bbmap/repair.sh in={fastq_1} in2={fastq_2} out={repaired_1} out2={repaired_2} outs={singletons}")
                print_substep(f"Ejecutando repair.sh (para descartar huérfanas antes de cutadapt): {repair_cmd}")
                subprocess.run(repair_cmd, shell=True, check=True)

                # A partir de aquí usaremos los repaired_1 y repaired_2
                actual_fastq_1 = repaired_1
                actual_fastq_2 = repaired_2

            else:
                # Si no hay que filtrar/descartar lecturas huérfanas, usamos los fastq originales
                actual_fastq_1 = fastq_1
                actual_fastq_2 = fastq_2

            # ------------------------------------------------------------------------
            # (2) Ejecutar cutadapt en Paired-end
            # ------------------------------------------------------------------------
            cutadapt_cmd = (f"cutadapt -a file:/mnt/netapp2/Store_uni/home/uvi/cm/lgv/DeNovo/inputs/adapters.fa -A file:/mnt/netapp2/Store_uni/home/uvi/cm/lgv/DeNovo/inputs/adapters.fa -o {read1_trimmed_temp} -p {read2_trimmed_temp} {actual_fastq_1} {actual_fastq_2}")

            print_substep(f"Ejecutando cutadapt (Paired-end): {cutadapt_cmd}")
            if not retry_command(MAX_PIPELINE_ATTEMPTS, cutadapt_cmd):
                raise RuntimeError("Error durante el trimming Paired-end.")

            # Renombrar los .temp a definitivos
            mv_trim_cmd_1 = f"mv {read1_trimmed_temp} {read1_trimmed}"
            mv_trim_cmd_2 = f"mv {read2_trimmed_temp} {read2_trimmed}"
            subprocess.run(mv_trim_cmd_1, shell=True, check=True, executable='/bin/bash')
            subprocess.run(mv_trim_cmd_2, shell=True, check=True, executable='/bin/bash')

        # ======================================
        # SINGLE-END
        # ======================================
        else:
            MAX_PIPELINE_ATTEMPTS=3
            # SINGLE-END => no hay reparación previa
            cutadapt_cmd = (f"cutadapt -a file:/mnt/netapp2/Store_uni/home/uvi/cm/lgv/DeNovo/inputs/adapters.fa -o {trimmed_reads_temp} {fastq_1}")
            print_substep(f"Ejecutando cutadapt (Single-end): {cutadapt_cmd}")
            if not retry_command(MAX_PIPELINE_ATTEMPTS, cutadapt_cmd):
                raise RuntimeError("Error durante el trimming Single-end.")

            mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
            subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')

        print_substep("[INFO] Trimming completado correctamente.")
        return True

    except Exception as e:
        print_error(f"[ERROR] Trimming failed: {e}")
        return False

def count_true_superranges_in(pA, pB, comp_global):
    """
    Devuelve cuántos rangos (xA, xB) con res=True (independientemente de anotar)
    contienen o coinciden con (pA, pB). Es decir: xA <= pA y xB >= pB.
    """
    count = 0
    for (xA, xB), info in comp_global.items():
        if info.get("res") is True:  # solo importa que sea True
            # si xA <= pA y xB >= pB => (xA, xB) contiene (pA, pB)
            if xA <= pA and xB >= pB:
                count += 1
    return count

def filtrar_rangos_mas_grandes(comparaciones_hechas_global):
    """
    Recorre todos los pares (pA, pB) con res == True y desanota (anotar=False)
    aquellos rangos "grandes" que contengan (pA' >= pA && pB' <= pB) algún subrango
    más pequeño que también sea True.
    """
    # 1) Recolectar todos los pares donde res == True (independientemente de anotar)
    true_pairs = []
    for (pA, pB), info in comparaciones_hechas_global.items():
        if info["res"] is True:
            true_pairs.append((pA, pB))

    # 2) Ordenarlos del más grande al más pequeño
    #    (p.ej. un rango de 100,200 tiene tamaño 100, uno de 100,105 tiene 5).
    true_pairs.sort(key=lambda x: x[1] - x[0], reverse=True)

    # 3) Para cada rango grande, vemos si contiene a otro más pequeño que sea True
    for (pA, pB) in true_pairs:
        # Si ya se desanotó en iteraciones anteriores, no hacemos nada
        if not comparaciones_hechas_global[(pA, pB)]["anotar"]:
            continue

        # Buscamos un subrango (pA2, pB2) con pA2>=pA, pB2<=pB y res==True
        for (pA2, pB2) in true_pairs:
            # Evitar compararse consigo mismo
            if (pA2, pB2) == (pA, pB):
                continue

            # Ver si es subrango
            if pA2 >= pA and pB2 <= pB:
                # Y es True
                if comparaciones_hechas_global[(pA2, pB2)]["res"] is True:
                    # => Marcamos el rango grande con anotar=False y salimos
                    comparaciones_hechas_global[(pA, pB)]["anotar"] = False
                    break

# =========================
# Preparación mínima de referencia (solo cuando --use_shiver_assembly)
# =========================
def prepare_reference_minimal(ref_fa: str):
    """
    Prepara la referencia SOLO con lo necesario, sin chequeos:
      1) samtools faidx ref.fasta
      2) picard CreateSequenceDictionary R=ref.fasta O=ref.dict
    Si algo falla, que falle (propaga excepción).
    """
    # 1) FAI
    subprocess.run(["samtools", "faidx", ref_fa], check=True)
    # 2) .dict (Picard)
    dict_out = str(Path(ref_fa).with_suffix(".dict"))
    subprocess.run(["picard", "CreateSequenceDictionary", f"R={ref_fa}", f"O={dict_out}"], check=True)

# =========================
# Helpers para SHIVER
# =========================
def _open_text_auto(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")

def _add_pair_suffix(in_path, out_path=None, pair=None):
    """
    Crea una copia FASTQ(.gz) añadiendo /1 o /2 al PRIMER token del header (@ID).
    Limpia /1 o /2 existente antes de añadir el correcto.
    """
    in_path = str(in_path)
    if pair not in ("1", "2"):
        raise ValueError("pair debe ser '1' o '2'")
    if out_path is None:
        if in_path.endswith(".fastq.gz"):
            out_path = in_path[:-9] + "_with12.fastq.gz"
        elif in_path.endswith(".fastq"):
            out_path = in_path[:-6] + "_with12.fastq"
        else:
            raise ValueError("La entrada debe ser .fastq o .fastq.gz")
    clean_suffix = re.compile(r'(\/[12])?$')
    with _open_text_auto(in_path, "rt") as fin, _open_text_auto(out_path, "wt") as fout:
        for i, line in enumerate(fin, 1):
            if i % 4 == 1 and line.startswith("@"):
                parts = line.rstrip("\n").split(maxsplit=1)
                first = parts[0]
                rest = (" " + parts[1]) if len(parts) > 1 else ""
                core = first[1:]
                core = clean_suffix.sub("", core)
                line = "@" + core + f"/{pair}" + rest + "\n"
            fout.write(line)
    return out_path

def _ensure_slash_suffix(path, pair):
    """Devuelve el mismo path si ya termina en /1 o /2; si no, crea *_with12.fastq(.gz) y devuelve esa ruta."""
    def _peek_first_id(p):
        with _open_text_auto(p, "rt") as fh:
            for i, line in enumerate(fh, 1):
                if i == 1 and line.startswith("@"):
                    return line.strip().split()[0]
                if i > 1:
                    break
        return ""
    first = _peek_first_id(str(path))
    if first.endswith("/1") or first.endswith("/2"):
        return str(path)
    return _add_pair_suffix(str(path), pair=pair)

def run_shiver_assembly(
    srr,
    fastq_r1_trimmed,
    fastq_r2_trimmed,
    work_dir,
    base_dir,
    shiver_dir,
    rnaviralspades,
    ref_alignment_fa,
    adapters_fa,
    primers_fa,
    shiver_config,
    spades_threads,
    spades_mem_gb,
    shiver_threads,
    spades_conda_env=None,
    spades_env_prefix=None,
):
    """
    Ejecuta el flujo:
      1) rnaviralspades.py -1 R1 -2 R2 -o HIV_rnaviral -t X -m Y
      2) shiver_init.sh MyInitDir config.sh MyRefAligment.fasta  MyAdapters.fasta MyPrimers.fasta
      3) shiver_align_contigs.sh MyInitDir config.sh contigs.fasta SRR
      4) Asegurar sufijo /1 y /2 en FASTQ trimeados (crea *_with12.fastq si hace falta)
      5) shiver_map_reads.sh MyInitDir config.sh contigs.fasta SRR SRR.blast SRR_cut_wRefs.fasta R1 R2
    Devuelve la ruta al BAM resultante: {work_dir}/{srr}.bam
    """
    # Normaliza rutas respecto a base_dir (no a CWD)
    from pathlib import Path
    root = Path(base_dir).expanduser().resolve()
    def _abs(p):
        p = Path(p).expanduser()
        return (p if p.is_absolute() else (root / p)).resolve()

    shiver_dir_abs   = _abs(shiver_dir)
    ref_alignment_fa = _abs(ref_alignment_fa)
    adapters_fa      = _abs(adapters_fa)
    primers_fa       = _abs(primers_fa)
    shiver_bin = shiver_dir_abs / "bin"
    shiver_init = shiver_bin / "shiver_init.sh"
    shiver_align = shiver_bin / "shiver_align_contigs.sh"
    shiver_map = shiver_bin / "shiver_map_reads.sh"
    config_sh = shiver_bin / "config.sh"
    work_dir  = Path(work_dir).expanduser().resolve()
    out_spades = work_dir / "HIV_rnaviral"
    if not (out_spades / "contigs.fasta").exists():
        # Ejecutar rnaviralSPAdes; si hay entorno conda indicado, usar conda run
        exe = str(rnaviralspades)  # usa la ruta completa; no uses .name

        cmd = [
            exe, "-1", str(fastq_r1_trimmed), "-2", str(fastq_r2_trimmed),
            "-o", str(out_spades), "-t", str(spades_threads), "-m", str(spades_mem_gb)
        ]
        print_substep(f"[SHIVER] Ejecutando rnaviralSPAdes: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    contigs = out_spades / "contigs.fasta"
    if not contigs.exists():
        raise FileNotFoundError(f"No se encontró contigs.fasta en {out_spades}")

    # 2) shiver_init
    init_dir = work_dir / "MyInitDir"
    if not init_dir.exists():
        init_dir.mkdir(parents=True, exist_ok=True)
    cmd_init = [
        str(shiver_init), str(init_dir), str(config_sh),
        str(ref_alignment_fa), str(adapters_fa), str(primers_fa)
    ]
    print_substep(f"[SHIVER] shiver_init: {' '.join(cmd_init)}")
    subprocess.run(cmd_init, check=True)

    # 3) shiver_align_contigs
    cmd_align = [
        str(shiver_align), str(init_dir), str(config_sh),
        str(out_spades / 'contigs.fasta'), str(srr)
    ]
    print_substep(f"[SHIVER] shiver_align_contigs: {' '.join(cmd_align)}")
    subprocess.run(cmd_align, check=True)

    # 4) Asegurar /1 y /2 en headers de los FASTQ trimeados
    r1_for_shiver = _ensure_slash_suffix(fastq_r1_trimmed, "1")
    r2_for_shiver = _ensure_slash_suffix(fastq_r2_trimmed, "2")
    print_substep(f"[SHIVER] FASTQ para mapear: {r1_for_shiver}  {r2_for_shiver}")

    # 5) shiver_map_reads
    blast_out = work_dir / f"{srr}.blast"
    cut_refs = work_dir / f"{srr}_cut_wRefs.fasta"
    # === SIEMPRE usar los *_with12.fastq ===
    r1_with12 = Path(str(fastq_r1_trimmed).replace(".fastq", "_with12.fastq"))
    r2_with12 = Path(str(fastq_r2_trimmed).replace(".fastq", "_with12.fastq"))

    # Crea symlinks con nombres específicos para evitar colisiones en el work_dir
    tmpdir = work_dir / "__shiver_inputs__"
    tmpdir.mkdir(exist_ok=True)
    r1_src = tmpdir / f"{srr}_1_to_shivermap.fastq"
    r2_src = tmpdir / f"{srr}_2_to_shivermap.fastq"
    for src, link in [(r1_with12, r1_src), (r2_with12, r2_src)]:
        if link.exists() or link.is_symlink():
            link.unlink()
        os.symlink(src, link)
    print_substep(f"[SHIVER] FASTQ para mapear (symlinks): {r1_src}  {r2_src}")

    cmd_map = [
        str(shiver_map), str(init_dir), str(config_sh),
        str(out_spades / 'contigs.fasta'),
        str(srr), f"{srr}.blast", f"{srr}_cut_wRefs.fasta",
        str(r1_src), str(r2_src)
    ]
    print_substep(f"[SHIVER] shiver_map_reads: {' '.join(cmd_map)}")
    subprocess.run(cmd_map, check=True)

    bam_out = work_dir / f"{srr}.bam"
    if not bam_out.exists():
        # Algunas versiones nombran el bam de otra forma; intenta localizar *.bam reciente
        cands = sorted(work_dir.glob("*.bam"), key=lambda p: p.stat().st_mtime, reverse=True)
        if cands:
            return str(cands[0])
        raise FileNotFoundError("No se encontró el BAM resultante de SHIVER en el directorio de trabajo.")
    return str(bam_out)

def slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False):
    # Paso 4: Trimming
    # Intentamos trimming primero sin filtrar huérfanas:
    if layout == "PAIRED":
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)
    if layout == "SINGLE":
        trim_ok = run_trimming(fastq_single, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)

    # Si falla y es PAIRED, reintentar con filter_unpaired=True
    if not trim_ok and layout == "PAIRED":
        print("[WARN] Trimming normal falló. Reintentando con filter_unpaired=True ...")
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=True)

    if not trim_ok:
        print_error("[ERROR] Complete trimming failure, aborting pipeline.")
        sys.exit(1)

    # Paso 5: Mapeo de los reads trimeados
    print_step("Step 5: Mapping trimmed reads and sorting BAM file...")
    slow_bam_temp = work_dir / f"{srr}_slow_mapped_temp.bam"
    slow_bam = work_dir / f"{srr}_slow_mapped.bam"
    slow_sam = work_dir / f"{srr}_slow_mapped.sam"
    read1_trimmed = work_dir / "reads_trimmed_1.fastq"
    read2_trimmed = work_dir / "reads_trimmed_2.fastq"
    trimmed_reads = work_dir / "reads_trimmed.fastq"
    if not slow_bam.exists():
        if layout == "PAIRED":
            run_bwa_mem(refseq, read1_trimmed, read2_trimmed, slow_sam, slow_bam_temp, platform, min_fraction_align_reads)
        else:
            run_bwa_mem(refseq, trimmed_reads, None, slow_sam, slow_bam_temp, platform, min_fraction_align_reads)
        mv_bam_cmd = f"mv {slow_bam_temp} {slow_bam}"
        subprocess.run(mv_bam_cmd, shell=True, check=True, executable='/bin/bash')

    else:
        print_substep("[INFO] BAM file already exists. Skipping mapping...")
    return pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, slow_bam, modo, noindel)

def fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False):
    """
    Pipeline "rápido":
      1) Descargar SRA
      2) Convertir a FASTQ
      3) Mapeo [HACER A PARTIR DE AQUÍ]
      5) Ordenar y eliminar duplicados
      6) Selección de posiciones con polimorfismos [intentar hacer una función rápida]
      7) Generar fastq con los reads de interés
      8) Trimming
      9) Mapeo
      8) Análisis True/False
    """

    # Paso 4: Mapeo de los reads sin trimear
    print_step("Paso 4: Mapeo de lecturas NO trimeadas y ordenado del bam...")
    fast_bam_temp = work_dir / f"{srr}_fast_mapped_temp.bam"
    fast_bam = work_dir / f"{srr}_fast_mapped.bam"
    fast_sam = work_dir / f"{srr}_fast_mapped.sam"
    if not fast_bam.exists():
        if layout == "PAIRED":
            run_bwa_mem(refseq, fastq_1, fastq_2, fast_sam, fast_bam_temp, platform, min_fraction_align_reads)
        else:
            run_bwa_mem(refseq, fastq_single, None, fast_sam, fast_bam_temp, platform, min_fraction_align_reads)
        mv_bam_cmd = f"mv {fast_bam_temp} {fast_bam}"
        subprocess.run(mv_bam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] El archivo BAM ya existe. Saltando mapeo...")

    # Paso 5: Eliminar duplicados
    print_step("Paso 5: Eliminando duplicados y creando índice...")

    # Eliminar duplicados
    dedup_bam = work_dir / f"{srr}_fast_dedup.bam"
    dedup_bam_temp = work_dir / f"{srr}_fast_dedup_temp.bam"
    if not dedup_bam.exists():
        dedup_cmd = f"samtools markdup -r {fast_bam} {dedup_bam_temp}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, dedup_cmd):
            print_error("[ERROR] Error removing duplicates.")
            return False
        mv_deputbam_cmd = f"mv {dedup_bam_temp} {dedup_bam}"
        subprocess.run(mv_deputbam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] Deduplicated file already exists. Skipping...")

    #Creando un índice
    dedup_bam_index = work_dir / f"{srr}_fast_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(dedup_bam))
        except Exception as e:
            print_error("[ERROR] Failed to create BAM index file:", e)
            return False
    else:
        print_substep("[INFO] BAM index file already exists. Skipping...")



    # Paso 6: Generar FASTQ con los reads de interés (con >1 polimorfismo)
    print_step("Paso 6: Extrayendo FASTQ con reads de interés (con >1 polimorfismo)...")
    poly_fastq = work_dir / "reads_poly.fastq"
    poly_fastq_temp = work_dir / "reads_poly_temp.fastq"

    if not poly_fastq.exists():
        # Para obtener el nombre de la referencia se abre el BAM deduplicado:
        with pysam.AlignmentFile(str(dedup_bam), "rb") as bam_temp:
            ref_names = bam_temp.references

        pileup_file_1 = work_dir / f"{srr}_busqueda_inicial.pileup"
        print_step("Generando archivo Pileup inicial primera búsqueda...")
        if pileup_file_1.exists():    
            print_substep("[INFO] El archivo pileup inicial existe. Continuando...")
        else:
            pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file_1}") # NO GUARDAMOS LOS NOMBRES DE LOS READS POR LO QUE VIMOS EN ANTERIORES VERSIONES (NO HAY CORRESPONDENCIA EXACTA CON LAS BASES). Si quisiésemos guardarlo habría que añadir esta coletilla "--output-QNAME". Creo que la razón de esa no correspondencia puede ser que solo se filtra la calidad en ciertas columnas (se puede comprobar).
            if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
                print("[ERROR] Failed to generate initial pileup file.")
                return False

        posiciones = parse_pileup_and_extract_polimorfismos(pileup_file_1, 5)
        print(f"La longitud de la selección con pileup inicial es: {len(posiciones)}")
        
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos_1,max_len_reads, reads_in_fin = scan_region_for_polymorphisms(str(dedup_bam), ref_names[0], posiciones, 5, ref_fa, min_base_qual, noindel)
        ref_fa.close()

        # Contar cuántas posiciones polimórficas aporta cada read
        read_poly_counts = {}
        for pos, allele_dict in dic_pol_pos_1.items():
            for allele, read_set in allele_dict.items():
                for rname in read_set:
                    read_poly_counts[rname] = read_poly_counts.get(rname, 0) + 1

        # Seleccionar los reads que aparecen en más de una posición (es decir, >1 polimorfismo)
        selected_reads = {rname for rname, count in read_poly_counts.items() if count > 1}

        if not selected_reads:
            print("[INFO] No se encontraron reads con >1 polimorfismo en la región.")
        else:
            # Guardar la lista de reads en un archivo temporal (uno por línea)
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
                tmp_list = tmp.name
                for rname in selected_reads:
                    tmp.write(rname + "\n")
            cmd = (
                f"samtools view -h -F 4 {dedup_bam} | "
                f"awk 'BEGIN{{while((getline line < \"{tmp_list}\")>0) r[line]=1}} "
                f"/^@/{{print;next}} ($1 in r)' | "
                f"samtools fastq - > {poly_fastq_temp}"
            )
            print_substep(f"Ejecutando: {cmd}")
            try:
                subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
                print(f"[INFO] Archivo FASTQ con reads de interés generado: {poly_fastq}")
            except subprocess.CalledProcessError as e:
                ###print_verbose(f"[ERROR] Could not generate FASTQ file of reads of interest: {e}")
                return False
            finally:
                os.remove(tmp_list)
            mv_poly_fastq_cmd = f"mv {poly_fastq_temp} {poly_fastq}"
            subprocess.run(mv_poly_fastq_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] FASTQ file of selected reads already exists. Skipping...")


    # Paso 7: Trimming
    print_step("Step 7: Trimming reads...")
    trimmed_reads = work_dir / "reads_trimmed.fastq"
    if trimmed_reads.exists():
        print_substep("[INFO] El archivo 'reads_trimmed.fastq' ya existe. Saltando trimming.")
    else:
        trimmed_reads_temp = work_dir / "reads_trimmed_temp.fastq"
        cutadapt_cmd = f"cutadapt -a file:/mnt/netapp2/Store_uni/home/uvi/cm/lgv/DeNovo/inputs/adapters.fa -o {trimmed_reads_temp} {poly_fastq}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cutadapt_cmd):
            print_error("Error durante el trimming de los reads (single).")
            return False
        mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
        subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')    

    # Paso 8: Mapeo de los reads trimeados seleccionados
    print_step("Paso 8: Mapeo de lecturas seleccionadas trimeadas y ordenado del bam...")
    fast_bam_selected = work_dir / f"{srr}_fast_mapped_selected.bam"
    fast_bam_selected_temp = work_dir / f"{srr}_fast_mapped_selected_temp.bam"
    fast_sam_selected = work_dir / f"{srr}_fast_mapped_selected.sam"
    if not fast_bam_selected.exists():
        run_bwa_mem(refseq, trimmed_reads, None, fast_sam_selected, fast_bam_selected_temp, platform,min_fraction_align_reads)
        mv_fast_bam_selected_cmd = f"mv {fast_bam_selected_temp} {fast_bam_selected}"
        subprocess.run(mv_fast_bam_selected_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] El archivo BAM ya existe. Saltando mapeo...")

    return pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, fast_bam_selected, modo, noindel)

def maybe_move_ref_and_csv(srr: str, work_dir: Path, base_dir: Path):
    """
    Mueve {srr}_ref.fasta y {srr}_BaseFreqs_WithHXB2.csv desde work_dir a:
      - base_dir/ref_fasta/
      - base_dir/csv_with_HXB2/
    Solo actúa si existe work_dir/{srr}_ref.fasta (i.e., flujo RAVI tras assembly).
    """
    try:
        ref_src = work_dir / f"{srr}_ref.fasta"
        csv_src = work_dir / f"{srr}_BaseFreqs_WithHXB2.csv"
        moved = False
        if ref_src.exists():
            ref_dst_dir = base_dir / "ref_fasta"
            ref_dst_dir.mkdir(parents=True, exist_ok=True)
            ref_dst = ref_dst_dir / ref_src.name
            if ref_dst.exists():
                ref_dst.unlink()
            shutil.move(str(ref_src), str(ref_dst))
            print_substep(f"[MOVE] {ref_src.name} → {ref_dst_dir}")
            moved = True
        if csv_src.exists():
            csv_dst_dir = base_dir / "csv_with_HXB2"
            csv_dst_dir.mkdir(parents=True, exist_ok=True)
            csv_dst = csv_dst_dir / csv_src.name
            if csv_dst.exists():
                csv_dst.unlink()
            shutil.move(str(csv_src), str(csv_dst))
            print_substep(f"[MOVE] {csv_src.name} → {csv_dst_dir}")
            moved = True
        if not moved:
            print_verbose("[MOVE] No hay archivos de assembly para mover (saltando).")
    except Exception as e:
        print_error(f"[MOVE] Error moviendo ficheros previos al VCF: {e}")

def pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, bam, modo, noindel=False):

    print_step("\nEntrando en pasos finales ... (comunes en slow y fast)\n")

    # Paso 6: Eliminar duplicados
    print_step("Paso final 1: Eliminando duplicados y creando índice...")

    # Eliminar duplicados
    dedup_bam = work_dir / f"{srr}_slow_dedup.bam"
    dedup_bam_temp = work_dir / f"{srr}_slow_dedup_temp.bam"

    if not dedup_bam.exists():
        dedup_cmd = f"samtools markdup -r {bam} {dedup_bam_temp}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, dedup_cmd):
            print("Error al eliminar duplicados.")
            return False
        mv_deputbam_cmd = f"mv {dedup_bam_temp} {dedup_bam}"
        subprocess.run(mv_deputbam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] El archivo deduplicado ya existe. Saltando...")

    #Creando un índice
    dedup_bam_index = work_dir / f"{srr}_slow_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            print(f"Estoy aquii: {dedup_bam}")
            pysam.index(str(dedup_bam))
        except Exception as e:
            print("Error al crear el índice del BAM deduplicado:", e)
            return False
    else:
        print_substep("[INFO] El archivo indexado del bam deduplicado ya existe. Saltando...")

    # Paso 7: Leyendo posiciones

    ############################################################################
    # PRUEBAS, CREACIÓN DE PILEUP POR ARCHIVO (RECUPERADO DE VERSIÓN ANTERIOR) #
    ############################################################################

    pileup_file = work_dir / f"{srr}_slow.pileup"
    print_step("Paso final 2: Generando archivo Pileup inicial...")
    if pileup_file.exists():    
        print_substep("[INFO] Pileup file exists. Continuing...")
    else:
        pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file}") # NO GUARDAMOS LOS NOMBRES DE LOS READS POR LO QUE VIMOS EN ANTERIORES VERSIONES (NO HAY CORRESPONDENCIA EXACTA CON LAS BASES). Si quisiésemos guardarlo habría que añadir esta coletilla "--output-QNAME". Creo que la razón de esa no correspondencia puede ser que solo se filtra la calidad en ciertas columnas (se puede comprobar).
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
            print("Error al generar el archivo Pileup inicial.")
            return False

    # Definir las combinaciones de parámetros que se evaluarán:
    from itertools import product
    all_param_combinations = list(product(MIN_READS_LIST, READS_GAMETOS_MINORITARIOS_LIST, THRESHOLD_GAMETO_FREQ_LIST))

    vcf_generado = False
    coinfection_detected = True #SOLO PARA LO DE RAVI, SI NO FALSE
    for (mr, rgm, freq_val) in tqdm(all_param_combinations):
        # ================================
        # Repetir PASOS 7 y 8 con estos parámetros
        # ================================

        print_substep(
            f"\n[ANALYSIS] Parámetros => MIN_READS={mr}, "
            f"reads_gametos_minoritarios={rgm}, THRESHOLD_GAMETO_FREQ={freq_val}"
        )

        global MIN_READS #se cambian esos parámetros a nivel global
        global reads_gametos_minoritarios
        global THRESHOLD_GAMETO_FREQ

        MIN_READS = mr
        reads_gametos_minoritarios = rgm
        THRESHOLD_GAMETO_FREQ = freq_val

        posiciones_pileup = parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS)
        print_substep(f"Length of selection using pileup is: {len(posiciones_pileup)}")

        # Si no hay posiciones, escribimos la línea con campos vacíos y pasamos a la siguiente iteración
        if len(posiciones_pileup) < 2:
            result_filename = base_dir / f"results_{modo}.txt"
            with open(result_filename, "a") as rf:
                rf.write(f"Sra:{srr}\tCoinf_detected:False\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:\tNo_Rec:\tBegFin_Reads:\n")
            continue  # Se omite todo lo demás de esta iteración

        ###print_verbose(posiciones_pileup)

        ############################################################################
        # PRUEBAS, CREACIÓN DE PILEUP POR ARCHIVO (RECUPERADO DE VERSIÓN ANTERIOR) #
        ############################################################################

        print_step("Final step 3: Reading mapped positions...")   
        bam = pysam.AlignmentFile(dedup_bam, "rb") # Abrir el archivo BAM para extraer el nombre de la referencia
        ref_names = bam.references  # Obtiene una tupla de nombres de referencia
        bam.close()
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos,max_len_reads, reads_in_fin=scan_region_for_polymorphisms(dedup_bam, ref_names[0], posiciones_pileup, MIN_READS, ref_fa, min_base_qual, noindel)
        print_substep(f"Maximum read length is: {max_len_reads}")
        # Crear una instancia de PrettyPrinter con un ancho de línea mayor si es necesario
        ###print_verbose(dic_pol_pos)
        ref_fa.close()

        # Paso 8: Análisis (True/False)
        print_step("Paso final 4: Análisis de posiciones consecutivas ...")
        results_file = base_dir / "results.txt"
        results,comparaciones_hechas_global,dist_por_pos = pipeline_agrupacion_por_distribuciones(dic_pol_pos,max_len=max_len_reads,threshold=0.1,min_size=2) # este treshold es para la agrupación por grupos distintos, no para el filtrado de gametos, ES DISTINTO
        filtrar_rangos_mas_grandes(comparaciones_hechas_global)
        #resultado=comparar_posiciones_consecutivas(dic_pol_pos)
        ###print_verbose(f"\n\n[DEBUG] comparaciones_hechas_global:\n {comparaciones_hechas_global}")
        coinfection_detected = True #SOLO PARA LO DE RAVI, SI NO FALSE
        if not vcf_generado:
            for pos, freq_dict in dist_por_pos.items():
                # freq_dict es p.ej. {'A':0.3, 'T':0.7}
                count_over_01 = sum(1 for f in freq_dict.values() if f > 0.1)
                if count_over_01 >= 2:
                    coinfection_detected = True
                    break
        # Filtramos las parejas True / False, pero solo si "anotar" == True
        true_list = []
        false_list = []

        for (pA, pB), info in comparaciones_hechas_global.items():
            if info["anotar"]:  # Solo guardamos si anotar == True
                if info["res"] is True:
                    # (1) Contar cuántos rangos True (anotar sea True O False) 
                    #     contienen (pA, pB)
                    bigger_count = count_true_superranges_in(pA, pB, comparaciones_hechas_global)
                    
                    # (2) Guardar en la tupla
                    #     Si antes guardabas (pA, pB, detail),
                    #     ahora añadimos bigger_count como campo extra:
                    true_list.append((pA, pB, info["detail"], bigger_count))

                elif info["res"] is False:
                    false_list.append((pA, pB, info["detail"]))

        # Formateamos las parejas True y False
        true_pairs = ";".join([f"{pA},{pB},{bigger_count},{detail}" for (pA, pB, detail, bigger_count) in true_list])
        false_pairs = ";".join([f"{a},{b},{detail}" for (a, b, detail) in false_list])

        # (1) Formatear reads_in_fin en forma "110,300;301,450;..."
        reads_in_fin_str = ";".join(f"{st},{en}" for (st, en) in reads_in_fin)

        # Nombramos el archivo de resultados (versión según tu ejemplo)
        result_filename = base_dir / f"results_{modo}.txt"

        # Abrimos o creamos el archivo en modo append y escribimos
        with open(result_filename, "a") as rf:
            rf.write(f"Sra:{srr}\tCoinf_detected:{coinfection_detected}\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\tBegFin_Reads:{reads_in_fin_str}\n")

    if coinfection_detected and not vcf_generado:
        # 1) Creamos la carpeta si no existe
        coinfection_dir = base_dir / "vcf_coinfecciones"
        coinfection_dir.mkdir(parents=True, exist_ok=True) #IMPORTANTE: DESCOMENTAR ESTO SI SE QUIERE GENERAR EL VCF

        # 2) Mover ref/csv SOLO si venimos de assembly (detectado por presencia de {srr}_ref.fasta)
        maybe_move_ref_and_csv(srr, work_dir, base_dir)

        # 3) Nombre del VCF
        vcf_path = coinfection_dir / f"{srr}.vcf"

        # 3) Comando para generar el VCF -> IMPORTANTE: DESCOMENTAR EL SUBPROCESS SI SE QUIERE GENERAR EL VCF
        #    Ejemplo con bcftools (requiere que lo tengas instalado):
        #    samtools mpileup -B -f ref.fasta dedup.bam | bcftools call -mv -Ov -o out.vcf
        cmd_vcf = (f"bcftools mpileup -d 100000 -f {refseq} {dedup_bam} | bcftools call -mv -Ov  --ploidy 1 -o {vcf_path}")
        print_substep(f"Ejecutando: {cmd_vcf}")
        subprocess.run(cmd_vcf, shell=True, check=True)

        print_substep(f"Se ha generado el VCF para potencial coinfección: {vcf_path}")
        vcf_generado = True
        
    #print_step("Pipeline completado exitosamente.")
    return True

def pasos_ravi(srr, base_dir, refseq, work_dir, platform, dedup_bam, modo, noindel=False):

    print(f"Se ha entrado con {noindel} en el modo Ravi")

    output_bam = str(dedup_bam).replace(".bam", "_output.bam")
    intermediate_bam_0 = str(dedup_bam).replace(".bam", "_intermediate_0.bam")
    intermediate_bam_1 = str(dedup_bam).replace(".bam", "_intermediate_1.bam")
    sort_paired_cmd = f"samtools sort -n -o {intermediate_bam_1} {dedup_bam}"
    retry_command(MAX_DOWNLOAD_ATTEMPTS, sort_paired_cmd)
    #rm_intermediate_bam_cmd=f"rm {intermediate_bam_0}"
    #subprocess.run(rm_intermediate_bam_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    fixmate_cmd = f"samtools fixmate -m {intermediate_bam_1} {intermediate_bam_0}"
    retry_command(MAX_DOWNLOAD_ATTEMPTS, fixmate_cmd)
    rm_intermediate_bam_cmd_2=f"rm {intermediate_bam_1}"
    subprocess.run(rm_intermediate_bam_cmd_2, shell=True, stderr=subprocess.PIPE, text=True)
    samtools_sort_cmd = f"samtools sort -o {output_bam} {intermediate_bam_0}"
    print_substep(f"Ejecutando: {samtools_sort_cmd}")
    result_sort = subprocess.run(samtools_sort_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if result_sort.returncode != 0:
        print_output(f"[STDERR]: {result_sort.stderr}")
        raise subprocess.CalledProcessError(result_sort.returncode, samtools_sort_cmd)
    rm_intermediate_bam_cmd_3=f"rm {intermediate_bam_0}"
    subprocess.run(rm_intermediate_bam_cmd_3, shell=True, stderr=subprocess.PIPE, text=True)

    print_substep(f"[INFO] Sorted BAM file created: {output_bam}")


    #Creando un índice
    dedup_bam_index = work_dir / f"{srr}_output.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(output_bam))
        except Exception as e:
            print("Error al crear el índice del BAM deduplicado:", e)
            return False
    else:
        print_substep("[INFO] El archivo indexado del bam deduplicado ya existe. Saltando...")

    # Paso 7: Leyendo posiciones

    ############################################################################
    # PRUEBAS, CREACIÓN DE PILEUP POR ARCHIVO (RECUPERADO DE VERSIÓN ANTERIOR) #
    ############################################################################

    pileup_file = work_dir / f"{srr}_ravi.pileup"
    output_bam_dir = work_dir / f"{output_bam}"
    print_step("Paso final 2: Generando archivo Pileup inicial...")
    if pileup_file.exists():    
        print_substep("[INFO] Pileup file exists. Continuing...")
    else:
        pileup_cmd = (f"samtools mpileup -A -B -Q {min_base_qual} -f {refseq} {output_bam_dir} > {pileup_file}") # NO GUARDAMOS LOS NOMBRES DE LOS READS POR LO QUE VIMOS EN ANTERIORES VERSIONES (NO HAY CORRESPONDENCIA EXACTA CON LAS BASES). Si quisiésemos guardarlo habría que añadir esta coletilla "--output-QNAME". Creo que la razón de esa no correspondencia puede ser que solo se filtra la calidad en ciertas columnas (se puede comprobar).
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
            print("Error al generar el archivo Pileup inicial.")
            return False

    # Definir las combinaciones de parámetros que se evaluarán:
    from itertools import product
    all_param_combinations = list(product(MIN_READS_LIST, READS_GAMETOS_MINORITARIOS_LIST, THRESHOLD_GAMETO_FREQ_LIST))

    vcf_generado = False
    coinfection_detected = True #SOLO PARA LO DE RAVI, SI NO FALSE
    for (mr, rgm, freq_val) in tqdm(all_param_combinations):
        # ================================
        # Repetir PASOS 7 y 8 con estos parámetros
        # ================================

        print_substep(
            f"\n[ANALYSIS] Parámetros => MIN_READS={mr}, "
            f"reads_gametos_minoritarios={rgm}, THRESHOLD_GAMETO_FREQ={freq_val}"
        )

        global MIN_READS #se cambian esos parámetros a nivel global
        global reads_gametos_minoritarios
        global THRESHOLD_GAMETO_FREQ

        MIN_READS = mr
        reads_gametos_minoritarios = rgm
        THRESHOLD_GAMETO_FREQ = freq_val

        posiciones_pileup = parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS)
        print_substep(f"Length of selection using pileup is: {len(posiciones_pileup)}")

        # Si no hay posiciones, escribimos la línea con campos vacíos y pasamos a la siguiente iteración
        if len(posiciones_pileup) < 2:
            result_filename = base_dir / f"results_{modo}.txt"
            with open(result_filename, "a") as rf:
                rf.write(f"Sra:{srr}\tCoinf_detected:False\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:\tNo_Rec:\tBegFin_Reads:\n")
            continue  # Se omite todo lo demás de esta iteración

        ###print_verbose(posiciones_pileup)

        ############################################################################
        # PRUEBAS, CREACIÓN DE PILEUP POR ARCHIVO (RECUPERADO DE VERSIÓN ANTERIOR) #
        ############################################################################

        print_step("Final step 3: Reading mapped positions...")   
        bam = pysam.AlignmentFile(output_bam, "rb") # Abrir el archivo BAM para extraer el nombre de la referencia
        ref_names = bam.references  # Obtiene una tupla de nombres de referencia
        bam.close()
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos,max_len_reads, reads_in_fin=scan_region_for_polymorphisms(output_bam, ref_names[0], posiciones_pileup, MIN_READS, ref_fa, min_base_qual, noindel)
        print_substep(f"Maximum read length is: {max_len_reads}")
        # Crear una instancia de PrettyPrinter con un ancho de línea mayor si es necesario
        ###print_verbose(dic_pol_pos)
        ref_fa.close()

        # Paso 8: Análisis (True/False)
        print_step("Paso final 4: Análisis de posiciones consecutivas ...")
        results_file = base_dir / "results.txt"
        results,comparaciones_hechas_global,dist_por_pos = pipeline_agrupacion_por_distribuciones(dic_pol_pos,max_len=max_len_reads,threshold=0.1,min_size=2) # este treshold es para la agrupación por grupos distintos, no para el filtrado de gametos, ES DISTINTO
        filtrar_rangos_mas_grandes(comparaciones_hechas_global)
        #resultado=comparar_posiciones_consecutivas(dic_pol_pos)
        ###print_verbose(f"\n\n[DEBUG] comparaciones_hechas_global:\n {comparaciones_hechas_global}")
        coinfection_detected = True #SOLO PARA LO DE RAVI, SI NO FALSE
        if not vcf_generado:
            for pos, freq_dict in dist_por_pos.items():
                # freq_dict es p.ej. {'A':0.3, 'T':0.7}
                count_over_01 = sum(1 for f in freq_dict.values() if f > 0.1)
                if count_over_01 >= 2:
                    coinfection_detected = True
                    break
        # Filtramos las parejas True / False, pero solo si "anotar" == True
        true_list = []
        false_list = []

        for (pA, pB), info in comparaciones_hechas_global.items():
            if info["anotar"]:  # Solo guardamos si anotar == True
                if info["res"] is True:
                    # (1) Contar cuántos rangos True (anotar sea True O False) 
                    #     contienen (pA, pB)
                    bigger_count = count_true_superranges_in(pA, pB, comparaciones_hechas_global)
                    
                    # (2) Guardar en la tupla
                    #     Si antes guardabas (pA, pB, detail),
                    #     ahora añadimos bigger_count como campo extra:
                    true_list.append((pA, pB, info["detail"], bigger_count))

                elif info["res"] is False:
                    false_list.append((pA, pB, info["detail"]))

        # Formateamos las parejas True y False
        true_pairs = ";".join([f"{pA},{pB},{bigger_count},{detail}" for (pA, pB, detail, bigger_count) in true_list])
        false_pairs = ";".join([f"{a},{b},{detail}" for (a, b, detail) in false_list])

        # (1) Formatear reads_in_fin en forma "110,300;301,450;..."
        reads_in_fin_str = ";".join(f"{st},{en}" for (st, en) in reads_in_fin)

        # Nombramos el archivo de resultados (versión según tu ejemplo)
        result_filename = base_dir / f"results_{modo}.txt"

        # Abrimos o creamos el archivo en modo append y escribimos
        with open(result_filename, "a") as rf:
            rf.write(f"Sra:{srr}\tCoinf_detected:{coinfection_detected}\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\tBegFin_Reads:{reads_in_fin_str}\n")

    if coinfection_detected and not vcf_generado:
        # 1) Creamos la carpeta si no existe
        coinfection_dir = base_dir / "vcf_coinfecciones"
        coinfection_dir.mkdir(parents=True, exist_ok=True) #IMPORTANTE: DESCOMENTAR ESTO SI SE QUIERE GENERAR EL VCF

        # 2) Mover ref/csv SOLO si venimos de assembly (detectado por presencia de {srr}_ref.fasta)
        maybe_move_ref_and_csv(srr, work_dir, base_dir)

        # 3) Nombre del VCF
        vcf_path = coinfection_dir / f"{srr}.vcf"

        # 3) Comando para generar el VCF -> IMPORTANTE: DESCOMENTAR EL SUBPROCESS SI SE QUIERE GENERAR EL VCF
        #    Ejemplo con bcftools (requiere que lo tengas instalado):
        #    samtools mpileup -B -f ref.fasta dedup.bam | bcftools call -mv -Ov -o out.vcf
        cmd_vcf = (f"bcftools mpileup -d 100000 -f {refseq} {output_bam} | bcftools call -mv -Ov  --ploidy 1 -o {vcf_path}")
        print_substep(f"Ejecutando: {cmd_vcf}")
        subprocess.run(cmd_vcf, shell=True, check=True)

        print_substep(f"Se ha generado el VCF para potencial coinfección: {vcf_path}")
        vcf_generado = True
        
    #print_step("Pipeline completado exitosamente.")
    return True


def inicio_comun(srr, base_dir, refseq, modo, noindel=False):
    """
    Pipeline "lento":
      1) Descargar SRA
      2) Convertir a FASTQ
    """
    if modo == "slow":
        print("[INFO] Iniciando pipeline SLOW...")
    if modo == "fast":
        print("[INFO] Iniciando pipeline FAST...")

    work_dir = base_dir / srr
    if not work_dir.exists():
        work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    print_step("Paso 1: Identificando plataforma...")
    attempt1 = 1
    result_platform_and_layout = None

    while attempt1 <= 10:
        try:
            result_platform_and_layout = identify_platform_and_layout(srr)
            if result_platform_and_layout:
                break  # Se obtuvo la información correctamente
        except Exception as e:
            print_substep(f"Intento {attempt1} fallido al obtener información para {srr}: {e}")
        attempt1 += 1
        time.sleep(random.randint(1, 20))

    if not result_platform_and_layout:
        print("Error al obtener información del archivo SRA tras varios intentos.")
        sys.exit(1)
        return False

    platform = result_platform_and_layout["platform"]
    layout = result_platform_and_layout["layout"]

    # Paso 2: Descargar SRA (si no existe)
    print_step("Step 2: Downloading data...")
    sra_file = work_dir / f"{srr}.sra"
    if not sra_file.exists():
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, f"prefetch --output-directory {base_dir} {srr}"):
            print("Error al descargar el archivo SRA.")
            return False
    else:
        print_substep("[INFO] El archivo SRA ya existe. Saltando descarga.")

    # Paso 3: Convertir SRA a FASTQ
    print_step("Step 3: Converting SRA to FASTQ...")
    fastq_1 = work_dir / f"{srr}_1_done.fastq"
    fastq_1_temp = work_dir / f"{srr}_1.fastq"
    fastq_2 = work_dir / f"{srr}_2_done.fastq"
    fastq_2_temp = work_dir / f"{srr}_2.fastq"
    fastq_single = work_dir / f"{srr}_done.fastq"
    fastq_single_temp = work_dir / f"{srr}.fastq"

    mv_fastq_cmd_1 = f"mv {fastq_1_temp} {fastq_1}"
    mv_fastq_cmd_2 = f"mv {fastq_2_temp} {fastq_2}"
    mv_fastq_cmd = f"mv {fastq_single_temp} {fastq_single}"

    # Lógica para ILLUMINA / ONT / PACBIO (similar a run_pipeline_fast)
    if platform in ["ILLUMINA","OXFORD NANOPORE", "ONT", "OXFORD_NANOPORE"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')


        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    elif platform in ["PACBIO","PACBIO_SMRT"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ _1 and _2 files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')

        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Archivo FASTQ single ya existe. Saltando conversión.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    if modo == "slow":     
        return slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel)
    if modo == "fast":     
        return fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel)


def _sra_to_fastq(srr: str, base_dir: Path, out_dir: Path, threads: int):
    work_dir = base_dir / srr
    if not work_dir.exists():
        work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    print_step("Paso 1: Identificando plataforma...")
    attempt1 = 1
    result_platform_and_layout = None

    while attempt1 <= 10:
        try:
            result_platform_and_layout = identify_platform_and_layout(srr)
            if result_platform_and_layout:
                break  # Se obtuvo la información correctamente
        except Exception as e:
            print_substep(f"Intento {attempt1} fallido al obtener información para {srr}: {e}")
        attempt1 += 1
        time.sleep(random.randint(1, 20))

    if not result_platform_and_layout:
        print("Error al obtener información del archivo SRA tras varios intentos.")
        sys.exit(1)
        return False

    platform = result_platform_and_layout["platform"]
    layout = result_platform_and_layout["layout"]

    # Paso 2: Descargar SRA (si no existe)
    print_step("Step 2: Downloading data...")
    sra_file = work_dir / f"{srr}.sra"
    if not sra_file.exists():
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, f"prefetch --output-directory {base_dir} {srr}"):
            print("Error al descargar el archivo SRA.")
            return False
    else:
        print_substep("[INFO] El archivo SRA ya existe. Saltando descarga.")

    # Paso 3: Convertir SRA a FASTQ
    print_step("Step 3: Converting SRA to FASTQ...")
    fastq_1 = work_dir / f"{srr}_1_done.fastq"
    fastq_1_temp = work_dir / f"{srr}_1.fastq"
    fastq_2 = work_dir / f"{srr}_2_done.fastq"
    fastq_2_temp = work_dir / f"{srr}_2.fastq"
    fastq_single = work_dir / f"{srr}_done.fastq"
    fastq_single_temp = work_dir / f"{srr}.fastq"

    mv_fastq_cmd_1 = f"mv {fastq_1_temp} {fastq_1}"
    mv_fastq_cmd_2 = f"mv {fastq_2_temp} {fastq_2}"
    mv_fastq_cmd = f"mv {fastq_single_temp} {fastq_single}"

    # Lógica para ILLUMINA / ONT / PACBIO (similar a run_pipeline_fast)
    if platform in ["ILLUMINA","OXFORD NANOPORE", "ONT", "OXFORD_NANOPORE"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')


        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    elif platform in ["PACBIO","PACBIO_SMRT"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ _1 and _2 files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')

        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Archivo FASTQ single ya existe. Saltando conversión.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error al convertir SRA a FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')
    return fastq_1, fastq_2

def inicio_comun_simulacion(fastq_1, fastq_2, simulacion_folder, base_dir, refseq, modo, platform_arg="ILLUMINA", noindel=False):
    """
    Versión simulada del pipeline slow:
      - Recibe dos archivos FASTQ ya existentes.
      - Omite la descarga y conversión de SRA.
      - Asume plataforma ILLUMINA y layout PAIRED.
      - Los nombres de los archivos de salida incluirán también el nombre del directorio de simulación.
    """
    # Define un nombre de directorio para la simulación
    sim_dir = simulacion_folder  # Este nombre puede modificarse según lo que necesites
    work_dir = base_dir / simulacion_folder
    os.chdir(work_dir)

    print(f"\nYOU HAVE ENTERED SIMULATION MODE OR ARE INDICATING SPECIFIC FASTQS\n")
    if modo == "slow":
        print("[INFO] Starting SLOW pipeline (simulation/fastq)...")
    if modo == "fast":
        print("[INFO] Iniciando pipeline FAST (simulación/fastq)...")

    
    # Obtener el nombre del directorio de simulación para usarlo en los nombres de archivo
    sim_dir_name = work_dir.name

    # Forzar plataforma y layout
    platform = platform_arg.upper()
    if fastq_2:
        layout = "PAIRED"
    else:
        layout = "SINGLE"
    print(f"[INFO] Platform identified: {platform}")
    print(f"[INFO] Layout: {layout}")
    print(f"\nFastq indicated -> Jumping to step 4 ...\n")
    srr=str(simulacion_folder).split("/")[-1]
    slow_bam = work_dir / f"{srr}_slow_mapped.bam"
    
    if modo == "slow":
        return slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_1, modo, noindel)
    if modo == "fast":
        return fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_1, modo, noindel)


##########################################################
# 3) PARSE DE ARGUMENTOS
##########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Pipeline que procesa SRRs o FASTQ en modo simulación, con parámetros configurables."
    )
    # Parámetros (listas) que se pueden pasar como CSV
    parser.add_argument(
        "--MIN_READS_LIST",
        default="5",
        help="Lista separada por comas con el número mínimo de lecturas. Ej: 5,10,20"
    )
    parser.add_argument(
        "--READS_GAMETOS_MINORITARIOS_LIST",
        default="1",
        help="Lista separada por comas con el número mínimo de lecturas en gametos minoritarios. Ej: 1,10,15"
    )
    parser.add_argument(
        "--THRESHOLD_GAMETO_FREQ_LIST",
        default="0.00",
        help="Lista separada por comas con umbrales de frecuencia de gametos. Ej: 0.05,0.1,0.20"
    )

    # Flags de modo
    parser.add_argument(
        "--slow",
        action="store_true",
        help="Usar pipeline lento en lugar del rápido."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Activa el modo detallado (verbose)."
    )
    parser.add_argument(
        "--noindel",
        action="store_true",
        help="Ignora inserciones y deleciones a la hora de escanear polimorfismos."
    )

    # Flag para modo simulación
    parser.add_argument(
        "--simulacion",
        action="store_true",
        help="Activa el modo simulación con 3 posicionales: <fastq1> [fastq2] <carpeta>."
    )

    # Argumentos posicionales (uso flexible)
    parser.add_argument(
        "positional_args",
        nargs="*",
        help="Argumentos posicionales. En modo normal: <SRR>. En modo simulación: <fastq1> [fastq2] <carpeta>."
    )

    parser.add_argument(
    "--virus",
    type=str,
    default=None,
    help="Ruta al archivo refseq del virus. Si no se especifica, se usará el refseq por defecto (refseq/NC_045512.2.fna)."
    )

    parser.add_argument(
    "--platform",
    type=str,
    default="ILLUMINA",
    help="Plataforma de secuenciación para modo simulación: ILLUMINA, OXFORD_NANOPORE o PACBIO."
    )

    parser.add_argument("--bam", help="BAM deduplicado e indexado de entrada")

    # ===== Integración SPAdes + SHIVER =====
    parser.add_argument(
        "--use_shiver_assembly",
        action="store_true",
        help="Usa rnaviralSPAdes + SHIVER para generar un BAM a partir de los FASTQ TRIMEADOS."
    )
    parser.add_argument(
        "--shiver_dir",
        default="software/shiver",
        help="Ruta a la carpeta de SHIVER (contiene bin/)."
    )
    parser.add_argument(
        "--rnaviralspades",
        default="./software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py",
        help="Ruta al script rnaviralspades.py."
    )
    parser.add_argument(
        "--spades_conda_env",
        default="none",
        help="Nombre del entorno conda para ejecutar rnaviralspades.py (usa 'none' para desactivar)"
    )
    parser.add_argument(
        "--shiver_config",
        default=None,
        help="Ruta a config.sh de SHIVER (por defecto: <shiver_dir>/bin/config.sh)."
    )
    parser.add_argument(
        "--shiver_threads",
        type=int,
        default=6,
        help="Hilos para el mapeo en SHIVER (smalt/bwa)."
    )
    parser.add_argument(
        "--spades_threads",
        type=int,
        default=8,
        help="Hilos para rnaviralSPAdes."
    )
    parser.add_argument(
        "--spades_mem_gb",
        type=int,
        default=12,
        help="Memoria (GB) para rnaviralSPAdes."
    )
    parser.add_argument(
        "--shiver_ref_alignment",
        default="MyRefAligment.fasta",
        help="Alineamiento de referencia para shiver_init."
    )
    parser.add_argument(
        "--shiver_adapters",
        default="MyAdapters.fasta",
        help="FASTA de adaptadores para shiver_init."
    )
    parser.add_argument(
        "--shiver_primers",
        default="MyPrimers.fasta",
        help="FASTA de primers para shiver_init."
    )

    return parser.parse_args()

##########################################################
# 4) MAIN
##########################################################
def main():
    print_header()
    global VERBOSE
    global MIN_READS_LIST
    global READS_GAMETOS_MINORITARIOS_LIST
    global THRESHOLD_GAMETO_FREQ_LIST
    global MAX_PIPELINE_ATTEMPTS

    # 1) Parse de argumentos
    args = parse_args()

    # 2) Ajustar listas globales (si el usuario pasa algo distinto)
    MIN_READS_LIST = [int(x) for x in args.MIN_READS_LIST.split(",")]
    READS_GAMETOS_MINORITARIOS_LIST = [int(x) for x in args.READS_GAMETOS_MINORITARIOS_LIST.split(",")]
    THRESHOLD_GAMETO_FREQ_LIST = [float(x) for x in args.THRESHOLD_GAMETO_FREQ_LIST.split(",")]

    # 3) Ajustar verbose y noindel
    VERBOSE = args.verbose
    noindel_option = args.noindel

    # 4) Rutas base
    base_dir = Path(__file__).resolve().parent

    if args.virus:
        refseq = Path(args.virus).resolve()
    else:
        refseq = (base_dir / "SARS_CoV-2" / "NC_045512.2.fna").resolve()

    # 5) Lógica de simulación vs. no simulación
    if args.simulacion:
        # Se esperan 3 posicionales: fastq1, fastq2, carpeta
        # if len(args.positional_args) < 2:
        #     print_error("[ERROR] In --simulacion (local) mode, 2-3 positional arguments required: <fastq1> [fastq2] <simulation_folder>")
        #     sys.exit(1)

        # fastq1 = args.positional_args[0]
        # if len(args.positional_args) == 3:
        #     fastq2 = args.positional_args[1]
        #     simulacion_folder = args.positional_args[2]

        # if len(args.positional_args) == 2:
        #     print(f"\nSolo se ha reconocido un fastq {fastq1}, se entiende modo single")
        #     fastq2 = False
        #     simulacion_folder = args.positional_args[1]

        if args.bam:                     # arrancamos directamente del BAM
            simulacion_folder = args.positional_args[0]
            modo = "ravi"
            dedup_bam = f"{simulacion_folder}.bam"
            print(f"dedup_bam: {dedup_bam}")
            work_dir = base_dir / simulacion_folder
            srr=str(simulacion_folder).split("/")[-1]
            os.chdir(work_dir)
            # IMPORTANTE: En modo --bam SIN --use_shiver_assembly
            # mantenemos el comportamiento anterior, sin tocar la referencia.
            ok = pasos_ravi(srr, base_dir, refseq, work_dir, args.platform, dedup_bam, modo, noindel=noindel_option)
        elif args.use_shiver_assembly:
            # === Flujo: trimming -> SPAdes+SHIVER -> RAVI ===
            # A) 3 posicionales: <fastq1> <fastq2> <carpeta> (como hasta ahora)
            # B) 1 posicional SRR/ERR/DRR: descarga con prefetch+fasterq-dump en <carpeta>=<SRR>
            modo = "ravi"
            print("Entré en modo shiver")
            if len(args.positional_args) >= 3:
                fastq_1 = Path(args.positional_args[0])
                fastq_2 = Path(args.positional_args[1])
                simulacion_folder = args.positional_args[2]
                srr = str(Path(simulacion_folder).name)
                work_dir = base_dir / simulacion_folder
                work_dir.mkdir(parents=True, exist_ok=True)
            elif len(args.positional_args) == 1 and re.match(r'^[ESD]RR[0-9]+$', args.positional_args[0].strip()):
                srr = args.positional_args[0].strip()
                simulacion_folder = srr
                work_dir = base_dir / simulacion_folder
                # === Copiamos la lógica EXACTA del modo que ya funciona ===
                # 1) prefetch al base_dir  -> base_dir/<SRR>/<SRR>.sra
                # 2) fasterq-dump con la RUTA al .sra -> work_dir/<SRR>_1.fastq/_2.fastq
                fastq_1, fastq_2 = _sra_to_fastq(
                    srr=srr,
                    base_dir=base_dir,
                    out_dir=work_dir,
                    threads=getattr(args, "shiver_threads", 4)
                )
            else:
                print_error("[ERROR] --use_shiver_assembly requiere <fastq1> <fastq2> <carpeta> o bien un único SRR/ERR/DRR.")
                sys.exit(1)

            # Asumimos PAIRED cuando existen R1/R2
            layout = "PAIRED"
            os.chdir(work_dir)

            # Trimming (usa tu run_trimming habitual)
            r1_trim = work_dir / "reads_trimmed_1.fastq"
            r2_trim = work_dir / "reads_trimmed_2.fastq"
            if not (r1_trim.exists() and r2_trim.exists()):
                _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=False)
                if not _ok:
                    print("[WARN] Trimming falló. Reintentando con filter_unpaired=True ...")
                    _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=True)
                if not _ok:
                    print_error("[ERROR] Trimming falló antes de SHIVER.")
                    sys.exit(1)

            bam_out = work_dir / f"{simulacion_folder}.bam"
            print(f"Existe {bam_out}, por lo que saltamos Shiver")
            if not bam_out.exists():
                # SPAdes + SHIVER
                bam_out = run_shiver_assembly(
                    srr=srr,
                    fastq_r1_trimmed=str(r1_trim),
                    fastq_r2_trimmed=str(r2_trim),
                    work_dir=str(work_dir),
                    base_dir=str(base_dir),
                    shiver_dir=args.shiver_dir,
                    rnaviralspades=args.rnaviralspades,
                    ref_alignment_fa=args.shiver_ref_alignment,
                    adapters_fa=args.shiver_adapters,
                    primers_fa=args.shiver_primers,
                    shiver_config=args.shiver_config,
                    spades_threads=args.spades_threads,
                    spades_mem_gb=args.spades_mem_gb,
                    shiver_threads=args.shiver_threads,
                    spades_conda_env=args.spades_conda_env,
                )

            refseq_new = work_dir / f"{simulacion_folder}_ref.fasta"
            print_substep(f"[SHIVER] BAM generado: {bam_out}, entrando en ravi y mapeando {bam_out} con {refseq_new}")
            ok = pasos_ravi(srr, base_dir, refseq_new, work_dir, args.platform, bam_out, modo, noindel=noindel_option)

        # if args.slow:
        #     # MODO SIMULACIÓN SLOW
        #     modo = "slow"
        #     ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, args.platform, noindel_option)
        # else:
        #     # MODO SIMULACIÓN FAST
        #     modo = "fast"
        #     ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, args.platform, noindel_option)

        if ok:
            sys.exit(0)
            print_step("Pipeline completado exitosamente.")
        else:
            print_error("[ERROR] Simulation pipeline failed.")
            sys.exit(1)

    # Si NO es simulación => interpretamos que el primer posicional es SRR
    if len(args.positional_args) < 1:
        print_error("[ERROR] Missing SRR. Or use --simulacion <fastq1> [fastq2] <simulation_folder> (for local data).")
        sys.exit(1)

    srr = args.positional_args[0]

    # 6) Elegir pipeline según los argumentos
    pipeline_attempt = 1
    while pipeline_attempt <= MAX_PIPELINE_ATTEMPTS:
        if args.bam:
            # flujo BAM
            simulacion_folder = args.positional_args[0] if len(args.positional_args) > 0 else srr
            modo = "ravi"
            dedup_bam = f"{simulacion_folder}.bam"
            print(f"dedup_bam: {dedup_bam}")
            work_dir = base_dir / simulacion_folder
            srr = str(simulacion_folder).split("/")[-1]
            os.chdir(work_dir)
            # IMPORTANTE: En modo --bam SIN --use_shiver_assembly
            # mantenemos el comportamiento anterior, sin tocar la referencia.
            ok = pasos_ravi(srr, base_dir, refseq, work_dir, args.platform, dedup_bam, modo, noindel=noindel_option)
        elif args.use_shiver_assembly:
            # === trimming -> SPAdes + SHIVER -> RAVI ===
            print("Entré en modo shiver")
            # (tu bloque existente: con <fastq1> <fastq2> <carpeta> o SRR)
            modo = "ravi"
            if len(args.positional_args) >= 3:
                fastq_1 = Path(args.positional_args[0])
                fastq_2 = Path(args.positional_args[1])
                simulacion_folder = args.positional_args[2]
                srr = str(Path(simulacion_folder).name)
                work_dir = base_dir / simulacion_folder
                work_dir.mkdir(parents=True, exist_ok=True)
            elif len(args.positional_args) == 1 and re.match(r'^[ESD]RR[0-9]+$', args.positional_args[0].strip()):
                srr = args.positional_args[0].strip()
                simulacion_folder = srr
                work_dir = base_dir / simulacion_folder
                work_dir.mkdir(parents=True, exist_ok=True)
                fastq_1, fastq_2 = _sra_to_fastq(
                    srr=srr,
                    base_dir=base_dir,
                    out_dir=work_dir,
                    threads=getattr(args, "shiver_threads", 4)
                )
            else:
                print_error("[ERROR] --use_shiver_assembly requiere <fastq1> <fastq2> <carpeta> o bien un único SRR/ERR/DRR.")
                sys.exit(1)

            # Asumimos PAIRED cuando existen R1/R2
            layout = "PAIRED"
            os.chdir(work_dir)

            # Trimming (usa tu run_trimming habitual)
            r1_trim = work_dir / "reads_trimmed_1.fastq"
            r2_trim = work_dir / "reads_trimmed_2.fastq"
            if not (r1_trim.exists() and r2_trim.exists()):
                _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=False)
                if not _ok:
                    print("[WARN] Trimming falló. Reintentando con filter_unpaired=True ...")
                    _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=True)
                if not _ok:
                    print_error("[ERROR] Trimming falló antes de SHIVER.")
                    sys.exit(1)

            # SPAdes + SHIVER
            bam_out = work_dir / f"{simulacion_folder}.bam"
            print(f"Existe {bam_out}, por lo que saltamos Shiver")
            if not bam_out.exists():
                bam_out = run_shiver_assembly(
                    srr=srr,
                    fastq_r1_trimmed=str(r1_trim),
                    fastq_r2_trimmed=str(r2_trim),
                    work_dir=str(work_dir),
                    base_dir=str(base_dir),
                    shiver_dir=args.shiver_dir,
                    rnaviralspades=args.rnaviralspades,
                    ref_alignment_fa=args.shiver_ref_alignment,
                    adapters_fa=args.shiver_adapters,
                    primers_fa=args.shiver_primers,
                    shiver_config=args.shiver_config,
                    spades_threads=args.spades_threads,
                    spades_mem_gb=args.spades_mem_gb,
                    shiver_threads=args.shiver_threads,
                    spades_conda_env=args.spades_conda_env,
                )

            refseq_new = work_dir / f"{simulacion_folder}_ref.fasta"
            print_substep(f"[SHIVER] BAM generado: {bam_out}, entrando en ravi y mapeando {bam_out} con {refseq_new}")
            ok = pasos_ravi(srr, base_dir, refseq_new, work_dir, args.platform, bam_out, modo, noindel=noindel_option)
        elif args.slow:
            # flujo SLOW clásico -> slow_comun(...)
            modo = "slow"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)
        else:
            # flujo FAST clásico -> fast_comun(...)
            modo = "fast"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)

        if ok:
            print_step("Pipeline completado exitosamente.")
            sys.exit(0)
        else:
            print_error(f"Error en el pipeline. Reintentando... (Intento {pipeline_attempt}/{MAX_PIPELINE_ATTEMPTS})")
            pipeline_attempt += 1

    print_error("[ERROR] Pipeline failed after multiple attempts.")
    sys.exit(1)

if __name__ == "__main__":
    main()
