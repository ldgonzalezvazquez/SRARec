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

############################################
# Global configuration
############################################
min_base_qual = 20
min_fraction_align_reads=0.9
#############
#Execution variables
#############
MAX_PIPELINE_ATTEMPTS = 3
MAX_DOWNLOAD_ATTEMPTS = 3
VERBOSE = False
#module load cesga/2020 gcc/system sra-toolkit/2.10.9 edirect/20240904 cutadapt/3.5 porechop/0.2.4 bwa/0.7.17 minimap2/2.20 htslib/1.19 samtools/1.19 bcftools/1.19
############################################
# Auxiliary functions
############################################

def print_header():
    # ANSI color codes
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

    {BOLD}{YELLOW}Programa:{RESET} {CYAN}SRA{GREEN}Rec{RESET}
    {BOLD}{YELLOW}Autor:{RESET} {MAGENTA}Luis Daniel González-Vázquez{RESET}
    {BOLD}{YELLOW}Correo:{RESET} {BLUE}luisdaniel.gonzalez@uvigo.es{RESET}
    {BOLD}{YELLOW}Año:{RESET} {RED}2025{RESET}
    {BOLD}{YELLOW}Versión:{RESET} {CYAN}1.{GREEN}0{RESET}
    """

    print(header)

def print_step(message):
    """Print a main step in white."""
    print(f"\033[1;37m{message}\033[0m")  # Bright white

def print_substep(message):
    """Print a substep or executed command in blue."""
    print(f"\033[1;34m{message}\033[0m")  # Blue

def print_output(message):
    """Print the output of the executed program in green."""
    print(f"\033[1;32m{message}\033[0m")  # Green

def print_error(message):
    """Print an error in red."""
    print(f"\033[1;31m[ERROR] {message}\033[0m")

def print_verbose(message):
    """Print extra messages only if VERBOSE = True."""
    if VERBOSE:
        # Another color can be chosen, for example magenta (35)
        print(f"\033[1;35m[VERBOSE] {message}\033[0m")

def identify_platform_and_layout(srr):
    """
    Identify the sequencing platform and read layout (single-end or paired-end) for a given SRR.
    Return a dictionary containing the platform name and layout.
    """
    try:
        cmd = f"esearch -db sra -query {srr} | efetch -format runinfo"
        print_verbose(f"Identify the sequencing platform: {cmd}")
        runinfo = subprocess.check_output(cmd, shell=True, text=True).strip()

        lines = runinfo.split('\n')
        if len(lines) < 2:
            raise ValueError("No information was found in RunInfo.")

        headers = lines[0].split(',')
        data = lines[1].split(',')

        platform_index = headers.index("Platform") if "Platform" in headers else -1
        layout_index = headers.index("LibraryLayout") if "LibraryLayout" in headers else -1

        if platform_index == -1 or layout_index == -1:
            raise ValueError("The 'Platform' or 'LibraryLayout' columns were not found in RunInfo.")

        platform = data[platform_index].strip()
        layout = data[layout_index].strip()

        print(f"[INFO] Platform: {platform}")
        print(f"[INFO] Read: {layout}")

        return {"platform": platform, "layout": layout}
    except Exception as e:
        print_error(f"Error in platform: {e}")
        sys.exit(1)

def retry_command(max_attempts, cmd, shell=True):
    """
    Run a command with retries in case of failure.
    """
    attempt = 1
    while attempt <= max_attempts:
        try:
            print_substep(f"Ejecutando: {cmd}")
            process = subprocess.Popen(
                cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = process.communicate() #waits for the process to finish and captures its output (stdout) and errors (stderr)
            print_output(stdout)
            if process.returncode != 0: #process.returncode returns the command exit code
                # Also show stderr for debugging
                print_output(f"[STDERR]: {stderr}") 
                raise subprocess.CalledProcessError(process.returncode, cmd) # subprocess is the module that runs the subprocess
            return True
        except subprocess.CalledProcessError:
            print_substep(f"Attemp {attempt} of {max_attempts} failed. Retrying...")
            attempt += 1
    print_error(f"CError after {max_attempts} attemps: {cmd}")
    return False

def filter_bam_by_alignment_fraction(input_bam, output_bam, min_fraction=0.7):
    """
    Filter the input BAM and write a new BAM (output_bam) containing only
    reads whose percentage of aligned bases (query_alignment_length / query_length)
    is at least min_fraction (default: 0.5 or 50%).

    Parameters:
      input_bam (str): Path to the input BAM.
      output_bam (str): Path to the filtered output BAM.
      min_fraction (float): Minimum mapping fraction required to keep the read.
    """
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    for read in in_bam.fetch():
        # If the read is unmapped, skip it.
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        
        # read.query_alignment_length: length of the aligned portion (excluding insertions and clips)
        # read.query_length: total read length
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
    Run bwa mem or minimap2, depending on the platform, to generate a SAM file
    and then convert it to BAM. 
    """
    try:
        # Define the read group correctly
        read_group = f"@RG\\tID:{Path(fastq1).stem}\\tSM:{Path(fastq1).stem}\\tLB:Lib\\tPL:ILLUMINA" #Path(fastq1).stem is a Python expression that uses pathlib to extract the base filename without extension

        # Adjustment according to the platform
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
            raise ValueError(f"Not supported platform: {platform}") #The Python raise statement manually raises an exception and indicates that an error prevents the program from continuing

        print_substep(f"Ejecutando: {map_cmd}")

        # Run the mapping step to generate the .sam file
        result = subprocess.run(map_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            # Also show stderr if it fails
            print_output(f"[STDERR]: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, map_cmd)

        print(f"[INFO] SAM file generated: {output_sam}")

        # Convert the .sam file to .bam
        # Step 1: Convert SAM to BAM
        intermediate_bam = str(output_bam).replace(".bam", "_intermediate.bam")
        samtools_view_cmd = f"samtools view -q 30 -bS {output_sam} -o {intermediate_bam}"
        print_substep(f"Ejecutando: {samtools_view_cmd}")
        result_view = subprocess.run(samtools_view_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result_view.returncode != 0:
            print_output(f"[STDERR]: {result_view.stderr}")
            raise subprocess.CalledProcessError(result_view.returncode, samtools_view_cmd)

        #Step 2: Filter reads by alignment quality (read coverage against the reference)
        intermediate_bam_00 = str(output_bam).replace(".bam", "_intermediate_00.bam")
        sort_cmd = f"samtools sort -o {intermediate_bam_00} {intermediate_bam} "
        retry_command(MAX_DOWNLOAD_ATTEMPTS, sort_cmd)
        subprocess.run(f"samtools index {intermediate_bam_00}", shell=True, check=True)
        intermediate_bam_0 = str(output_bam).replace(".bam", "_intermediate_0.bam")
        filter_bam_by_alignment_fraction(intermediate_bam_00,intermediate_bam_0,min_fraction_align_reads)


        # Step 2: Sort the BAM by coordinate (important for mpileup)
        # <-- ADAPTED (removed "-n")
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

        print_substep(f"[INFO] BAM file created: {output_bam}")

    except subprocess.CalledProcessError as e: #as e → Capture the exception in variable e to access its attributes
        print(f"[ERROR] Error : {e.cmd}\nCódigo de salida: {e.returncode}") #e.cmd → The command that caused the error
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Error : {e}")
        sys.exit(1)

def _get_base_or_indel_at_pos_indel(read, pos, ref_fasta=None, ref_name=None, noindel=False):
    """
    Return (event, quality) for reference position 'pos' (0-based).
    
    - If noindel=True, insertions and deletions are IGNORED:
        * Only check whether a base exists in the read (qPos is not None).
        * If present, return that base and its quality.
        * If there is no base or there is a deletion, return (None, None).
    
    - If noindel=False, use the standard logic:
        * First look for an insertion "between pos-1 and pos" ..
        * Then look for a base/deletion at pos.
        * Combine them if both exist: p.ej. "A+3ACT", "-1A+3ACT".
    """

    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    # ----------------------------------------------------------------
    # "noindel" MODE => BASES ONLY
    # ----------------------------------------------------------------
    if noindel:
        for qPos, rPos in aligned_pairs:
            if rPos == pos:
                # If the read contributes a base at this position (qPos is not None),
                # return (base, quality). Ignore insertions/deletions.
                if qPos is not None:
                    base_str  = read.query_sequence[qPos].upper()
                    base_qual = read.query_qualities[qPos]
                    return (base_str, base_qual)
                else:
                    # There is a deletion -> ignore it
                    pass
        # If we reach this point, there was no base at the position => return None
        return (None, None)

    # ----------------------------------------------------------------
    # NORMAL MODE => DETECT INSERTIONS + BASE/DELETION
    # ----------------------------------------------------------------
    insertion_str  = None
    insertion_qual = None

    # (1) Look for an INSERTION between pos-1 and pos
    for i in range(len(aligned_pairs) - 1):
        q1, r1 = aligned_pairs[i]
        q2, r2 = aligned_pairs[i + 1]

        # "Entre pos-1 y pos": r1 == pos y r2 is None...
        if r1 is not None and r1 == pos and r2 is None and q2 is not None:
            # Detect insertion
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

            # INSTEAD OF USING THE FIRST ONE, TAKE THE MINIMUM QUALITY ACROSS ALL BASES
            ins_quals = [read.query_qualities[q] for q in ins_qpos]
            insertion_qual = min(ins_quals)

            break  # Take only the first insertion found

    # (2) Look for a base or deletion at 'pos'
    base_str  = None
    base_qual = None

    for i, (qPos, rPos) in enumerate(aligned_pairs):
        if rPos == pos:
            if qPos is not None:
                # => there is a base
                base_str  = read.query_sequence[qPos].upper()
                base_qual = read.query_qualities[qPos]
            else:
                # => deletion
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
                base_qual = 0  # Quality 0 (by convention for deletions)

            break

    # (3) Combine results according to what was found
    if not base_str and not insertion_str:
        return (None, None)

    if base_str and insertion_str:
        # Ej: "A+3ACT" o "-1A+3ACT"
        combined_str = base_str + insertion_str
        # Choose quality
        if base_qual == 0 and insertion_qual is not None:
            # E.g. if the base is a deletion => 0, 
            # then assign the insertion quality
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
    Return True if 'deletion_current' is a continuation of 'deletion_prev'
    considering that the positional difference is delta.
    """
    n_prev, seq_prev = parse_deletion(deletion_prev) #returns deletion size and sequence
    n_cur, seq_cur = parse_deletion(deletion_current)
    if n_prev is None or n_cur is None:
        return False
    if n_prev != n_cur + delta: #The size of the current deletion plus delta must match the size of the previously annotated deletion
        return False
    if seq_prev[delta:] == seq_cur: #returns the sequence from delta onward, i.e. what that sequence should be without the difference (the non-common part)
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

        # 1) Collect reads
        for read in bam.fetch(ref_name, pos, pos + 1):
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            # Update min/max
            r_start = read.reference_start
            r_end   = read.reference_end
            if read.query_name not in read_minmax:
                read_minmax[read.query_name] = [r_start, r_end]
            else:
                # (Optional) Adjust if the read spans more
                if r_start < read_minmax[read.query_name][0]:
                    read_minmax[read.query_name][0] = r_start
                if r_end > read_minmax[read.query_name][1]:
                    read_minmax[read.query_name][1] = r_end

            base_str, base_qual = _get_base_or_indel_at_pos_indel(
                read, pos, ref_fasta, ref_name, noindel
            )
            if not base_str:
                continue

            # Filter by quality if it is not a deletion
            if not base_str.startswith('-') and base_qual < min_base_qual:
                continue

            base_reads[base_str].add(read.query_name)
            len_reads.append(len(read.query_sequence))

        # 2) Remove conflicting reads (those assigned to >1 variant)
        read_to_variants = defaultdict(list)
        for variant, rset in base_reads.items():
            for read_name in rset:
                read_to_variants[read_name].append(variant)

        for read_name, variants_list in read_to_variants.items():
            if len(variants_list) > 1:
                # Remove that read from all variants
                for v in variants_list:
                    base_reads[v].discard(read_name)

        # 3) Filter by min_reads
        passing_bases = {
            variant: rnames
            for variant, rnames in base_reads.items()
            if len(rnames) >= min_reads
        }

        # 4) Filter continuous deletions
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

        # 5) Save only if there are >= 2 variants
        if len(passing_bases) >= 2:
            results[pos] = passing_bases

    bam.close()

    # Convert read_minmax into a list
    reads_in_fin = []
    for read_name, (st, en) in read_minmax.items():
        reads_in_fin.append((st, en))

    # Return the usual output
    return results, max(len_reads), reads_in_fin


#END FUNCTION BLOCK SCAN_REGION_FOR_POLYMORPHISMS
#START FUNCTION BLOCK COMPARAR_POSICIONES_CONSECUTIVAS:

def is_indel(variant):
    # Function to check whether a variant is an indel (insertion or deletion)
    return variant.startswith('-') or '+' in variant
    
    # This function determines how to format a variant in the "detail" field.
    #  - If it is not an insertion, return it unchanged.
    #  - If it is an insertion, check whether there is an "X+Nseq" where X is one nucleotide.
    #    and X appears as a standalone variant in variants_dict => remove X, leaving "+Nseq".

def format_variant(variants_dict, variant):
    # 1) If it is not an insertion, return it unchanged
    #    (deletion => keep it, base => keep it).
    if not '+' in variant:
        return variant
        
    # In "A+2CT"...
    #   - the "A" is the previous part
    #   - check whether "A" exists as a standalone variant in variants_dict
    #   - if yes, return "+2CT"; otherwise, return "A+2CT"

    
    idx_plus = variant.find('+')   
    # base_previa is everything before '+'
    base_previa = variant[:idx_plus]
    # The rest is the part including '+'
    insercion = variant[idx_plus:]  # p.ej. "+2CT"

    # and it exists as a standalone variant => return only the insertion part
    if base_previa in variants_dict:
        # there is a standalone "A" variant => remove the previous base
        return insercion
    else:
        # keep the previous base
        return variant

def check_comparacion_custom(posA, polyA, posB, polyB):
    """
    Perform the four-gamete test between posA and posB.
    
    Return (result, detail), where:
      - result = True/False/None
      - detail:
         * If True:  "vA1.vA2,vB1.vB2 / vA1'.vA2',vB1'.vB2' / ..." 
                    (for each quadruple that formed four gametes)
         * If False: 
             - if ALL intersections are indel-based => "posAIndels,posBIndels"
             - otherwise => ""
         * If None:  "" (no intersections)
         
    EXTRA:
      - For insertions of type "A+2CT", if "A" exists alone at that position, 
        it is not written ("+2CT" instead). If it does not exist alone, keep "A+2CT".
    """
    # 1) Coverage at each position (sum of all variant sets)
    coverageA = sum(len(rset) for rset in polyA.values())
    coverageB = sum(len(rset) for rset in polyB.values())
    coverage_min = min(coverageA, coverageB)    


    # Control variables
    any_intersection = False       # To check whether there was at least one intersection
    only_indels = True            # Will remain True if all intersections contain an indel in A or B
    posA_indels = set()           # Indels found in posA
    posB_indels = set()           # Indels found in posB

    # Store the quadruples that DO form the four gametes
    quadruplas_ok = []

    # Store information for reruns with more restrictive parameters
    gametos_data = {} 

    variantesA = list(polyA.keys())  # p.ej. ["A", "G", "A+2CT", "-3C", ...]
    variantesB = list(polyB.keys())

    for (vA1, vA2) in itertools.combinations(variantesA, 2):
        for (vB1, vB2) in itertools.combinations(variantesB, 2):
            gametos = {
                (vA1, vB1): False,
                (vA1, vB2): False,
                (vA2, vB1): False,
                (vA2, vB2): False
            }

            # Check intersections
            for (gA, gB) in gametos:
                readsA = polyA[gA]
                readsB = polyB[gB]
                inter = readsA.intersection(readsB)
                
                inter_count=len(inter)

                # Filter 1: minimum read threshold
                enough_reads = (inter_count >= reads_gametos_minoritarios)
                
                # Filter 2: frequency threshold (relative to coverage_min)
                # coverage_min may be 0; avoid division by zero
                if coverage_min > 0:
                    freq = inter_count / coverage_min
                else:
                    freq = 0.0
                
                # Less restrictive frequency filter (ej. >= 0.1 => 10%)
                enough_freq = (freq >= THRESHOLD_GAMETO_FREQ)

                # Store the information for this gamete
                # Use (gA,gB) as the key in gametos_data
                gametos_data[(gA, gB)] = {"intersect_count": inter_count, "freq": freq}

                if enough_reads and enough_freq:
                    gametos[(gA, gB)] = True
                    any_intersection = True
                    # Check whether there is an indel in gA or gB
                    if not (is_indel(gA) or is_indel(gB)):
                        only_indels = False
                    else:
                        # if gA is an indel => store it in posA_indels
                        if is_indel(gA):
                            posA_indels.add(gA)
                        # if gB is an indel => store it in posB_indels
                        if is_indel(gB):
                            posB_indels.add(gB)

            # If all four gametes were found
            if all(gametos.values()):
                # Build the output "vA1.vA2,vB1.vB2"
                # but format each one with format_variant
                fmtA1 = format_variant(polyA, vA1)
                fmtA2 = format_variant(polyA, vA2)
                fmtB1 = format_variant(polyB, vB1)
                fmtB2 = format_variant(polyB, vB2)

                combos_str = f"{fmtA1}.{fmtA2},{fmtB1}.{fmtB2}"
                quadruplas_ok.append(combos_str)

    # Determine the result
    if quadruplas_ok:
        # => we have at least one quadruple with four gametes
        detalle = "/".join(quadruplas_ok)
        return (True, detalle)
    else:
        # => four gametes were not formed
        if any_intersection:
            # => there were intersections
            if only_indels:
                # => all intersections are indel-based
                # format the indel list
                # and apply format_variant to each one
                if posA_indels:
                    posA_str = "/".join(sorted(format_variant(polyA, ind) for ind in posA_indels))
                else:
                    posA_str = ""

                if posB_indels:
                    posB_str = "/".join(sorted(format_variant(polyB, ind) for ind in posB_indels))
                else:
                    posB_str = ""

                return (False, f"{posA_str},{posB_str}")
            else:
                # => there was at least one base-base intersection => no detail
                return (False, "")
        else:
            # => no intersection => None
            return (None, "")

def comparar_posiciones_consecutivas(dic_pol_pos):
    """
    Compare pairs of consecutive positions in dic_pol_pos using 
    check_comparacion_custom, which returns (result, detail).
    
    Parameters:
      dic_pol_pos (dict): {pos -> {variante: set(read_names), ...}}
      
    Return:
      dict: {
          'True':  [(posA, posB, detalle), ...],
          'False': [(posA, posB, detalle), ...],
          'None':  [(posA, posB, detalle), ...]
      }
      where:
        - detalle is the string returned by check_comparacion_custom (p.ej. "A.G,B.T/...").
    """

    comparaciones = {
        'True': [],
        'False': [],
        'None': []
    }

    # Get the sorted positions
    sorted_positions = sorted(dic_pol_pos.keys())

    # Iterate over pairs of consecutive positions
    for i in range(len(sorted_positions) - 1):
        posA = sorted_positions[i]
        posB = sorted_positions[i + 1]

        polyA = dic_pol_pos[posA]
        polyB = dic_pol_pos[posB]

        # Call check_comparacion_custom, which returns (result, detail)
        resultado, detalle = check_comparacion_custom(posA, polyA, posB, polyB)

        # Store according to the value of 'result'
        if resultado is True:
            comparaciones['True'].append((posA, posB, detalle))
        elif resultado is False:
            comparaciones['False'].append((posA, posB, detalle))
        elif resultado is None:
            comparaciones['None'].append((posA, posB))

    return comparaciones

#END FUNCTION BLOCK COMPARAR_POSICIONES_CONSECUTIVAS
#FAST IMPLEMENTATION-SPECIFIC FUNCTION:
def extract_polymorphic_fastq(bam_path, ref_name, posiciones_a_analizar, output_fastq, min_reads=5, ref_fasta_path=None,noindel=False):
    """
    Extract reads mapping to polymorphic positions within the region [start, end)
    from the BAM file and generate a FASTQ file with those reads.
    
    Parameters:
      bam_path: Path to the BAM file.
      ref_name: Reference name (ej: "NC_045512.2").
      start: Start position (0-based).
      end: End position (0-based, no inclusiva).
      output_fastq: Path to the resulting FASTQ file.
      min_reads: Minimum number of reads required to consider a base/indel.
      ref_fasta_path: (Optional) Path to the reference FASTA file, used to extract the real sequence
                      in case of deletions.
    """
    # Open the reference FASTA if the path is provided.
    ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None

    # Obtain polymorphic positions using the previously defined function.
    polymorphic_positions,max_len_reads, reads_in_fin = scan_region_for_polymorphisms(bam_path, ref_name, posiciones_a_analizar, min_reads, ref_fasta, min_base_qual, noindel)
    
    # Close the FASTA if it was opened.
    if ref_fasta:
        ref_fasta.close()
    
    # Collect the names of all reads contributing a base/indel at polymorphic positions.
    reads_of_interest = set()
    for pos, bases in polymorphic_positions.items():
        for rnames in bases.values():
            reads_of_interest.update(rnames)
    
    if not reads_of_interest:
        print("[INFO] No reads of interest were found in the specified region.")
        return

    # Save the read list to a temporary file (one per line)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp_list = tmp.name
        for rname in reads_of_interest:
            tmp.write(rname + "\n")
    
    # Build the command to extract reads of interest and generate a FASTQ file.
    #  1) samtools view -h -F 4: extract mapped reads (with header)
    #  2) awk: filter lines whose first field (read name) is in the temporary file.
    #  3) samtools fastq -: convert the filtered SAM to FASTQ and write it to the output file.
    cmd = (
        f"samtools view -h -F 4 {bam_path} | "
        f"awk 'BEGIN{{while((getline line < \"{tmp_list}\")>0) r[line]=1}} "
        f"/^@/{{print;next}} ($1 in r)' | "
        f"samtools fastq - > {output_fastq}"
    )
    
    print(f": {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        print(f"[INFO] FASTQ file generated: {output_fastq}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR]  {e}")
        sys.exit(1)
    finally:
        # Remove the temporary read-list file.
        os.remove(tmp_list)

#END FASTA IMPLEMENTATION-SPECIFIC FUNCTION

#START PILEUP-RELATED FUNCTIONS

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
            # skip '^' + mapQuality
            i += 2

        elif c == '$':
            i += 1

        elif c in ['+', '-']:
            sign = c
            i += 1
            num_str = ''
            # read the length
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
                # no number, ignore
                continue
        
        elif c.upper() in ('A', 'C', 'G', 'T'):
            b = c.upper()
            base_counts[b] = base_counts.get(b, 0) + 1
            i += 1

        else:
            # Unsupported character -> ignore it
            i += 1

    bases_cumplen = [b for b,cnt in base_counts.items() if cnt >= MIN_READS]
    return (len(bases_cumplen) >= 2)


def parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS): 
    with open(pileup_file, 'r') as f:
        posiciones_polimorficas = []

        for line in f:
            fields = line.strip().split('\t')
            print_verbose(fields)
            if len(fields) < 6: 
                continue
            pos = int(fields[1])
            #print(f"Position: {pos}")
            bases = str(fields[4])
            if len(bases) < (MIN_READS*2):
                continue
            polimorfismo_si_o_no = quick_check_col5(bases,MIN_READS)
            if polimorfismo_si_o_no:
                print_verbose(f"Posición analizada: {pos}, bases analizadas: {bases}")
                posiciones_polimorficas.append(pos-1) #Adjust to 0-based coordinates for reading with pysam

    return posiciones_polimorficas

#END PILEUP-RELATED FUNCTIONS
#START SEARCH-ALGORITHM FUNCTIONS
def construir_distribuciones(dic_pol_pos):
    """
    Convert dic_pol_pos (pos -> {base -> set(reads), ...})
    into a dict pos -> { base -> percentage }.

    Verify that no read appears in multiple variants at the same position,
    calculate (#reads_variant / total_reads). 
    Print the result for debugging.
    """
    distributions_por_pos = {}

    print_substep("Building polymorphism distributions ...")
    for pos, variant_dict in tqdm(dic_pol_pos.items()):
        all_reads = set()
        for variant, rset in variant_dict.items():
            for r in rset:
                if r in all_reads:
                    raise ValueError(
                        f" '{r}' {pos}."
                    )
                all_reads.add(r)

        total_reads = sum(len(rset) for rset in variant_dict.values())
        if total_reads == 0:
            continue

        dist = {}
        for variant, rset in variant_dict.items():
            dist[variant] = len(rset) / total_reads

        distributions_por_pos[pos] = dist

    print_verbose(f"\n[DEBUG] Distributions by position: {distributions_por_pos}")
    return distributions_por_pos


def generar_todos_grupos_por_ventana(positions, max_read_length, min_size=2):
    """
    Generate ALL 'windows' (subsets) where 
    (pos_final - pos_inicial) <= max_read_length.
    
    Then remove windows that are a 
    strict subset of a larger window.
    """

    # == PART 1: Generate windows as before ==
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

    # == PART 2: Remove subsets ==
    # 1) Sort by size in descending order
    todos_grupos_ordenados = sorted(todos_grupos, key=len, reverse=True)

    # 2) Iterate from largest to smallest
    ventanas_finales = []
    for v in todos_grupos_ordenados:
        # If 'v' is a strict subset of any element in 'ventanas_finales'
        # discard it
        # Otherwise, add it
        if not any((v < x) for x in ventanas_finales):
            # 'v < x' significa v.issubset(x) y v != x
            ventanas_finales.append(v)

    return ventanas_finales

##############################################################################
# 1.A - Ignore the base: compare fractions
##############################################################################

def agrupar_por_distribucion_individual(distributions_por_pos, threshold=0.05):
    """
    Group positions based on the individual fractions in their distribution.
    
    For each position (e.g. 100, 101), each base (e.g. 'A', 'G') and its percentage,
    try to add the position to an existing group whose average (for that sub-distribution)
    is within the threshold (threshold) of the position fraction.
    
    Each group is represented as a dictionary with:
       - 'promedio': numeric value (float) that is the mean of the fractions in the group.
       - 'posiciones': a set with the positions that contributed.
       
    When a new position is added, the average is updated:
       new_average = (promedio * N + current_percentage) / (N+1),
       where N is the number of positions already in the group.
    
    If the position does not fit any group for that fraction, a new group is created.
    
    At the end, groups that are exactly equal to or subsets of others are removed.
    
    Return a list of groups.
    """
    grupos = []  # Each group is a dictionary: {'promedio': float, 'posiciones': set(), 'base': base}
    
    print_verbose(f"\n------------------------Distritutions by pos:\n {distributions_por_pos}\n")

    # Process positions in increasing order (although the order may not matter)
    for pos in sorted(distributions_por_pos.keys()):
        print_verbose(f"\nPOSICIÓN: {pos}")
        dist = list(set(distributions_por_pos[pos].values())) #Ejemplo [0.6, 0.4]
        print_verbose(f"DIST: {dist}")
        # For each base and its percentage at the position:
        for porcentaje in dist:
            print_verbose(f"Porcentaje: {porcentaje}, grupos len: {len(grupos)}")
            if len(grupos) > 0 and any(abs(porcentaje - grupo['promedio']) <= threshold for grupo in grupos):
                for grupo in grupos:
                    if abs(porcentaje - grupo['promedio']) <= threshold:
                        grupo['posiciones'].add(pos)
                        print_verbose(f"Last media: {grupo['promedio']}: {porcentaje}")
                        N = len(grupo['posiciones'])
                        grupo['promedio'] = (grupo['promedio'] * (N - 1) + porcentaje) / N
                        print_verbose(f"New media: {grupo['promedio']}")
            else:
                print_verbose("The group does not exist: creating ...")
                nuevo_grupo = {'promedio': porcentaje, 'posiciones': {pos}}
                print_verbose(f"Nuevo grupo: {nuevo_grupo}")
                grupos.append(nuevo_grupo)               
    
    print_verbose(f"Disorder groups: {grupos}")

    grupos = sorted(grupos, key=lambda grupo: len(grupo['posiciones']), reverse=True)

    print_verbose(f"Grupos: {grupos}")

    grupos_finales = []
    for grupo in grupos:
        pos_set = grupo['posiciones']
        if not any(pos_set == otro['posiciones'] for otro in grupos_finales): #use <= to also remove subgroups
            grupos_finales.append(grupo)
    print_verbose(f"Grupos finales: {grupos_finales}") 

    return grupos_finales

##############################################################################
# LOGIC 
##############################################################################

def existe_subtrue_mas_pequeno(p_first, p_last, comp_global):
    """
    Return True if comp_global contains a pair (pA, pB):
      1) comp_global[(pA, pB)]["res"] == True
      2) pA >= p_first y pB <= p_last
      3) (pA, pB) != (p_first, p_last)
      4) anotar == True (optional, if this should also require anotar=True)
    """
    for (pA, pB), info in comp_global.items():
        if info.get("res") is True and info.get("anotar", True) is True:
            if pA >= p_first and pB <= p_last and (pA, pB) != (p_first, p_last):
                return True
    return False

def desanotar_rangos_mas_grandes(p_first, p_last, comp_global):
    """
    Mark with anotar=False all ranges (pA, pB) en 'comp_global' que:
      - have res=True, anotar=True,
      - are different from (p_first, p_last),
      - contain (p_first, p_last), es decir pA <= p_first y pB >= p_last.

    This is used to 'discard' larger ranges when we prefer to keep
    the new range (smaller or more 'specific').
    """
    for (pA, pB), info in comp_global.items():
        # Must be True and marked for annotation
        if info.get("res") is True and info.get("anotar", True) is True:
            # Must be a different range de (p_first, p_last)
            if (pA, pB) != (p_first, p_last):
                # And contain it (inclusive)
                if pA <= p_first and pB >= p_last:
                    print_verbose(f"{(pA, pB)}")
                    comp_global[(pA, pB)]["anotar"] = False

def refine_compare(pos_list, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global): #CHECK WHETHER THIS CAN BE MADE INTO AN INDEPENDENT FUNCTION
    
    resultados = []

    if len(pos_list) < 2:
        print_verbose("JUST ONE POSITION")
        return

    p_first = pos_list[0]
    p_last  = pos_list[-1]
    print_verbose(f"Comparing: p_first = {p_first}, p_last = {p_last}")
    if p_first == p_last:
        return

    key = (p_first, p_last)


    if key in comparaciones_hechas_global:
        info = comparaciones_hechas_global[key]
        res = info["res"]
        detail = info["detail"]
        print_verbose(f"Compared: {key}: {res}, {detail} (intersected={info['intersected']})")
        # If it was already marked in groups (intersected=True), do not refine further
        if info["intersected"]:
            return
        else:
            # If it appears for the first time in groups, update the flag to True
            comparaciones_hechas_global[key]["intersected"] = True
            print_verbose(f" {key}")
    else:
        res, detail = check_comparacion_custom(p_first, dic_pol_pos[p_first], p_last, dic_pol_pos[p_last])
        # Here, in the refined-group context, the flag is set to True
        comparaciones_hechas_global[key] = {"res": res, "detail": detail, "intersected": True, "anotar": True}
        print_verbose(f"Resultado de {key} (refinamiento de grupo): {res}, {detail}, intersected: True")

    print_verbose(f"Key={key} – Res={res}")


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

        # Subdivide according to parity
        if n % 2 == 1: #odd
            mid_index = n // 2
            refine_compare(pos_list[:mid_index+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) #The +1 is correct because this is accessed like a range, so the last element is not included
            refine_compare(pos_list[mid_index:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) #Here, mid_index refers to the first element.
        else:
            half = n // 2
            refine_compare(pos_list[:half+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)
            refine_compare(pos_list[half-1:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

        # ---------------------------------------------
        # NEW ANNOTATION LOGIC (SIN before_count):
        # ---------------------------------------------
        # Store this result (res=True) in comparaciones_hechas_global,
        # but do not decide whether to annotate until checking whether a smaller sub-True already exists:
        comparaciones_hechas_global[key] = {
            "res": True,
            "detail": detail,
            "intersected": True,
            # anotar will be defined depending on whether subTrue exists
            "anotar": True  # provisional value
        }

        # If a smaller sub-True already exists, set anotar=False
        if existe_subtrue_mas_pequeno(p_first, p_last, comparaciones_hechas_global):
            print_verbose(f" {key}")
            comparaciones_hechas_global[key]["anotar"] = False
        else:
            print_verbose(f" {key}.")

        desanotar_rangos_mas_grandes(p_first, p_last, comparaciones_hechas_global)
        return

    else:  # res is None
        if len(pos_list) == 2:
            resultado_a_anotar=(p_first, p_last, None, detail)
            if resultado_a_anotar not in resultados:    
                resultados.append(resultado_a_anotar)
            return

        # Reduce p_last one by one (from right to left)
        for idx in range(len(pos_list) - 2, 0, -1): #For each iteration, this returns an idx that is one position less each time from the last position, skipping the first; for example, if it were 3 it would be 3,2,1 (len would be 4).
            list_candidate=pos_list[:idx+1]
            p_first = list_candidate[0]
            p_last  = list_candidate[-1]
            resultado_a_anotar=(p_first, p_last, None, detail)
            if resultado_a_anotar not in resultados:    
                resultados.append(resultado_a_anotar)
            refine_compare(list_candidate, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)
        return

##############################################################################
# compare_fn and comparar_grupos_interseccion
##############################################################################

def comparar_grupos_interseccion(grupos, check_comparacion_custom, dic_pol_pos, comparaciones_hechas_global):
    print_verbose("\nENTERING THE comparar_grupos_interseccion FUNCTION\n")
    
    resultados = {'True': [], 'False': [], 'None': []}
    etiquetados = {}
    print_verbose(f"grupos: {grupos}")
    
    # Build a dictionary associating each position with its label(s), i.e. the averages of the groups it belongs to
    for grupo in grupos:
        for p in grupo["posiciones"]:
            if p in etiquetados:
                etiquetados[p] = f"{etiquetados[p]},{grupo['promedio']}" #This concatenates the average values found at that position: it takes the averages already present and appends the new one separated by a comma
            else:
                etiquetados[p] = str(grupo["promedio"])
    
    print_verbose(f"etiquetados: {etiquetados}")
    posiciones_ordenadas = sorted(etiquetados.keys())
    print_verbose(f"posiciones_ordenadas: {posiciones_ordenadas}")
    
    # Iterate over all pairs (pA, pB) with pA < pB to avoid repetitions
    for i in range(len(posiciones_ordenadas)):
        for j in range(i + 1, len(posiciones_ordenadas)):
            pA = posiciones_ordenadas[i]
            pB = posiciones_ordenadas[j]
            # Convert the labels to float lists for comparison
            gA = sorted(map(float, etiquetados[pA].split(',')))
            gB = sorted(map(float, etiquetados[pB].split(',')))
            
            print_verbose(f"Comparando: pA: {pA}, gA: {gA} | pB: {pB}, gB: {gB}")
            
            # 
            if gA == gB:
                print_verbose(f"{pA} y {pB} are equal; .")
                continue
            
            # Use the sorted tuple (minimum, maximum) to avoid duplicates
            pair = (min(pA, pB), max(pA, pB))
            if pair not in comparaciones_hechas_global:
                res, detail = check_comparacion_custom(pA, dic_pol_pos[pA], pB, dic_pol_pos[pB])
                print_verbose(f"Key={pair} – Res={res}")
                comparaciones_hechas_global[pair] = {"res": res, "detail": detail, "intersected": False, "anotar": True}
                if existe_subtrue_mas_pequeno(pair[0], pair[1], comparaciones_hechas_global):
                    print_verbose(f"{pair}")
                    comparaciones_hechas_global[pair]["anotar"] = False
                else:
                    print_verbose(f" {pair}")

                desanotar_rangos_mas_grandes(pair[0], pair[1], comparaciones_hechas_global)
                print_verbose(f": {comparaciones_hechas_global}")
                polyA = dic_pol_pos[pA]
                polyB = dic_pol_pos[pB]
                print_verbose(f"Comparing position {pA} (reads: {polyA}) with {pB} (reads: {polyB})")
                res, detail = check_comparacion_custom(pA, polyA, pB, polyB)
                if res is True:
                    resultados['True'].append((pA, pB, detail))
                elif res is False:
                    resultados['False'].append((pA, pB, detail))
                else:
                    resultados['None'].append((pA, pB, detail))
            else:
                print_verbose(f" {pair}")
    
    return resultados

##############################################################################
# 2. MAIN FUNCTION
##############################################################################

def pipeline_agrupacion_por_distribuciones(dic_pol_pos, max_len=5, threshold=0.05, min_size=2):
    """
    1) Build distributions
    2) Generate windows
    3) Use a global set 'comparaciones_hechas_global' to avoid repeating 
       comparisons across ALL windows
    4) Para cada ventana:
       - Filter positions => sub_dist
       - Group while ignoring the base
       - Compare BETWEEN groups (adjacent)
       - Compare WITHIN each group using the refined function
    5) Return a dict with the information
    """

    dist_por_pos = construir_distribuciones(dic_pol_pos)

    print_verbose(f"dist_por_pos: {dist_por_pos}")

    ventanas = generar_todos_grupos_por_ventana(dist_por_pos.keys(), max_len, min_size)

    print_verbose(f"ventanas: {ventanas}")

    results_all_windows = {}
    comparaciones_hechas_global = {} 

    print_substep("Analyzing positions ...")
    for i, window in enumerate(ventanas):
        sub_dist = {pos: dist_por_pos[pos] for pos in window}
        grupos_subdist = agrupar_por_distribucion_individual(sub_dist, threshold)

        print_verbose(f"grupos_subdist: {grupos_subdist}")

        comp_entre = comparar_grupos_interseccion(grupos_subdist, check_comparacion_custom, dic_pol_pos, comparaciones_hechas_global)

        print_verbose(f"comp_entre: {comp_entre}")

        comp_en_grupos = {}
        for j, g in enumerate(grupos_subdist):
            print_verbose(f"G:{g}\n")
            # Use the refined version
            pos_sorted = sorted(g["posiciones"])
            comp_en_grupos[j] = refine_compare(pos_sorted, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

            print_verbose(f"comp_en_grupos[j]: {comp_en_grupos[j]}")

        results_all_windows[i] = {"window_positions": sorted(window),"grupos": grupos_subdist,"comparaciones_entre_grupos": comp_entre,"comparaciones_en_grupos": comp_en_grupos}

        print_verbose(f"results_all_windows[i]: {results_all_windows[i]}")

    return results_all_windows, comparaciones_hechas_global, dist_por_pos
# END SEARCH-ALGORITHM FUNCTIONS

def run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False):
    """
    Perform read trimming. If 'filter_unpaired=True' and 'layout' is PAIRED,
    first run repair.sh to discard/repair orphan pairs,
    then run the standard trimming with cutadapt.
    """

    print_step("Step 4: Read trimming...")

    trimmed_reads = work_dir / "reads_trimmed.fastq"
    read1_trimmed = work_dir / "reads_trimmed_1.fastq"
    read2_trimmed = work_dir / "reads_trimmed_2.fastq"
    read1_trimmed_temp = work_dir / "reads_trimmed_1_temp.fastq"
    read2_trimmed_temp = work_dir / "reads_trimmed_2_temp.fastq"
    trimmed_reads_temp = work_dir / "reads_trimmed.fastq.temp"

    # If the output files already exist, skip this step
    if (read1_trimmed.exists() and read2_trimmed.exists()) or trimmed_reads.exists():
        print_substep("[INFO] Trimming files already exist. Skipping trimming.")
        return True

    try:
        # ======================================
        # PAIRED-END
        # ======================================
        if layout == "PAIRED":
            # ------------------------------------------------------------------------
            # (1) If filter_unpaired=True -> run repair.sh BEFORE trimming
            # ------------------------------------------------------------------------
            MAX_PIPELINE_ATTEMPTS=1
            if filter_unpaired:
                MAX_PIPELINE_ATTEMPTS=3
                # repair output paths
                repaired_1 = work_dir / "reads_repaired_1.fastq"
                repaired_2 = work_dir / "reads_repaired_2.fastq"
                singletons = work_dir / "singletons.fastq"

                repair_cmd = (
                    f"/opt/cesga/2022/software/Compiler/gcccore/system/bbmap/39.19/repair.sh "
                    f"in={fastq_1} in2={fastq_2} "
                    f"out={repaired_1} out2={repaired_2} outs={singletons}"
                    )
                print_substep(f"Running repair.sh: {repair_cmd}")
                subprocess.run(repair_cmd, shell=True, check=True)

                # From here onward, use repaired_1 and repaired_2
                actual_fastq_1 = repaired_1
                actual_fastq_2 = repaired_2

            else:
                # If orphan reads do not need to be filtered/discarded, use the original FASTQ files
                actual_fastq_1 = fastq_1
                actual_fastq_2 = fastq_2

            # ------------------------------------------------------------------------
            # (2) Run cutadapt in paired-end mode
            # ------------------------------------------------------------------------
            cutadapt_cmd = (f"cutadapt -a file:{base_dir}/adapters.fa -A file:{base_dir}/adapters.fa -o {read1_trimmed_temp} -p {read2_trimmed_temp} {actual_fastq_1} {actual_fastq_2}")

            print_substep(f"Ejecutando cutadapt (Paired-end): {cutadapt_cmd}")
            if not retry_command(MAX_PIPELINE_ATTEMPTS, cutadapt_cmd):
                raise RuntimeError("Error during paired-end trimming.")

            # Rename .temp files to final names
            mv_trim_cmd_1 = f"mv {read1_trimmed_temp} {read1_trimmed}"
            mv_trim_cmd_2 = f"mv {read2_trimmed_temp} {read2_trimmed}"
            subprocess.run(mv_trim_cmd_1, shell=True, check=True, executable='/bin/bash')
            subprocess.run(mv_trim_cmd_2, shell=True, check=True, executable='/bin/bash')

        # ======================================
        # SINGLE-END
        # ======================================
        else:
            MAX_PIPELINE_ATTEMPTS=3
            # SINGLE-END => there is no previous repair step
            cutadapt_cmd = (f"cutadapt -a file:{base_dir}/adapters.fa -o {trimmed_reads_temp} {fastq_1}")
            print_substep(f"Ejecutando cutadapt (Single-end): {cutadapt_cmd}")
            if not retry_command(MAX_PIPELINE_ATTEMPTS, cutadapt_cmd):
                raise RuntimeError("Error during single-end trimming.")

            mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
            subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')

        print_substep("[INFO] Trimming completed successfully.")
        return True

    except Exception as e:
        print_error(f"Trimming falló: {e}")
        return False

def count_true_superranges_in(pA, pB, comp_global):
    """
    Return how many ranges (xA, xB) con res=True (regardless of anotar)
    contain or match (pA, pB). Es decir: xA <= pA y xB >= pB.
    """
    count = 0
    for (xA, xB), info in comp_global.items():
        if info.get("res") is True:  # only True matters
            # si xA <= pA y xB >= pB => (xA, xB) contiene (pA, pB)
            if xA <= pA and xB >= pB:
                count += 1
    return count

def slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False):
    # Step 4: Trimming
    # First try trimming without filtering orphan reads:
    if layout == "PAIRED":
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)
    if layout == "SINGLE":
        trim_ok = run_trimming(fastq_single, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)

    # If it fails and layout is PAIRED, retry with filter_unpaired=True
    if not trim_ok and layout == "PAIRED":
        print("[WARN] Standard trimming failed. Retrying with filter_unpaired=True ...")
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=True)

    if not trim_ok:
        print_error("Complete trimming failure, aborting pipeline.")
        sys.exit(1)

    # Step 5: Mapping trimmed reads
    print_step("Step 5: Mapping trimmed reads and sorting the BAM...")
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
        print_substep("[INFO] The BAM file already exists. Skipping mapping...")
    return pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, slow_bam, modo, noindel=False)

def fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False):
    """
    "fast" pipeline:
      1) Download SRA
      2) Convert to FASTQ
      3) Mapping [CONTINUE FROM HERE]
      5) Sort and remove duplicates
      6) Selection of polymorphic positions [try to implement a fast function]
      7) Generate FASTQ with reads of interest
      8) Trimming
      9) Mapping
      8) True/False analysis
    """

    # Step 4: Mapping untrimmed reads
    print_step("Step 4: Mapping UNTRIMMED reads and sorting the BAM...")
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
        print_substep("[INFO] The BAM file already exists. Skipping mapping...")

    # Step 5: Remove duplicates
    print_step("Step 5: Removing duplicates and creating index...")

    dedup_bam = work_dir / f"{srr}_fast_dedup.bam"
    dedup_bam_temp = work_dir / f"{srr}_fast_dedup_temp.bam"
    if not dedup_bam.exists():
        dedup_cmd = f"samtools markdup -r {fast_bam} {dedup_bam_temp}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, dedup_cmd):
            print_error("Error removing duplicates.")
            return False
        mv_deputbam_cmd = f"mv {dedup_bam_temp} {dedup_bam}"
        subprocess.run(mv_deputbam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] The deduplicated file already exists. Skipping...")

    #Creating an index
    dedup_bam_index = work_dir / f"{srr}_fast_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(dedup_bam))
        except Exception as e:
            print_error("Error creating the deduplicated BAM index:", e)
            return False
    else:
        print_substep("[INFO] The indexed deduplicated BAM file already exists. Skipping...")



    # Step 6: Generate FASTQ with reads of interest (with >1 polymorphism)
    print_step("Step 6: Extracting FASTQ with reads of interest (with >1 polymorphism)...")
    poly_fastq = work_dir / "reads_poly.fastq"
    poly_fastq_temp = work_dir / "reads_poly_temp.fastq"

    if not poly_fastq.exists():
        # Open the deduplicated BAM to obtain the reference name:
        with pysam.AlignmentFile(str(dedup_bam), "rb") as bam_temp:
            ref_names = bam_temp.references

        pileup_file_1 = work_dir / f"{srr}_busqueda_inicial.pileup"
        print_step("Generating initial pileup file for the first search...")
        if pileup_file_1.exists():    
            print_substep("[INFO] The initial pileup file already exists. Continuing...")
        else:
            pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file_1}") # READ NAMES ARE NOT SAVED BECAUSE, IN PREVIOUS VERSIONS, THERE WAS NO EXACT CORRESPONDENCE WITH THE BASES. To save them, the "--output-QNAME" option would need to be added. The reason for this mismatch may be that quality filtering is only applied to certain columns (this can be checked).
            if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
                print("Error generating the initial pileup file.")
                return False

        posiciones = parse_pileup_and_extract_polimorfismos(pileup_file_1, 5)
        print(f": {len(posiciones)}")
        
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos_1,max_len_reads, reads_in_fin = scan_region_for_polymorphisms(str(dedup_bam), ref_names[0], posiciones, 5, ref_fa, min_base_qual, noindel)
        ref_fa.close()

        # Count how many polymorphic positions each read contributes
        read_poly_counts = {}
        for pos, allele_dict in dic_pol_pos_1.items():
            for allele, read_set in allele_dict.items():
                for rname in read_set:
                    read_poly_counts[rname] = read_poly_counts.get(rname, 0) + 1

        # Select reads that appear in more than one position (i.e. >1 polymorphism)
        selected_reads = {rname for rname, count in read_poly_counts.items() if count > 1}

        if not selected_reads:
            print("[INFO] No reads with >1 polymorphism were found in the region.")
        else:
            # Save the read list to a temporary file (one per line)
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
                print(f"[INFO]  FASTQ file generated: {poly_fastq}")
            except subprocess.CalledProcessError as e:
                print_verbose(f"[ERROR] FASTQ file could not be generated: {e}")
                return False
            finally:
                os.remove(tmp_list)
            mv_poly_fastq_cmd = f"mv {poly_fastq_temp} {poly_fastq}"
            subprocess.run(mv_poly_fastq_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] The FASTQ file with selected reads already exists. Skipping...")


    # Step 7: Trimming
    print_step("Step 7: Read trimming...")
    trimmed_reads = work_dir / "reads_trimmed.fastq"
    if trimmed_reads.exists():
        print_substep("[INFO] The 'reads_trimmed.fastq' file already exists. Skipping trimming.")
    else:
        trimmed_reads_temp = work_dir / "reads_trimmed_temp.fastq"
        cutadapt_cmd = f"cutadapt -a file:{base_dir}/adapters.fa -o {trimmed_reads_temp} {poly_fastq}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cutadapt_cmd):
            print_error("Error during read trimming (single-end).")
            return False
        mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
        subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')    

    # Step 8: Mapping selected trimmed reads
    print_step("Step 8: Mapping selected trimmed reads and sorting the BAM...")
    fast_bam_selected = work_dir / f"{srr}_fast_mapped_selected.bam"
    fast_bam_selected_temp = work_dir / f"{srr}_fast_mapped_selected_temp.bam"
    fast_sam_selected = work_dir / f"{srr}_fast_mapped_selected.sam"
    if not fast_bam_selected.exists():
        run_bwa_mem(refseq, trimmed_reads, None, fast_sam_selected, fast_bam_selected_temp, platform,min_fraction_align_reads)
        mv_fast_bam_selected_cmd = f"mv {fast_bam_selected_temp} {fast_bam_selected}"
        subprocess.run(mv_fast_bam_selected_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] The BAM file already exists. Skipping mapping...")

    return pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, fast_bam_selected, modo, noindel=False)


def pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, bam, modo, noindel=False):

    print_step("\nEntering final steps ... (common to slow and fast)\n")

    # Step 6: Remove duplicates
    print_step("Final step 1: Removing duplicates and creating index...")

    # Eliminar duplicados
    dedup_bam = work_dir / f"{srr}_slow_dedup.bam"
    dedup_bam_temp = work_dir / f"{srr}_slow_dedup_temp.bam"

    if not dedup_bam.exists():
        dedup_cmd = f"samtools markdup -r {bam} {dedup_bam_temp}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, dedup_cmd):
            print("Error removing duplicates.")
            return False
        mv_deputbam_cmd = f"mv {dedup_bam_temp} {dedup_bam}"
        subprocess.run(mv_deputbam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] The deduplicated file already exists. Skipping...")

    #Creating an index
    dedup_bam_index = work_dir / f"{srr}_slow_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(dedup_bam))
        except Exception as e:
            print("Error creating the deduplicated BAM index:", e)
            return False
    else:
        print_substep("[INFO] The indexed deduplicated BAM file already exists. Skipping...")

    # Step 7: Reading positions

    ############################################################################
    # TESTS, PILEUP FILE CREATION (RECOVERED FROM A PREVIOUS VERSION) #
    ############################################################################

    pileup_file = work_dir / f"{srr}_slow.pileup"
    print_step("Final step 2: Generating initial pileup file...")
    if pileup_file.exists():    
        print_substep("[INFO] The pileup file already exists. Continuing...")
    else:
        pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file}") # READ NAMES ARE NOT SAVED BECAUSE, IN PREVIOUS VERSIONS, THERE WAS NO EXACT CORRESPONDENCE WITH THE BASES. To save them, the "--output-QNAME" option would need to be added. The reason for this mismatch may be that quality filtering is only applied to certain columns (this can be checked).
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
            print("Error generating the initial pileup file.")
            return False

    # Define the parameter combinations to be evaluated:
    from itertools import product
    all_param_combinations = list(product(MIN_READS_LIST, READS_GAMETOS_MINORITARIOS_LIST, THRESHOLD_GAMETO_FREQ_LIST))

    for (mr, rgm, freq_val) in tqdm(all_param_combinations):
        # ================================
        # Repeat STEPS 7 and 8 with these parameters
        # ================================

        print_substep(
            f"\n[ANALYSIS] Parámetros => MIN_READS={mr}, "
            f"reads_gametos_minoritarios={rgm}, THRESHOLD_GAMETO_FREQ={freq_val}"
        )

        global MIN_READS #these parameters are changed globally
        global reads_gametos_minoritarios
        global THRESHOLD_GAMETO_FREQ

        MIN_READS = mr
        reads_gametos_minoritarios = rgm
        THRESHOLD_GAMETO_FREQ = freq_val

        posiciones_pileup = parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS)
        print_substep(f"La longitud de la selección con pileup es: {len(posiciones_pileup)}")

        # If there are no positions, write the line with empty fields and move to the next iteration
        if len(posiciones_pileup) < 2:
            result_filename = base_dir / f"results_{modo}.txt"
            with open(result_filename, "a") as rf:
                rf.write(
                    f"Min_read_num:{mr}\t"
                    f"Min_gameto_num:{rgm}\t"
                    f"Min_gameto_freq:{freq_val}\t"
                    f"Sra:{srr}\tPlatform:{platform}\t"
                    f"Rec:\tNo_Rec:\tBegFin_Reads:\n"  # Empty fields
                )
            continue  # Skip the rest of this iteration

        print_verbose(posiciones_pileup)

        ############################################################################
        # TESTS, PILEUP FILE CREATION (RECOVERED FROM A PREVIOUS VERSION) #
        ############################################################################

        print_step("Final step 3: Reading mapped positions...")   
        bam = pysam.AlignmentFile(dedup_bam, "rb") # Open the BAM file to extract the reference name
        ref_names = bam.references  # Gets a tuple of reference names
        bam.close()
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos,max_len_reads, reads_in_fin=scan_region_for_polymorphisms(dedup_bam, ref_names[0], posiciones_pileup, MIN_READS, ref_fa, min_base_qual, noindel)
        print_substep(f"EL TAMAÑO MÁXIMO DE READ ES: {max_len_reads}")
        # Create a PrettyPrinter instance with a larger line width if needed
        print_verbose(dic_pol_pos)
        ref_fa.close()

        # Step 8: Analysis (True/False)
        print_step("Final step 4: Consecutive-position analysis ...")
        results_file = base_dir / "results.txt"
        results,comparaciones_hechas_global,dist_por_pos = pipeline_agrupacion_por_distribuciones(dic_pol_pos,max_len=max_len_reads,threshold=0.1,min_size=2) # this threshold is for grouping distinct groups, not for gamete filtering; IT IS DIFFERENT
        #resultado=comparar_posiciones_consecutivas(dic_pol_pos)
        print_verbose(f"\n\n[DEBUG] comparaciones_hechas_global:\n {comparaciones_hechas_global}")
        coinfection_detected = False
        # Filter the True / False pairs, but only if "anotar" == True
        true_list = []
        false_list = []

        for (pA, pB), info in comparaciones_hechas_global.items():
            if info["anotar"]:  # Only save if anotar == True
                if info["res"] is True:
                    bigger_count = count_true_superranges_in(pA, pB, comparaciones_hechas_global)
                    
                    # (2) Save in the tuple
                    #     If previously saving (pA, pB, detail),
                    #     now add bigger_count as an extra field:
                    true_list.append((pA, pB, info["detail"], bigger_count))

                elif info["res"] is False:
                    false_list.append((pA, pB, info["detail"]))

        # Format the True and False pairs
        true_pairs = ";".join([f"{pA},{pB},{bigger_count},{detail}" for (pA, pB, detail, bigger_count) in true_list])
        false_pairs = ";".join([f"{a},{b},{detail}" for (a, b, detail) in false_list])

        # (1) Format reads_in_fin as "110,300;301,450;..."
        reads_in_fin_str = ";".join(f"{st},{en}" for (st, en) in reads_in_fin)

        # Name the results file (version according to your example)
        result_filename = base_dir / f"results_{modo}.txt"

        # Open or create the file in append mode and write to it
        with open(result_filename, "a") as rf:
            rf.write(f"Min_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tSra:{srr}\tPlatform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\tBegFin_Reads:{reads_in_fin_str}\n")
        
    print_step("Pipeline completed successfully.")
    return True

def inicio_comun(srr, base_dir, refseq, modo, noindel=False):
    """
    "slow" pipeline:
      1) Download SRA
      2) Convert to FASTQ
    """
    if modo == "slow":
        print("[INFO] Starting SLOW pipeline...")
    if modo == "fast":
        print("[INFO] Starting FAST pipeline...")

    work_dir = base_dir / srr
    if not work_dir.exists():
        work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    print_step("Step 1: Identifying platform...")
    attempt1 = 1
    result_platform_and_layout = None

    while attempt1 <= 10:
        try:
            result_platform_and_layout = identify_platform_and_layout(srr)
            if result_platform_and_layout:
                break  # The information was obtained successfully
        except Exception as e:
            print_substep(f"Attemp {attempt} failed {srr}: {e}")
        attempt1 += 1
        time.sleep(random.randint(1, 20))

    if not result_platform_and_layout:
        print("Error obtaining information from the SRA file after several attempts.")
        sys.exit(1)
        return False

    platform = result_platform_and_layout["platform"]
    layout = result_platform_and_layout["layout"]

    # Step 2: Download SRA (if it does not exist)
    print_step("Step 2: Downloading data...")
    sra_file = work_dir / f"{srr}.sra"
    if not sra_file.exists():
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, f"prefetch --output-directory {base_dir} {srr}"):
            print("Error downloading the SRA file.")
            return False
    else:
        print_substep("[INFO] The SRA file already exists. Skipping download.")

    # Step 3: Convert SRA to FASTQ
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

    # Logic for ILLUMINA / ONT / PACBIO (similar to run_pipeline_fast)
    if platform in ["ILLUMINA","OXFORD NANOPORE", "ONT", "OXFORD_NANOPORE"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error converting SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')


        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error converting SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    elif platform in ["PACBIO","PACBIO_SMRT"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ _1 and _2 files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error converting SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')

        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Error converting SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    if modo == "slow":     
        return slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False)
    if modo == "fast":     
        return fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False)

def inicio_comun_simulacion(fastq_1, fastq_2, simulacion_folder, base_dir, refseq, modo, noindel=False):
    """
    Simulated version of the slow pipeline:
      - Receives two existing FASTQ files.
      - Skips SRA download and conversion.
      - Assumes ILLUMINA platform and PAIRED layout.
      - Output filenames will also include the simulation directory name.
    """
    # Define a directory name for the simulation
    sim_dir = simulacion_folder  # This name can be modified as needed
    work_dir = base_dir / simulacion_folder
    os.chdir(work_dir)

    print(f"\nSIMULATION MODE\n")
    if modo == "slow":
        print("[INFO] Starting SLOW pipeline (simulation/FASTQ)...")
    if modo == "fast":
        print("[INFO] Starting FAST pipeline (simulation/FASTQ)...")

    
    # Obtain the simulation directory name to use in filenames
    sim_dir_name = work_dir.name

    # Force platform and layout
    platform = "ILLUMINA"
    if fastq_2:
        layout = "PAIRED"
    else:
        layout = "SINGLE"
    print(f"[INFO] Platform: {platform}")
    print(f"[INFO] Read: {layout}")
    print(f"\nFastq -> going to step 4 ...\n")
    srr=str(simulacion_folder).split("/")[-1]
    slow_bam = work_dir / f"{srr}_slow_mapped.bam"
    
    if modo == "slow":
        return slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_1, modo, noindel=False)
    if modo == "fast":
        return fast_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_1, modo, noindel=False)


##########################################################
# 3) ARGUMENT PARSING
##########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Pipeline that processes SRRs or FASTQ files in simulation mode, with configurable parameters."
    )
    # Parameters (listas) que se pueden pasar como CSV
    parser.add_argument(
        "--MIN_READS_LIST",
        default="10",
        help="Comma-separated list with the minimum number of reads. Example: 5,10,20"
    )
    parser.add_argument(
        "--READS_GAMETOS_MINORITARIOS_LIST",
        default="10",
        help="Comma-separated list with the minimum number of reads in minor gametes. Example: 1,10,15"
    )
    parser.add_argument(
        "--THRESHOLD_GAMETO_FREQ_LIST",
        default="0.05",
        help="Comma-separated list with gamete-frequency thresholds. Example: 0.05,0.1,0.20"
    )

    # Mode flags
    parser.add_argument(
        "--slow",
        action="store_true",
        help="Use the slow pipeline instead of the fast one."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed mode (verbose)."
    )
    parser.add_argument(
        "--noindel",
        action="store_true",
        help="Ignore insertions and deletions when scanning polymorphisms."
    )

    # Flag for simulation mode
    parser.add_argument(
        "--simulacion",
        action="store_true",
        help="Enable simulation mode with 3 positional arguments: <fastq1> [fastq2] <folder>."
    )

    # Positional arguments (uso flexible)
    parser.add_argument(
        "positional_args",
        nargs="*",
        help="Positional arguments. In normal mode: <SRR>. In simulation mode: <fastq1> [fastq2] <folder>."
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

    # 1) Parse arguments
    args = parse_args()

    # 2) Set global lists (if the user provides different values)
    MIN_READS_LIST = [int(x) for x in args.MIN_READS_LIST.split(",")]
    READS_GAMETOS_MINORITARIOS_LIST = [int(x) for x in args.READS_GAMETOS_MINORITARIOS_LIST.split(",")]
    THRESHOLD_GAMETO_FREQ_LIST = [float(x) for x in args.THRESHOLD_GAMETO_FREQ_LIST.split(",")]

    # 3) Set verbose and noindel
    VERBOSE = args.verbose
    noindel_option = args.noindel

    # 4) Base paths
    base_dir = Path(__file__).resolve().parent
    refseq = base_dir / "refseq" / "NC_045512.2.fna"

    # 5) Simulation vs. non-simulation logic
    if args.simulacion:
        # Expecting 3 positional arguments: fastq1, fastq2, folder
        if len(args.positional_args) < 2:
            print_error("In --simulacion mode, 2-3 positional arguments are required: <fastq1> [fastq2] <folderSimulacion>")
            sys.exit(1)

        fastq1 = args.positional_args[0]
        if len(args.positional_args) == 3:
            fastq2 = args.positional_args[1]
            simulacion_folder = args.positional_args[2]

        if len(args.positional_args) == 2:
            print(f"\n Just one fastq {fastq1}, mode single")
            fastq2 = False
            simulacion_folder = args.positional_args[1]

        if args.slow:
            # SLOW SIMULATION MODE
            modo = "slow"
            ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, noindel_option)
        else:
            # FAST SIMULATION MODE
            modo = "fast"
            ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, noindel_option)

        if ok:
            sys.exit(0)
        else:
            print_error("Simulation pipeline failed.")
            sys.exit(1)

    # If this is NOT a simulation => interpret the first positional argument as SRR
    if len(args.positional_args) < 1:
        print_error("Missing SRR (or use --simulacion with <fastq1> <fastq2> <simulationFolder>).")
        sys.exit(1)

    srr = args.positional_args[0]

    # 6) Choose slow or fast pipeline
    pipeline_attempt = 1
    while pipeline_attempt <= MAX_PIPELINE_ATTEMPTS:
        if args.slow:
            modo = "slow"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)
        else:
            modo = "fast"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)

        if ok:
            print_step("Pipeline completed successfully.")
            sys.exit(0)
        else:
            print_error(f"Error in the pipeline. Retrying... (Attemp {pipeline_attempt}/{MAX_PIPELINE_ATTEMPTS})")
            pipeline_attempt += 1

    print_error("Pipeline failed after multiple attempts.")
    sys.exit(1)

if __name__ == "__main__":
    main()
