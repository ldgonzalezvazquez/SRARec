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
# Global Configuration
############################################
min_base_qual = 20
min_fraction_align_reads=0.9
############################
# Execution variables
#########################
MAX_PIPELINE_ATTEMPTS = 3
MAX_DOWNLOAD_ATTEMPTS = 3
VERBOSE = False


def _to_one_based(pos) -> int:
    """Convert internal 0-based coordinates to 1-based for results output."""
    return int(pos) + 1



# Do not store per-read start/end ranges (saves RAM/time)
SAVE_READS_MINMAX = True
if sys.getrecursionlimit() < 20000:
    sys.setrecursionlimit(20000)


############################################
# Auxiliary Functions
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

    {BOLD}{YELLOW}Program:{RESET} {CYAN}SRA{GREEN}Rec{RESET}
    {BOLD}{YELLOW}Authors:{RESET} {MAGENTA}Luis Daniel González-Vázquez & Darren P. Martin{RESET}
    {BOLD}{YELLOW}Email:{RESET} {BLUE}luisdaniel.gonzalez@uvigo.es{RESET}
    {BOLD}{YELLOW}Year:{RESET} {RED}2026{RESET}
    {BOLD}{YELLOW}Version:{RESET} {CYAN}1.{GREEN}3{RESET}
    """

    print(header)

def print_step(message):
    """Print a main step in white."""
    print(f"\033[1;37m{message}\033[0m")  # Bright white

def print_substep(message):
    """Print a substep (executed command) in blue."""
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
        # You can choose another color, e.g. magenta (35)
        print(f"\033[1;35m[VERBOSE] {message}\033[0m")

def identify_platform_and_layout(srr):
    """
    Identifies the sequencing platform and read type (Single End or Paired End) for a given SRR.
    Returns a dictionary with the platform name and layout.
    """
    try:
        cmd = f"esearch -db sra -query {srr} | efetch -format runinfo"
        ###print_verbose(f"Identifying platform with command: {cmd}")
        runinfo = subprocess.check_output(cmd, shell=True, text=True).strip()

        lines = runinfo.split('\n')
        if len(lines) < 2:
            raise ValueError("No information found in RunInfo.")

        headers = lines[0].split(',')
        data = lines[1].split(',')

        platform_index = headers.index("Platform") if "Platform" in headers else -1
        layout_index = headers.index("LibraryLayout") if "LibraryLayout" in headers else -1

        if platform_index == -1 or layout_index == -1:
            raise ValueError("Could not find 'Platform' or 'LibraryLayout' columns in RunInfo.")

        platform = data[platform_index].strip()
        layout = data[layout_index].strip()

        print(f"[INFO] Platform identified: {platform}")
        print(f"[INFO] Layout: {layout}")

        return {"platform": platform, "layout": layout}
    except Exception as e:
        print_error(f"[ERROR] Could not identify sequencing platform and layout: {e}")
        sys.exit(1)

def retry_command(max_attempts, cmd, shell=True):
    """
    Executes a command with retries on failure.
    """
    attempt = 1
    while attempt <= max_attempts:
        try:
            print_substep(f"Running: {cmd}")
            process = subprocess.Popen(
                cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = process.communicate() # waits for the process to finish and captures stdout and stderr
            print_output(stdout)
            if process.returncode != 0: # process.returncode returns the command exit code
                # Also show stderr for debugging
                print_output(f"[STDERR]: {stderr}") 
                raise subprocess.CalledProcessError(process.returncode, cmd) # subprocess is the module that runs the subprocess
            return True
        except subprocess.CalledProcessError:
            print_substep(f"Attempt {attempt} of {max_attempts} failed. Retrying...")
            attempt += 1
    print_error(f"Command failed after {max_attempts} attempts: {cmd}")
    return False

def filter_bam_by_alignment_fraction(input_bam, output_bam, min_fraction=0.7):
    """
    Filters the input BAM and writes a new BAM (output_bam) containing only
    reads whose fraction of aligned bases (query_alignment_length / query_length)
    is at least min_fraction (default 0.5 or 50%).

    Parameters:
      input_bam (str): Path to the input BAM.
      output_bam (str): Path to the filtered output BAM.
      min_fraction (float): Minimum mapping fraction to keep the read.
    """
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    for read in in_bam.fetch():
        # If the read is unmapped, skip it.
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        
        # read.query_alignment_length: length of the aligned portion (excludes insertions and clips)
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
    Runs bwa mem (or minimap2, depending on platform) to generate a SAM file
    and then converts it to BAM.
    """
    try:
        # Define the read group correctly
        read_group = f"@RG\\tID:{Path(fastq1).stem}\\tSM:{Path(fastq1).stem}\\tLB:Lib\\tPL:ILLUMINA" # Path(fastq1).stem uses pathlib to extract the base filename without extension

        # Adjust according to platform
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
            raise ValueError(f"Unsupported platform: {platform}") # raise manually raises an exception when an error prevents the program from continuing

        print_substep(f"Running: {map_cmd}")

        # Run mapping to generate the .sam file
        result = subprocess.run(map_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            # Also show stderr on failure
            print_output(f"[STDERR]: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, map_cmd)

        print(f"[INFO] SAM file generated: {output_sam}")

        # Convert the .sam file to .bam
        # Step 1: Convert SAM to BAM
        intermediate_bam = str(output_bam).replace(".bam", "_intermediate.bam")
        samtools_view_cmd = f"samtools view -q 30 -bS {output_sam} -o {intermediate_bam}"
        print_substep(f"Running: {samtools_view_cmd}")
        result_view = subprocess.run(samtools_view_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result_view.returncode != 0:
            print_output(f"[STDERR]: {result_view.stderr}")
            raise subprocess.CalledProcessError(result_view.returncode, samtools_view_cmd)

        # Step 2: Filter reads by alignment quality (read coverage against the reference)
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
        print_substep(f"Running: {samtools_sort_cmd}")
        result_sort = subprocess.run(samtools_sort_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        if result_sort.returncode != 0:
            print_output(f"[STDERR]: {result_sort.stderr}")
            raise subprocess.CalledProcessError(result_sort.returncode, samtools_sort_cmd)

        print_substep(f"[INFO] Sorted BAM file created: {output_bam}")

    except subprocess.CalledProcessError as e: # as e -> capture the exception in variable e to access its attributes
        print(f"[ERROR] Error executing: {e.cmd}\nExit code: {e.returncode}") # e.cmd -> the command that caused the error
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        sys.exit(1)

def _get_base_or_indel_at_pos_indel(read, pos, ref_fasta=None, ref_name=None, noindel=False):
    """
    Returns (event, quality) for reference position 'pos' (0-based).

    - If noindel=True, insertions and deletions are IGNORED:
        * Only check whether a base exists in the read (qPos is not None).
        * If so, return that base and its quality.
        * If not, or if there is a deletion, return (None, None).

    - If noindel=False, use the usual logic:
        * First search for an insertion "between pos-1 and pos" (same code as before).
        * Then search for base/deletion at pos.
        * Combine if both exist: e.g. "A+3ACT", "-1A+3ACT".
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
                    # Deletion present -> ignore it
                    pass
        # If we reach here, there was no base at the position -> return None
        return (None, None)

    # ----------------------------------------------------------------
    # NORMAL MODE => DETECT INSERTIONS + BASE/DELETION
    # ----------------------------------------------------------------
    insertion_str  = None
    insertion_qual = None

    # (1) Search for INSERTION between pos-1 and pos
    for i in range(len(aligned_pairs) - 1):
        q1, r1 = aligned_pairs[i]
        q2, r2 = aligned_pairs[i + 1]

        # "Between pos-1 and pos": r1 == pos and r2 is None...
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

            break  # Take only the first insertion we find

    # (2) Search for base or deletion at 'pos'
    base_str  = None
    base_qual = None

    for i, (qPos, rPos) in enumerate(aligned_pairs):
        if rPos == pos:
            if qPos is not None:
                # => base present
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
                base_qual = 0  # Quality 0 (convention for deletions)

            break

    # (3) Combine results according to what was found
    if not base_str and not insertion_str:
        return (None, None)

    if base_str and insertion_str:
        # E.g. "A+3ACT" or "-1A+3ACT"
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
    Returns True if 'deletion_current' is a continuation of 'deletion_prev'
    given that the position difference is delta.
    For example, if deletion_prev = "-5CTTCT" and delta = 1, deletion_current is expected to be "-4TTCT".
    There is no need to handle deletions with insertions afterward because these are single-base events and,
    even if repeated with two bases at the previous position, we need that for naming the insertion.
    """
    n_prev, seq_prev = parse_deletion(deletion_prev) # returns deletion size and sequence
    n_cur, seq_cur = parse_deletion(deletion_current)
    if n_prev is None or n_cur is None:
        return False
    if n_prev != n_cur + delta: # current deletion size plus delta must equal the previously annotated deletion size
        return False
    if seq_prev[delta:] == seq_cur: # returns the sequence from delta onward, i.e. what that sequence should be without the difference (the non-shared part)
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

            # Update min/max only if we want to store it
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

            # Filter by quality if it is not a deletion
            if not base_str.startswith('-') and base_qual < min_base_qual:
                continue

            base_reads[base_str].add(read.query_name)
            len_reads.append(len(read.query_sequence))

        # 2) Remove conflicting reads (those present in >1 variant)
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

        # 5) Store only if >= 2 variants
        if len(passing_bases) >= 2:
            results[pos] = passing_bases

    bam.close()

    # Convert read_minmax to a list (or empty if not stored)
    if SAVE_READS_MINMAX:
        reads_in_fin = [(st, en) for (st, en) in read_minmax.values()]
    else:
        reads_in_fin = []

    # Return the usual values
    return results, max(len_reads), reads_in_fin


# END SCAN_REGION_FOR_POLYMORPHISMS FUNCTION BLOCK
# START COMPARAR_POSICIONES_CONSECUTIVAS FUNCTION BLOCK:

def is_indel(variant):
    # Function to determine whether a variant is an indel (insertion or deletion)
    return variant.startswith('-') or '+' in variant
    
    # This function determines how to format a variant in "detail".
    #  - If it is not an insertion, return it unchanged.
    #  - If it is an insertion, check for "X+Nseq" where X is 1 nucleotide.
    #      and X appears as a standalone variant in variants_dict => remove X, leaving "+Nseq".

def format_variant(variants_dict, variant):
    # 1) If it is not an insertion, return it as-is
    #    (deletion => unchanged, base => unchanged).
    if not '+' in variant:
        return variant
        
    # In "A+2CT"...
    #   - "A" is the preceding part
    #   - check whether "A" exists as a standalone variant in variants_dict
    #   - if yes, return "+2CT"; otherwise return "A+2CT"
    # LOGIC:
    # Find the first '+' and split
    
    idx_plus = variant.find('+')   
    # base_previa is everything before '+'
    base_previa = variant[:idx_plus]
    # Remainder is the part including '+'
    insercion = variant[idx_plus:]  # e.g. "+2CT"

    # and it exists as a standalone variant => return only the insertion part
    if base_previa in variants_dict:
        # there is a standalone "A" variant => remove the preceding base
        return insercion
    else:
        # keep the preceding base
        return variant

def check_comparacion_custom(posA, polyA, posB, polyB):
    """
    Performs the 4-gamete test between posA and posB.

    Returns (result, detail), where:
      - result = True/False/None
      - detail:
         * If True:  "vA1.vA2,vB1.vB2 / vA1'.vA2',vB1'.vB2' / ..."
                    (for each quadruplet that formed 4 gametes)
         * If False:
             - if ALL intersections are indel-based => "posAIndels,posBIndels"
             - otherwise => ""
         * If None:  "" (no intersections)

    EXTRA:
      - For insertions like "A+2CT", if "A" exists alone at that position,
        it is not written ("+2CT" instead). If it does not exist alone, keep "A+2CT".
    """
    # 1) Coverage at each position (sum of all variant sets)
    coverageA = sum(len(rset) for rset in polyA.values())
    coverageB = sum(len(rset) for rset in polyB.values())
    coverage_min = min(coverageA, coverageB)    


    # Control variables
    any_intersection = False       # Track whether there was at least one intersection
    only_indels = True            # Remains True if ALL intersections have an indel in A or B
    posA_indels = set()           # Indels found at posA
    posB_indels = set()           # Indels found at posB

    # Store quadruplets that DO form the 4 gametes
    quadruplas_ok = []

    variantesA = list(polyA.keys())  # e.g. ["A", "G", "A+2CT", "-3C", ...]
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

            # Check intersections
            for (gA, gB) in gametos:
                readsA = polyA[gA]
                readsB = polyB[gB]
                inter = readsA.intersection(readsB)
                
                inter_count=len(inter)
                if (inter_count != 0) and (inter_count < min_inter): # Can always become more restrictive, but not less, otherwise False results would also change
                    min_inter=inter_count

                # Filter 1: minimum read threshold
                enough_reads = (inter_count >= reads_gametos_minoritarios)
                
                # Filter 2: frequency threshold (relative to coverage_min)
                # coverage_min may be 0; avoid division by zero
                if coverage_min > 0:
                    freq = inter_count / coverage_min
                else:
                    freq = 0.0

                if freq < min_freq:
                    min_freq = freq
                
                # Less restrictive frequency filter (e.g. >= 0.1 => 10%)
                enough_freq = (freq >= THRESHOLD_GAMETO_FREQ)

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

            # If all 4 gametes were satisfied
            if all(gametos.values()):
                # Build output "vA1.vA2,vB1.vB2"
                # formatting each one with format_variant
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
                    

    # Determine the result
    if quadruplas_ok:
        # => at least one quadruplet with 4 gametes
        detalle = "/".join(quadruplas_ok)
        return (True, detalle)
    else:
        # => 4 gametes were not formed
        if any_intersection:
            # => there were overlaps
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

                return (False, f"{n_gametos_true},{posA_str},{posB_str}")
            else:
                # => there was a base-base intersection => no detail
                return (False, f"{n_gametos_true}")
        else:
            # => no intersection => None
            return (None, "")

def comparar_posiciones_consecutivas(dic_pol_pos):
    """
    Compares consecutive position pairs in dic_pol_pos using
    check_comparacion_custom, which returns (result, detail).

    Parameters:
      dic_pol_pos (dict): {pos -> {variant: set(read_names), ...}}

    Returns:
      dict: {
          'True':  [(posA, posB, detail), ...],
          'False': [(posA, posB, detail), ...],
          'None':  [(posA, posB, detail), ...]
      }
      where:
        - detail is the string returned by check_comparacion_custom (e.g. "A.G,B.T/...").
    """

    comparaciones = {
        'True': [],
        'False': [],
        'None': []
    }

    # Get sorted positions
    sorted_positions = sorted(dic_pol_pos.keys())

    # Iterate over consecutive position pairs
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

# END COMPARAR_POSICIONES_CONSECUTIVAS FUNCTION BLOCK
# FUNCTION EXCLUSIVE TO FAST IMPLEMENTATION:
def extract_polymorphic_fastq(bam_path, ref_name, posiciones_a_analizar, output_fastq, min_reads=5, ref_fasta_path=None,noindel=False):
    """
    Extracts reads mapping to polymorphic positions within region [start, end)
    from the BAM file and generates a FASTQ file with those reads.

    Parameters:
      bam_path: Path to the BAM file.
      ref_name: Reference name (e.g. "NC_045512.2").
      start: Start position (0-based).
      end: End position (0-based, non-inclusive).
      output_fastq: Path to the resulting FASTQ file.
      min_reads: Minimum number of reads required to consider a base/indel.
      ref_fasta_path: (Optional) Path to the reference FASTA file, to extract the actual sequence
                      in case of deletions.
    """
    # Open the reference FASTA if a path is provided.
    ref_fasta = pysam.FastaFile(ref_fasta_path) if ref_fasta_path else None

    # Get polymorphic positions using the function already defined.
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
        print("[INFO] No reads of interest found in the specified region.")
        return

    # Save the read list in a temporary file (one per line)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp_list = tmp.name
        for rname in reads_of_interest:
            tmp.write(rname + "\n")
    
    # Build the command to extract reads of interest and generate a FASTQ.
    # The command does the following:
    #  1) samtools view -h -F 4: extract mapped reads (with header)
    #  2) awk: filter lines whose first field (read name) is in the temporary file.
    #  3) samtools fastq -: convert filtered SAM to FASTQ, writing to the output file.
    cmd = (
        f"samtools view -h -F 4 {bam_path} | "
        f"awk 'BEGIN{{while((getline line < \"{tmp_list}\")>0) r[line]=1}} "
        f"/^@/{{print;next}} ($1 in r)' | "
        f"samtools fastq - > {output_fastq}"
    )
    
    print(f"Running: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        print(f"[INFO] FASTQ file with reads of interest generated: {output_fastq}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Could not generate FASTQ file with reads of interest: {e}")
        sys.exit(1)
    finally:
        # Remove the temporary read-list file.
        os.remove(tmp_list)

# END FUNCTION EXCLUSIVE TO FAST IMPLEMENTATION

# START PILEUP USAGE FUNCTIONS

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
            # Character we do not handle -> ignore it
            i += 1

    bases_cumplen = [b for b,cnt in base_counts.items() if cnt >= MIN_READS]
    return (len(bases_cumplen) >= 2)


def parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS): # CHECK WHETHER THIS IS 0-BASED OR NOT? WITH THIS METHOD WE SELECT POSITION 701 AS THE FIRST; WITH THE OTHER WE SELECTED 700. -> NOW ADJUSTED TO 1-BASED -> I COUNTED POSITIONS ON THE REFERENCE AND 39 (1-BASED) WAS ALSO 39 IN THE PILEUP
    with open(pileup_file, 'r') as f:
        posiciones_polimorficas = []

        for line in f:
            fields = line.strip().split('\t')
            ####print_verbose(fields)
            if len(fields) < 6: # CHANGE THIS IF READ NAMES ARE NOT EXTRACTED LATER
                continue
            pos = int(fields[1])
            #print(f"Position: {pos}")
            bases = str(fields[4])
            if len(bases) < (MIN_READS*2):
                continue
            polimorfismo_si_o_no = quick_check_col5(bases,MIN_READS)
            if polimorfismo_si_o_no:
                ###print_verbose(f"Analyzed position: {pos}, analyzed bases: {bases}")
                posiciones_polimorficas.append(pos-1) # Adjust to 0-based for reading with pysam

    return posiciones_polimorficas

# END PILEUP USAGE FUNCTIONS
# START SEARCH ALGORITHM FUNCTIONS
def construir_distribuciones(dic_pol_pos):
    """
    Converts dic_pol_pos (pos -> {base -> set(reads), ...})
    into a dict pos -> { base -> percentage }.

    Verifies that no read appears in multiple variants at the same position,
    computes (#reads_variant / total_reads).
    Prints the result for debugging.
    """
    distributions_por_pos = {}

    print_substep("Building polymorphism distributions ...")
    for pos, variant_dict in tqdm(dic_pol_pos.items()):
        all_reads = set()
        for variant, rset in variant_dict.items():
            for r in rset:
                if r in all_reads: # UNCOMMENT THIS WHEN CODE ADAPTATION IS DONE AND TO INVESTIGATE WHY READS ARE REPEATED
                    raise ValueError(
                        f"Read '{r}' appears in multiple variants at position {pos}."
                    )
                all_reads.add(r)

        total_reads = sum(len(rset) for rset in variant_dict.values())
        if total_reads == 0:
            continue

        dist = {}
        for variant, rset in variant_dict.items():
            dist[variant] = len(rset) / total_reads

        distributions_por_pos[pos] = dist

    ###print_verbose(f"\n[DEBUG] Distributions per position: {distributions_por_pos}")
    return distributions_por_pos


def generar_todos_grupos_por_ventana(positions, max_read_length, min_size=2):
    """
    Generates ALL 'windows' (subsets) where
    (pos_final - pos_inicial) <= max_read_length.

    Then removes windows that are a strict subset of a larger window.
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
    # 1) Sort descending by size
    todos_grupos_ordenados = sorted(todos_grupos, key=len, reverse=True)

    # 2) Iterate large->small
    ventanas_finales = []
    for v in todos_grupos_ordenados:
        # If 'v' is a strict subset of any in 'ventanas_finales', REVIEW AND UNDERSTAND THIS CODE WELL -> IT WORKS, BUT NEEDS TO BE UNDERSTOOD
        # discard it
        # Otherwise, add it
        if not any((v < x) for x in ventanas_finales):
            # 'v < x' means v.issubset(x) and v != x
            ventanas_finales.append(v)

    # 3) Optionally, we could reorder 'ventanas_finales' 
    #    according to preference (by size asc/desc or 
    #    return to the order in which they were generated).
    #    Here we leave them in descending order by size:
    return ventanas_finales

##############################################################################
# 1.A - Ignore the base: compare fractions
##############################################################################

def agrupar_por_distribucion_individual(distributions_por_pos, threshold=0.05):
    """
    Groups positions based on the individual fractions of their distribution.

    For each position (e.g. 100, 101, etc.), for each base (e.g. 'A', 'G', etc.) and its percentage,
    try to add the position to an existing group whose average (for that sub-distribution)
    is within the threshold of the position fraction.

    Each group is represented as a dictionary with:
       - 'promedio': numeric value (float) that is the mean of fractions in the group.
       - 'posiciones': a set with the positions that contributed.

    When a new position is added, the average is updated:
       new_average = (average * N + current_percentage) / (N+1),
       where N is the number of positions already in the group.

    If the position does not fit any group for that fraction, create a new group.

    At the end, remove groups that are strictly equal to or subsets of others.

    Returns a list of groups.
    """
    grupos = []  # Each group is a dictionary: {'promedio': float, 'posiciones': set(), 'base': base}
    
    ###print_verbose(f"\n------------------------NEW WINDOW: Distributions per pos:\n {distributions_por_pos}\n")

    # Process positions in ascending order (although order may not matter)
    for pos in sorted(distributions_por_pos.keys()):
        ###print_verbose(f"\nPOSITION: {pos}")
        dist = list(set(distributions_por_pos[pos].values())) # Example [0.6, 0.4]
        ###print_verbose(f"DIST: {dist}")
        # For each base and its percentage at the position:
        for porcentaje in dist:
            ####print_verbose(f"Percentage: {porcentaje}, grupos len: {len(grupos)}")
            if len(grupos) > 0 and any(abs(porcentaje - grupo['promedio']) <= threshold for grupo in grupos):
                for grupo in grupos:
                    if abs(porcentaje - grupo['promedio']) <= threshold:
                        grupo['posiciones'].add(pos)
                        ###print_verbose(f"Previous mean: {grupo['promedio']}, adjusting with percentage: {porcentaje}")
                        N = len(grupo['posiciones'])
                        grupo['promedio'] = (grupo['promedio'] * (N - 1) + porcentaje) / N
                        ###print_verbose(f"New mean: {grupo['promedio']}")
            else:
                ###print_verbose("Group does not exist: Creating ...")
                nuevo_grupo = {'promedio': porcentaje, 'posiciones': {pos}}
                ###print_verbose(f"New group: {nuevo_grupo}")
                grupos.append(nuevo_grupo)               
    
    ###print_verbose(f"Unsorted groups: {grupos}")

    grupos = sorted(grupos, key=lambda grupo: len(grupo['posiciones']), reverse=True)

    ###print_verbose(f"Groups: {grupos}")

    grupos_finales = []
    for grupo in grupos:
        pos_set = grupo['posiciones']
        # If another group in grupos_finales has a larger position set that contains pos_set,
        # then discard this group.
        if not any(pos_set == otro['posiciones'] for otro in grupos_finales): # use <= to also remove subgroups
            grupos_finales.append(grupo)
    ###print_verbose(f"Final groups: {grupos_finales}") 

    return grupos_finales

##############################################################################
# REFINED LOGIC WITHIN GROUP
##############################################################################

def existe_subtrue_mas_pequeno(p_first, p_last, comp_global):
    """
    Returns True if within comp_global there exists a pair (pA, pB)
    such that:
      1) comp_global[(pA, pB)]["res"] == True
      2) pA >= p_first and pB <= p_last
      3) (pA, pB) != (p_first, p_last)
      4) anotar == True (optional, if you also require anotar=True)
    """
    for (pA, pB), info in comp_global.items():
        if info.get("res") is True and info.get("anotar", True) is True:
            if pA >= p_first and pB <= p_last and (pA, pB) != (p_first, p_last):
                return True
    return False

def desanotar_rangos_mas_grandes(p_first, p_last, comp_global):
    """
    Marks anotar=False for all ranges (pA, pB) in 'comp_global' that:
      - have res=True, anotar=True,
      - are distinct from (p_first, p_last),
      - contain (p_first, p_last), i.e. pA <= p_first and pB >= p_last.

    This is used to 'discard' larger ranges when we prefer to keep
    the new range (smaller or more specific).
    """
    for (pA, pB), info in comp_global.items():
        # Must be True and marked for annotation
        if info.get("res") is True and info.get("anotar", True) is True:
            # Must be a range distinct from (p_first, p_last)
            if (pA, pB) != (p_first, p_last):
                # And contain it (inclusive)
                if pA <= p_first and pB >= p_last:
                    ###print_verbose(f"   -> Marking anotar=False for large range {(pA, pB)}")
                    comp_global[(pA, pB)]["anotar"] = False

def refine_compare(pos_list, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global): # SEE IF THIS CAN BE A STANDALONE FUNCTION
    
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
        ###print_verbose(f"Already compared {key}: {res}, {detail} (intersected={info['intersected']})")
        # If already marked in groups (intersected=True), do not refine further
        if info["intersected"]:
            return
        else:
            # If it appears for the first time in groups, update the flag to True
            comparaciones_hechas_global[key]["intersected"] = True
            ###print_verbose(f"Updating flag to True for {key} (already done in intersection)")
    else:
        res, detail = check_comparacion_custom(p_first, dic_pol_pos[p_first], p_last, dic_pol_pos[p_last])
        # Here, in the refined-group context, the flag is set to True
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

        # Subdivide according to parity
        if n % 2 == 1: # odd
            mid_index = n // 2
            refine_compare(pos_list[:mid_index+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) # The +1 is correct because it is accessed like a range so the last element is not included
            refine_compare(pos_list[mid_index:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global) # Here mid_index refers to the first element.
        else:
            half = n // 2
            refine_compare(pos_list[:half+1], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)
            refine_compare(pos_list[half-1:], dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

        # ---------------------------------------------
        # NEW ANNOTATION LOGIC (WITHOUT before_count):
        # ---------------------------------------------
        # Store this result (res=True) in comparaciones_hechas_global,
        # but do not decide anotar until checking whether a smaller sub-True already exists:
        comparaciones_hechas_global[key] = {
            "res": True,
            "detail": detail,
            "intersected": True,
            # anotar will be defined depending on whether subTrue exists
            "anotar": True  # provisional value
        }

        # If a smaller sub-True already exists, set anotar=False
        if existe_subtrue_mas_pequeno(p_first, p_last, comparaciones_hechas_global):
            ###print_verbose(f"A more specific True already exists within {key}; annotate with anotar=False.")
            comparaciones_hechas_global[key]["anotar"] = False
        else:
            ###print_verbose(f"No smaller sub-True exists for {key}; annotate with anotar=True.")
            pass

        desanotar_rangos_mas_grandes(p_first, p_last, comparaciones_hechas_global)
        return

    else:  # res is None
        if len(pos_list) == 2:
            resultado_a_anotar=(p_first, p_last, None, detail)
            if resultado_a_anotar not in resultados:    
                resultados.append(resultado_a_anotar)
            return

        # Reduce p_last one by one (right to left)
        for idx in range(len(pos_list) - 2, 0, -1): # Each iteration returns idx one position earlier from the last, skipping the first; e.g. if len were 4, idx would be 3,2,1.
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
    ###print_verbose("\nENTERING comparar_grupos_interseccion\n")
    
    resultados = {'True': [], 'False': [], 'None': []}
    etiquetados = {}
    ###print_verbose(f"grupos: {grupos}")
    
    # Build a dictionary mapping each position to its label(s) (group averages it belongs to)
    for grupo in grupos:
        for p in grupo["posiciones"]:
            if p in etiquetados:
                etiquetados[p] = f"{etiquetados[p]},{grupo['promedio']}" # Concatenate average values found at this position; take existing averages and append the new one separated by comma
            else:
                etiquetados[p] = str(grupo["promedio"])
    
    ###print_verbose(f"etiquetados: {etiquetados}")
    posiciones_ordenadas = sorted(etiquetados.keys())
    ###print_verbose(f"posiciones_ordenadas: {posiciones_ordenadas}")
    
    # Iterate over all pairs (pA, pB) with pA < pB to avoid duplicates
    for i in range(len(posiciones_ordenadas)):
        for j in range(i + 1, len(posiciones_ordenadas)):
            pA = posiciones_ordenadas[i]
            pB = posiciones_ordenadas[j]
            # Convert labels to float lists for comparison
            gA = sorted(map(float, etiquetados[pA].split(',')))
            gB = sorted(map(float, etiquetados[pB].split(',')))
            
            ###print_verbose(f"Comparing: pA: {pA}, gA: {gA} | pB: {pB}, gB: {gB}")
            
            # Compare only if groups are different
            if gA == gB:
                ###print_verbose(f"Groups for {pA} and {pB} are equal; skipping comparison.")
                continue
            
            # Use the ordered tuple (minimum, maximum) to avoid duplicates
            pair = (min(pA, pB), max(pA, pB))
            if pair not in comparaciones_hechas_global:
                res, detail = check_comparacion_custom(pA, dic_pol_pos[pA], pB, dic_pol_pos[pB])
                ####print_verbose(f"Key={pair} – Res={res}")
                comparaciones_hechas_global[pair] = {"res": res, "detail": detail, "intersected": False, "anotar": True}
                if existe_subtrue_mas_pequeno(pair[0], pair[1], comparaciones_hechas_global):
                    ###print_verbose(f"A more specific True already exists within {pair}; annotate with anotar=False.")
                    comparaciones_hechas_global[pair]["anotar"] = False
                else:
                    ###print_verbose(f"No smaller sub-True exists for {pair}; annotate with anotar=True.")
                    pass

                desanotar_rangos_mas_grandes(pair[0], pair[1], comparaciones_hechas_global)
                ###print_verbose(f"Global comparisons made: {comparaciones_hechas_global}")
                polyA = dic_pol_pos[pA]
                polyB = dic_pol_pos[pB]
                ###print_verbose(f"Comparing position {pA} (reads: {polyA}) with {pB} (reads: {polyB})")
                res, detail = check_comparacion_custom(pA, polyA, pB, polyB)
                if res is True:
                    resultados['True'].append((pA, pB, detail))
                elif res is False:
                    resultados['False'].append((pA, pB, detail))
                else:
                    resultados['None'].append((pA, pB, detail))
            else:
                ###print_verbose(f"Pair {pair} was already compared; skipping.")
                pass
    
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
    4) For each window:
       - Filter positions => sub_dist
       - Group ignoring the base
       - Compare BETWEEN groups (adjacent)
       - Compare WITHIN each group with the refined function
    5) Return a dict with the info
    """

    dist_por_pos = construir_distribuciones(dic_pol_pos)

    ###print_verbose(f"dist_por_pos: {dist_por_pos}")

    ventanas = generar_todos_grupos_por_ventana(dist_por_pos.keys(), max_len, min_size)

    ###print_verbose(f"ventanas: {ventanas}")

    results_all_windows = {}
    comparaciones_hechas_global = {} 

    print_substep("Analyzing positions ...")
    for i, window in enumerate(ventanas):
        sub_dist = {pos: dist_por_pos[pos] for pos in window}
        grupos_subdist = agrupar_por_distribucion_individual(sub_dist, threshold)

        ###print_verbose(f"grupos_subdist: {grupos_subdist}")

        comp_entre = comparar_grupos_interseccion(grupos_subdist, check_comparacion_custom, dic_pol_pos, comparaciones_hechas_global)

        ###print_verbose(f"comp_entre: {comp_entre}")

        comp_en_grupos = {}
        for j, g in enumerate(grupos_subdist):
            ###print_verbose(f"G:{g}\n")
            # Use the refined version
            pos_sorted = sorted(g["posiciones"])
            comp_en_grupos[j] = refine_compare(pos_sorted, dic_pol_pos, check_comparacion_custom, comparaciones_hechas_global)

            ###print_verbose(f"comp_en_grupos[j]: {comp_en_grupos[j]}")

        results_all_windows[i] = {"window_positions": sorted(window),"grupos": grupos_subdist,"comparaciones_entre_grupos": comp_entre,"comparaciones_en_grupos": comp_en_grupos}

        ###print_verbose(f"results_all_windows[i]: {results_all_windows[i]}")

    return results_all_windows, comparaciones_hechas_global, dist_por_pos
# END SEARCH ALGORITHM FUNCTIONS

def run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False):
    """
    Performs read trimming. If 'filter_unpaired=True' and 'layout' is PAIRED,
    first runs repair.sh to discard/adjust orphan pairs,
    then normal trimming with cutadapt.
    """

    print_step("Step 4: Trimming reads...")

    trimmed_reads = work_dir / "reads_trimmed.fastq"
    read1_trimmed = work_dir / "reads_trimmed_1.fastq"
    read2_trimmed = work_dir / "reads_trimmed_2.fastq"
    read1_trimmed_temp = work_dir / "reads_trimmed_1_temp.fastq"
    read2_trimmed_temp = work_dir / "reads_trimmed_2_temp.fastq"
    trimmed_reads_temp = work_dir / "reads_trimmed.fastq.temp"

    # If output files already exist, skip the process
    if (read1_trimmed.exists() and read2_trimmed.exists()) or trimmed_reads.exists():
        print_substep("[INFO] Trimmed read files already exist. Skipping trimming.")
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
                # Repair output paths
                repaired_1 = work_dir / "reads_repaired_1.fastq"
                repaired_2 = work_dir / "reads_repaired_2.fastq"
                singletons = work_dir / "singletons.fastq"

                repair_cmd = (f"repair.sh in={fastq_1} in2={fastq_2} out={repaired_1} out2={repaired_2} outs={singletons}")
                print_substep(f"Running repair.sh (to discard orphan reads before cutadapt): {repair_cmd}")
                subprocess.run(repair_cmd, shell=True, check=True)

                # From here on use repaired_1 and repaired_2
                actual_fastq_1 = repaired_1
                actual_fastq_2 = repaired_2

            else:
                # If orphan reads do not need filtering/discarding, use the original fastq files
                actual_fastq_1 = fastq_1
                actual_fastq_2 = fastq_2

            # ------------------------------------------------------------------------
            # (2) Run cutadapt in paired-end mode
            # ------------------------------------------------------------------------
            cutadapt_cmd = (f"cutadapt -a file:{base_dir}/complementary_data/adapters.fa -A file:{base_dir}/complementary_data/adapters.fa -o {read1_trimmed_temp} -p {read2_trimmed_temp} {actual_fastq_1} {actual_fastq_2}")
            print_substep(f"Running cutadapt (paired-end): {cutadapt_cmd}")
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
            # SINGLE-END => no prior repair step
            cutadapt_cmd = (f"cutadapt -a file:{base_dir}/complementary_data/adapters.fa -o {trimmed_reads_temp} {fastq_1}")
            print_substep(f"Running cutadapt (single-end): {cutadapt_cmd}")
            if not retry_command(MAX_PIPELINE_ATTEMPTS, cutadapt_cmd):
                raise RuntimeError("Error during single-end trimming.")

            mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
            subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')

        print_substep("[INFO] Trimming completed successfully.")
        return True

    except Exception as e:
        print_error(f"[ERROR] Trimming failed: {e}")
        return False

def count_true_superranges_in(pA, pB, comp_global):
    """
    Returns how many ranges (xA, xB) with res=True (regardless of anotar)
    contain or match (pA, pB). That is: xA <= pA and xB >= pB.
    """
    count = 0
    for (xA, xB), info in comp_global.items():
        if info.get("res") is True:  # only matters that it is True
            # if xA <= pA and xB >= pB => (xA, xB) contains (pA, pB)
            if xA <= pA and xB >= pB:
                count += 1
    return count

def filtrar_rangos_mas_grandes(comparaciones_hechas_global):
    """
    Iterates over all pairs (pA, pB) with res == True and un-annotates (anotar=False)
    those "large" ranges that contain (pA' >= pA && pB' <= pB) some smaller
    subrange that is also True.
    """
    # 1) Collect all pairs where res == True (regardless of anotar)
    true_pairs = []
    for (pA, pB), info in comparaciones_hechas_global.items():
        if info["res"] is True:
            true_pairs.append((pA, pB))

    # 2) Sort them from largest to smallest
    #    (e.g. a range 100,200 has size 100; one of 100,105 has size 5).
    true_pairs.sort(key=lambda x: x[1] - x[0], reverse=True)

    # 3) For each large range, check whether it contains another smaller True range
    for (pA, pB) in true_pairs:
        # If already un-annotated in previous iterations, do nothing
        if not comparaciones_hechas_global[(pA, pB)]["anotar"]:
            continue

        # Search for a subrange (pA2, pB2) with pA2>=pA, pB2<=pB and res==True
        for (pA2, pB2) in true_pairs:
            # Avoid comparing with itself
            if (pA2, pB2) == (pA, pB):
                continue

            # Check whether it is a subrange
            if pA2 >= pA and pB2 <= pB:
                # And it is True
                if comparaciones_hechas_global[(pA2, pB2)]["res"] is True:
                    # => Mark the large range with anotar=False and exit
                    comparaciones_hechas_global[(pA, pB)]["anotar"] = False
                    break

# =========================
# Minimal reference preparation (only when --use_shiver_assembly)
# =========================
def download_reference_by_accession(accession, out_fasta):
    """
    Downloads a reference from NCBI using efetch.
    Requires Entrez Direct installed, as the script already uses esearch/efetch.
    """

    out_fasta = Path(out_fasta).expanduser().resolve()
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    if out_fasta.exists() and out_fasta.stat().st_size > 0:
        print_substep(f"[REF] Reference already exists: {out_fasta}")
        return out_fasta

    cmd = (
        f"efetch -db nucleotide -id {accession} -format fasta "
        f"> {out_fasta}"
    )

    print_substep(f"[REF] Downloading reference {accession}: {cmd}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

    if not out_fasta.exists() or out_fasta.stat().st_size == 0:
        raise RuntimeError(f"Reference download failed or produced empty file: {out_fasta}")

    return out_fasta


def prepare_reference_indices(ref_fa):
    """
    Creates the required indices:
      - BWA: .amb, .ann, .bwt, .pac, .sa
      - Samtools: .fai
      - Picard: .dict
    """

    ref_fa = Path(ref_fa).expanduser().resolve()

    if not ref_fa.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fa}")

    bwa_index_file = Path(str(ref_fa) + ".amb")
    faidx_file = Path(str(ref_fa) + ".fai")
    dict_file = ref_fa.parent / f"{ref_fa.stem}.dict"

    if bwa_index_file.exists():
        print_substep("[REF] BWA index already exists. Skipping.")
    else:
        subprocess.run(["bwa", "index", str(ref_fa)], check=True)
        print_substep("[REF] BWA index created.")

    if faidx_file.exists():
        print_substep("[REF] Samtools faidx already exists. Skipping.")
    else:
        subprocess.run(["samtools", "faidx", str(ref_fa)], check=True)
        print_substep("[REF] Samtools faidx created.")

    if dict_file.exists():
        print_substep("[REF] Picard dictionary already exists. Skipping.")
    else:
        subprocess.run(
            ["picard", "CreateSequenceDictionary", f"R={str(ref_fa)}", f"O={str(dict_file)}"],
            check=True
        )
        print_substep("[REF] Picard dictionary created.")

    return ref_fa


def create_custom_reference_index(custom_ref, base_dir, virus_name):
    """
    Creates a custom index in reference_indexes/<virus_name>/.

    Each time it runs:
      --index_only --custom_ref reference.fasta --virus <folder_name>

    the previous folder in reference_indexes/<virus_name>/ is deleted and a new one is created.
    """

    custom_ref = Path(custom_ref).expanduser().resolve()
    folder_name = str(virus_name).lower().replace(" ", "_")

    if not custom_ref.exists():
        raise FileNotFoundError(f"Custom reference FASTA not found: {custom_ref}")

    custom_index_dir = base_dir / "reference_indexes" / folder_name

    if custom_index_dir.exists():
        shutil.rmtree(custom_index_dir)

    custom_index_dir.mkdir(parents=True, exist_ok=True)

    custom_ref_indexed = custom_index_dir / "reference.fna"

    shutil.copy2(str(custom_ref), str(custom_ref_indexed))

    print_step(f"[REF] Creating custom reference index in: {custom_index_dir}")

    prepare_reference_indices(custom_ref_indexed)

    return custom_ref_indexed


def get_custom_reference_from_index(base_dir, virus_name):
    """
    Uses the custom index in reference_indexes/<virus_name>/ created with:
      --index_only --custom_ref reference.fasta --virus <folder_name>
    """

    folder_name = str(virus_name).lower().replace(" ", "_")
    custom_index_dir = base_dir / "reference_indexes" / folder_name
    custom_ref = custom_index_dir / "reference.fna"

    if not custom_ref.exists():
        raise FileNotFoundError(
            f"No custom reference index found at {custom_index_dir}.\n"
            "First run:\n"
            "python3 SRARec1.3.py --index_only --custom_ref /path/to/reference.fasta --virus <folder_name>"
        )

    required_files = [
        Path(str(custom_ref) + ".amb"),
        Path(str(custom_ref) + ".ann"),
        Path(str(custom_ref) + ".bwt"),
        Path(str(custom_ref) + ".pac"),
        Path(str(custom_ref) + ".sa"),
        Path(str(custom_ref) + ".fai"),
        custom_ref.parent / f"{custom_ref.stem}.dict",
    ]

    missing = [str(f) for f in required_files if not f.exists()]

    if missing:
        raise FileNotFoundError(
            "Custom reference exists, but some index files are missing:\n"
            + "\n".join(missing)
            + "\nRun again:\n"
            + "python3 SRARec1.3.py --index_only --custom_ref /path/to/reference.fasta --virus "
            + folder_name
        )

    print_substep(f"[REF] Using custom indexed reference: {custom_ref}")

    return custom_ref


def get_reference_by_virus(base_dir, virus):
    """
    Selects and prepares the reference.

    Default:
      - SARS-CoV-2: NC_045512.2

    Optional:
      - HIV-1 / VIH: NC_001802.1
      - any other name: uses the latest custom index created.
    """

    if virus is None:
        virus = "sars-cov-2"

    virus = str(virus).lower().replace("_", "-")

    references_dir = base_dir / "reference_indexes"

    if virus in ["sars-cov-2", "sarscov2", "covid", "covid-19"]:
        accession = "NC_045512.2"
        ref = references_dir / "sars_cov_2" / f"{accession}.fna"
        ref = download_reference_by_accession(accession, ref)
        ref = prepare_reference_indices(ref)
        return ref

    elif virus in ["hiv", "hiv-1", "hiv1", "vih", "vih-1"]:
        accession = "NC_001802.1"
        ref = references_dir / "hiv_1" / f"{accession}.fna"
        ref = download_reference_by_accession(accession, ref)
        ref = prepare_reference_indices(ref)
        return ref

    else:
        folder_name = str(virus).lower().replace(" ", "_")
        print_substep(
            f"[REF] Virus '{virus}' is not SARS-CoV-2 or HIV/VIH. "
            f"Using custom reference index in reference_indexes/{folder_name}/."
        )
        return get_custom_reference_from_index(base_dir, folder_name)
# =========================
# Helpers for SHIVER
# =========================
def _open_text_auto(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")

def _add_pair_suffix(in_path, out_path=None, pair=None):
    """
    Creates a FASTQ(.gz) copy adding /1 or /2 to the FIRST header token (@ID).
    Cleans any existing /1 or /2 before adding the correct suffix.
    """
    in_path = str(in_path)
    if pair not in ("1", "2"):
        raise ValueError("pair must be '1' or '2'")
    if out_path is None:
        if in_path.endswith(".fastq.gz"):
            out_path = in_path[:-9] + "_with12.fastq.gz"
        elif in_path.endswith(".fastq"):
            out_path = in_path[:-6] + "_with12.fastq"
        else:
            raise ValueError("Input must be .fastq or .fastq.gz")
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
    """Return the same path if it already ends in /1 or /2; otherwise create *_with12.fastq(.gz) and return that path."""
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
    spades,
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
    Runs the workflow:
      1) rnaviralspades.py -1 R1 -2 R2 -o HIV_rnaviral -t X -m Y
      2) shiver_init.sh MyInitDir config.sh MyRefAligment.fasta  MyAdapters.fasta MyPrimers.fasta
      3) shiver_align_contigs.sh MyInitDir config.sh contigs.fasta SRR
      4) Ensure /1 and /2 suffixes on trimmed FASTQ (creates *_with12.fastq if needed)
      5) shiver_map_reads.sh MyInitDir config.sh contigs.fasta SRR SRR.blast SRR_cut_wRefs.fasta R1 R2
    Returns the path to the resulting BAM: {work_dir}/{srr}.bam
    """
    # Normalize paths relative to base_dir (not CWD)
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
      # Run SPAdes in RNA viral mode:
    # spades.py --rnaviral -1 R1 -2 R2 -o OUT -t THREADS -m MEM
        cmd = [
            str(spades),
            "--rnaviral",
            "-1", str(fastq_r1_trimmed),
            "-2", str(fastq_r2_trimmed),
            "-o", str(out_spades),
            "-t", str(spades_threads),
            "-m", str(spades_mem_gb)
        ]

        print_substep(f"[SHIVER] Running SPAdes RNA viral: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    contigs = out_spades / "contigs.fasta"
    if not contigs.exists():
        raise FileNotFoundError(f"contigs.fasta not found in {out_spades}")

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

    # 4) Ensure /1 and /2 in trimmed FASTQ headers
    r1_for_shiver = _ensure_slash_suffix(fastq_r1_trimmed, "1")
    r2_for_shiver = _ensure_slash_suffix(fastq_r2_trimmed, "2")
    print_substep(f"[SHIVER] FASTQ for mapping: {r1_for_shiver}  {r2_for_shiver}")

    # 5) shiver_map_reads
    blast_out = work_dir / f"{srr}.blast"
    cut_refs = work_dir / f"{srr}_cut_wRefs.fasta"
    # === ALWAYS use the *_with12.fastq files ===
    r1_with12 = Path(str(fastq_r1_trimmed).replace(".fastq", "_with12.fastq"))
    r2_with12 = Path(str(fastq_r2_trimmed).replace(".fastq", "_with12.fastq"))

    # Create symlinks with specific names to avoid collisions in work_dir
    tmpdir = work_dir / "__shiver_inputs__"
    tmpdir.mkdir(exist_ok=True)
    r1_src = tmpdir / f"{srr}_1_to_shivermap.fastq"
    r2_src = tmpdir / f"{srr}_2_to_shivermap.fastq"
    for src, link in [(r1_with12, r1_src), (r2_with12, r2_src)]:
        if link.exists() or link.is_symlink():
            link.unlink()
        os.symlink(src, link)
    print_substep(f"[SHIVER] FASTQ for mapping (symlinks): {r1_src}  {r2_src}")

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
        # Some versions name the BAM differently; try to locate the most recent *.bam
        cands = sorted(work_dir.glob("*.bam"), key=lambda p: p.stat().st_mtime, reverse=True)
        if cands:
            return str(cands[0])
        raise FileNotFoundError("SHIVER output BAM not found in the working directory.")
    return str(bam_out)

def slow_comun(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, modo, noindel=False):
    # Step 4: Trimming
    # Try trimming first without filtering orphans:
    if layout == "PAIRED":
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)
    if layout == "SINGLE":
        trim_ok = run_trimming(fastq_single, fastq_2, layout, base_dir, work_dir, filter_unpaired=False)

    # If it fails and layout is PAIRED, retry with filter_unpaired=True
    if not trim_ok and layout == "PAIRED":
        print("[WARN] Normal trimming failed. Retrying with filter_unpaired=True ...")
        trim_ok = run_trimming(fastq_1, fastq_2, layout, base_dir, work_dir, filter_unpaired=True)

    if not trim_ok:
        print_error("[ERROR] Complete trimming failure, aborting pipeline.")
        sys.exit(1)

    # Step 5: Map trimmed reads
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
    "Fast" pipeline:
      1) Download SRA
      2) Convert to FASTQ
      3) Mapping [START FROM HERE]
      5) Sort and remove duplicates
      6) Select positions with polymorphisms [try to make a fast function]
      7) Generate fastq with reads of interest
      8) Trimming
      9) Mapping
      8) True/False analysis
    """

    # Step 4: Map untrimmed reads
    print_step("Step 4: Mapping untrimmed reads and sorting BAM...")
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
        print_substep("[INFO] BAM file already exists. Skipping mapping...")

    # Step 5: Remove duplicates
    print_step("Step 5: Removing duplicates and creating index...")

    # Remove duplicates
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

    # Creating an index
    dedup_bam_index = work_dir / f"{srr}_fast_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(dedup_bam))
        except Exception as e:
            print_error("[ERROR] Failed to create BAM index file:", e)
            return False
    else:
        print_substep("[INFO] BAM index file already exists. Skipping...")



    # Step 6: Generate FASTQ with reads of interest (with >1 polymorphism)
    print_step("Step 6: Extracting FASTQ with reads of interest (>1 polymorphism)...")
    poly_fastq = work_dir / "reads_poly.fastq"
    poly_fastq_temp = work_dir / "reads_poly_temp.fastq"

    if not poly_fastq.exists():
        # Open the deduplicated BAM to obtain the reference name:
        with pysam.AlignmentFile(str(dedup_bam), "rb") as bam_temp:
            ref_names = bam_temp.references

        pileup_file_1 = work_dir / f"{srr}_busqueda_inicial.pileup"
        print_step("Generating initial pileup file (first search pass)...")
        if pileup_file_1.exists():    
            print_substep("[INFO] Initial pileup file already exists. Continuing...")
        else:
            pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file_1}") # WE DO NOT STORE READ NAMES AS SEEN IN PREVIOUS VERSIONS (NO EXACT CORRESPONDENCE WITH BASES). To store them, add "--output-QNAME". The mismatch may be because quality is filtered only in certain columns (can be checked).
            if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
                print("[ERROR] Failed to generate initial pileup file.")
                return False

        posiciones = parse_pileup_and_extract_polimorfismos(pileup_file_1, 5)
        print(f"Initial pileup selection length: {len(posiciones)}")
        
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos_1,max_len_reads, reads_in_fin = scan_region_for_polymorphisms(str(dedup_bam), ref_names[0], posiciones, 5, ref_fa, min_base_qual, noindel)
        ref_fa.close()

        # Count how many polymorphic positions each read contributes
        read_poly_counts = {}
        for pos, allele_dict in dic_pol_pos_1.items():
            for allele, read_set in allele_dict.items():
                for rname in read_set:
                    read_poly_counts[rname] = read_poly_counts.get(rname, 0) + 1

        # Select reads appearing at more than one position (i.e. >1 polymorphism)
        selected_reads = {rname for rname, count in read_poly_counts.items() if count > 1}

        if not selected_reads:
            print("[INFO] No reads with >1 polymorphism found in the region.")
        else:
            # Save the read list in a temporary file (one per line)
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
            print_substep(f"Running: {cmd}")
            try:
                subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
                print(f"[INFO] FASTQ file with reads of interest generated: {poly_fastq}")
            except subprocess.CalledProcessError as e:
                ###print_verbose(f"[ERROR] Could not generate FASTQ file of reads of interest: {e}")
                return False
            finally:
                os.remove(tmp_list)
            mv_poly_fastq_cmd = f"mv {poly_fastq_temp} {poly_fastq}"
            subprocess.run(mv_poly_fastq_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] FASTQ file of selected reads already exists. Skipping...")


    # Step 7: Trimming
    print_step("Step 7: Trimming reads...")
    trimmed_reads = work_dir / "reads_trimmed.fastq"
    if trimmed_reads.exists():
        print_substep("[INFO] reads_trimmed.fastq already exists. Skipping trimming.")
    else:
        trimmed_reads_temp = work_dir / "reads_trimmed_temp.fastq"
        cutadapt_cmd = f"cutadapt -a file:{base_dir}/complementary_data/adapters.fa -o {trimmed_reads_temp} {poly_fastq}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cutadapt_cmd):
            print_error("Error during single-end read trimming.")
            return False
        mv_trim_cmd = f"mv {trimmed_reads_temp} {trimmed_reads}"
        subprocess.run(mv_trim_cmd, shell=True, check=True, executable='/bin/bash')    

    # Step 8: Map selected trimmed reads
    print_step("Step 8: Mapping selected trimmed reads and sorting BAM...")
    fast_bam_selected = work_dir / f"{srr}_fast_mapped_selected.bam"
    fast_bam_selected_temp = work_dir / f"{srr}_fast_mapped_selected_temp.bam"
    fast_sam_selected = work_dir / f"{srr}_fast_mapped_selected.sam"
    if not fast_bam_selected.exists():
        run_bwa_mem(refseq, trimmed_reads, None, fast_sam_selected, fast_bam_selected_temp, platform,min_fraction_align_reads)
        mv_fast_bam_selected_cmd = f"mv {fast_bam_selected_temp} {fast_bam_selected}"
        subprocess.run(mv_fast_bam_selected_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] BAM file already exists. Skipping mapping...")

    return pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, fast_bam_selected, modo, noindel)

def maybe_move_ref_and_csv(srr: str, work_dir: Path, base_dir: Path):
    """
    Moves {srr}_ref.fasta and {srr}_BaseFreqs_WithHXB2.csv from work_dir to:
      - base_dir/ref_fasta/
      - base_dir/csv_with_HXB2/
    Only acts if work_dir/{srr}_ref.fasta exists (i.e. RAVI flow after assembly).
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
            print_verbose("[MOVE] No assembly files to move (skipping).")
    except Exception as e:
        print_error(f"[MOVE] Error moving files before VCF generation: {e}")

def pasos_comunes_slow_fast(srr, base_dir, refseq, work_dir, platform, layout, fastq_1, fastq_2, fastq_single, bam, modo, noindel=False):

    print_step("\nEntering final steps ... (shared by slow and fast pipelines)\n")

    # Step 6: Remove duplicates
    print_step("Final step 1: Removing duplicates and creating index...")

    # Remove duplicates
    dedup_bam = work_dir / f"{srr}_slow_dedup.bam"
    dedup_bam_temp = work_dir / f"{srr}_slow_dedup_temp.bam"

    if not dedup_bam.exists():
        dedup_cmd = f"samtools markdup -r {bam} {dedup_bam_temp}"
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, dedup_cmd):
            print("Failed to remove duplicates.")
            return False
        mv_deputbam_cmd = f"mv {dedup_bam_temp} {dedup_bam}"
        subprocess.run(mv_deputbam_cmd, shell=True, check=True, executable='/bin/bash')
    else:
        print_substep("[INFO] Deduplicated file already exists. Skipping...")

    # Creating an index
    dedup_bam_index = work_dir / f"{srr}_slow_dedup.bam.bai"
    if not dedup_bam_index.exists():
        try:
            print(f"Estoy aquii: {dedup_bam}")
            pysam.index(str(dedup_bam))
        except Exception as e:
            print("Failed to create deduplicated BAM index:", e)
            return False
    else:
        print_substep("[INFO] Deduplicated BAM index already exists. Skipping...")

    # Step 7: Reading positions

    ############################################################################
    # TESTS, PILEUP CREATION PER FILE (RECOVERED FROM PREVIOUS VERSION) #
    ############################################################################

    pileup_file = work_dir / f"{srr}_slow.pileup"
    print_step("Final step 2: Generating initial pileup file...")
    if pileup_file.exists():    
        print_substep("[INFO] Pileup file exists. Continuing...")
    else:
        pileup_cmd = (f"samtools mpileup -B -Q {min_base_qual} -f {refseq} {dedup_bam} > {pileup_file}") # WE DO NOT STORE READ NAMES AS SEEN IN PREVIOUS VERSIONS (NO EXACT CORRESPONDENCE WITH BASES). To store them, add "--output-QNAME". The mismatch may be because quality is filtered only in certain columns (can be checked).
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
            print("Failed to generate initial pileup file.")
            return False

    # Define parameter combinations to evaluate:
    from itertools import product
    all_param_combinations = list(product(MIN_READS_LIST, READS_GAMETOS_MINORITARIOS_LIST, THRESHOLD_GAMETO_FREQ_LIST))

    vcf_generado = False
    coinfection_detected = False
    for (mr, rgm, freq_val) in tqdm(all_param_combinations):
        # ================================
        # Repeat STEPS 7 and 8 with these parameters
        # ================================

        print_substep(
            f"\n[ANALYSIS] Parameters => MIN_READS={mr}, "
            f"reads_gametos_minoritarios={rgm}, THRESHOLD_GAMETO_FREQ={freq_val}"
        )

        global MIN_READS # these parameters are changed globally
        global reads_gametos_minoritarios
        global THRESHOLD_GAMETO_FREQ

        MIN_READS = mr
        reads_gametos_minoritarios = rgm
        THRESHOLD_GAMETO_FREQ = freq_val

        posiciones_pileup = parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS)
        print_substep(f"Length of selection using pileup is: {len(posiciones_pileup)}")

        # If there are no positions, write the line with empty fields and move to the next iteration
        if len(posiciones_pileup) < 2:
            result_filename = base_dir / f"results_{modo}.txt"
            with open(result_filename, "a") as rf:
                rf.write(f"Sra:{srr}\tCoinf_detected:False\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:\tNo_Rec:\tBegFin_Reads:\n")
            continue  # Skip the rest of this iteration

        ###print_verbose(posiciones_pileup)

        ############################################################################
        # TESTS, PILEUP CREATION PER FILE (RECOVERED FROM PREVIOUS VERSION) #
        ############################################################################

        print_step("Final step 3: Reading mapped positions...")   
        bam = pysam.AlignmentFile(dedup_bam, "rb") # Open the BAM file to extract the reference name
        ref_names = bam.references  # Get a tuple of reference names
        bam.close()
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos,max_len_reads, reads_in_fin=scan_region_for_polymorphisms(dedup_bam, ref_names[0], posiciones_pileup, MIN_READS, ref_fa, min_base_qual, noindel)
        print_substep(f"Maximum read length is: {max_len_reads}")
        # Create a PrettyPrinter instance with a larger line width if needed
        ###print_verbose(dic_pol_pos)
        ref_fa.close()

        # Step 8: Analysis (True/False)
        print_step("Final step 4: Analyzing consecutive positions...")
        results_file = base_dir / "results.txt"
        results,comparaciones_hechas_global,dist_por_pos = pipeline_agrupacion_por_distribuciones(dic_pol_pos,max_len=max_len_reads,threshold=0.1,min_size=2) # this threshold is for grouping distinct groups, not for gamete filtering; IT IS DIFFERENT
        filtrar_rangos_mas_grandes(comparaciones_hechas_global)
        # result=comparar_posiciones_consecutivas(dic_pol_pos)
        ###print_verbose(f"\n\n[DEBUG] comparaciones_hechas_global:\n {comparaciones_hechas_global}")
        coinfection_detected = False
        if not vcf_generado:
            for pos, freq_dict in dist_por_pos.items():
                # freq_dict is e.g. {'A':0.3, 'T':0.7}
                count_over_01 = sum(1 for f in freq_dict.values() if f > 0.1)
                if count_over_01 >= 2:
                    coinfection_detected = True
                    break
        # Filter True / False pairs, but only if "anotar" == True
        true_list = []
        false_list = []

        for (pA, pB), info in comparaciones_hechas_global.items():
            if info["anotar"]:  # Store only if anotar == True
                if info["res"] is True:
                    # (1) Count how many True ranges (anotar True OR False) 
                    #     contain (pA, pB)
                    bigger_count = count_true_superranges_in(pA, pB, comparaciones_hechas_global)
                    
                    # (2) Store in the tuple
                    #     If you previously stored (pA, pB, detail),
                    #     now add bigger_count as an extra field:
                    true_list.append((pA, pB, info["detail"], bigger_count))

                elif info["res"] is False:
                    false_list.append((pA, pB, info["detail"]))

        # Format True and False pairs
        true_pairs = ";".join(
            [f"{_to_one_based(pA)},{_to_one_based(pB)},{bigger_count},{detail}"
             for (pA, pB, detail, bigger_count) in true_list]
        )
        false_pairs = ";".join(
            [f"{_to_one_based(a)},{_to_one_based(b)},{detail}" for (a, b, detail) in false_list]
        )

        # (1) Format reads_in_fin as "110,300;301,450;..."
        reads_in_fin_str = ";".join(
            f"{_to_one_based(st)},{_to_one_based(en)}" for (st, en) in reads_in_fin
        )

        # Name the results file (version per your example)
        result_filename = base_dir / f"results_{modo}.txt"

        # Open or create the file in append mode and write
        with open(result_filename, "a") as rf:
            rf.write(f"Sra:{srr}\tCoinf_detected:{coinfection_detected}\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\tBegFin_Reads:{reads_in_fin_str}\n")

    if coinfection_detected and not vcf_generado:
        # 1) Create the folder if it does not exist
        coinfection_dir = base_dir / "vcf_coinfecciones"
        coinfection_dir.mkdir(parents=True, exist_ok=True) # IMPORTANT: UNCOMMENT THIS IF VCF GENERATION IS DESIRED

        # 2) Move ref/csv ONLY if coming from assembly (detected by presence of {srr}_ref.fasta)
        maybe_move_ref_and_csv(srr, work_dir, base_dir)

        # 3) VCF filename
        vcf_path = coinfection_dir / f"{srr}.vcf"

        # 3) Command to generate the VCF -> IMPORTANT: UNCOMMENT THE SUBPROCESS IF VCF GENERATION IS DESIRED
        #    Example with bcftools (requires bcftools installed):
        #    samtools mpileup -B -f ref.fasta dedup.bam | bcftools call -mv -Ov -o out.vcf
        cmd_vcf = (f"bcftools mpileup -d 100000 -f {refseq} {dedup_bam} | bcftools call -mv -Ov  --ploidy 1 -o {vcf_path}")
        print_substep(f"Running: {cmd_vcf}")
        subprocess.run(cmd_vcf, shell=True, check=True)

        print_substep(f"VCF generated for potential coinfection screening: {vcf_path}")
        vcf_generado = True
        
    #print_step("Pipeline completed successfully.")
    return True

def pasos_ravi(srr, base_dir, refseq, work_dir, platform, dedup_bam, modo, noindel=False):

    print(f"Entered Bam mode with noindel={noindel}")

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
    print_substep(f"Running: {samtools_sort_cmd}")
    result_sort = subprocess.run(samtools_sort_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if result_sort.returncode != 0:
        print_output(f"[STDERR]: {result_sort.stderr}")
        raise subprocess.CalledProcessError(result_sort.returncode, samtools_sort_cmd)
    rm_intermediate_bam_cmd_3=f"rm {intermediate_bam_0}"
    subprocess.run(rm_intermediate_bam_cmd_3, shell=True, stderr=subprocess.PIPE, text=True)

    print_substep(f"[INFO] Sorted BAM file created: {output_bam}")


    # Creating an index
    dedup_bam_index = work_dir / f"{srr}_output.bam.bai"
    if not dedup_bam_index.exists():
        try:
            pysam.index(str(output_bam))
        except Exception as e:
            print("Failed to create deduplicated BAM index:", e)
            return False
    else:
        print_substep("[INFO] Deduplicated BAM index already exists. Skipping...")

    # Step 7: Reading positions

    ############################################################################
    # TESTS, PILEUP CREATION PER FILE (RECOVERED FROM PREVIOUS VERSION) #
    ############################################################################

    pileup_file = work_dir / f"{srr}_ravi.pileup"
    output_bam_dir = work_dir / f"{output_bam}"
    print_step("Final step 2: Generating initial pileup file...")
    if pileup_file.exists():    
        print_substep("[INFO] Pileup file exists. Continuing...")
    else:
        pileup_cmd = (f"samtools mpileup -A -B -Q {min_base_qual} -f {refseq} {output_bam_dir} > {pileup_file}") # WE DO NOT STORE READ NAMES AS SEEN IN PREVIOUS VERSIONS (NO EXACT CORRESPONDENCE WITH BASES). To store them, add "--output-QNAME". The mismatch may be because quality is filtered only in certain columns (can be checked).
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, pileup_cmd):
            print("Failed to generate initial pileup file.")
            return False

    # Define parameter combinations to evaluate:
    from itertools import product
    all_param_combinations = list(product(MIN_READS_LIST, READS_GAMETOS_MINORITARIOS_LIST, THRESHOLD_GAMETO_FREQ_LIST))

    vcf_generado = False
    coinfection_detected = False
    for (mr, rgm, freq_val) in tqdm(all_param_combinations):
        # ================================
        # Repeat STEPS 7 and 8 with these parameters
        # ================================

        print_substep(
            f"\n[ANALYSIS] Parameters => MIN_READS={mr}, "
            f"reads_gametos_minoritarios={rgm}, THRESHOLD_GAMETO_FREQ={freq_val}"
        )

        global MIN_READS # these parameters are changed globally
        global reads_gametos_minoritarios
        global THRESHOLD_GAMETO_FREQ

        MIN_READS = mr
        reads_gametos_minoritarios = rgm
        THRESHOLD_GAMETO_FREQ = freq_val

        posiciones_pileup = parse_pileup_and_extract_polimorfismos(pileup_file, MIN_READS)
        print_substep(f"Length of selection using pileup is: {len(posiciones_pileup)}")

        # If there are no positions, write the line with empty fields and move to the next iteration
        if len(posiciones_pileup) < 2:
            result_filename = base_dir / f"results_{modo}.txt"
            with open(result_filename, "a") as rf:
                rf.write(f"Sra:{srr}\tCoinf_detected:False\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:\tNo_Rec:\tBegFin_Reads:\n")
            continue  # Skip the rest of this iteration

        ###print_verbose(posiciones_pileup)

        ############################################################################
        # TESTS, PILEUP CREATION PER FILE (RECOVERED FROM PREVIOUS VERSION) #
        ############################################################################

        print_step("Final step 3: Reading mapped positions...")   
        bam = pysam.AlignmentFile(output_bam, "rb") # Open the BAM file to extract the reference name
        ref_names = bam.references  # Get a tuple of reference names
        bam.close()
        ref_fa = pysam.FastaFile(str(refseq))
        dic_pol_pos,max_len_reads, reads_in_fin=scan_region_for_polymorphisms(output_bam, ref_names[0], posiciones_pileup, MIN_READS, ref_fa, min_base_qual, noindel)
        print_substep(f"Maximum read length is: {max_len_reads}")
        # Create a PrettyPrinter instance with a larger line width if needed
        ###print_verbose(dic_pol_pos)
        ref_fa.close()

        # Step 8: Analysis (True/False)
        print_step("Final step 4: Analyzing consecutive positions...")
        results_file = base_dir / "results.txt"
        results,comparaciones_hechas_global,dist_por_pos = pipeline_agrupacion_por_distribuciones(dic_pol_pos,max_len=max_len_reads,threshold=0.1,min_size=2) # this threshold is for grouping distinct groups, not for gamete filtering; IT IS DIFFERENT
        filtrar_rangos_mas_grandes(comparaciones_hechas_global)
        # result=comparar_posiciones_consecutivas(dic_pol_pos)
        ###print_verbose(f"\n\n[DEBUG] comparaciones_hechas_global:\n {comparaciones_hechas_global}")
        coinfection_detected = False
        if not vcf_generado:
            for pos, freq_dict in dist_por_pos.items():
                # freq_dict is e.g. {'A':0.3, 'T':0.7}
                count_over_01 = sum(1 for f in freq_dict.values() if f > 0.1)
                if count_over_01 >= 2:
                    coinfection_detected = True
                    break
        # Filter True / False pairs, but only if "anotar" == True
        true_list = []
        false_list = []

        for (pA, pB), info in comparaciones_hechas_global.items():
            if info["anotar"]:  # Store only if anotar == True
                if info["res"] is True:
                    # (1) Count how many True ranges (anotar True OR False) 
                    #     contain (pA, pB)
                    bigger_count = count_true_superranges_in(pA, pB, comparaciones_hechas_global)
                    
                    # (2) Store in the tuple
                    #     If you previously stored (pA, pB, detail),
                    #     now add bigger_count as an extra field:
                    true_list.append((pA, pB, info["detail"], bigger_count))

                elif info["res"] is False:
                    false_list.append((pA, pB, info["detail"]))

        # Format True and False pairs
        true_pairs = ";".join(
            [f"{_to_one_based(pA)},{_to_one_based(pB)},{bigger_count},{detail}"
             for (pA, pB, detail, bigger_count) in true_list]
        )
        false_pairs = ";".join(
            [f"{_to_one_based(a)},{_to_one_based(b)},{detail}" for (a, b, detail) in false_list]
        )

        # (1) Format reads_in_fin as "110,300;301,450;..."
        reads_in_fin_str = ";".join(
            f"{_to_one_based(st)},{_to_one_based(en)}" for (st, en) in reads_in_fin
        )

        # Name the results file (version per your example)
        result_filename = base_dir / f"results_{modo}.txt"

        # Open or create the file in append mode and write
        with open(result_filename, "a") as rf:
            rf.write(f"Sra:{srr}\tCoinf_detected:{coinfection_detected}\tMin_read_num:{mr}\tMin_gameto_num:{rgm}\tMin_gameto_freq:{freq_val}\tPlatform:{platform}\tRec:{true_pairs}\tNo_Rec:{false_pairs}\tBegFin_Reads:{reads_in_fin_str}\n")

    if coinfection_detected and not vcf_generado:
        # 1) Create the folder if it does not exist
        coinfection_dir = base_dir / "vcf_coinfecciones"
        coinfection_dir.mkdir(parents=True, exist_ok=True) # IMPORTANT: UNCOMMENT THIS IF VCF GENERATION IS DESIRED

        # 2) Move ref/csv ONLY if coming from assembly (detected by presence of {srr}_ref.fasta)
        maybe_move_ref_and_csv(srr, work_dir, base_dir)

        # 3) VCF filename
        vcf_path = coinfection_dir / f"{srr}.vcf"

        # 3) Command to generate the VCF -> IMPORTANT: UNCOMMENT THE SUBPROCESS IF VCF GENERATION IS DESIRED
        #    Example with bcftools (requires bcftools installed):
        #    samtools mpileup -B -f ref.fasta dedup.bam | bcftools call -mv -Ov -o out.vcf
        cmd_vcf = (f"bcftools mpileup -d 100000 -f {refseq} {output_bam} | bcftools call -mv -Ov  --ploidy 1 -o {vcf_path}")
        print_substep(f"Running: {cmd_vcf}")
        subprocess.run(cmd_vcf, shell=True, check=True)

        print_substep(f"VCF generated for potential coinfection screening: {vcf_path}")
        vcf_generado = True
        
    #print_step("Pipeline completed successfully.")
    return True


def inicio_comun(srr, base_dir, refseq, modo, noindel=False):
    """
    "Slow" pipeline:
      1) Download SRA
      2) Convert to FASTQ
    """
    if modo == "slow":
        print("[INFO] Iniciando pipeline SLOW...")
    if modo == "fast":
        print("[INFO] Iniciando pipeline FAST...")

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
                break  # Information retrieved successfully
        except Exception as e:
            print_substep(f"Attempt {attempt1} failed to retrieve information for {srr}: {e}")
        attempt1 += 1
        time.sleep(random.randint(1, 20))

    if not result_platform_and_layout:
        print("Failed to retrieve SRA file information after multiple attempts.")
        sys.exit(1)
        return False

    platform = result_platform_and_layout["platform"]
    layout = result_platform_and_layout["layout"]

    # Step 2: Download SRA (if it does not exist)
    print_step("Step 2: Downloading data...")
    sra_file = work_dir / f"{srr}.sra"
    if not sra_file.exists():
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, f"prefetch --output-directory {base_dir} {srr}"):
            print("Failed to download SRA file.")
            return False
    else:
        print_substep("[INFO] SRA file already exists. Skipping download.")

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
                print_substep("[INFO] FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')


        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    elif platform in ["PACBIO","PACBIO_SMRT"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ _1 and _2 files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')

        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
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

    print_step("Step 1: Identifying platform...")
    attempt1 = 1
    result_platform_and_layout = None

    while attempt1 <= 10:
        try:
            result_platform_and_layout = identify_platform_and_layout(srr)
            if result_platform_and_layout:
                break  # Information retrieved successfully
        except Exception as e:
            print_substep(f"Attempt {attempt1} failed to retrieve information for {srr}: {e}")
        attempt1 += 1
        time.sleep(random.randint(1, 20))

    if not result_platform_and_layout:
        print("Failed to retrieve SRA file information after multiple attempts.")
        sys.exit(1)
        return False

    platform = result_platform_and_layout["platform"]
    layout = result_platform_and_layout["layout"]

    # Step 2: Download SRA (if it does not exist)
    print_step("Step 2: Downloading data...")
    sra_file = work_dir / f"{srr}.sra"
    if not sra_file.exists():
        if not retry_command(MAX_DOWNLOAD_ATTEMPTS, f"prefetch --output-directory {base_dir} {srr}"):
            print("Failed to download SRA file.")
            return False
    else:
        print_substep("[INFO] SRA file already exists. Skipping download.")

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
                print_substep("[INFO] FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')


        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fasterq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')

    elif platform in ["PACBIO","PACBIO_SMRT"]:
        if layout == "PAIRED":
            if fastq_1.exists() and fastq_2.exists():
                print_substep("[INFO] FASTQ _1 and _2 files already exist. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} --split-files -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd_1, shell=True, check=True, executable='/bin/bash')
                subprocess.run(mv_fastq_cmd_2, shell=True, check=True, executable='/bin/bash')

        else:  # SINGLE
            if fastq_single.exists():
                print_substep("[INFO] Single-end FASTQ file already exists. Skipping conversion.")
            else:
                cmd_fastq = f"fastq-dump {sra_file} -O {work_dir}"
                if not retry_command(MAX_DOWNLOAD_ATTEMPTS, cmd_fastq):
                    print("Failed to convert SRA to FASTQ.")
                    return False
                subprocess.run(mv_fastq_cmd, shell=True, check=True, executable='/bin/bash')
    return fastq_1, fastq_2

def inicio_comun_simulacion(fastq_1, fastq_2, simulacion_folder, base_dir, refseq, modo, platform_arg="ILLUMINA", noindel=False):
    """
    Local FASTQ pipeline entry point (slow or fast):
      - Uses existing FASTQ file(s) directly.
      - Skips SRA download and conversion.
      - Platform is taken from platform_arg; layout is PAIRED or SINGLE.
      - Output file names include the output folder name.
    """
    # Define the output folder name for local FASTQ mode
    sim_dir = simulacion_folder  # This name can be changed as needed
    work_dir = base_dir / simulacion_folder
    os.chdir(work_dir)

    print(f"\nLocal FASTQ mode: using provided FASTQ files\n")
    if modo == "slow":
        print("[INFO] Starting SLOW pipeline (local FASTQ)...")
    if modo == "fast":
        print("[INFO] Starting FAST pipeline (local FASTQ)...")

    
    # Get the output folder name for use in output filenames
    sim_dir_name = work_dir.name

    # Force platform and layout
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
# 3) ARGUMENT PARSING
##########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="SRARec pipeline for processing SRA accessions, local FASTQ files, "
            "or pre-existing BAM files to detect recombination/coinfection signals."
    )
    # Threshold parameters (comma-separated lists)
    parser.add_argument(
        "--MIN_READS_LIST",
        default="10",
        help="Comma-separated list of minimum read counts. Example: 10,20"
    )
    parser.add_argument(
        "--READS_GAMETOS_MINORITARIOS_LIST",
        default="5",
        help="Comma-separated list of minimum read counts required for minor gametes. Example: 5,10"
    )
    parser.add_argument(
        "--THRESHOLD_GAMETO_FREQ_LIST",
        default="0.05",
        help="Comma-separated list of minimum gamete-frequency thresholds. Example: 0.05,0.10"
    )

    # Mode flags
    parser.add_argument(
        "--slow",
        action="store_true",
        help="Run the slow pipeline instead of the default fast pipeline."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed verbose output."
    )
    parser.add_argument(
        "--noindel",
        action="store_true",
        help="Ignore insertions and deletions when scanning polymorphic positions."
    )

    # Flag for local FASTQ mode
    parser.add_argument(
        "--fastq",
        action="store_true",
        help="Run local FASTQ mode. Expected positional arguments: "
            "<fastq1> [fastq2] <output_folder>."
    )

    # Positional arguments (flexible usage)
    parser.add_argument(
        "positional_args",
        nargs="*",
        help="Positional arguments. Normal mode: <SRR>. "
            "Local FASTQ mode: <fastq1> [fastq2] <output_folder>. "
            "SHIVER assembly mode: either <fastq1> <fastq2> <output_folder> or a single SRR/ERR/DRR accession."
    )

    parser.add_argument(
    "--virus",
    type=str,
    default=None,
    help=   "Virus/reference to use. Default: sars-cov-2. "
            "Use hiv/vih for HIV-1 NC_001802.1. "
            "Any other value uses reference_indexes/<virus_name>/ created with --index_only."
    )
    parser.add_argument(
    "--custom_ref",
    type=str,
    default=None,
    help="Custom FASTA reference. Used only with --index_only and --virus."
    )

    parser.add_argument(
    "--index_only",
    action="store_true",
    help="Only create/overwrite the custom reference index and exit."
    )

    parser.add_argument(
    "--platform",
    type=str,
    default="ILLUMINA",
    help=   "Sequencing platform label written to results. "
            "Required for mapping in --fastq mode. "
            "In --bam mode it is metadata only (default: ILLUMINA). "
            "In SRA mode the platform is detected automatically from RunInfo."
    )

    parser.add_argument("--bam", help="Input BAM file mode. The BAM is expected to be available in the selected working folder.")

    # ===== SPAdes + SHIVER integration =====
    parser.add_argument(
        "--use_shiver_assembly",
        action="store_true",
        help="Use rnaviral plus SHIVER to generate a BAM file from trimmed paired-end FASTQ files."
    )
    parser.add_argument(
        "--shiver_dir",
        default="$CONDA_PREFIX",
        help="Path to the SHIVER installation directory containing the bin/ folder."
    )
    parser.add_argument(
        "--spades",
        default="spades.py",
        help="SPAdes command or path. By default, spades.py is used from PATH. "
            "RNA viral mode is run internally as: spades.py --rnaviral."
    )
    parser.add_argument(
        "--spades_conda_env",
        default="none",
        help="Conda environment name for running SPAdes. Use 'none' to disable this option."
    )
    parser.add_argument(
        "--shiver_config",
        default=None,
        help="Path to the SHIVER config.sh file. If omitted, <shiver_dir>/bin/config.sh is used."
    )
    parser.add_argument(
        "--shiver_threads",
        type=int,
        default=6,
        help="Number of threads for SHIVER read mapping."
    )
    parser.add_argument(
        "--spades_threads",
        type=int,
        default=8,
        help="Number of threads for rnaviral (spades)."
    )
    parser.add_argument(
        "--spades_mem_gb",
        type=int,
        default=12,
        help="Memory in GB for rnaviral (spades)."
    )
    parser.add_argument(
        "--shiver_ref_alignment",
        default="MyRefAligment.fasta",
        help="Reference alignment FASTA file used by shiver_init."
    )
    parser.add_argument(
        "--shiver_adapters",
        default="MyAdapters.fasta",
        help="Adapter FASTA file used by shiver_init."
    )
    parser.add_argument(
        "--shiver_primers",
        default="MyPrimers.fasta",
        help="Primer FASTA file used by shiver_init. You can use adapters.fa"
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

    # 2) Adjust global lists (if the user passes different values)
    MIN_READS_LIST = [int(x) for x in args.MIN_READS_LIST.split(",")]
    READS_GAMETOS_MINORITARIOS_LIST = [int(x) for x in args.READS_GAMETOS_MINORITARIOS_LIST.split(",")]
    THRESHOLD_GAMETO_FREQ_LIST = [float(x) for x in args.THRESHOLD_GAMETO_FREQ_LIST.split(",")]

    # 3) Adjust verbose and noindel
    VERBOSE = args.verbose
    noindel_option = args.noindel

    # 4) Base paths
    base_dir = Path(__file__).resolve().parent

    try:
        if args.index_only:
            if not args.custom_ref:
                print_error("[ERROR] --index_only requires --custom_ref.")
                sys.exit(1)
            if not args.virus:
                print_error("[ERROR] --index_only requires --virus <folder_name> (subfolder under reference_indexes/).")
                sys.exit(1)

            refseq = create_custom_reference_index(args.custom_ref, base_dir, args.virus)

            print_step("[REF] Custom index created successfully.")
            print_substep(f"[REF] Index stored in: reference_indexes/{str(args.virus).lower().replace(' ', '_')}/")
            print_substep(f"[REF] To use it, run with: --virus {args.virus}")
            sys.exit(0)

        refseq = get_reference_by_virus(base_dir, args.virus)
        if refseq is None:
            raise ValueError("refseq is None after reference selection.")

        print_substep(f"[REF] Virus selected: {args.virus}")
        print_substep(f"[REF] Reference used: {refseq}")

    except Exception as e:
        print_error(f"[ERROR] Reference/index preparation failed: {e}")
        sys.exit(1)

    # 5) FASTQ mode vs SRA/BAM mode
    if args.fastq:
        # 3 positional args expected: fastq1, fastq2, folder
        # if len(args.positional_args) < 2:
        #     print_error("[ERROR] In --simulacion (local) mode, 2-3 positional arguments required: <fastq1> [fastq2] <simulation_folder>")
        #     sys.exit(1)

        # fastq1 = args.positional_args[0]
        # if len(args.positional_args) == 3:
        #     fastq2 = args.positional_args[1]
        #     simulacion_folder = args.positional_args[2]

        # if len(args.positional_args) == 2:
        #     print(f"\nOnly one fastq recognized {fastq1}, interpreting as single mode")
        #     fastq2 = False
        #     simulacion_folder = args.positional_args[1]

        if args.bam:                     # start directly from BAM
            simulacion_folder = args.positional_args[0]
            modo = "bam"
            dedup_bam = f"{simulacion_folder}.bam"
            print(f"dedup_bam: {dedup_bam}")
            work_dir = base_dir / simulacion_folder
            srr=str(simulacion_folder).split("/")[-1]
            os.chdir(work_dir)
            # IMPORTANT: In --bam mode WITHOUT --use_shiver_assembly
            # keep the previous behavior without changing the reference.
            ok = pasos_ravi(srr, base_dir, refseq, work_dir, args.platform, dedup_bam, modo, noindel=noindel_option)
        elif args.use_shiver_assembly:
            # === Flow: trimming -> SPAdes+SHIVER -> RAVI ===
            # A) 3 positional args: <fastq1> <fastq2> <folder> (as before)
            # B) 1 positional SRR/ERR/DRR: download with prefetch+fasterq-dump into <folder>=<SRR>
            modo = "bam"
            print("Entered SHIVER assembly mode")
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
                # === Copy the EXACT logic from the mode that already works ===
                # 1) prefetch to base_dir  -> base_dir/<SRR>/<SRR>.sra
                # 2) fasterq-dump with PATH to .sra -> work_dir/<SRR>_1.fastq/_2.fastq
                fastq_1, fastq_2 = _sra_to_fastq(
                    srr=srr,
                    base_dir=base_dir,
                    out_dir=work_dir,
                    threads=getattr(args, "shiver_threads", 4)
                )
            else:
                print_error("[ERROR] --use_shiver_assembly requires <fastq1> <fastq2> <output_folder> or a single SRR/ERR/DRR accession.")
                sys.exit(1)

            # Assume PAIRED when R1/R2 exist
            layout = "PAIRED"
            os.chdir(work_dir)

            # Trimming (uses your usual run_trimming)
            r1_trim = work_dir / "reads_trimmed_1.fastq"
            r2_trim = work_dir / "reads_trimmed_2.fastq"
            if not (r1_trim.exists() and r2_trim.exists()):
                _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=False)
                if not _ok:
                    print("[WARN] Trimming failed. Retrying with filter_unpaired=True ...")
                    _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=True)
                if not _ok:
                    print_error("[ERROR] Trimming failed before SHIVER.")
                    sys.exit(1)

            bam_out = work_dir / f"{simulacion_folder}.bam"
            print(f"{bam_out} already exists; skipping SHIVER")
            if not bam_out.exists():
                # SPAdes + SHIVER
                bam_out = run_shiver_assembly(
                    srr=srr,
                    fastq_r1_trimmed=str(r1_trim),
                    fastq_r2_trimmed=str(r2_trim),
                    work_dir=str(work_dir),
                    base_dir=str(base_dir),
                    shiver_dir=args.shiver_dir,
                    spades=args.spades,
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
            print_substep(f"[SHIVER] BAM generated: {bam_out}; entering Bam and mapping {bam_out} with {refseq_new}")
            ok = pasos_ravi(srr, base_dir, refseq_new, work_dir, args.platform, bam_out, modo, noindel=noindel_option)

        # if args.slow:
        #     # LOCAL SIMULATION SLOW MODE
        #     modo = "slow"
        #     ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, args.platform, noindel_option)
        # else:
        #     # LOCAL SIMULATION FAST MODE
        #     modo = "fast"
        #     ok = inicio_comun_simulacion(fastq1, fastq2, simulacion_folder, base_dir, refseq, modo, args.platform, noindel_option)

        if ok:
            sys.exit(0)
            print_step("Pipeline completed successfully.")
        else:
            print_error("[ERROR] FASTQ pipeline failed.")
            sys.exit(1)

    # If not FASTQ mode, the first positional argument is the SRR accession
    if len(args.positional_args) < 1:
        print_error("[ERROR] Missing SRR. Or use --fastq <fastq1> [fastq2] <output_folder> (for local data).")
        sys.exit(1)

    srr = args.positional_args[0]

    # 6) Choose pipeline according to arguments
    pipeline_attempt = 1
    while pipeline_attempt <= MAX_PIPELINE_ATTEMPTS:
        if args.bam:
            # BAM flow
            simulacion_folder = args.positional_args[0] if len(args.positional_args) > 0 else srr
            modo = "bam"
            dedup_bam = f"{simulacion_folder}.bam"
            print(f"dedup_bam: {dedup_bam}")
            work_dir = base_dir / simulacion_folder
            srr = str(simulacion_folder).split("/")[-1]
            os.chdir(work_dir)
            # IMPORTANT: In --bam mode WITHOUT --use_shiver_assembly
            # keep the previous behavior without changing the reference.
            ok = pasos_ravi(srr, base_dir, refseq, work_dir, args.platform, dedup_bam, modo, noindel=noindel_option)
        elif args.use_shiver_assembly:
            # === trimming -> SPAdes + SHIVER -> RAVI ===
            print("Entered SHIVER assembly mode")
            # (your existing block: with <fastq1> <fastq2> <folder> or SRR)
            modo = "bam"
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
                print_error("[ERROR] --use_shiver_assembly requires <fastq1> <fastq2> <output_folder> or a single SRR/ERR/DRR accession.")
                sys.exit(1)

            # Assume PAIRED when R1/R2 exist
            layout = "PAIRED"
            os.chdir(work_dir)

            # Trimming (uses your usual run_trimming)
            r1_trim = work_dir / "reads_trimmed_1.fastq"
            r2_trim = work_dir / "reads_trimmed_2.fastq"
            if not (r1_trim.exists() and r2_trim.exists()):
                _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=False)
                if not _ok:
                    print("[WARN] Trimming failed. Retrying with filter_unpaired=True ...")
                    _ok = run_trimming(str(fastq_1), str(fastq_2), layout, base_dir, work_dir, filter_unpaired=True)
                if not _ok:
                    print_error("[ERROR] Trimming failed before SHIVER.")
                    sys.exit(1)

            # SPAdes + SHIVER
            bam_out = work_dir / f"{simulacion_folder}.bam"
            print(f"{bam_out} already exists; skipping SHIVER")
            if not bam_out.exists():
                bam_out = run_shiver_assembly(
                    srr=srr,
                    fastq_r1_trimmed=str(r1_trim),
                    fastq_r2_trimmed=str(r2_trim),
                    work_dir=str(work_dir),
                    base_dir=str(base_dir),
                    shiver_dir=args.shiver_dir,
                    spades=args.spades,
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
            print_substep(f"[SHIVER] BAM generated: {bam_out}; entering Bam and mapping {bam_out} with {refseq_new}")
            ok = pasos_ravi(srr, base_dir, refseq_new, work_dir, args.platform, bam_out, modo, noindel=noindel_option)
        elif args.slow:
            # classic SLOW flow -> slow_comun(...)
            modo = "slow"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)
        else:
            # classic FAST flow -> fast_comun(...)
            modo = "fast"
            ok = inicio_comun(srr, base_dir, refseq, modo, noindel_option)

        if ok:
            print_step("Pipeline completed successfully.")
            sys.exit(0)
        else:
            print_error(f"Pipeline error. Retrying... (attempt {pipeline_attempt}/{MAX_PIPELINE_ATTEMPTS})")
            pipeline_attempt += 1

    print_error("[ERROR] Pipeline failed after multiple attempts.")
    sys.exit(1)

if __name__ == "__main__":
    main()
