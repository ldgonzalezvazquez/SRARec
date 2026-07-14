import argparse
from Bio import SeqIO
import random
import math
import os
import subprocess
import re
import csv

def mutate_sequence(seq, pi_diversity, insertion_rate, deletion_rate, max_indel_size=5):
    """
    Introduce mutations into the original sequence according to pi diversity (π).
    """
    nucleotides = ['A', 'T', 'C', 'G']
    seq_list = list(seq)
    num_sites = len(seq_list)

    num_insertions = int(num_sites * insertion_rate) if insertion_rate else 0
    num_deletions = int(num_sites * deletion_rate) if deletion_rate else 0
    num_substitutions = (pi_diversity - (num_insertions + num_deletions)) if pi_diversity else 0

    total_mutations = num_insertions + num_deletions + num_substitutions
    if total_mutations <= 0:
        return "".join(seq_list)

    mutation_positions = random.sample(range(num_sites), total_mutations)
    insert_positions = set(random.sample(mutation_positions, num_insertions)) if num_insertions > 0 else set()
    delete_positions = set(random.sample(mutation_positions, num_deletions)) if num_deletions > 0 else set()

    new_seq = []
    i = 0
    while i < len(seq_list):
        if i in mutation_positions:
            if i in insert_positions:
                insert_size = random.randint(1, max_indel_size)
                insertion = ''.join(random.choices(nucleotides, k=insert_size))
                new_seq.append(insertion + seq_list[i])
            elif i in delete_positions:
                delete_size = random.randint(1, max_indel_size)
                i += (delete_size - 1)
            else:
                # sustitución
                new_nuc = random.choice([n for n in nucleotides if n != seq_list[i]])
                new_seq.append(new_nuc)
        else:
            new_seq.append(seq_list[i])
        i += 1

    return "".join(new_seq)


def generate_recombinants(seq1, seq2):
    """
    Generate two recombinant sequences by exchanging sequence fragments at a random breakpoint.    """
    breakpoint = random.randint(1, len(seq1) - 1)
    recombinant1 = seq1[:breakpoint] + seq2[breakpoint:]
    recombinant2 = seq2[:breakpoint] + seq1[breakpoint:]
    return recombinant1, recombinant2, breakpoint


def run_ngsngs_in_directory(fasta_path, output_dir, coverage, read_length):
    """
    Temporarily change the working directory to output_dir,
    run NGSNGS so FASTQ files are generated inside output_dir,
    and then return to the original directory.
    """
    old_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        fasta_base = os.path.basename(fasta_path)
        command = [
            "/mnt/lustre/scratch/nlsas/home/uvi/cm/apc/NGSNGS-master/ngsngs",
            "-i", fasta_base,
            "-c", str(coverage),
            "-l", str(read_length),
            "-seq", "PE",
            "-qs", "100",
            "-f", "fq",
            "-o", "./{args.output_prefix}"
        ]
        subprocess.run(command, check=True)
    finally:
        os.chdir(old_cwd)  


def run_pipeline_simulaciones(output_prefix, min_reads, reads_gametos, threshold_gameto):
    """
    Run SRARec.py with '--simulacion'.
    """
    command = [
        "./SRARec.py",
        "*_R1.fq",
        "*_R2.fq",
        output_prefix,
        "--slow",
        "--simulacion",
        "--MIN_READS_LIST", min_reads,
        "--READS_GAMETOS_MINORITARIOS_LIST", reads_gametos,
        "--THRESHOLD_GAMETO_FREQ_LIST", threshold_gameto
    ]
    subprocess.run(command, check=True)


def verificar_breakpoint(archivo_results, output_prefix, breakpoint):
    """
   Check whether the breakpoint is present in the results file.

    The function searches for lines containing Sra:<output_prefix>,
    parses all Rec ranges, and returns:
    - 'true' if at least one Rec range contains the breakpoint
    - 'false' if no Rec range contains the breakpoint

    It also reports the number of false-positive Rec ranges, defined here as
    Rec ranges that do not include the simulated breakpoint.
    """
    pattern_sra = re.compile(rf"Sra:{output_prefix}")
    pattern_rec = re.compile(r"\tRec:(.*?)\tNo_Rec")

    found_any = False
    global false_positives
    false_positives = 0  

    with open(archivo_results, 'r') as f:
        for linea in f:
            # Only analyze lines matching Sra:<output_prefix>
            if not pattern_sra.search(linea):
                continue

            match_rec = pattern_rec.search(linea)
            if match_rec:
                rangos_texto = match_rec.group(1).strip().split(';')
                for rango in rangos_texto:
                    parts = rango.split(',')
                    if len(parts) >= 2:
                        try:
                            inicio, fin = int(parts[0]), int(parts[1])
                            if (inicio-100) <= breakpoint <= (fin+100):
                                found_any = True
                            else:
                                false_positives += 1
                        except ValueError:
                            pass

    if found_any:
        print(
            f"\033[92mAt least one range includes breakpoint {breakpoint} in {output_prefix}.\033[0m"
        )
        print(f"\033[93mFalse positives: {false_positives}\033[0m")
        return "true"
    else:
        print(
            f"\033[91mNo range includes breakpoint {breakpoint} in {output_prefix}.\033[0m"
        )
        print(f"\033[93mFalse positives: {false_positives}\033[0m")
        return "false"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Mutate genetic sequences, generate recombinants, and run NGS simulations."
    )
    parser.add_argument("--output_prefix", type=str, required=True, help="Output folder and prefix for generated files.")
    parser.add_argument("--results_file", type=str, required=True, help="Results file used to verify whether the breakpoint was detected.")
    parser.add_argument("--coverage", type=int, default=100, help="Coverage value for NGSNGS.")
    parser.add_argument("--read_length", type=int, default=300, help="Read length for NGSNGS.")
    parser.add_argument("--diversity", type=float, help="Pi diversity (π) factor to determine the number of mutations.")
    parser.add_argument("--insertion_rate", type=float, help="Insertion rate.")
    parser.add_argument("--deletion_rate", type=float, help="Deletion rate.")
    parser.add_argument("--mutacion_intrahost", action="store_true", help="Apply a second round of intrahost mutations")
    parser.add_argument("--intrahost_diversity", type=float, help=" Intrahost pi diversity (π).")
    parser.add_argument("--intrahost_insertion_rate", type=float, help="Intrahost insertion rate.")
    parser.add_argument("--intrahost_deletion_rate", type=float, help="Intrahost deletion rate.")
    parser.add_argument("--parental1", type=str, help="FASTA file containing the first parental sequence.")
    parser.add_argument("--parental2", type=str, help="FASTA file containing the second parental sequence.")
    parser.add_argument("--MIN_READS_LIST", type=str, default="5,12,30", help="MIN_READS_LIST for pipeline.")
    parser.add_argument("--READS_GAMETOS_MINORITARIOS_LIST", type=str, default="5,20", help="READS_GAMETOS for pipeline.")
    parser.add_argument("--THRESHOLD_GAMETO_FREQ_LIST", type=str, default="0.01,0.05,0.1", help="THRESHOLD_GAMETO for pipeline.")
    return parser.parse_args()


def main():
    args = parse_arguments()

    # 1) Create output folder if it does not exist
    os.makedirs(args.output_prefix, exist_ok=True)

    # 2) Obtain parental sequences
    if args.parental1 and args.parental2:
        record1 = next(SeqIO.parse(args.parental1, "fasta"))
        record2 = next(SeqIO.parse(args.parental2, "fasta"))
        sequence_1 = str(record1.seq)
        sequence_2 = str(record2.seq)
    else:
        reference_path = "./refseq/NC_045512.2.fna"
        record = next(SeqIO.parse(reference_path, "fasta"))
        original_sequence = str(record.seq)
        pi_diversity = math.ceil(args.diversity * len(original_sequence)) if args.diversity else 0
        sequence_1 = mutate_sequence(original_sequence, pi_diversity, args.insertion_rate or 0, args.deletion_rate or 0)
        sequence_2 = mutate_sequence(original_sequence, pi_diversity, args.insertion_rate or 0, args.deletion_rate or 0)

    # 3) Recombinants
    recombinant1, recombinant2, breakpoint = generate_recombinants(sequence_1, sequence_2)

    # 4) Aply intrahost mutations
    if args.mutacion_intrahost:
        pi_intrahost = 0
        if args.intrahost_diversity:
            pi_intrahost = math.ceil(args.intrahost_diversity * len(sequence_1))
        sequence_1 = mutate_sequence(sequence_1, pi_intrahost, args.intrahost_insertion_rate or 0, args.intrahost_deletion_rate or 0)
        sequence_2 = mutate_sequence(sequence_2, pi_intrahost, args.intrahost_insertion_rate or 0, args.intrahost_deletion_rate or 0)
        recombinant1 = mutate_sequence(recombinant1, pi_intrahost, args.intrahost_insertion_rate or 0, args.intrahost_deletion_rate or 0)
        recombinant2 = mutate_sequence(recombinant2, pi_intrahost, args.intrahost_insertion_rate or 0, args.intrahost_deletion_rate or 0)

    # 5) Save FASTA file inside the output folder
    fasta_name = os.path.join(args.output_prefix, f"{args.output_prefix}.fasta")
    with open(fasta_name, "w") as f:
        f.write(f">{args.output_prefix}_parental1\n{sequence_1}\n")
        f.write(f">{args.output_prefix}_parental2\n{sequence_2}\n")
        f.write(f">{args.output_prefix}_recombinante1\n{recombinant1}\n")
        f.write(f">{args.output_prefix}_recombinante2\n{recombinant2}\n")

    # 6) Run NGSNGS and generate FASTQ files inside the output folder
    run_ngsngs_in_directory(fasta_name, args.output_prefix, args.coverage, args.read_length)

    # 7) run simulation folder
    run_pipeline_simulaciones(
        args.output_prefix,
        args.MIN_READS_LIST,
        args.READS_GAMETOS_MINORITARIOS_LIST,
        args.THRESHOLD_GAMETO_FREQ_LIST
    )

    # 8) Verify breakpoint => 'true'/'false'
    result_status = verificar_breakpoint(args.results_file, args.output_prefix, breakpoint)

    # 9) Save parameters and result to CSV
    csv_filename = "results_simulations.csv"
    fieldnames = [
        "output_prefix", "results_file", "coverage", "read_length",
        "diversity", "insertion_rate", "deletion_rate",
        "mutacion_intrahost", "intrahost_diversity", "intrahost_insertion_rate", "intrahost_deletion_rate",
        "parental1", "parental2", "MIN_READS_LIST", "READS_GAMETOS_MINORITARIOS_LIST", "THRESHOLD_GAMETO_FREQ_LIST",
        "breakpoint", "result","false_positives"
    ]

    final_result = "true" if result_status == "true" else "false"

    row = {
        "output_prefix": args.output_prefix,
        "results_file": args.results_file,
        "coverage": args.coverage,
        "read_length": args.read_length,
        "diversity": args.diversity if args.diversity is not None else "",
        "insertion_rate": args.insertion_rate if args.insertion_rate is not None else "",
        "deletion_rate": args.deletion_rate if args.deletion_rate is not None else "",
        "mutacion_intrahost": "True" if args.mutacion_intrahost else "",
        "intrahost_diversity": args.intrahost_diversity if args.intrahost_diversity is not None else "",
        "intrahost_insertion_rate": args.intrahost_insertion_rate if args.intrahost_insertion_rate is not None else "",
        "intrahost_deletion_rate": args.intrahost_deletion_rate if args.intrahost_deletion_rate is not None else "",
        "parental1": args.parental1 if args.parental1 else "",
        "parental2": args.parental2 if args.parental2 else "",
        "MIN_READS_LIST": args.MIN_READS_LIST,
        "READS_GAMETOS_MINORITARIOS_LIST": args.READS_GAMETOS_MINORITARIOS_LIST,
        "THRESHOLD_GAMETO_FREQ_LIST": args.THRESHOLD_GAMETO_FREQ_LIST,
        "breakpoint": breakpoint,
        "result": final_result,
        "false_positives": false_positives
    }

    file_exists = os.path.isfile(csv_filename)
    with open(csv_filename, "a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"Process completed. FASTQ files were saved in folder '{args.output_prefix}'.")
    print("CSV results file: results_simulations.csv")

if __name__ == "__main__":
    main()
