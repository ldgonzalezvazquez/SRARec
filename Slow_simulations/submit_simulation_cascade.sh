#!/bin/bash

##############################################################################
# 1) Create or overwrite param_list.txt with all parameter combinations
##############################################################################
echo "[submit_simulation_cascade.sh] Creating param_list.txt with all parameter combinations..."

cat << 'END_PY' > make_param_list.py
import itertools

# Adjust these values according to your needs
coverages = [200,350,500]
read_lengths = [300,500,1000]
diversities = [0.0022]
insertion_rates = [0.00,0.00010033]
deletion_rates = [0.00,0.00053507]
min_reads_list = [5,10,20]
reads_gametos_minor_list = [5,10,20]
threshold_gameto_freq_list = [0.05,0.1]
intrahost_diversities = [0.00,0.0003]
intrahost_insertion_rates = [0.00,0.00005]
intrahost_deletion_rates = [0.00,0.00005]

REPLICATES = 2  # Number of replicates

lines = []
for cov, rl, d, ins, dele, mr, rm, th, id, iins, idele in itertools.product(
    coverages,
    read_lengths,
    diversities,
    insertion_rates,
    deletion_rates,
    min_reads_list,
    reads_gametos_minor_list,
    threshold_gameto_freq_list,
    intrahost_diversities,
    intrahost_insertion_rates,
    intrahost_deletion_rates
):
    for r in range(1, REPLICATES + 1):
        # Format: coverage read_length diversity insertion_rate deletion_rate
        #          MIN_READS READS_MINOR THRESHOLD intrahost_diversities intrahost_insertion_rates intrahost_deletion_rates replicate
        lines.append(f"{cov} {rl} {d} {ins} {dele} {mr} {rm} {th} {id} {iins} {idele} {r}")

with open("param_list.txt", "w") as fh:
    for line in lines:
        fh.write(line + "\n")
END_PY

python3 make_param_list.py
RC=$?
if [ $RC -ne 0 ]; then
echo "[submit_simulation_cascade.sh] param_list.txt was created with $N lines."
    exit 1
fi

N=$(wc -l < param_list.txt)
echo "[submit_master.sh] Se creó param_list.txt con $N líneas."
if [ "$N" -eq 0 ]; then
    echo "[ERROR] param_list.txt is empty."
    exit 1
fi

##############################################################################
# 2) Create pipeline_param.sh with the following logic:
#    - load modules
#    - read INDEX from param_list.txt
#    - submit the next job (INDEX + 200) at the beginning
#    - run the Python simulation script with the selected parameters
##############################################################################
cat << 'EOF' > pipeline_param.sh
#!/bin/bash
#SBATCH -C clk
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time 01:00:00
#SBATCH --output $STORE/outputs_sbatch/

## ------------------------------------------------------------------------------
## MÓDULOS
## ------------------------------------------------------------------------------
module load cesga/2020 gcc/system sra-toolkit/2.10.9 edirect/20240904 cutadapt/3.5 porechop/0.2.4 bwa/0.7.17 minimap2/2.20 htslib/1.19 samtools/1.19

echo "[pipeline_param.sh] Starting on host=$(hostname), jobID=$SLURM_JOB_ID"

INDEX="$1"
if [ -z "$INDEX" ]; then
    echo "[ERROR] Usage: $0 <INDEX>"
    exit 1
fi

# ------------------------------------------------------------------------------
# 1) Calculate the next line to process (INDEX + 200)
#    If it exists in param_list.txt, submit the next job with afterany.
#    This is done at the beginning so the next job is submitted regardless of
#    whether the current Python run succeeds or fails.
# ------------------------------------------------------------------------------
NEXT=$((INDEX + 200))
TOTAL=$(wc -l < param_list.txt)
if [ "$NEXT" -lt "$TOTAL" ]; then
    echo "[INFO] Submitting sbatch 'pipeline_param.sh $NEXT' with afterany:$SLURM_JOB_ID"
    sbatch --job-name="param_${NEXT}" \
           --dependency=afterany:${SLURM_JOB_ID} \
           --output="outputs_sbatch/slurm_%j.out" \
           pipeline_param.sh "$NEXT"
else
    echo "[INFO] No further job submitted (NEXT=$NEXT >= TOTAL=$TOTAL)."
fi

# ------------------------------------------------------------------------------
# 2) Read line INDEX from param_list.txt and parse its parameters
# ------------------------------------------------------------------------------
LINE=$(sed -n "$((INDEX+1))p" param_list.txt)
if [ -z "$LINE" ]; then
    echo "[ERROR] No data found at line $INDEX in param_list.txt"
    exit 1
fi

# Format: coverage read_length diversity insertion_rate deletion_rate MIN_READS READS_MINOR THRESHOLD replicate
COV=$(echo "$LINE" | awk '{print $1}')
RL=$(echo  "$LINE" | awk '{print $2}')
D=$(echo   "$LINE" | awk '{print $3}')
I=$(echo   "$LINE" | awk '{print $4}')
DEL=$(echo "$LINE" | awk '{print $5}')
MR=$(echo  "$LINE" | awk '{print $6}')
RM=$(echo  "$LINE" | awk '{print $7}')
TH=$(echo  "$LINE" | awk '{print $8}')
ID=$(echo   "$LINE" | awk '{print $9}')
II=$(echo   "$LINE" | awk '{print $10}')
IDEL=$(echo "$LINE" | awk '{print $11}')
REP=$(echo "$LINE" | awk '{print $12}')


echo "[INFO] INDEX=$INDEX => $LINE"

# ------------------------------------------------------------------------------
# 3) Run simulations.py with the selected parameters
# ------------------------------------------------------------------------------
OUT_PREFIX="sim_${COV}_${RL}_${D}_${I}_${DEL}_${ID}_${II}_${IDEL}_${MR}_${RM}_${TH}_rep${REP}"

echo "Current directory: $(pwd)"
ls -l



echo "[INFO] Python version in this job:"
python --version
which python

echo "[INFO] Launching simulations.py"

python3 -u simulations.py --results_file results_slow.txt --coverage "$COV" --read_length "$RL" --diversity "$D" --insertion_rate "$I" --deletion_rate "$DEL" --MIN_READS_LIST "$MR" --READS_GAMETOS_MINORITARIOS_LIST "$RM" --THRESHOLD_GAMETO_FREQ_LIST "$TH" --output_prefix "$OUT_PREFIX" --mutacion_intrahost --intrahost_diversity "$ID" --intrahost_insertion_rate "$II" --intrahost_deletion_rate "$IDEL" 

RC_PY=$?
echo "[pipeline_param.sh] Python finished with exit code=$RC_PY"

# Remove the temporary simulation directory
rm -rf "$OUT_PREFIX"

echo "[pipeline_param.sh] Finished jobID=$SLURM_JOB_ID"
exit 0
EOF

chmod +x pipeline_param.sh
echo "[submit_simulation_cascade.sh] Created pipeline_param.sh"

##############################################################################
# Submit the first 200 jobs, corresponding to indices 0..199
##############################################################################
mkdir -p outputs_sbatch

MAX_CHAINS=5
if [ "$N" -lt "$MAX_CHAINS" ]; then
    MAX_CHAINS="$N"
fi

echo "[submit_simulation_cascade.sh] Submitting the first $MAX_CHAINS jobs..."

for (( i=0; i<MAX_CHAINS; i++ )); do
    sbatch --job-name="param_${i}" \
           --output="outputs_sbatch/slurm_%j.out" \
           --error="outputs_sbatch/slurm_%j.err" \
           pipeline_param.sh "$i"
done

echo "[submit_simulation_cascade.sh] Job submission process completed."
echo "[submit_simulation_cascade.sh] Use squeue or sacct to monitor the chained jobs."
