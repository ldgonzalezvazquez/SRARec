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
