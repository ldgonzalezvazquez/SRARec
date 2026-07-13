#!/bin/bash

##############################################################################
##############################################################################

PARAM_FILE="param_list_new.txt"

if [ ! -f "$PARAM_FILE" ]; then
    echo "[ERROR] No existe $PARAM_FILE"
    exit 1
fi

N=$(wc -l < "$PARAM_FILE")
echo "[submit_missing.sh] $PARAM_FILE has $N missing lines."

if [ "$N" -eq 0 ]; then
    echo "[INFO] Finished."
    exit 0
fi

##############################################################################
##############################################################################

cat << 'EOF' > pipeline_param_missing.sh
#!/bin/bash
#SBATCH -C clk
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time 01:00:00

module load cesga/2020 gcc/system sra-toolkit/2.10.9 edirect/20240904 cutadapt/3.5 porechop/0.2.4 bwa/0.7.17 minimap2/2.20 htslib/1.19 samtools/1.19

echo "[pipeline_param_missing.sh] Starting in host=$(hostname), jobID=$SLURM_JOB_ID"

INDEX="$1"
PARAM_FILE="$2"

if [ -z "$INDEX" ] || [ -z "$PARAM_FILE" ]; then
    echo "[ERROR] Usage: $0 <INDEX> <PARAM_FILE>"
    exit 1
fi

if [ ! -f "$PARAM_FILE" ]; then
    echo "[ERROR] Does not exist $PARAM_FILE"
    exit 1
fi

##############################################################################
##############################################################################

NEXT=$((INDEX + 100))
TOTAL=$(wc -l < "$PARAM_FILE")

if [ "$NEXT" -lt "$TOTAL" ]; then
    echo "[INFO] Launching next job: INDEX=$NEXT"
    sbatch --job-name="missing_${NEXT}" \
           --dependency=afterany:${SLURM_JOB_ID} \
           --output="outputs_sbatch_missing/slurm_%j.out" \
           --error="outputs_sbatch_missing/slurm_%j.err" \
           pipeline_param_missing.sh "$NEXT" "$PARAM_FILE"
else
    echo "[INFO] Finished: NEXT=$NEXT >= TOTAL=$TOTAL"
fi

##############################################################################
##############################################################################

LINE=$(sed -n "$((INDEX+1))p" "$PARAM_FILE")

if [ -z "$LINE" ]; then
    echo "[ERROR] No data in line $INDEX de $PARAM_FILE"
    exit 1
fi

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

OUT_PREFIX="sim_${COV}_${RL}_${D}_${I}_${DEL}_${ID}_${II}_${IDEL}_${MR}_${RM}_${TH}_rep${REP}"

echo "[INFO] INDEX=$INDEX"
echo "[INFO] LINE=$LINE"
echo "[INFO] OUT_PREFIX=$OUT_PREFIX"

echo "[INFO] Python:"
which python3
python3 --version

##############################################################################
##############################################################################

python3 -u simulations.py \
    --results_file results_fast.txt \
    --coverage "$COV" \
    --read_length "$RL" \
    --diversity "$D" \
    --insertion_rate "$I" \
    --deletion_rate "$DEL" \
    --MIN_READS_LIST "$MR" \
    --READS_GAMETOS_MINORITARIOS_LIST "$RM" \
    --THRESHOLD_GAMETO_FREQ_LIST "$TH" \
    --output_prefix "$OUT_PREFIX" \
    --mutacion_intrahost \
    --intrahost_diversity "$ID" \
    --intrahost_insertion_rate "$II" \
    --intrahost_deletion_rate "$IDEL"

RC_PY=$?
echo "[pipeline_param_missing.sh] Python ended with code=$RC_PY"

rm -rf "$OUT_PREFIX"

echo "[pipeline_param_missing.sh] END jobID=$SLURM_JOB_ID"
exit 0
EOF

chmod +x pipeline_param_missing.sh

##############################################################################
##############################################################################

mkdir -p outputs_sbatch_missing

MAX_CHAINS=100

if [ "$N" -lt "$MAX_CHAINS" ]; then
    MAX_CHAINS="$N"
fi

echo "[submit_missing.sh] Sending the first $MAX_CHAINS jobs..."

for (( i=0; i<MAX_CHAINS; i++ )); do
    sbatch --job-name="missing_${i}" \
           --output="outputs_sbatch_missing/slurm_%j.out" \
           --error="outputs_sbatch_missing/slurm_%j.err" \
           pipeline_param_missing.sh "$i" "$PARAM_FILE"
done

echo "[submit_missing.sh] Completed."
