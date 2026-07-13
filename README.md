## Running the simulations

The simulations are designed to run on CESGA using Slurm.

### 1. Enter the simulation directory

Choose the workflow you want to run:

```bash
cd Fast_simulations
```

or:

```bash
cd Slow_simulations
```

### 2. Check the dependencies

Before launching the simulations, make sure that the required CESGA modules, Python packages and NGSNGS installation are available.

```bash
module load cesga/2020 gcc/system cutadapt/3.5 \
  bwa/0.7.17 minimap2/2.20 htslib/1.19 samtools/1.19

python3 -c "import Bio, pysam, tqdm"
```

The path to the NGSNGS executable must also be correctly defined in `simulations.py`.

### 3. Give execution permissions

```bash
chmod +x SRARec.py
chmod +x submit_simulation_cascade.sh
chmod +x submit_missing*.sh
```

### 4. Set the simulation parameters

The simulation parameters are defined at the beginning of `submit_simulation_cascade.sh`.

For example:

```python
coverages = [200, 350, 500]
read_lengths = [300, 500, 1000]
diversities = [0.0022]

insertion_rates = [0.00, 0.00010033]
deletion_rates = [0.00, 0.00053507]

min_reads_list = [5, 10, 20]
reads_gametos_minor_list = [5, 10, 20]
threshold_gameto_freq_list = [0.05, 0.1]

intrahost_diversities = [0.00, 0.0003]
intrahost_insertion_rates = [0.00, 0.00005]
intrahost_deletion_rates = [0.00, 0.00005]

REPLICATES = 2
```

Every possible combination of these values will be simulated.

### 5. Submit the simulations

Once the parameters are ready, run:

```bash
./submit_simulation_cascade.sh
```

This script creates the parameter list and submits all simulations to the CESGA queue.

You can check the running jobs with:

```bash
squeue -u "$USER"
```

### 6. Check for missing simulations

After all jobs have finished, check whether any configurations failed or are missing:

```bash
python3 error_check.py
```

This creates a new file with the missing configurations:

```text
param_list_new.txt
```

If `param_list_new.txt` is not empty, rerun the missing simulations with:

```bash
./submit_missing.sh
```

Repeat this step until no missing configurations remain:

```bash
python3 error_check.py
wc -l param_list_new.txt
```

### 7. Generate the final results

Once all simulations have finished successfully, run:

```bash
python3 generate_figures_csv.py
```

This script uses `results_simulations.csv` to generate the final simulation summaries and figures.
