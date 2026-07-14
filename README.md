## Running the simulations

The simulation workflow was originally prepared for [CESGA (Galician Supercomputing Technology Centre)](https://www.cesga.es/) using Slurm, but it can be adapted to any Linux or HPC environment with the required dependencies installed. The current submission scripts use Slurm commands (`sbatch`, `squeue`, `sacct`) and CESGA-specific module loading, so these parts should be modified if another queue system is used.

### 1. Enter the simulation directory

If the repository has not yet been downloaded, clone the branch containing the simulation workflows:

```bash
git clone --branch SRARec_simulations --single-branch \
  https://github.com/ldgonzalezvazquez/SRARec.git

cd SRARec
```

Then enter the directory corresponding to the workflow you want to run (fast or slow simulations)

If the repository has already been cloned, simply navigate to the corresponding directory.

```bash
cd SRARec/Fast_simulations
```
```bash
cd SRARec/Fast_simulations
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

Before running the workflow for the first time, give execution permissions to the shell scripts used to generate the parameter combinations, submit the jobs and execute each simulation.

To give execution permission to all shell scripts in the current simulation directory:

```bash
chmod +x *.sh
```
Alternatively, permissions can be assigned individually:

```bash
chmod +x <script_name>
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
