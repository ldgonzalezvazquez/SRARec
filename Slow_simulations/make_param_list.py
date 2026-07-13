import itertools

# Adjust these values according to your needs
coverages = [200,350,500]
read_lengths = [300,500,1000]
diversities = [0.001,0.005,0.01]
insertion_rates = [0.00,0.00010032]
deletion_rates = [0.00,0.00010032]
min_reads_list = [5,10,20]
reads_gametos_minor_list = [5,10,20]
threshold_gameto_freq_list = [0.05,0.1,0.2]
intrahost_diversities = [0.00,0.0003]
intrahost_insertion_rates = [0.00,0.00003344]
intrahost_deletion_rates = [0.00,0.00003344]


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
