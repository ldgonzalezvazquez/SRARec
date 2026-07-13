import csv

# 1. Read the results already generated (output_prefix from the CSV)
existing_results = set()

with open("results_simulations.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        existing_results.add(row['output_prefix'])

# 2. Check param_list.txt line by line
missing_params = []

with open("param_list.txt") as f:
    for line in f:
        values = line.strip().split()
        if len(values) != 12:
            continue  # Malformed line

        # Reorder the values according to the OUT_PREFIX format
        COV  = values[0]
        RL   = values[1]
        D    = values[2]
        I    = values[3]
        DEL  = values[4]
        MR   = values[5]
        RM   = values[6]
        TH   = values[7]
        ID   = values[8]
        II   = values[9]
        IDEL = values[10]
        REP  = values[11]

        out_prefix = f"sim_{COV}_{RL}_{D}_{I}_{DEL}_{ID}_{II}_{IDEL}_{MR}_{RM}_{TH}_rep{REP}"

        if out_prefix not in existing_results:
            missing_params.append(line)

# 3. Write the lines that are not present in results_simulaciones.csv
with open("param_list_new.txt", "w") as f_out:
    for line in missing_params:
        f_out.write(line)

# 4. Print summary information
print(f"[INFO] Total number of parameters in param_list.txt: {len(open('param_list.txt').readlines())}")
print(f"[INFO] Parameters NOT found in results_simulaciones.csv: {len(missing_params)}")
print(f"[INFO] The file 'param_list_new.txt' has been created")
