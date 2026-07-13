#!/usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt

# CSV file to read
csv_file = "results_simulations.csv"  # Must be in the same directory

# Check if the file exists
if not os.path.exists(csv_file):
    print(f"Error: The file {csv_file} was not found.")
    exit()

# Load data
df = pd.read_csv(csv_file)

# Parameters to analyze (excluding 'breakpoint')
param_columns = [
    'coverage', 'read_length', 'diversity', 'insertion_rate', 'deletion_rate',
    'mutacion_intrahost', 'intrahost_diversity', 'intrahost_insertion_rate',
    'intrahost_deletion_rate', 'MIN_READS_LIST', 'READS_GAMETOS_MINORITARIOS_LIST',
    'THRESHOLD_GAMETO_FREQ_LIST'
]

# Create output directory
output_dir = "simulation_results_graphs"
os.makedirs(output_dir, exist_ok=True)

# Determine subplot arrangement
num_params = len(param_columns)
cols = 3
rows = (num_params + cols - 1) // cols

# Create the figure
fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
axes = axes.flatten()

# Explanation above the plots
fig.suptitle(
    "\nSegmented Line = 100% of (Detections + False Positives):\n"
    "   - Green: Real events\n"
    "   - Red: False positives\n\n"
    "Blue Bar: Detection rate\n"
    "Orange Bar: Lines with FP\n"
    "Black Bar: Lines with True + FP",
    fontsize=12, fontweight='bold', color='black', ha='center'
)

def compute_metrics(subdf):
    """
    Returns a Series with the relevant counts and computed rates
    for one value of the parameter.
    """
    total_simulations = len(subdf)
    
    # 'result' -> True = detection
    detection_count = subdf['result'].sum()
    
    # Count of FP for that group
    false_positive_count = subdf['falsos_positivos'].sum()
    # # lines that have >=1 FP
    false_positive_cases = (subdf['falsos_positivos'] > 0).sum()
    # lines that have True & FP
    cases_true_fp = len(subdf[(subdf['result'] == True) & (subdf['falsos_positivos'] > 0)])
    
    return pd.Series({
        'total_simulations': total_simulations,
        'detection_count': detection_count,
        'false_positive_count': false_positive_count,
        'false_positive_cases': false_positive_cases,
        'cases_true_fp': cases_true_fp
    })

for i, param in enumerate(param_columns):
    if param in df.columns and df[param].nunique() > 1:  # Avoid constant parameters
        grouped = df.groupby(param).apply(compute_metrics).reset_index()

        # Calculate bars (heights in %)
        grouped['detection_rate'] = (grouped['detection_count'] / grouped['total_simulations']) * 100
        grouped['fp_event_rate'] = (grouped['false_positive_cases'] / grouped['total_simulations']) * 100
        grouped['true_fp_rate'] = (grouped['cases_true_fp'] / grouped['total_simulations']) * 100
        
        ax = axes[i]
        width = 0.25  # width for each bar
        x_positions = range(len(grouped[param]))
        
        # BLUE bar: detection rate
        bars_blue = ax.bar(
            [x for x in x_positions],
            grouped['detection_rate'],
            width=width,
            color='blue', alpha=0.6
        )
        
        # ORANGE bar: lines with FP
        bars_orange = ax.bar(
            [x + width for x in x_positions],
            grouped['fp_event_rate'],
            width=width,
            color='orange', alpha=0.6
        )
        
        # BLACK bar: lines with True + FP
        bars_black = ax.bar(
            [x + 2*width for x in x_positions],
            grouped['true_fp_rate'],
            width=width,
            color='black', alpha=0.6
        )
        
        # Place absolute counts above each bar
        for bar, val in zip(bars_blue, grouped['detection_count']):
            ax.text(
                bar.get_x() + bar.get_width()/2,
                bar.get_height() + 2,
                str(int(val)),
                ha='center', fontsize=10, color='black'
            )
        for bar, val in zip(bars_orange, grouped['false_positive_cases']):
            ax.text(
                bar.get_x() + bar.get_width()/2,
                bar.get_height() + 2,
                str(int(val)),
                ha='center', fontsize=10, color='black'
            )
        for bar, val in zip(bars_black, grouped['cases_true_fp']):
            ax.text(
                bar.get_x() + bar.get_width()/2,
                bar.get_height() + 2,
                str(int(val)),
                ha='center', fontsize=10, color='black'
            )

        # =========== SEGMENTED LINE (Green+Red) =========== 
        # This represents 100% of "detections + false positives" for each param.
        # We'll draw it as a horizontal line at y=105 (5% above the top).
        
        line_height = 105
        for idx, row in grouped.iterrows():
            # total "events" = detection_count + false_positive_count
            total_events = row['detection_count'] + row['false_positive_count']
            if total_events <= 0:
                # no segment to show
                continue
            
            # fraction that is detection (green) vs FP (red)
            green_fraction = row['detection_count'] / total_events
            red_fraction = row['false_positive_count'] / total_events
            
            # We'll center the line near the same X as the bar group
            x_center = idx + width  # shift by ~width to align with middle group
            line_width = 0.5        # total width for the line segment
            start_x = x_center - line_width/2
            
            green_length = green_fraction * line_width
            red_length = red_fraction * line_width
            
            # Green part
            ax.hlines(
                y=line_height,
                xmin=start_x,
                xmax=start_x + green_length,
                color='green',
                linewidth=5
            )
            # Red part
            ax.hlines(
                y=line_height,
                xmin=start_x + green_length,
                xmax=start_x + green_length + red_length,
                color='red',
                linewidth=5
            )
            
            # optional text indicating real vs FP absolute counts
            # place them slightly above or below the line if needed
            # e.g. ax.text(...)

        # Format axis
        ax.set_ylabel("Percentage (%)")
        ax.set_xlabel(param)
        ax.set_xticks([x + width for x in x_positions])
        ax.set_xticklabels(grouped[param].astype(str), rotation=45, ha='right')
        ax.set_ylim(0, 110)  # 0% to 110% to show the line at ~105

# Hide empty subplots if fewer parameters
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Explanation text below
fig.text(
    0.5, 0.01,
    "Calculation Details:\n"
    "Blue bar = detection_rate = (detections / total_simulations) * 100; label = detection_count\n"
    "Orange bar = % lines with >=1 FP = (FP_cases / total_simulations) * 100; label = lines with FP\n"
    "Black bar = % lines that are True + FP = (cases_true_fp / total_simulations) * 100; label = lines True & FP\n"
    "Segmented line = 100% of (detection_count + false_positive_count): green = real events, red = false positives.\n",
    ha='center', fontsize=10, color='black', bbox=dict(facecolor='white', alpha=0.8)
)

plt.tight_layout(rect=[0, 0.05, 1, 0.92])
plt.savefig(os.path.join(output_dir, "all_simulation_graphs_miguel.png"))
plt.show()

print(f"All graphs have been saved in: {output_dir}")