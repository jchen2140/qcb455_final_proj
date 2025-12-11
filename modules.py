import numpy as np, pandas as pd, subprocess
from statsmodels.stats.multitest import fdrcorrection
import subprocess
import time

# YUBI: modified script to run on local machine, returns smaller graph

print("Loading data...")
# Load GLS results

genes = pd.read_csv('genes.txt', header=None).squeeze()
GLS_p = pd.DataFrame(np.load('GLS_p.npy'), columns=genes, index=genes)

# Compute and save weights for ClusterONE

# Threshold for filtering edges. We keep only edges where the weight is strong.
# Weights are 1 - FDR, so a weight >= 0.95 corresponds to an FDR <= 0.05.
# Test with strict threshold of 0.99 and weaker thresholds
WEIGHT_THRESHOLD = 0.6

print("Computing weights...")
stacked_p = GLS_p.stack()
stacked_p = stacked_p[stacked_p.index.get_level_values(0) <
                      stacked_p.index.get_level_values(1)]
fdr = pd.Series(fdrcorrection(stacked_p)[1], index=stacked_p.index)
weights = 1 - fdr

# YUBI ADDED: CRITICAL FILTERING STEP
# Reduce the graph density by keeping only edges that meet the weight threshold.
initial_edge_count = len(weights)
weights = weights[weights >= WEIGHT_THRESHOLD]
final_edge_count = len(weights)

if final_edge_count == 0:
    print("\n--- ERROR ---")
    print("Filtering removed all edges. Check your WEIGHT_THRESHOLD or input data.")
    print("-------------------")
    exit(1)
    
print(f"Initial edges: {initial_edge_count:,}. Filtered edges (weight >= {WEIGHT_THRESHOLD}): {final_edge_count:,}")
print(f"Reduction ratio: {(1 - (final_edge_count / initial_edge_count)) * 100:.2f}%")
print(f"Saving filtered weights to disk...")

weight_file = 'ClusterONE_weights.txt'
weights.to_csv(weight_file, sep=' ', header=None)

# Run ClusterONE at each d, and save results
print("Running ClusterONE...")
# for d in 0.2, 0.5, 0.9:
    # subprocess.check_call(f'java -jar cluster_one-1.0.jar {weight_file} '
                          # f'-d {d} -F csv > modules_d_{d}.csv', shell=True)

# Set the memory allocation flag. Adjust the number (12G) based on your Mac's total RAM (16G). This is pushing it.
JAVA_MEMORY_FLAG = '-Xmx12G' 

for d in 0.2, 0.5, 0.9:
    subprocess.check_call(f'java {JAVA_MEMORY_FLAG} -jar cluster_one-1.0.jar {weight_file} '
                          f'-d {d} -F csv > modules_d_{d}.csv', shell=True)
    
print("Done.")