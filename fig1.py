# Written by Yubi

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

print("Loading data...")

# 1. Define the file path
file_path = 'vizdf.tsv'

# 2. Load the TSV file
df = pd.read_csv(file_path, sep='\t')

print("Cleaning data...")

# 3. Clean and combine cluster columns
# The columns are cleaned by converting to lowercase string and stripping whitespace to handle
# 'unknown', 'Unknown', 'NaN', and other null representations consistently.
df['roarke_clusters'] = df['roarke_clusters'].fillna('').astype(str).str.lower().str.strip()
df['top_enriched_clusters'] = df['top_enriched_clusters'].fillna('').astype(str).str.lower().str.strip()

# Function to prioritize cluster annotations
def get_cluster_label(row):
    """Prioritizes roarke_clusters, then top_enriched_clusters, otherwise labels as 'Unknown'."""
    # Check for values that are likely just null/unknown markers
    null_markers = ['nan', 'unknown', '']

    roarke = row['roarke_clusters']
    enriched = row['top_enriched_clusters']

    if roarke not in null_markers:
        return roarke
    elif enriched not in null_markers:
        return enriched
    else:
        return 'Unknown'

print("Annotating clusters...")
# Apply the function to create a single Cluster column
df['Cluster'] = df.apply(get_cluster_label, axis=1)

# Separate data into Known and Unknown for plotting
df_known = df[df['Cluster'] != 'Unknown'].copy()
df_unknown = df[df['Cluster'] == 'Unknown'].copy()

known_clusters = df_known['Cluster'].unique()
print(f"Identified {len(known_clusters)} known clusters for annotation.")

print("Defining color palette...")
# 4. Define Color Palette
# Light gray for the large 'Unknown' background, and a distinct palette for the known clusters.
num_known_clusters = len(known_clusters)
color_map = {}
color_map['Unknown'] = '#CCCCCC' # Light grey

if num_known_clusters > 0:
    # Use a colorblind-friendly, high-contrast palette for the known clusters
    known_colors = sns.color_palette("bright", n_colors=num_known_clusters).as_hex()
    for cluster, color in zip(known_clusters, known_colors):
        color_map[cluster] = color

print("Creating scatter plot...")
# 5. Create the Scatter Plot
plt.figure(figsize=(12, 10))

# 6. Plot the 'Unknown' genes first (as background)
plt.scatter(
    df_unknown['hUMAP_x'],
    df_unknown['hUMAP_y'],
    color=color_map['Unknown'],
    s=1,          # Small size for background
    alpha=0.2,    # Highly translucent
    label=f'Unknown ({len(df_unknown)} genes)'
)

print("Plotting known clusters...")
# 7. Plot the known clusters on top
for cluster in known_clusters:
    subset = df_known[df_known['Cluster'] == cluster]
    plt.scatter(
        subset['hUMAP_x'],
        subset['hUMAP_y'],
        color=color_map[cluster],
        s=20,         # Larger size to make annotated points stand out
        alpha=0.8,
        edgecolors='black', # Optional: black border for better visibility
        linewidths=0.5,
        label=f'{cluster} ({len(subset)} genes)'
    )

# 8. Add plot details and legend
plt.title('Gene Co-essentiality UMAP with Annotated Clusters', fontsize=16)
plt.xlabel('hUMAP_x', fontsize=12)
plt.ylabel('hUMAP_y', fontsize=12)

# Create a clean legend
plt.legend(
    loc='center left',
    bbox_to_anchor=(1.01, 0.5), # Place legend outside the plot
    markerscale=5,             # Increase marker size in legend for clarity
    title="Gene Clusters",
    fontsize=10,
    title_fontsize=12
)
plt.grid(False)
plt.gca().set_aspect('equal', adjustable='box') # Ensure equal aspect ratio for UMAP
plt.tight_layout()

# Save the figure
output_filename = 'gene_cluster_umap_plot.png'
plt.savefig(output_filename)

print(f"Visualization saved to '{output_filename}'")
