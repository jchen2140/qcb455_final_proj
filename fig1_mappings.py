# Written by Yubi

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

# --- 1. DEFINE THE GO TERM MAPPING (RAW INPUT) ---
# REMOVED GENE COUNTS FROM THE VALUES (DISPLAY NAMES)
GO_TERM_MAPPING = {
    'biological process:regulation of defense response to virus by virus': 'Interferon signaling',
    'cellular component:u4/u6 x u5 tri-snrnp complex': 'Spliceosome',
    'cellular component:peroxisomal membrane': 'Peroxisome',
    'molecular function:telomeric DNA binding': 'Telomere',
    'biological process:tRNA processing': 'rRNA/ncRNA processing',
    'cellular component:nucleosome': 'Nucleosome',
    'biological process:canonical glycolysis': 'Glycolysis',
    'biological process:transcription by rna polymerase iii': 'RNA polymerase III',
    'biological process:rna 3\'-end processing': 'rRNA/ncRNA processing',
    'biological process:clathrin-dependent endocytosis': 'Endocytosis',
    'biological process:atp hydrolysis coupled proton transport': 'Endosome',
    'biological process:intra-golgi vesicle-mediated transport': 'ER-to-Golgi body transport',
    'biological process:atp-dependent chromatin remodeling': 'Chromatin remodeling',
    'biological process:regulation of transcription from rna polymerase ii promoter in response to hypoxia': 'Mediator',
    'biological process:gpi anchor biosynthetic process': 'GPI synthesis',
    'biological process:positive regulation of dna repair': 'DNA repair',
    'cellular component:cenp-a containing nucleosome assembly': 'Centromere',
    'biological process:interstrand cross-link repair': 'DNA repair',
    'cellular component:mediator complex': 'Mediator',
    'biological process:viral budding via host escrt complex': 'Endocytosis',
    'biological process:dna replication initiation': 'DNA replication',
    'molecular function:snorna binding': 'rRNA/ncRNA processing',
    'biological process:glycogen metabolic process': 'Glycolysis',
    'molecular function:telomerase rna binding': 'Telomere',
    'biological process:ceramide biosynthetic process': 'Lipid biosynthesis',
    # YUBI ADDED MORE annotations
    'biological process:positive regulation of tor signaling': 'Regulation of mTORC1',
    'biological process:phosphatidylinositol 3-kinase signaling': 'mTORC2/PI3K/AKT signaling',
    'biological process:phagophore assembly site' : 'Autophagy'
}

# YUBI: GO TERM MAPPING including gene counts
'''
GO_TERM_MAPPING_RAW = {
    # Key: GO Term from input (must be lowercase and stripped)
    # Value: Simpler Biological Name for the legend
    'biological process:regulation of defense response to virus by virus (15 genes)': 'Interferon signaling (15 genes)',
    'cellular component:U4/U6 x U5 tri-snRNP complex (13 genes)': 'Spliceosome (13 genes)',
    'cellular component:peroxisomal membrane (28 genes)': 'Peroxisome (28 genes)',
    'molecular function:telomeric DNA binding (10 genes)': 'Telomere (10 genes)',
    'biological process:tRNA processing (17 genes)': 'rRNA/ncRNA processing (17 genes)',
    'cellular component:nucleosome (63 genes)': 'Nucleosome (63 genes)',
    'biological process:canonical glycolysis (15 genes)': 'Glycolysis (15 genes)',
    'biological process:transcription by rna polymerase iii (21 genes)': 'RNA polymerase III (21 genes)',
    'biological process:rna 3\'-end processing (10 genes)': 'rRNA/ncRNA processing (10 genes)',
    'biological process:clathrin-dependent endocytosis (15 genes)': 'Endocytosis (15 genes)',
    'biological process:atp hydrolysis coupled proton transport (22 genes)': 'Endosome (22 genes)',
    'biological process:intra-golgi vesicle-mediated transport (19 genes)': 'ER-to-Golgi body transport (19 genes)',
    'biological process:atp-dependent chromatin remodeling (25 genes)': 'Chromatin remodeling (25 genes)',
    'biological process:regulation of transcription from rna polymerase ii promoter in response to hypoxia (27 genes)': 'Mediator (27 genes)',
    'biological process:gpi anchor biosynthetic process (32 genes)': 'GPI synthesis (32 genes)',
    'biological process:positive regulation of dna repair (5 genes)': 'DNA repair (5 genes)',
    'cellular component:cenp-a containing nucleosome assembly (20 genes)': 'Centromere (20 genes)',
    'biological process:interstrand cross-link repair (22 genes)': 'DNA repair (22 genes)',
    'cellular component:mediator complex (30 genes)': 'Mediator (30 genes)',
    'biological process:viral budding via host escrt complex (4 genes)': 'Endocytosis (4 genes)',
    'biological process:dna replication initiation (11 genes)': 'DNA replication (11 genes)',
    'molecular function:snorna binding (7 genes)': 'rRNA/ncRNA processing (7 genes)',
    'biological process:glycogen metabolic process (7 genes)': 'Glycolysis (7 genes)',
    'molecular function:telomerase rna binding (9 genes)': 'Telomere (9 genes)',
    'biological process:ceramide biosynthetic process (15 genes)': 'Lipid biosynthesis (15 genes)'
}
'''

# Helper function to remove the gene count parenthesis for robust lookup
def strip_gene_count(term: str) -> str:
    """Removes the '(XX genes)' or '(XX gene)' pattern from a string."""
    # Regex finds '(123 genes)' or '(1 gene)' and replaces it with an empty string
    return re.sub(r'\s*\(\d+\s*genes?\)$', '', term).strip()

# --- 2. PRE-PROCESS MAPPING FOR ROBUST LOOKUP ---
# Keys are now GO terms without counts, values are the final desired display labels.
'''
GO_TERM_MAPPING = {
    strip_gene_count(k): v
    for k, v in GO_TERM_MAPPING_RAW.items()
}
'''


# 3. Define the file path
file_path = 'vizdf.tsv'

print("Loading data...")

# 4. Load the TSV file
df = pd.read_csv(file_path, sep='\t')


print("Cleaning data...")

# 5. Clean and combine cluster columns
# The columns are cleaned by converting to lowercase string and stripping whitespace
df['roarke_clusters'] = df['roarke_clusters'].fillna('').astype(str).str.lower().str.strip()
df['top_enriched_clusters'] = df['top_enriched_clusters'].fillna('').astype(str).str.lower().str.strip()

# Function to prioritize cluster annotations AND apply the relabeling
def get_final_display_label(row, mapping):
    """
    Prioritizes roarke_clusters, then top_enriched_clusters.
    If a valid term is found, the gene count is stripped, and the term is checked
    against the robust mapping.
    """
    # Markers indicating null or explicit 'unknown' status
    null_markers = ['nan', 'unknown', '']

    roarke = row['roarke_clusters']
    enriched = row['top_enriched_clusters']

    # 1. Determine the best available GO term (already lowercase/stripped)
    best_term = None
    if roarke not in null_markers:
        best_term = roarke
    elif enriched not in null_markers:
        best_term = enriched

    # 2. Check for mapping or default to 'Unknown'
    if best_term:
        # Strip the dynamic gene count before checking against the robust mapping keys
        clean_key = strip_gene_count(best_term)
        if clean_key in mapping:
            # Use the simplified biological name from the mapping value
            return mapping[clean_key]
        else:
            # This term was a valid GO term but not in the mapping dictionary
            return 'Unknown'
    else:
        # This covers cases where both clusters are null/unknown
        return 'Unknown'

print("Annotating and relabeling clusters...")
# Apply the function to create the final display label column
df['Display_Label'] = df.apply(lambda row: get_final_display_label(row, GO_TERM_MAPPING), axis=1)

# Separate data into Known (Relabeled) and Unknown (Final 'Unknown' status) for plotting
df_known = df[df['Display_Label'] != 'Unknown'].copy()
df_unknown = df[df['Display_Label'] == 'Unknown'].copy()

known_labels = df_known['Display_Label'].unique()
print(f"Identified {len(known_labels)} relabeled biological processes for annotation.")
print(f"Genes labeled as 'Unknown': {len(df_unknown)}")


print("Defining color palette...")
# 6. Define Color Palette
num_known_labels = len(known_labels)
color_map = {}
color_map['Unknown'] = '#CCCCCC'

if num_known_labels > 0:
    palette_name = "tab10" if num_known_labels <= 10 else "husl"
    known_colors = sns.color_palette(palette_name, n_colors=num_known_labels).as_hex()
    for label, color in zip(known_labels, known_colors):
        color_map[label] = color

print("Creating scatter plot...")
# 7. Create the Scatter Plot
# Increased figure width to accommodate the legend
plt.figure(figsize=(15, 10)) # Increased width from 12 to 15

# 8. Plot the 'Unknown' genes first (as background)
plt.scatter(
    df_unknown['hUMAP_x'],
    df_unknown['hUMAP_y'],
    color=color_map['Unknown'],
    s=1,
    alpha=0.2,
    label=f'Unknown ({len(df_unknown)} genes)'
)

# YUBI NEW CODE: WRITES CLUSTER LABELS NEXT TO CLUSTER ON THE PLOT

# MANUAL OVERRIDE FOR LABEL POSITIONS
# Add labels here if they overlap. Coordinates are (x, y) for the anchor point of the text.
# ha='right', va='bottom' means the text box will expand to the top-left from this (x,y).
MANUAL_LABEL_POSITIONS = {
    'Mediator': (-6, 5) # YUBI: adjust new position for Mediator to prevent overlap
}

print("Plotting and annotating known clusters...")

# 9. Plot and Annotate the relabeled clusters
for label in known_labels:
    subset = df_known[df_known['Display_Label'] == label]
    color = color_map[label]
    
    # 10a. Plot the points
    plt.scatter(
        subset['hUMAP_x'],
        subset['hUMAP_y'],
        color=color,
        s=20,
        alpha=0.8,
        edgecolors='black',
        linewidths=0.5,
    )
    
    # 10b. Determine annotation position: Use manual override if available, else calculate.
    if label in MANUAL_LABEL_POSITIONS:
        annotate_x, annotate_y = MANUAL_LABEL_POSITIONS[label]
        # For manually placed labels, ha/va might need to be adjusted with the position
        # For consistency, we'll keep the default ha='right', va='bottom' for now,
        # and assume the (x,y) in MANUAL_LABEL_POSITIONS is the bottom-right corner of the label.
    else:
        # Calculate median for the center
        annotate_x = subset['hUMAP_x'].median() 
        annotate_y = subset['hUMAP_y'].median() 
        
        # Apply offset to the top-left of the center
        annotate_x -= 0.5 # Negative offset to push LEFT
        annotate_y += 0.5 # Positive offset to push UP
    
    # 10c. Place the label
    plt.text(
        annotate_x, 
        annotate_y, 
        label,
        fontsize=6,
        fontweight='bold',
        color=color,
        bbox=dict(facecolor='white', alpha=0.9, edgecolor=color, linewidth=1, boxstyle="round,pad=0.15"), 
        ha='right',        
        va='bottom'        
    )

# 10. Add plot details
plt.title('Gene Co-essentiality UMAP with Direct Cluster Annotation', fontsize=16)
plt.xlabel('hUMAP_x', fontsize=12)
plt.ylabel('hUMAP_y', fontsize=12)


plt.grid(False)
plt.gca().set_aspect('equal', adjustable='box') 
plt.tight_layout() 

# Save the figure
output_filename = 'gene_cluster_umap_edge_annotation_more.png'
plt.savefig(output_filename)

print(f"Visualization saved to '{output_filename}'")

# OLD CODE: WRITES CLUSTER LABELS IN LEGEND

'''
print("Plotting known clusters...")
# 9. Plot the relabeled clusters on top
for label in known_labels:
    subset = df_known[df_known['Display_Label'] == label]
    plt.scatter(
        subset['hUMAP_x'],
        subset['hUMAP_y'],
        color=color_map[label],
        s=20,
        alpha=0.8,
        edgecolors='black',
        linewidths=0.5,
        label=f'{label}'
    )

# 10. Add plot details and legend
plt.title('Gene Co-essentiality UMAP with Relabeled Biological Processes', fontsize=16)
plt.xlabel('hUMAP_x', fontsize=12)
plt.ylabel('hUMAP_y', fontsize=12)

# Create a clean legend
plt.legend(
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),
    markerscale=5,
    title="Biological Process Clusters",
    fontsize=10,
    title_fontsize=12
)
plt.grid(False)
plt.gca().set_aspect('equal', adjustable='box')

# Adjust layout: Manually set the right margin to prevent legend cutoff
# This line replaces plt.tight_layout() for more control
plt.subplots_adjust(right=0.8) # Adjust this value as needed. 0.8 leaves 20% of figure width for legend

# Save the figure
output_filename = 'gene_cluster_umap_mapping_fixed_layout.png'
plt.savefig(output_filename)

print(f"Visualization saved to '{output_filename}'")
'''
