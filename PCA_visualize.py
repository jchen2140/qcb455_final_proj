import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
import warnings
from load_screens import load_screens
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# YUBI: this file creates a UMAP plot of all of the genes using the PCA reduced dimensionality representation of the genes
# It highlights the tightest modules found by using the Spearman correlation analysis I ran previously

# Suppress UMAP Numba warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# --- CONFIGURATION ---
N_COMPONENTS = 2
N_NEIGHBORS = 30  # Adjust based on data size (often 15 to 50 is a good range)
MIN_DIST = 0.1     # Adjust for tightness of clusters (0.0 for tight, 0.99 for loose)
# YUBI: set this to desired number of visualized clusters
N_VISUAL_CLUSTERS = 25 # Set to the number of tightest clusters to visualize
BACKGROUND_CATEGORY = 'Other Genes' # YUBI: make this a more descripive name

# --- 1. DEFINE THE GO TERM MAPPING (RAW INPUT) ---
# REMOVED GENE COUNTS FROM THE VALUES (DISPLAY NAMES)
biological_mapping = {
    144: "Peroxisome",
    2: "Mitochondrial Translational Machinery",
    8: "Complex V",
    89: "SAGA Complex",
    4: "Complex I",
    197: "Ubiquitin-proteasome",
    182: "Antimicrobial Peptides",
    146: "N-linked glycosylation",
    11: "DNA Damage Response",
    29: "tRNA Modification",
    156: "ER Membrane Protein",
    145: "UFM1 conjugation",
    5: "Complexes III, IV",
    122: "BAF Complex",
    298: "Nucleotide Metabolism",
    3: "Mitochondrial Translational Machinery",
    91: "Focal Adhesion",
    38: "Unknown",
    98: "HOPS",
    90: "Actin Cytoskeleton",
    6: "Pyruvate processing",
    117: "DNA Damage Response",
    97: "mTORC1 signaling",
    92: "Selenocysteine Incorporation",
    125: "TGF-beta signal"
}


# ----------------------------------------------------------------------------------
# 1. LOAD DATA AND PRE-PROCESSING
# NOTE: The script assumes your gene names are in the index/first column of 
# loadings_matrix.csv and that the cluster file is formatted as described.
# Load batch-corrected screens

print("loading screens...")

screens = load_screens()

# transpose data so that the genes are the columns
df_pca_input = screens.T

df_shape = df_pca_input.shape
print(f"Shape of the screens DataFrame: {df_shape}")


print("normalize data...")
# Normalize the data
# YUBI: I'm not sure how important this is
scaler = StandardScaler()
# Note: We fit and transform the data
scaled_data = scaler.fit_transform(df_pca_input)

# Run targeted PCA and calculate spearman correlation

# PCA with only top 232 principal components
N_PCS = 232
print(f"\nRunning PCA to extract the top {N_PCS} components...")

# Run PCA to get the specified number of components (232)
pca_reduced = PCA(n_components=N_PCS)
# Fit and transform is not needed here, we just need the components/loadings
pca_reduced.fit(scaled_data)

# Explained Variance for the top 232 components
total_explained_variance = np.sum(pca_reduced.explained_variance_ratio_)
print(f"The first {N_PCS} PCs explain {total_explained_variance:.2%} of the total variance.")
# Should be 80%


# EXTRACT LOADINGS (THE REDUCED REPRESENTATION OF GENES)
print("extract loadings...")

# The 'components_' attribute is the Loadings matrix (PCs x Genes).
# We transpose it to make the Genes the index (rows) and the PCs the columns.
# Shape is now (17634 Genes x 232 PCs). This is the reduced dimension
# representation of the *genes*.
loadings_matrix_df = pd.DataFrame(
    pca_reduced.components_.T,
    index=df_pca_input.columns,
    columns=[f'PC {i+1}' for i in range(N_PCS)]
)

print(f"Loadings matrix dataframe shape (Genes x PCs): {loadings_matrix_df.shape}")

print("Loading clustered data...")

try:
    # Load the cluster assignments (Module_ID and Gene_Names)
    clusters_df = pd.read_csv('TOM_clusters_0.99.csv')
    # YUBI: Load the Spearman Correlation Matrix for making umap projection
    # using this means I have to run the code on Adroit because that is where the csv file is
    spearman_corr_df = pd.read_csv('gene_spearman_corr_matrix.csv', index_col=0)
    
except FileNotFoundError as e:
    print(f"Error: Could not find file {e.filename}. Please check your file paths.")
    

# 2. CLEAN AND MERGE DATA

print("clean and merge data...")

# Unpivot the cluster dataframe from wide format to long format
# This creates a simple mapping: Gene_Name -> Module_ID
cluster_map_df = clusters_df.melt(
    id_vars=['Module_ID'],
    value_vars=clusters_df.columns[1:],
    value_name='Gene_Name'
).dropna(subset=['Gene_Name'])
cluster_map_df = cluster_map_df[['Gene_Name', 'Module_ID']].set_index('Gene_Name')

# Merge the loadings matrix with the cluster mapping
# This aligns the 232D vector with its cluster ID
# We use an inner join to only keep genes that are present in both files (clustered genes)
merged_df = loadings_matrix_df.merge(
    cluster_map_df,
    left_index=True,
    right_index=True,
    how='inner'
)

# Explicitly convert Module_ID to string to prevent mismatch errors during filtering
merged_df['Module_ID'] = merged_df['Module_ID'].astype(str)

# YUBI NEW
# Map the Module_ID (number) to the Biological Function (name)
merged_df['Biological_Function'] = merged_df['Module_ID'].astype(int).map(biological_mapping)

# Define feature columns once for reuse
feature_cols = [col for col in merged_df.columns if col not in ['Module_ID', 'Biological_Function']]

# ----------------------------------------------------------------------------------
# 5. APPLY UMAP DIMENSIONALITY REDUCTION (Visualizing Gene Similarity)
# ----------------------------------------------------------------------------------

# Create a copy of the entire merged dataset
full_umap_df = merged_df.copy()

# Run UMAP on ALL genes
X = full_umap_df[feature_cols].values

print(f"\nRunning UMAP on all {len(X)} genes ({len(loadings_matrix_df)} total) from 232D to {N_COMPONENTS}D...")

# YUBI: using spearman correlation matrix as the pre-computed metric in the umap projection

# Use the full Spearman correlation matrix (N x N)
if spearman_corr_df is not None:
    # 1. Convert Correlation to Distance: 1 - Correlation <--- ADDED
    D = 1 - spearman_corr_df.values
    # 2. UMAP initialization: Set metric to 'precomputed'
    reducer = umap.UMAP(
        n_components=N_COMPONENTS,
        n_neighbors=N_NEIGHBORS,
        min_dist=MIN_DIST,
        metric='precomputed',
        random_state=42
    )
    # 3. Fit and Transform using the distance matrix (D)
    umap_embedding = reducer.fit_transform(D) 

else:
    # Fallback to default Euclidean distance on the 232D vectors
    print("Falling back to Euclidean distance on 232D vectors.")
    X = full_umap_df[feature_cols].values
    reducer = umap.UMAP(
        n_components=N_COMPONENTS,
        n_neighbors=N_NEIGHBORS,
        min_dist=MIN_DIST,
        random_state=42
    )
    umap_embedding = reducer.fit_transform(X)

# Add the 2D coordinates back to the plotting DataFrame
full_umap_df['UMAP_1'] = umap_embedding[:, 0]
full_umap_df['UMAP_2'] = umap_embedding[:, 1]

# ----------------------------------------------------------------------------------
# 6. WITHIN-MODULE TIGHTNESS ANALYSIS IN UMAP SPACE (Find the 5 tightest modules)
# ----------------------------------------------------------------------------------

print("computing within module tightness in UMAP 2D space...")

# define all_modules
all_modules = full_umap_df['Biological_Function'].unique().tolist()

# Columns containing the UMAP coordinates
umap_cols = ['UMAP_1', 'UMAP_2']

# Calculate UMAP 2D Centroids (The average 2D vector for all genes in a module)
umap_centroid_df = full_umap_df.groupby('Biological_Function')[umap_cols].mean()

print(f"\nCalculating 2D tightness scores for all {len(all_modules)} modules...")

umap_module_tightness = {}
for module_name in all_modules:
    # 1. Get genes belonging to this module
    module_genes = full_umap_df[full_umap_df['Biological_Function'] == module_name]
    
    if module_genes.empty:
        continue
        
    # 2. Get the features for these genes (2D UMAP vectors)
    X_umap = module_genes[umap_cols].values
    
    # 3. Get the 2D centroid for this module
    C_umap = umap_centroid_df.loc[module_name].values.reshape(1, -1)
    
    # 4. Calculate Euclidean distances from all genes in the module to the 2D centroid
    # np.linalg.norm is used for efficient Euclidean distance calculation
    # This correctly measures the visual spread in the 2D plot.
    distances_umap = np.linalg.norm(X_umap - C_umap, axis=1)
    
    # 5. Store the mean distance (lower mean distance = tighter module)
    umap_module_tightness[module_name] = np.mean(distances_umap)

# Convert to Series and sort ascending (lowest distance is tightest)
tightness_series_umap = pd.Series(umap_module_tightness)
# Re-assign the tightest modules based on the NEW 2D tightness scores
tightest_modules_series = tightness_series_umap.sort_values(ascending=True).head(N_VISUAL_CLUSTERS)
tightest_module_names = tightest_modules_series.index.tolist()

print(f"\n--- Top {N_VISUAL_CLUSTERS} Tightest Modules (Lowest Mean **UMAP 2D Euclidean Distance** to Centroid) ---")
# The distance value indicates the average spread of genes in the 2D plot.
print(tightest_modules_series.to_string())

# ----------------------------------------------------------------------------------
# 7. PREPARE DATA FOR FULL VISUALIZATION (Highlighting Tightest Clusters)
# ----------------------------------------------------------------------------------

print("prepating data for full visualization...")
# full_umap_df already contains UMAP_1 and UMAP_2 columns

# Create a new categorical column for plotting, using the NEWLY CALCULATED tightest_module_names
full_umap_df['Plot_Category'] = np.where(
    full_umap_df['Biological_Function'].isin(tightest_module_names),
    full_umap_df['Biological_Function'],
    BACKGROUND_CATEGORY
)

# Define the order of categories for plotting. Placing 'Other Genes' last ensures it's 
# plotted underneath the colored clusters when plotted in layers.
category_order = tightest_module_names + [BACKGROUND_CATEGORY]
full_umap_df['Plot_Category'] = pd.Categorical(full_umap_df['Plot_Category'], categories=category_order)

# --- DATA VALIDATION CHECK ---
highlighted_gene_count = full_umap_df[full_umap_df['Plot_Category'] != BACKGROUND_CATEGORY].shape[0]
print(f"\n[Validation] Total genes in the {N_VISUAL_CLUSTERS} highlighted modules: {highlighted_gene_count}")
# -----------------------------

# Define the order of categories for plotting. Placing 'Other Genes' last ensures it's 
# plotted underneath the colored clusters when plotted in layers.
category_order = tightest_module_names + [BACKGROUND_CATEGORY]
full_umap_df['Plot_Category'] = pd.Categorical(full_umap_df['Plot_Category'], categories=category_order)

# YUBI: new code to calculate centroids
annotation_centroids = full_umap_df[full_umap_df['Plot_Category'] != BACKGROUND_CATEGORY].groupby('Biological_Function')[umap_cols].mean()

# Manual Override Dictionary for Label Placement
# Define specific (UMAP_1, UMAP_2) coordinates for labels that overlap.
# Coordinates are chosen based on the provided plot to avoid overlap with clusters/labels.
MANUAL_LABEL_POSITIONS = {
    'Complex I': (4.5, 5.35), 
    'ER Membrane Protein': (5.2, 4.2),
    'N-linked glycosylation': (5.5, 4.3), 
    'Actin Cytoskeleton': (6.3, 4.25), 
    'HOPS': (5.9, 4.4),
    'mTORC1 signaling': (6.1, 4.0),
    'UFM1 conjugation': (5.3, 3.4),
    'tRNA Modification': (4.7, 3.3),
}

# ----------------------------------------------------------------------------------
# 6. CREATE THE FULL 2D VISUALIZATION
# ----------------------------------------------------------------------------------

print("Generating 2D plot...")

# Define a custom color palette: 20 colors for the tightest modules, 
# and a light gray for the 'Other Genes' background.
# Using 'tab20' for 20 distinct, aesthetically pleasing colors.
# Using 'husl' for more than 20 clusters
highlight_colors = sns.color_palette("husl", N_VISUAL_CLUSTERS).as_hex()
background_color = '#CCCCCC' # Light Gray

custom_palette = dict(zip(tightest_module_names, highlight_colors))
custom_palette[BACKGROUND_CATEGORY] = background_color


sns.set_style("whitegrid")
plt.figure(figsize=(14, 12))

# FIX: Create filtered background data and convert Plot_Category to object type.
# This strips the categorical metadata (the 20 module IDs) so Seaborn only sees
# the 'Other Genes' level, resolving the ValueError.
background_data = full_umap_df[full_umap_df['Plot_Category'] == BACKGROUND_CATEGORY].copy()
background_data['Plot_Category'] = background_data['Plot_Category'].astype(object)

# Plot the 'Other Genes' first to form the background
sns.scatterplot(
    x='UMAP_1', 
    y='UMAP_2', 
    hue='Plot_Category', 
    data=background_data, # Using the fixed background_data
    palette={BACKGROUND_CATEGORY: background_color},
    s=10,    # <<< DECREASED SIZE for less background clutter
    alpha=0.3, # <<< YUBI: make background points darker
    linewidth=0,
    legend=False # Do not add 'Other Genes' to the legend yet
)

# Plot the highlighted modules over the background
scatter = sns.scatterplot(
    x='UMAP_1', 
    y='UMAP_2', 
    hue='Plot_Category', 
    hue_order=category_order[:-1], # Only plot the highlighted modules in this layer
    data=full_umap_df[full_umap_df['Plot_Category'] != BACKGROUND_CATEGORY], 
    palette=custom_palette, 
    s=30,    # <<< INCREASED SIZE for visibility
    alpha=1.0, # <<< INCREASED ALPHA (fully opaque)
    linewidth=0,
    # YUBI: remove legend once I get it working
    legend=False
)

# 3. Add Direct Cluster Annotation 
# We iterate over the tightest module names and place the text labels.
for module_name in tightest_module_names:
    
    # Exclude plottingUnknown modules for now
    if module_name == 'Unknown':
        continue
    
    # 1. Check for manual override position
    if module_name in MANUAL_LABEL_POSITIONS:
        annotate_x, annotate_y = MANUAL_LABEL_POSITIONS[module_name]
        # For manual positions, use standard alignment to control placement precisely
        ha_align = 'center'
        va_align = 'center'
    else:
        # 2. Use calculated centroid position with an offset
        # Get the calculated centroid coordinates
        annotate_x = annotation_centroids.loc[module_name, 'UMAP_1']
        annotate_y = annotation_centroids.loc[module_name, 'UMAP_2']
        
        # Apply a small offset from the centroid to the top-right
        offset_x = 0.05
        offset_y = 0.05
        annotate_x += offset_x
        annotate_y += offset_y
        
        # Use left/bottom alignment for automatic placements
        ha_align = 'left'
        va_align = 'bottom'
    
    # Get the color for the label
    color = custom_palette[module_name]

    # Adjust the position slightly to offset the label from the center point
    # We will offset it to the top-right of the centroid for readability
    offset_x = 0.05
    offset_y = 0.05
    
    # Place the label
    plt.text(
        annotate_x + offset_x,  # Offset X position
        annotate_y + offset_y,  # Offset Y position
        # Use line breaks for long function names to prevent overlap
        module_name.replace(' ', '\n'), 
        fontsize=7,
        fontweight='bold',
        color='black', # Use black text for better contrast against the background
        # Add a white box around the text for maximum readability over cluster points
        bbox=dict(
            facecolor='white', 
            alpha=0.8, 
            edgecolor=color, # Edge color matches the cluster color
            linewidth=1, 
            boxstyle="round,pad=0.2"
        ), 
        ha='left',        
        va='bottom'        
    )

# Configure the plot appearance
plt.title(
    f'Top {N_VISUAL_CLUSTERS} Tightest Modules from Hierarchical Clustering', 
    fontsize=20, 
    fontweight='bold'
)
plt.xlabel('UMAP Component 1', fontsize=14)
plt.ylabel('UMAP Component 2', fontsize=14)


# ----------------------------------------------------------------------------------
# MODIFICATION: SET AXIS LIMITS TO ZOOM IN
# ----------------------------------------------------------------------------------

# Set the X-axis limit (UMAP Component 1) up to 8
plt.xlim(
    left=3,
    right=8
)

# Set the Y-axis limit (UMAP Component 2) up to 7
plt.ylim(
    bottom=2, 
    top=7  
)

'''
# YUBI: remove legend once I get it working
# Customize the legend
# The legend will now show all highlighted modules
legend = scatter.legend(
    title='Module Names', 
    fontsize=8, # Smaller font size to accommodate 20 entries
    title_fontsize=12, 
    bbox_to_anchor=(1.05, 1), 
    loc='upper left',
    markerscale=2,
    ncol=1 # Ensure the legend is in a single column
)
legend.get_title().set_fontweight('bold')
'''

# Ensure layout is tight to prevent legend from being cut off
plt.tight_layout()

# Save the plot
# NOTE YUBI: save plot once I confirm results
plt.savefig('HC_plot_top25_biological.png', dpi=300)

plt.show()

print("\nScript execution complete. Check 'HC_plot_top25_biological.png' for the visualization.")