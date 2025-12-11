import numpy as np
from load_screens import load_screens
from scipy.special import stdtr
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram # Added for hierarchical clustering
from scipy.cluster.hierarchy import fcluster
from typing import List, Tuple 
import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.stats import linregress

# This file computes the reduced dimensionality representation of the genes using PCA
# It computes a Spearman correlation matrix over all genes to create modules of genes
# I am clustering with TOM-based hierarichal clustering and WGCNA core functions


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
loadings_matrix = pd.DataFrame(
    pca_reduced.components_.T,
    index=df_pca_input.columns,
    columns=[f'PC {i+1}' for i in range(N_PCS)]
)

print(f"Loadings matrix shape (Genes x PCs): {loadings_matrix.shape}")


'''
# CALCULATE SPEARMAN CORRELATION BETWEEN GENE LOADINGS

# We now calculate the correlation between the rows of the loadings matrix.
# This finds the Spearman correlation between the 232-dimensional loading
# vectors for every pair of the 17,634 genes.

print("\nCalculating the full 17634 x 17634 Spearman correlation matrix...")
# Note: This operation can be time-consuming due to the size (17634 * 17634)
# Ensure your environment has sufficient memory for this calculation.
gene_correlation_matrix = loadings_matrix.T.corr(method='spearman')

print(f"Gene-Gene Correlation Matrix Shape: {gene_correlation_matrix.shape}")
'''


# WGCNA Core Functions
def calculate_signed_adjacency(corr_matrix: pd.DataFrame, beta: int) -> pd.DataFrame:
    """
    Calculates the Signed Adjacency Matrix A based on correlation rho and soft-threshold beta.
    A_ij = ((1 + rho_ij) / 2) ^ beta
    This mapping ensures 0 < A < 1 and only strongly POSITIVELY correlated genes have high adjacency.
    """
    print(f"Calculating Adjacency Matrix using soft-threshold power beta={beta}...")
    # Convert correlation [-1, 1] to Similarity [0, 1]
    similarity_matrix = (corr_matrix + 1) / 2
    # Apply soft-thresholding
    adjacency_matrix = similarity_matrix ** beta
    # Set diagonal to 1
    np.fill_diagonal(adjacency_matrix.values, 1)
    return adjacency_matrix

def calculate_tom_distance(adjacency_matrix: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the Topological Overlap Matrix (TOM) and TOM-based distance (1 - TOM).
    """
    print("Calculating Topological Overlap Matrix (TOM)...")
    A = adjacency_matrix.values
    
    # Calculate k (Node connectivity)
    k = A.sum(axis=1) - 1 # sum of all connections excluding the self-connection (A_ii=1)
    
    # Calculate TOM numerator (L) and denominator
    # L_ij = sum_u(A_iu * A_uj) + A_ij
    L = A @ A + A # Matrix multiplication A@A is sum_u(A_iu * A_uj)

    # Denominator: min(k_i, k_j) + 1 - A_ij
    # Create a matrix where M_ij = min(k_i, k_j)
    k_min_matrix = np.minimum.outer(k, k)
    
    # TOM_ij = L_ij / (min(k_i, k_j) + 1 - A_ij)
    # The term +1 in the denominator is critical for TOM to be between 0 and 1
    # Note: WGCNA documentation sometimes uses slightly different denominator forms, 
    # but L_ij / (k_min_matrix - A + 1) is a common robust formulation.
    TOM = L / (k_min_matrix + 1 - A)
    np.fill_diagonal(TOM, 1) # Ensure diagonal is exactly 1

    print("Converting TOM to condensed distance vector (1 - TOM)...")
    # Convert the full TOM matrix into a condensed (1D) distance vector
    tom_distance_vector = 1 - TOM[np.triu_indices(TOM.shape[0], k=1)]
    
    return tom_distance_vector, TOM

def calculate_mean_connectivity(adjacency_matrix: pd.DataFrame) -> float:
    """
    Calculates the mean network connectivity (k_bar) from the Adjacency Matrix.
    k_i = sum_j (A_ij) (excluding the self-connection, which is 1)
    k_bar = mean(k_i)
    """
    # Calculate the degree (connectivity) for each node (gene)
    # The sum of each row (or column) of A is k_i
    # We subtract 1 from the sum because the diagonal A_ii is set to 1
    k_i = adjacency_matrix.values.sum(axis=1) - 1 
    
    # Calculate the mean connectivity
    mean_k = np.mean(k_i)
    return mean_k

def calculate_scale_free_fit(adjacency_matrix: pd.DataFrame) -> float:
    """
    Calculates the R^2 of the linear fit for the scale-free topology.
    This involves fitting a line to log(P(k)) vs log(k).
    """
    # 1. Calculate connectivity (k_i)
    k_i = adjacency_matrix.values.sum(axis=1) - 1
    
    # 2. Get the distribution P(k) (Frequency of each k value)
    # np.histogram counts how many genes have a certain connectivity k
    hist, bin_edges = np.histogram(k_i, bins=50) # Use a reasonable number of bins
    
    # 3. Use the mid-points of the bins as k values for fitting
    k_mid = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Filter out bins with zero frequency (log(0) is undefined)
    non_zero_freq = hist > 0
    k_mid = k_mid[non_zero_freq]
    P_k = hist[non_zero_freq] / len(k_i) # Normalize to get P(k)
    
    # 4. Perform the linear fit on log-log plot: log(P(k)) = -gamma * log(k) + constant
    log_k = np.log(k_mid)
    log_P_k = np.log(P_k)
    
    # Perform linear regression: linregress returns (slope, intercept, r_value, p_value, stderr)
    slope, intercept, r_value, p_value, stderr = linregress(log_k, log_P_k)
    
    # R^2 is the square of the correlation coefficient (r_value)
    # Note: WGCNA documentation specifically uses the R^2 of the linear model fit, 
    # but uses the sign of the slope (negative) to indicate a valid fit.
    R_squared = r_value**2
    
    # A true scale-free fit requires a NEGATIVE slope (gamma > 0)
    # WGCNA R^2 is typically the square of the correlation coefficient.
    return R_squared

# Crucial Step for Gene-Gene Correlation:
# We transpose the matrix so that Genes are COLUMNS and PCs are ROWS, 
# allowing pd.corr() to calculate gene-gene correlation.
data_df = loadings_matrix.T 

# NOTE YUBI: edit parameters to set sample size
# Adroit cluster runs full data
N_CLUSTER_SAMPLE = 17634
sample_data = data_df.iloc[:, :N_CLUSTER_SAMPLE]
print("\nCalculating Sample Spearman correlation matrix...")
# Correlation matrix (Genes x Genes)
sample_corr_matrix = sample_data.corr(method='spearman')
print(f"Sample Correlation Matrix Shape: {sample_corr_matrix.shape}")

# 2. WGCNA Steps: Adjacency and TOM Distance
# NOTE: In a real WGCNA analysis, you would systematically choose beta (e.g., between 6 and 20)
# to satisfy the scale-free topology criterion.
SOFT_THRESHOLD_BETA = 6 
ADJ_matrix = calculate_signed_adjacency(sample_corr_matrix, SOFT_THRESHOLD_BETA)
# Calculate metrics
R2 = calculate_scale_free_fit(ADJ_matrix)
mean_k = calculate_mean_connectivity(ADJ_matrix)
    
print(f" R^2={R2:.4f}, Mean k={mean_k:.4f}")

TOM_distance_vector, TOM_matrix = calculate_tom_distance(ADJ_matrix)

# 3. Hierarchical Clustering (using TOM Distance)
# Linkage method remains 'ward' as requested
print("Performing Hierarchical Clustering with 'ward' linkage on TOM Distance...")
linkage_matrix = linkage(TOM_distance_vector, method='ward')

# 4. Dynamic Tree Cut (Module Assignment)
# A common starting point for TOM-based clustering is a cut height between 0.9 and 0.99.
# YUBI NOTE: raise cut threshold to 0.99 to get larger clusters in the full dataset
CUT_THRESHOLD = 0.99
print(f"Assigning modules using Dynamic Tree Cut simulation (fcluster) at height {CUT_THRESHOLD}...")
# fcluster returns an array of cluster labels (1, 2, 3, ...)
module_assignments = fcluster(linkage_matrix, CUT_THRESHOLD, criterion='distance')
num_modules = len(np.unique(module_assignments))
print(f"Detected {num_modules} modules in the sample of {N_CLUSTER_SAMPLE} genes.")


# --- Visualization: Dendrogram with Module Colors ---

plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(18, 8))

# Define colors for the modules
# Use a color map and cycle through colors for modules
# NOTE YUBI: this part of code is currently not working, dendrogram is not being colored
module_colors = plt.colormaps.get_cmap('Spectral')
color_list = [module_colors(i-1) for i in module_assignments]

# Plot the dendrogram with colors based on fcluster assignment
dendrogram(
    linkage_matrix,
    orientation='top',
    labels=sample_corr_matrix.index.tolist(),
    leaf_rotation=90,
    leaf_font_size=10,
    color_threshold=CUT_THRESHOLD, # Draw a line at the cut height
    above_threshold_color='grey',
    link_color_func=lambda x: color_list[x] if x < len(color_list) else 'grey', # Color the leaves
    show_leaf_counts=True,
)

plt.title(f'WGCNA-like Module Detection (TOM Distance, Ward Linkage, $\\beta={SOFT_THRESHOLD_BETA}$)', fontsize=18)
plt.ylabel('TOM Distance (1 - TOM)', fontsize=14)
plt.xlabel(f'Genes (Sample of {N_CLUSTER_SAMPLE})', fontsize=14)
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.axhline(y=CUT_THRESHOLD, color='r', linestyle='--', label=f'Module Cut Threshold ({CUT_THRESHOLD})')
plt.legend()
plt.tight_layout()

# SAVE DENDROGRAM TO PNG
DENDROGRAM_FILENAME = f'dendrogram_TOM_{CUT_THRESHOLD}.png'
plt.savefig(DENDROGRAM_FILENAME)


# NOTE YUBI: do not show plot on adroit cluster
# plt.show()
print(f"\nThe visualization now shows modules derived from the TOM distance and a clear cut point.")
print(f"Dendrogram saved to {DENDROGRAM_FILENAME}.")

# --- Save Results (Updated to include Module ID and save module lists) ---

# Create a DataFrame of the module assignment for the sample
module_df = pd.DataFrame({
    'Gene': sample_corr_matrix.index,
    'Module_ID': module_assignments
})


# SAVE CLUSTERED GENES FILE (TOM_clusters.csv)
print("Saving TOM_clusters.csv...")

# Group by Module_ID and aggregate gene names into a list for each module
clustered_genes = module_df.groupby('Module_ID')['Gene'].apply(list).reset_index()

# Rename the columns for clarity
clustered_genes.columns = ['Module_ID', 'Genes']

# Convert the list of genes into a comma-separated string for easy CSV saving/reading
# This makes it easier to save to CSV where each cell is a single value (string).
clustered_genes['Genes'] = clustered_genes['Genes'].apply(lambda x: ', '.join(x))

# Find the maximum number of genes in any module to determine the output columns
max_genes_in_module = clustered_genes['Genes'].apply(lambda x: len(x.split(', ')))
max_genes_count = max_genes_in_module.max()

# Split the comma-separated string back into multiple columns for the requested format 
# (where each row is a cluster and values are genes).
# This creates a sparse DataFrame where each column is a "Gene Slot"
module_lists_df = clustered_genes['Genes'].str.split(', ', expand=True)

# Prepend the Module ID back to the front
module_lists_df.insert(0, 'Module_ID', clustered_genes['Module_ID'])

# Save the final file
CLUSTERS_FILENAME = f'TOM_clusters_{CUT_THRESHOLD}.csv'

# Create the full header list: ['Module_ID', 'Gene_1', 'Gene_2', ...]
full_header = ['Module_ID'] + [f'Gene_{i+1}' for i in range(max_genes_count)]

module_lists_df.to_csv(CLUSTERS_FILENAME, index=False, header=full_header)
print(f"Successfully saved module definitions to {CLUSTERS_FILENAME}. (Total modules: {num_modules})")






'''

# HIERARCHICAL CLUSTERING (MODULE DETECTION)

print("HIERARCHICAL CLUSTERING (MODULE DETECTION)...")

# OPTIMIZATION NOTE: For 17635 genes, the N^2 distance matrix is ~310 million elements.
# This may still exceed available memory and runtime limits, 
# even with optimized scipy functions.
N_CLUSTER_SAMPLE = 50 
print(f"Clustering a sample of the first {N_CLUSTER_SAMPLE} genes for visualization...")

# Extract the sample data (crucial for pdist)
# The data must be in (observations/genes x features/samples) format
sample_data = loadings_matrix.T.iloc[:N_CLUSTER_SAMPLE]

# 1. Calculate the distance directly from the original data (more efficient)
# Use 'correlation' metric in pdist, then convert to the desired 1 - |r| distance
print("calculating condensed distance matrix...")
# pdist with 'correlation' calculates (1 - Pearson correlation).
# Spearman requires the data to be rank-transformed first, then Pearson calculated.
# Since we already have the Spearman correlation matrix (gene_correlation_matrix),
# we must calculate the distance from that matrix efficiently.

# Extract the sample correlation matrix
sample_corr_matrix = gene_correlation_matrix.iloc[:N_CLUSTER_SAMPLE, :N_CLUSTER_SAMPLE]

# 1. Convert correlation to condensed distance: 1D array of 1 - abs(r)
# We use np.triu_indices to get the upper triangle of the correlation matrix 
# and convert it to the required condensed format (1D array)
# This avoids creating the full N x N distance matrix explicitly
upper_tri_indices = np.triu_indices(sample_corr_matrix.shape[0], k=1)
condensed_distance_vector = 1 - np.abs(sample_corr_matrix.values[upper_tri_indices])

# 2. Perform Hierarchical/Linkage Clustering
# The 'linkage' function accepts the condensed distance vector directly (BEST FOR PERFORMANCE)
print("performing hierarchical clustering with condensed distance...")
# 'ward' method is highly recommended for minimizing within-cluster variance.
linkage_matrix = linkage(condensed_distance_vector, method='ward') 


print("plotting dendrogram...")
# 3. Plot the Dendrogram for Module Identification
plt.figure(figsize=(15, 6))
# The truncate_mode and p parameters are useful for large samples
dendrogram(
    linkage_matrix,
    orientation='top',
    labels=sample_corr_matrix.index.tolist(),
    leaf_rotation=90,
    leaf_font_size=8,
    # Uncomment the following lines for visualizing the full 17k genes (if you manage to cluster them)
    truncate_mode='lastp',  # Show only the last 'p' merged clusters
    p=50,                   # Only display the top 50 merges
    show_leaf_counts=True,
    color_threshold=0.7 # Example threshold line for cutting the tree
)
# The y-axis shows the distance at which the clusters were merged.
# A horizontal line (e.g., at a distance of 0.7) represents the 'cut' to define modules.
plt.title(f'Hierarchical Clustering Dendrogram (Ward Linkage, 1 - |Spearman ρ| Distance)', fontsize=16)
plt.ylabel('Distance (1 - |Spearman ρ|)', fontsize=12)
plt.xlabel(f'Genes (Sample of {N_CLUSTER_SAMPLE})', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.axhline(y=0.7, color='r', linestyle='--', label='Module Cut Threshold')
plt.legend()
plt.tight_layout()
plt.show()

# 4. (Optional) Clustermap with Explicit Cluster Assignments
# Once you have a 'cut' (e.g., at distance=0.7), you can use the fcluster function
# to assign cluster labels, and then reorder the heatmap based on those labels.
# This makes the modules visually obvious on the heatmap.


max_d = 0.7 # Distance threshold for cutting
clusters = fcluster(linkage_matrix, max_d, criterion='distance')
# Sort the genes by their assigned cluster
cluster_order = sample_corr_matrix.index[np.argsort(clusters)]

print("plotting clustered heatmap with cluster order...")
sns.clustermap(
    sample_corr_matrix.loc[cluster_order, cluster_order], # Reorder both axes
    row_cluster=False, # Already ordered by cluster
    col_cluster=False, # Already ordered by cluster
    cmap='coolwarm',
    vmin=-1, vmax=1,
    linewidths=0.5,
    figsize=(12, 12)
)
plt.suptitle(f'Heatmap Reordered by Modules (Cut at Distance {max_d})', y=1.02, fontsize=16)
plt.show()

print(f"The Dendrogram is the key visualization for setting the module distance threshold.")
'''

'''

# EXTRACT STRONGLY CORRELATED GENE PAIRS using vectorized numpy approach
def extract_strong_correlations(corr_matrix: pd.DataFrame, threshold: float, top_n: int = 10, return_cond: bool = False):
    """
    Converts a square correlation matrix into a filtered list of strongly correlated pairs.
    
    Args:
        corr_matrix: The square N x N correlation matrix.
        threshold: The minimum absolute correlation (|r|) to be considered "strong".
        top_n: The maximum number of results to display.
        
    Returns:
        A DataFrame of the top N strongest pairs.
    """
    corr_array = corr_matrix.values
    gene_names = corr_matrix.index.tolist()
    
    # 1. Get indices of the upper triangle (excludes diagonal and duplicates)
    upper_tri_indices = np.triu_indices(corr_array.shape[0], k=1)
    
    # 2. Get correlation values from the upper triangle
    correlations = corr_array[upper_tri_indices]
    
    # 3. Filter indices where absolute correlation meets the threshold
    strong_indices_mask = np.abs(correlations) >= threshold
    
    # 4. Filter the correlation values and their original row/column indices
    final_correlations = correlations[strong_indices_mask]
    row_indices = upper_tri_indices[0][strong_indices_mask]
    col_indices = upper_tri_indices[1][strong_indices_mask]
    
    # 5. Create the final results DataFrame
    results_df = pd.DataFrame({
        'Gene A': [gene_names[i] for i in row_indices],
        'Gene B': [gene_names[i] for i in col_indices],
        'Correlation': final_correlations
    })

    # 6. Sort by absolute correlation and take the top N or return all
    # 6. Sort by absolute correlation and take the top N using the 'key' argument
    # This avoids creating and dropping the temporary 'Abs Correlation' column.
    results_df = results_df.sort_values(
        by='Correlation', 
        key=np.abs, 
        ascending=False
    )
    
    # if return_cond is True, return the full dataframe
    if return_cond:
        return results_df
    else:
        return results_df.head(top_n)


# Define your desired threshold (e.g., 0.75 for strong co-regulation)
CORRELATION_THRESHOLD = 0.75 
# YUBI ASK: how many pairs do we want to compare in the results?
# This was only used for preliminary analysis, I want all of the results that are significant
TOP_N_PAIRS = 20
# set to False if I only want top N pairs
return_full_results = True

print("extracting strongly correlated gene pairs...")
top_correlated_pairs = extract_strong_correlations(
    gene_correlation_matrix, 
    CORRELATION_THRESHOLD, 
    TOP_N_PAIRS,
    return_full_results
)

# HIERARCHICAL CLUSTERING AND HEATMAP VISUALIZATION

print("HIERARCHICAL CLUSTERING (MODULE DETECTION)...")
# Due to the large size (17634 genes), we will only cluster the first 50 genes for quick visualization.
# Clustering 17634 genes takes a significant amount of time and resources.

# YUBI: test with sample of 50, but change to full sample when running on Adroit
N_CLUSTER_SAMPLE = 50 
print(f"Clustering a sample of the first {N_CLUSTER_SAMPLE} genes for visualization...")

# Extract the sample correlation matrix
sample_corr_matrix = gene_correlation_matrix.iloc[:N_CLUSTER_SAMPLE, :N_CLUSTER_SAMPLE]

# 1. Convert correlation to distance: distance = 1 - abs(correlation)
print("converting correlation to distance...")
distance_matrix = 1 - np.abs(sample_corr_matrix)

# 2. Perform Hierarchical/Linkage Clustering
# 'ward' method minimizes the variance within each cluster
# YUBI: single or average may have faster performance but less accurate results
print("performing hierarchical clustering...")
linkage_matrix = linkage(distance_matrix, method='ward')

# 3. Plot Clustered Heatmap (Heatmap + Dendrogram)
print("plotting clustered heatmap...")
plt.figure(figsize=(12, 12))
sns.clustermap(
    sample_corr_matrix, 
    row_linkage=linkage_matrix, 
    col_linkage=linkage_matrix,
    cmap='coolwarm', # Colormap for correlation (Red for high positive, Blue for high negative)
    vmin=-1, vmax=1,
    linewidths=0.5,
    figsize=(12, 12)
)
plt.suptitle(f'Hierarchical Clustering (Co-expression) of {N_CLUSTER_SAMPLE} Genes', y=1.02, fontsize=16)
plt.show()
print(f"The Clustermap visualization shows groups (modules) of co-expressed genes.")
print(f"Genes close together on the dendrogram are highly correlated (either positively or negatively).")
print("\nAnalysis Complete. Check your file system for 'strong_gene_correlations.csv'.")
'''

# DISPLAY RESULTS: for displaying top N results only
'''
print(f"\nTOP {TOP_N_PAIRS} GENE PAIRS with |Spearman Correlation| >= {CORRELATION_THRESHOLD}")

if top_correlated_pairs.empty:
    print(f"No gene pairs found with an absolute correlation above the {CORRELATION_THRESHOLD} threshold.")
    print("Consider lowering the CORRELATION_THRESHOLD value for more results.")
else:
    # Use to_markdown for clean table printing
    print(top_correlated_pairs.to_markdown(index=False, numalign="left", stralign="left", floatfmt=".4f"))

print("\nSAMPLE CORRELATION RESULTS (First 10 Genes)")
sample_correlation = gene_correlation_matrix.iloc[:10, :10]
print(sample_correlation.to_markdown(numalign="left", stralign="left", floatfmt=".4f"))


print("saving significant gene pairs as strong_gene_correlations.csv...")

# Save the file containing ALL pairs above the threshold
OUTPUT_FILENAME = 'strong_gene_correlations.csv'

if not top_correlated_pairs.empty:
    # Save the DataFrame to a CSV file without the index
    top_correlated_pairs.to_csv(OUTPUT_FILENAME, index=False, float_format='%.6f')
    print(f"Successfully saved {len(top_correlated_pairs):,} gene pairs with |r| >= {CORRELATION_THRESHOLD} to {OUTPUT_FILENAME}")
    
else:
    print(f"No gene pairs found with an absolute correlation above the {CORRELATION_THRESHOLD} threshold. File not saved.")
'''


'''
# DISPLAY RESULTS (Sample Subset)

print("\nSAMPLE CORRELATION RESULTS (First 10 Genes)")
# The full matrix is too large to display, so print a small sample.
sample_correlation = gene_correlation_matrix.iloc[:10, :10]

print(sample_correlation.to_markdown(numalign="left", stralign="left", floatfmt=".4f"))

# Visualize the correlation of a small subset
plt.figure(figsize=(8, 7))
sns.heatmap(sample_correlation, annot=True, fmt=".2f", cmap='coolwarm',
            cbar_kws={'label': 'Spearman Correlation (in PC space)'},
            linewidths=.5, linecolor='black')
plt.title('Spearman Correlation of Top 10 Genes based on 255 PCA Loadings', fontsize=12)
plt.yticks(rotation=0)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

print("\nAnalysis Complete. The 'gene_correlation_matrix' contains the 17634x17634 results.")
'''

# Run full PCA for initial analysis
'''
print("run PCA...")

# run PCA
pca_full = PCA()
pca_full.fit(scaled_data)


# Calculate cumulative explained variance
explained_variance_ratio_cumsum = np.cumsum(pca_full.explained_variance_ratio_)

print("generating plots...")


# Plot the explained variance to determine the optimal number of components (Scree Plot)
plt.figure(figsize=(10, 5))
plt.plot(explained_variance_ratio_cumsum, marker='o', linestyle='-', color='indigo')
plt.title('Cumulative Explained Variance by Principal Component')
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Explained Variance Ratio')
plt.grid(True, linestyle='--', alpha=0.6)

# Highlight where 80% or 90% of variance is explained (common cutoff points)
cutoff_80 = np.argmax(explained_variance_ratio_cumsum >= 0.8) + 1
cutoff_90 = np.argmax(explained_variance_ratio_cumsum >= 0.9) + 1
plt.axvline(x=cutoff_80, color='r', linestyle='--', label=f'80% Variance at PC {cutoff_80}')
plt.axvline(x=cutoff_90, color='g', linestyle='--', label=f'90% Variance at PC {cutoff_90}')
plt.legend()
plt.tight_layout()
plt.show()

print(f"\nSummary of Variance Explained:")
print(f"80% of variance explained by the first {cutoff_80} PCs.")
print(f"90% of variance explained by the first {cutoff_90} PCs.")

print("run PCA with top 2 components only...")
'''

# Exploratory data analsysis
'''
# B. Run PCA for Visualization (Top 2 Components)
# For visualization, we typically select the top 2 (or 3) components.
N_COMPONENTS = 2
pca_2d = PCA(n_components=N_COMPONENTS)
principal_components = pca_2d.fit_transform(scaled_data)

# Create a DataFrame for the 2D plot
pc_df = pd.DataFrame(data = principal_components,
                     columns = ['PC 1', 'PC 2'],
                     index = df_pca_input.index)

print(f"Reduced dimension data shape: {pc_df.shape}")
print(f"Explained Variance for PC 1: {pca_2d.explained_variance_ratio_[0]:.2%}")
print(f"Explained Variance for PC 2: {pca_2d.explained_variance_ratio_[1]:.2%}")

print("visualize top 2 principal components...")


# --- 4. VISUALIZATION ---
plt.figure(figsize=(12, 8))
sns.scatterplot(
    x='PC 1',
    y='PC 2',
    data=pc_df,
    s=70,
    alpha=0.8,
    color='darkblue'
)

plt.title('PCA of Cell Line Gene Expression Data (PC 1 vs PC 2)', fontsize=16)
plt.xlabel(f'Principal Component 1 ({pca_2d.explained_variance_ratio_[0]*100:.2f}%)', fontsize=12)
plt.ylabel(f'Principal Component 2 ({pca_2d.explained_variance_ratio_[1]*100:.2f}%)', fontsize=12)
plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
plt.axvline(0, color='grey', linestyle='--', linewidth=0.5)
plt.legend(title='Hypothetical Group')
plt.grid(True, linestyle=':', alpha=0.5)
plt.show()

print("done")
'''

# YUBI: add code to save results once I'm sure that they are good