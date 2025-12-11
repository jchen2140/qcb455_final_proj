# Written by Yubi

import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO
import warnings
from scipy.spatial.distance import euclidean
from load_screens import load_screens
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# YUBI: this file creates a UMAP plot of all of the genes using the PCA reduced dimensionality representation of the genes
# It highlights all of the GO terms associated with the genes

# --- CONFIGURATION ---
N_COMPONENTS = 2
N_NEIGHBORS = 15     # <--- DECREASED to prioritize local structure (tighter clusters)
MIN_DIST = 0.05      # <--- DECREASED to allow points to pack more densely 
BACKGROUND_CATEGORY = 'Unknown' # The value used in the GO term columns for unlabeled genes
N_TOP_GO_TERMS = 10  # **The number of tightest GO terms to highlight**


# Suppress UMAP Numba warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ----------------------------------------------------------------------------------
# 1. LOAD ALL DATA AND PERFORM INITIAL MERGE
# ----------------------------------------------------------------------------------

print("Loading and merging gene data and cluster assignments...")

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
    # Load the cluster assignments (Module_ID and gene names)
    clusters_df = pd.read_csv('TOM_clusters_0.99.csv')
    # Load the new TSV file
    go_term_df = pd.read_csv('vizdf.tsv', sep='\t')
    # NEW: Load the Spearman Correlation Matrix
    spearman_corr_df = pd.read_csv('gene_spearman_corr_matrix.csv', index_col=0)
    
except FileNotFoundError as e:
    print(f"Error: Could not find file {e.filename}. Please check your file paths.")
    
    
# YUBI ADDED: Helper function to clean GO term strings
def clean_go_term(go_term):
    """
    Removes the GO domain prefix (e.g., 'biological process:', 'cellular component:')
    from a GO term string for a cleaner legend.
    """
    if ':' in go_term:
        # Split by the first colon and return the part after it, stripping leading/trailing whitespace
        return go_term.split(':', 1)[1].strip()
    return go_term
    
# ----------------------------------------------------------------------------------
# 2. CLEAN AND MERGE DATA (Consolidating two GO columns)
# ----------------------------------------------------------------------------------

print("\nConsolidating GO term data...")
# Prepare GO term data for merge
# YUBI: not a fan of renaming gene_names column to Gene_Name, but I will keep it for now
go_term_df = go_term_df.rename(columns={'gene_names': 'Gene_Name'}).set_index('Gene_Name')

# Ensure all 'Unknown' or missing values are treated as the BACKGROUND_CATEGORY string
go_term_df['roarke_clusters'] = go_term_df['roarke_clusters'].fillna(BACKGROUND_CATEGORY)
go_term_df['top_enriched_clusters'] = go_term_df['top_enriched_clusters'].fillna(BACKGROUND_CATEGORY)

# Coalesce the two columns into a single 'GO_Term' column.
# np.where checks the first condition: if 'roarke_clusters' is NOT Unknown, use it.
# Otherwise, it checks the second column: if it's NOT Unknown, use it.
# If both are 'Unknown', the result is 'Unknown'.
go_term_df['GO_Term'] = np.where(
    go_term_df['roarke_clusters'] != BACKGROUND_CATEGORY, 
    go_term_df['roarke_clusters'], 
    go_term_df['top_enriched_clusters']
)

# Final merge: Add the consolidated GO terms to the feature data
full_analysis_df = loadings_matrix_df.merge(
    go_term_df[['GO_Term']],
    left_index=True,
    right_index=True,
    how='inner' # Use inner join to only keep genes present in both files
)

# YUBI: should I add .fillna({'GO_Term': BACKGROUND_CATEGORY}) here?

# Align the Spearman matrix to the genes in the analysis dataframe
gene_list = full_analysis_df.index.tolist()
try:
    spearman_corr_df = spearman_corr_df.loc[gene_list, gene_list]
except KeyError:
    print("\nError: Spearman matrix gene list does not match the analysis gene list. Reverting to Euclidean distance.")
    spearman_corr_df = None # Fallback to default UMAP metric
    

# Define feature columns (only the 232 dimensions)
feature_cols = [col for col in full_analysis_df.columns if col != 'GO_Term']

# ----------------------------------------------------------------------------------
# 2. IDENTIFY TOP 10 TIGHTEST GO TERMS
# ----------------------------------------------------------------------------------

print("\nCalculating GO Term tightness (mean distance to centroid in 232D space)...")

# Filter out the 'Unknown' genes for centroid calculation
known_go_df = full_analysis_df[full_analysis_df['GO_Term'] != BACKGROUND_CATEGORY].copy()
all_go_terms = known_go_df['GO_Term'].unique().tolist()

# 1. Calculate the Centroid (Mean Vector) for every GO Term
centroid_df = known_go_df.groupby('GO_Term')[feature_cols].mean()

# 2. Calculate Tightness Score (Mean Euclidean Distance to Centroid)
go_term_tightness = {}
MIN_GENES_FOR_TIGHTNESS = 5 # Prevent calculation on groups of 1 or 2 genes

for go_term in all_go_terms:
    go_genes = known_go_df[known_go_df['GO_Term'] == go_term]
    
    if go_genes.shape[0] < MIN_GENES_FOR_TIGHTNESS:
        # Skip very small groups, as tightness is unreliable
        continue 
        
    X_go = go_genes[feature_cols].values
    C_go = centroid_df.loc[go_term].values.reshape(1, -1)
    
    # Calculate Euclidean distance for all genes in the group to the centroid
    distances = np.linalg.norm(X_go - C_go, axis=1)
    go_term_tightness[go_term] = np.mean(distances)

# 3. Select the Top N_TOP_GO_TERMS (Lowest distance = highest tightness)
tightness_series = pd.Series(go_term_tightness)
top_tightest_go_terms = tightness_series.sort_values(ascending=True).head(N_TOP_GO_TERMS).index.tolist()

print(f"\n[Result] The top {N_TOP_GO_TERMS} tightest GO terms have been identified.")
# Print the identified GO terms for validation
print("Top 10 Tightest GO Terms:", "\n".join(top_tightest_go_terms))

# ----------------------------------------------------------------------------------
# 3. PREPARE PLOTTING DATA
# ----------------------------------------------------------------------------------

full_umap_df = full_analysis_df.copy()

# YUBI ADDED
# 1. Create the cleaned list of top GO terms for the legend/plotting categories
highlight_categories_full = sorted(top_tightest_go_terms)
highlight_categories_cleaned = [clean_go_term(term) for term in highlight_categories_full]

# Create a mapping dictionary: Full GO Term -> Cleaned GO Term
go_term_map = dict(zip(highlight_categories_full, highlight_categories_cleaned))

# Rename this to 'Plot_Category_Full' so the next step works
# Set Plot_Category: Use the actual GO term for the top 10, otherwise use 'Unknown'
full_umap_df['Plot_Category_Full'] = np.where(
    full_umap_df['GO_Term'].isin(top_tightest_go_terms),
    full_umap_df['GO_Term'],
    BACKGROUND_CATEGORY
)

# YUBI ADDED
# Apply the mapping to get the final 'Plot_Category' with short names
def map_or_keep_unknown(term):
    """Map a full GO term to its cleaned name, or keep 'Unknown' unchanged."""
    return go_term_map.get(term, term) # If term is 'Unknown', .get() returns the default 'Unknown'

full_umap_df['Plot_Category'] = full_umap_df['Plot_Category_Full'].apply(map_or_keep_unknown)

# Identify the categories to highlight (the top 10 GO terms)
# NOW USING THE CLEANED GO TERMS
highlight_categories = highlight_categories_cleaned # Use the cleaned list!
N_HIGHLIGHT_CATEGORIES = len(highlight_categories)


# Define the order of categories for plotting.
category_order = highlight_categories + [BACKGROUND_CATEGORY]
full_umap_df['Plot_Category'] = pd.Categorical(full_umap_df['Plot_Category'], categories=category_order)



# ----------------------------------------------------------------------------------
# 4. APPLY UMAP DIMENSIONALITY REDUCTION (Visualizing Gene Similarity)
# ----------------------------------------------------------------------------------

X = full_umap_df[feature_cols].values
print(f"\nRunning UMAP on all {len(X)} genes from 232D to {N_COMPONENTS}D...")

print("\nRunning UMAP using PRECOMPUTED Spearman Distance...") # <--- CHANGED/ADDED

# Use the full Spearman correlation matrix (N x N)
if spearman_corr_df is not None:
    # 1. Convert Correlation to Distance: 1 - Correlation <--- ADDED
    D = 1 - spearman_corr_df.values # <--- ADDED
    # 2. UMAP initialization: Set metric to 'precomputed' <--- CHANGED/ADDED
    reducer = umap.UMAP(
        n_components=N_COMPONENTS,
        n_neighbors=N_NEIGHBORS,
        min_dist=MIN_DIST,
        metric='precomputed', # <--- KEY CHANGE: Use the custom distance matrix
        random_state=42
    )
    # 3. Fit and Transform using the distance matrix (D) <--- CHANGED
    umap_embedding = reducer.fit_transform(D) # <--- CHANGED: Passes the distance matrix D instead of feature matrix X

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

umap_embedding = reducer.fit_transform(X)
full_umap_df['UMAP_1'] = umap_embedding[:, 0]
full_umap_df['UMAP_2'] = umap_embedding[:, 1]


# ----------------------------------------------------------------------------------
# 5. CREATE THE 2D VISUALIZATION
# ----------------------------------------------------------------------------------

print("Generating 2D plot highlighting the top 10 tightest GO terms...")

# Define color palette: Use a distinct palette for the 10 terms, gray for the background.
highlight_colors = sns.color_palette("tab10", N_HIGHLIGHT_CATEGORIES).as_hex() 
background_color = '#CCCCCC' 

custom_palette = dict(zip(highlight_categories, highlight_colors))
custom_palette[BACKGROUND_CATEGORY] = background_color


sns.set_style("whitegrid")
plt.figure(figsize=(14, 12))

# 1. Plot the 'Unknown' background layer first (Fixes seaborn categorical error)
background_data = full_umap_df[full_umap_df['Plot_Category'] == BACKGROUND_CATEGORY].copy()
background_data['Plot_Category'] = background_data['Plot_Category'].astype(object)

sns.scatterplot(
    x='UMAP_1', 
    y='UMAP_2', 
    hue='Plot_Category', 
    data=background_data, 
    palette={BACKGROUND_CATEGORY: background_color},
    s=10,    
    alpha=0.05, 
    linewidth=0,
    legend=False 
)

# 2. Plot the highlighted GO terms over the background
scatter = sns.scatterplot(
    x='UMAP_1', 
    y='UMAP_2', 
    hue='Plot_Category', 
    hue_order=highlight_categories, # hue_order uses the cleaned list
    data=full_umap_df[full_umap_df['Plot_Category'] != BACKGROUND_CATEGORY], 
    palette=custom_palette, 
    s=40,    
    alpha=1.0, 
    linewidth=0 
)

# Configure the plot appearance
plt.title(
    f'UMAP Projection Highlighting Top {N_TOP_GO_TERMS}', 
    fontsize=20, 
    fontweight='bold'
)
plt.xlabel('UMAP Component 1', fontsize=14)
plt.ylabel('UMAP Component 2', fontsize=14)

# Customize the legend
legend = scatter.legend(
    title=f'Top {N_TOP_GO_TERMS} Tightest GO Terms', 
    fontsize=10, 
    title_fontsize=12, 
    bbox_to_anchor=(1.05, 1), 
    loc='upper left',
    markerscale=2,
    ncol=1 
)
legend.get_title().set_fontweight('bold')

plt.tight_layout()

# YUBI NOTE: save figure once I confirm results
plt.savefig('umap_top10_go_spearman.png', dpi=300)

# do not show plot on adroit
# plt.show()

print("\nScript execution complete. Check 'umap_top10_tightest_go_terms.png' for the visualization.")
