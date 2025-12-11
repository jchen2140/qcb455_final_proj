# Written by Yubi

import numpy as np
from load_screens import load_screens
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# This file computes the reduced dimensionality representation of the genes using PCA
# It computes a Spearman correlation matrix over all genes to create modules of genes only
# It saves the resulting matrix in my repo as gene_spearman_corr_matrix.csv


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

# Save Results

# Define the directory and filename for saving the correlation matrix
FILE_NAME = 'gene_spearman_corr_matrix.csv'

print(f"\nSaving correlation matrix as {FILE_NAME}...")

# Save the DataFrame to a CSV file.
# Since the matrix is large (17634x17634), CSV is generally a good, portable format.
sample_corr_matrix.to_csv(FILE_NAME, index=True) # index=True saves the gene names as the first column

print("Correlation matrix saved successfully.")
