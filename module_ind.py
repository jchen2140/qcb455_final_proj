import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
COORDINATE_FILE = 'vizdf.tsv'       # Placeholder name for the file with hUMAP coordinates
MODULE_FILE = 'Supplementary_Data_2.xlsx' # Name of the Excel file
# YUBI: adjust parameter to target specific module
TARGET_MODULE = 2067

print("loading coordinate data...")
# 1. Load the coordinate data
try:
    coord_df = pd.read_csv(COORDINATE_FILE, sep='\t')
    # Assuming the column with gene names is the first column, we rename it for consistency.
    # coord_df = coord_df.rename(columns={coord_df.columns[0]: 'Gene'})
    
except FileNotFoundError:
    print(f"Error: Coordinate file '{COORDINATE_FILE}' not found. Please ensure the TSV file is uploaded.")
    exit()

print("loading module data...")
# 2. Load the module data from the 'Co-essential Modules' tab
module_df = pd.read_excel(MODULE_FILE, sheet_name="Co-essential Modules", header=2)

print("extracting genes from target module...")
# 3. Extract the genes for the target module (Module #1)
module_ind_data = module_df[module_df['Module #'] == TARGET_MODULE]

if module_ind_data.empty:
    print(f"Error: Module #{TARGET_MODULE} not found in the 'Co-essential Modules' sheet.")
    exit()


all_cols = module_ind_data.columns.tolist()
try:
    start_index = all_cols.index('Genes')
    # Select all columns from 'Genes' up to the last column
    gene_columns = all_cols[start_index:]
except ValueError:
    print("Error: Could not find a column named 'Genes'. Please check the column name.")
    exit()

# Extract all gene names from the relevant columns and flatten them into a single, clean list
gene_list_raw = module_ind_data[gene_columns].values.flatten()

# Clean the list: remove NaNs, empty strings, and convert to a unique set
module_genes = set(
    str(g).strip() for g in gene_list_raw if pd.notna(g) and str(g).strip() != ''
)

print(f"Found {len(module_genes)} unique genes for Module #{TARGET_MODULE}.")


print("merging gene list with coordinates...")
# 4. Merge gene list with coordinates
# Filter the coordinate dataframe to only include genes in the target module
df_to_plot = coord_df[coord_df['gene_names'].isin(module_genes)].copy()

if df_to_plot.empty:
    print(f"Warning: None of the genes in Module #{TARGET_MODULE} were found in the coordinate file.")
    exit()

print("plotting all genes as background, and target Module highlighted...")
# 5. Plotting (All genes as background, Target Module highlighted)
plt.figure(figsize=(10, 8))

# Plot all genes as faint background
plt.scatter(
    coord_df['hUMAP_x'],
    coord_df['hUMAP_y'],
    color='#CCCCCC', # Light grey background
    s=1,
    alpha=0.3,
    label='All Genes (Background)'
)

# Plot Target Module genes, which is the primary focus
plt.scatter(
    df_to_plot['hUMAP_x'],
    df_to_plot['hUMAP_y'],
    color='#e31a1c', # Bright Red for highlight
    s=25,        # Larger size to stand out
    alpha=0.9,
    edgecolors='black',
    linewidths=0.5,
    label=f'Module #{TARGET_MODULE} Genes ({len(df_to_plot)} genes)'
)

plt.title(f'hUMAP Coordinates for Gene Module #{TARGET_MODULE}', fontsize=16)
plt.xlabel('hUMAP_x', fontsize=12)
plt.ylabel('hUMAP_y', fontsize=12)
plt.legend(loc='best')
plt.grid(False)
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()

# Save the figure
output_filename = f'module_{TARGET_MODULE}_umap_plot.png'
plt.savefig(output_filename)

print(f"Visualization saved to '{output_filename}'")