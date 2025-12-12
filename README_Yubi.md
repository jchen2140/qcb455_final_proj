I wrote or significantly modified scripts from Wainberg et al.'s shared code to create Figure 1a and Figure 1b.

bias_corrected_data.py returned the normalized and bias corrected gene expression data as a csv file to share the dataset with team members for downstream data analysis. I wrote this file completely.
fig1_mappings.py created the full Figure 1a using the diffusion map and UMAP reduction, as well as added the manual annotations of module functions to the plot. This script was based on fig1.py with all modifications added by Yubi.
fig1.py created Figure 1a without the module labels. This script was based on the authors' script (generate_layout.py).
gene_pairs.py conducted the Generalized Least Squares algorithm on the bias corrected and normalized gene expression data and saved the results in npy and txt files for downstream analysis. I modified this script to run on the 3.14 python environment.
load_screens.py normalized and bias corrected the gene expression data. I modified this script to run on the 3.14 python environment.
module_analysis.rmd analyzed the quality and statistical significance of the modules found through ClusterONE analysis. I wrote this script entirely.
module_ind.py plotted individual modules in 2D space to analyze the statistical relationships and structure for elucidating underlying biological connectiveness. This script was based on the authors' script (generate_layout.py) with significant modifications by me.
modules.py ran the ClusterONE algorithm on the GLS matrix to generate co-essential modules. I modified this script to test different thresholds for graph construction.
qcb_env.yml set up the local python environment for running all code. I wrote this entirely.
view_GLS.py viewed the npy results of the GLS algorithm for sharing with team members for downstream analysis. I wrote this entirely.
PCA.py ran the full PCA dimensionality reduction of the normalized and bias corrected gene expression data, computed a full spearman correlation matrix over all genes, and conducted hierarchical clustering of the genes. This script included testing different parameters. I wrote this script entirely.
PCA_GO_visualize.py plotted the results of the hierarchical clustering analysis in a 2D space using UMAP. This created Figure 1c, which was ultimately scratched from the project. I wrote this entirely.
PCA_only.py computed the full PCA dimensionality reduction of the normalized and bias corrected gene expression data and computed a full spearman correlation matrix over all genes, which was saved on the Adroit cluster because it was too large for the local machine. I wrote this entirely.
PCA_visualize.py used the computed results from PCA_only.py to create a UMAP reduction and plot the gene modules in a 2D space to create Figure 1b. It also added the manual annotations as module labels onto the plot. I wrote this script entirely.
job.slurm was the script that was modified and used to run all jobs on Adroit cluster. I wrote this entirely.
