### RMarkdown Files Overview
- **FinalProjectFig2.Rmd**  
  - Generates **Figure 2** (corresponds to **Figure 1b** in the original paper) 
  - Compares GLS vs. Pearson p-value distributions

- **FinalProjectFig3.Rmd**  
  - Generates **Figures 3a** and **3b** in our project:  
    - Figure 3a corresponds to **Figure 3a** in the original paper  
    - Figure 3b corresponds to **Figure 3c** in the original paper

### Outputs
- **fig3a.jpg** — mTORC1-like module (our Figure 3a / original Figure 3a)  
- **fig3b.jpg** — Autophagy-like module (our Figure 3b / original Figure 3c)  
- **gls_fig2.png** — GLS p-value histogram (our Figure 2 / original Figure 1b)  
- **pearson_fig2.png** — Pearson p-value histogram (our Figure 2 / original Figure 1b)

### How to reproduce figures
1. Run gene_pairs.py
   - This script is not in this repo as I did not write it, but can be found [here](https://github.com/kundajelab/coessentiality/blob/master/gene_pairs.py)
   - This script was provided in the paper to compute Pearson correlation p-values used for our Figure 2, outputting the files `Pearson_p.npy` and `Pearson_sign.npy`
   - To compute the Pearson p-values, I replaced the covariance matrix with an identity matrix
     - Keeping the covariance matrix computes the GLS p-values, which are also used for Figure 2 
2. Run FinalProjectFig2.Rmd
   - **Make sure `GLS_p.npy` and `Pearson_p.npy` are in the working directory**
     - `GLS_p.npy` is the matrix of GLS p-values between every pair of genes. An earlier step in our project outputted this, but it can be recreated using `gene_pairs.py` if the identity matrix is replaced with a covariance matrix
     - `Pearson_p.npy` is the matrix of Pearson p-values between every pair of genes. This is outputted by `gene_pairs.py`
   - Extracts the **strict lower triangle** of each matrix
   - Filters to finite p-values between 0 and 1
   - Creates separate histograms for GLS and Pearson p-values
   - Uses 100 bins over [0, 1], draws a red vertical line at the median p-value, and annotates each plot with the median
   - Finally, saves **`gls_fig2.jpg`** and **`pearson_fig2.jpg`**.
     - These together form **Figure 2** in our project (corresponding to **Figure 1b** in the original paper)
3. Run FinalProjectFig3.Rmd
   - **Make sure `GLS_p.npy`, `GLS_sign.npy`, and `modules_d_0.5.csv` are in the working directory**
     - `GLS_sign.npy` is the matrix of GLS sign values for every pair of genes
     - `modules_d_0.5.csv` contains a "Members" column with a space-separated list of genes detected at density parameter d = 0.5 for each module
   - Loads ClusterONE modules for clusters 512 (mTORC1-like) and 1434 (autophagy-like)
   - Maps module genes using `genes.txt` to their corresponding p-value and sign
   - For each gene pair, calculates the signed edge weight $w_{ij} = -\log_{10}(p_{ij}) \cdot s_{ij}$
   - Only edges with p < 1e-5 are kept
   - Top neighbors per gene
     - For mTORC1 (our Fig 3a): **4 strongest neighbors per gene**  
     - For autophagy (our Fig 3b): **8 strongest neighbors per gene**
   - Builds network
     - positions nodes using a **force-directed layout** weighted by ($|w_{ij}|$)
     - draws edges with thickness proportional to co-essentiality strength
     - colors nodes by **protein complex** membership
     - uses shapes (circle vs square) to distinguish known vs newly identified pathway members
     - adds gene labels
  - Finally, the figure is saved as **`fig3a.jpg`** or **`fig3b.jpg`**.
