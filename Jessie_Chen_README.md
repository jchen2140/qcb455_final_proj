## `cancer_type_dependencies.r`

This script identifies **cancer-type–specific dependencies** of predefined **gene modules** using gene essentiality screening data.

### Inputs

* **Gene essentiality matrix** (loaded via `load_screens()`):
  Rows represent genes and columns represent cell lines. Cancer type is inferred from cell-line column names.
* **Module file:** `modules_d_{d}.csv`
  Contains a `"Members"` column with space-separated gene lists defining each module.

### Workflow

1. **Load and preprocess data**

   * Load essentiality screens and gene modules.
   * Remove cancer types represented by only one cell line.
   * One-hot encode remaining cancer types.

2. **Covariance adjustment**

   * Apply a Cholesky transform of the inverse covariance matrix to decorrelate cell-line measurements.

3. **Gene-level regression**

   * Fit an OLS model for each gene:
     *essentiality ~ cancer-type indicators*.
   * Store p-values for each cancer type.

4. **Module-level meta-analysis (ACAT)**

   * For each gene module × cancer-type pair, combine gene-level p-values using ACAT.

5. **Multiple-testing correction**

   * Apply FDR correction separately within each cancer type.
   * Keep results with FDR < 0.5.

### Output

* **`cancer_type_dependencies.tsv`**
  A ranked table of significant module–cancer-type associations with columns:
  **Rank, Module genes, Cancer type, p, FDR**.

-------------------------------------------------------------------
# `load_screens.r`

## Overview

`load_screens.r` loads and preprocesses gene essentiality screening data. It harmonizes cell-line identifiers, corrects technical bias using olfactory genes, and outputs a cleaned gene-by-cell-line matrix for downstream analysis.

---

## Inputs

### 1. `gene_effect.csv`

* Main gene essentiality matrix.
* First column contains gene identifiers (converted to rownames).
* Remaining columns are Broad Institute cell-line IDs.

### 2. `sample_info.csv`

Used to map:

* **`Broad_ID` → `CCLE_name`**
  Columns required:
* `Broad_ID`
* `CCLE_name`

### 3. `olfactory_genes.txt`

* A list of gene names (one per line).
* Used for PCA-based bias correction.

---

## Workflow

### 1. Load and Format Gene Essentiality Data

* Read `gene_effect.csv`.
* Convert the first column to rownames.
* Transpose so matrix is **genes × cell lines**.
* Trim rownames to the first word (to clean gene symbols).

### 2. Map Cell-Line Identifiers

* Read `sample_info.csv`.
* Match Broad IDs in `gene_effect.csv` to CCLE names.
* Retain only cell lines present in both files.
* Rename columns to CCLE names.

### 3. PCA-Based Bias Correction Using Olfactory Genes

1. Load olfactory gene list.
2. Subset to genes present in the screen.
3. If ≥4 olfactory genes available:

   * Run PCA on **cell lines × olfactory genes**.
   * Reconstruct the signal explained by up to the top 4 PCs.
   * Subtract these PC effects from the corresponding genes in the main matrix.
   * This removes common technical variation associated with olfactory genes.

### 4. Final Formatting

* Drop the last 4 columns to match the behavior of the Python preprocessing.
* Return a cleaned data frame of:
  **rows = genes, columns = cell lines (CCLE names).**

---

## Output

### Returned Object

A **preprocessed gene essentiality matrix** with:

* Harmonized cell-line names
* PCA-corrected values (using olfactory genes)
* Final dimension after removing last 4 columns

This matrix is used by downstream scripts such as `cancer_type_dependencies.r`.

-------------------------------------------------------------------
# figure6a.r

## Overview

`figure6a.r` visualizes the differential essentiality of gene modules across cancer types using the results file `cancer_type_dependencies.tsv`. It computes significance metrics, selects top tissue types, and generates a publication-ready plot.

## Inputs

* **cancer_type_dependencies.tsv**: tab-separated file with columns:

  * Rank
  * Module genes
  * Cancer type
  * p (ACAT p-value)
  * FDR

## Workflow

1. **Load and clean data**

   * Read TSV and clean column names using `janitor::clean_names()`
   * Rename key variables to `ModuleGenes`, `CancerType`, `Pvalue`, `FDR`
   * Compute `NegLogP = -log10(Pvalue)` and `NegLogFDR = -log10(FDR)`

2. **Compute FDR thresholds per cancer type**

   * Calculate mean `NegLogFDR` per cancer type for reference lines in the plot

3. **Select top 20 cancer types**

   * Rank by median `NegLogP`
   * Filter both plot data and FDR thresholds to top 20
   * Set factor levels for consistent ordering

4. **Visualization**

   * Jittered points of `-log10(p)` per module–cancer-type pair
   * Red horizontal segments for FDR thresholds
   * Minimal theme, rotated x-axis labels, title: "Differential essentiality of modules by tissue type"

5. **Save figure**

   * Output: `figure6a_alphabetized.png`

## Output

* **figure6a_and_fdr.png**: plot showing differential essentiality of modules across top cancer types and count of significant modules per tissue type.
