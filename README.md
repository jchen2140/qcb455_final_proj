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


