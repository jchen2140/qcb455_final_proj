# Written by Jessie Chen

# cancer_type_dependencies.R
library(tidyverse)
library(MASS)
library(matrixStats)
library(data.table)
library(stats)

# Load helper functions for reading/cleaning screen data
source("load_screens.R")

# Load screen matrix (genes x cell lines)
screens <- load_screens()
genes <- rownames(screens)

# Load module definitions from the three thresholds (note: only one was
# loaded in the python version)
modules <- list()
for (d in c(0.2, 0.5, 0.9)) {
  df <- read_csv(paste0("modules_d_", d, ".csv"), show_col_types = FALSE) %>%
    select(Members)
  
  # Extract gene lists for each module
  for (members in df$Members) {
    gene_list <- str_split(gsub('"', '', members), "\\s+")[[1]]
    modules <- append(modules, list(gene_list))
  }
}

# Annotate each cell line with its cancer type
cell_line_cancer_types <- colnames(screens) %>%
  str_split("_") %>%
  map_chr(~ str_to_title(paste(.[-1], collapse=" ")))

# Remove cancer types that appear only once
cancer_type_frequencies <- table(cell_line_cancer_types)
singletons <- names(cancer_type_frequencies[cancer_type_frequencies == 1])
mask <- !(cell_line_cancer_types %in% singletons)

# Filter screens and labels to remove singletons
screens <- screens[, mask]
cell_line_cancer_types <- cell_line_cancer_types[mask]

# One-hot encode cancer types
one_hot_cancer_types <- model.matrix(~ cell_line_cancer_types - 1)
colnames(one_hot_cancer_types) <- gsub("cell_line_cancer_types", "", colnames(one_hot_cancer_types))
cancer_types <- colnames(one_hot_cancer_types)

# Compute covariance-adjusted design matrix
cov_mat <- cov(screens)
chol_inv <- t(chol(ginv(cov_mat)))
inputs <- chol_inv %*% cbind(1, one_hot_cancer_types)

# Run gene-level regressions
gene_ps <- matrix(NA, nrow = nrow(screens), ncol = ncol(one_hot_cancer_types))
rownames(gene_ps) <- genes
colnames(gene_ps) <- cancer_types

for (i in seq_len(nrow(screens))) {
  output <- chol_inv %*% as.numeric(screens[i, ])
  model <- lm(output ~ inputs[, -1])
  gene_ps[i, ] <- summary(model)$coefficients[-1, 4]
}
gene_ps <- as.data.frame(gene_ps)

# ACAT function for combining p-values
ACAT <- function(ps) {
  p <- mean(tan((0.5 - ps) * pi))
  1 - pcauchy(p)
}

# Aggregate gene-level p-values to module-level p-values
results <- list()
for (gene_set in modules) {
  gene_set_ps <- gene_ps[rownames(gene_ps) %in% gene_set, , drop = FALSE]
  for (cancer_type in cancer_types) {
    p_meta <- ACAT(as.numeric(gene_set_ps[[cancer_type]]))
    results <- append(results, list(list(
      Module_genes = paste(gene_set, collapse = ", "),
      Cancer_type = cancer_type,
      p = p_meta
    )))
  }
}

# Convert list to a data frame
results_df <- rbindlist(results)

# Add FDR, rank, and clean final table
results_df <- results_df %>%
  group_by(Cancer_type) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  arrange(p) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

fwrite(results_df, "cancer_type_dependencies.tsv", sep = "\t", quote = TRUE, scipen = 2)
