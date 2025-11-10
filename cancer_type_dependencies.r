# cancer_type_dependencies.R
library(tidyverse)
library(MASS)
library(matrixStats)
library(data.table)
library(stats)


source("load_screens.R")


screens <- load_screens()
genes <- rownames(screens)

modules <- list()
for (d in c(0.2, 0.5, 0.9)) {
  df <- read_csv(paste0("modules_d_", d, ".csv"), show_col_types = FALSE) %>%
    select(Members)
  
  for (members in df$Members) {
    gene_list <- str_split(gsub('"', '', members), "\\s+")[[1]]
    modules <- append(modules, list(gene_list))
  }
}
cell_line_cancer_types <- colnames(screens) %>%
  str_split("_") %>%
  map_chr(~ str_to_title(paste(.[-1], collapse=" ")))

cancer_type_frequencies <- table(cell_line_cancer_types)
singletons <- names(cancer_type_frequencies[cancer_type_frequencies == 1])
mask <- !(cell_line_cancer_types %in% singletons)
screens <- screens[, mask]
cell_line_cancer_types <- cell_line_cancer_types[mask]

one_hot_cancer_types <- model.matrix(~ cell_line_cancer_types - 1)
colnames(one_hot_cancer_types) <- gsub("cell_line_cancer_types", "", colnames(one_hot_cancer_types))
cancer_types <- colnames(one_hot_cancer_types)

cov_mat <- cov(screens)
chol_inv <- t(chol(ginv(cov_mat)))
inputs <- chol_inv %*% cbind(1, one_hot_cancer_types)

gene_ps <- matrix(NA, nrow = nrow(screens), ncol = ncol(one_hot_cancer_types))
rownames(gene_ps) <- genes
colnames(gene_ps) <- cancer_types

for (i in seq_len(nrow(screens))) {
  output <- chol_inv %*% as.numeric(screens[i, ])
  model <- lm(output ~ inputs[, -1])
  gene_ps[i, ] <- summary(model)$coefficients[-1, 4]
}
gene_ps <- as.data.frame(gene_ps)

ACAT <- function(ps) {
  p <- mean(tan((0.5 - ps) * pi))
  1 - pcauchy(p)
}

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

results_df <- rbindlist(results)

results_df <- results_df %>%
  group_by(Cancer_type) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  arrange(p) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

fwrite(results_df, "cancer_type_dependencies.tsv", sep = "\t", quote = TRUE, scipen = 2)
