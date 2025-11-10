library(tidyverse)

load_screens <- function() {
  # Load screens: first column as rownames
  screens <- readr::read_csv("gene_effect.csv", show_col_types = FALSE) %>%
    tibble::column_to_rownames(var = colnames(.)[1]) %>%
    t() %>%
    as.data.frame()
  
  # Keep only first word of rownames
  rownames(screens) <- sapply(strsplit(rownames(screens), " "), `[`, 1)
  
  # Map Broad ID to CCLE name
  cell_lines <- read.csv("sample_info.csv", check.names = FALSE)
  cell_lines <- cell_lines[, c("Broad_ID", "CCLE_name")]
  rownames(cell_lines) <- cell_lines$Broad_ID
  
  # Reorder columns to match cell_lines
  # Use intersect to avoid errors if cell_lines has IDs not in screens
  common_broad_ids <- intersect(rownames(cell_lines), colnames(screens))
  screens <- screens[, common_broad_ids, drop = FALSE]
  colnames(screens) <- cell_lines[common_broad_ids, "CCLE_name"]
  
  # Bias-correct using olfactory genes
  olfactory_genes <- trimws(readLines("olfactory_genes.txt"))
  olfactory_genes <- intersect(olfactory_genes, rownames(screens))
  olfactory_data <- screens[olfactory_genes, , drop = FALSE]
  olfactory_data <- na.omit(olfactory_data)
  
  if (nrow(olfactory_data) >= 4) {
    # PCA on transpose: rows = cell lines, columns = olfactory genes
    pca <- prcomp(t(olfactory_data), center = TRUE, scale. = FALSE)
    
    # Determine number of PCs to use (max 4)
    n_pcs <- min(4, ncol(pca$rotation), ncol(pca$x))
    
    if (n_pcs > 0) {
      # Reconstruct t(olfactory_data) using top PCs: scores %*% t(loadings)
      reconstructed_t_olfactory <- pca$x[, 1:n_pcs, drop = FALSE] %*% 
                                     t(pca$rotation[, 1:n_pcs, drop = FALSE])
      
      # Add the mean (center) back
      reconstructed_t_olfactory <- reconstructed_t_olfactory + 
                                     matrix(pca$center, 
                                            nrow = nrow(pca$x), 
                                            ncol = length(pca$center), 
                                            byrow = TRUE)
      
      # Transpose to get [genes x cell lines]
      top_PC_effects <- t(reconstructed_t_olfactory)
      
      # Assign rownames to ensure correct subtraction
      rownames(top_PC_effects) <- rownames(pca$rotation)
      
      # Subtract top PCs from screens (olfactory genes only)
      # Make sure to align by name
      common_genes <- intersect(rownames(screens), rownames(top_PC_effects))
      screens[common_genes, ] <- screens[common_genes, ] - top_PC_effects[common_genes, ]
    }
  }
  
  # Remove last 4 columns (to match Python behavior)
  screens <- screens[, 1:(ncol(screens) - 4), drop = FALSE]
  
  return(screens)
}