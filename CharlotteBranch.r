# put in terminal "git pull" every time I write code BEFORE I START

# install the necessary packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("umap")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("umap")



######################### FIGURE 6a, 6b ######################################

# read in the file containing the bias-corrected CERES scores for genes
bias_corrected <- read.csv("data/bias_corrected_data.csv") |>
  column_to_rownames("X")


### plotting co-essentiality of genes in bone cancer

# filter the total list of genes to just include bone cancer lines
bias_corrected_bone <- bias_corrected |>
  select(contains("BONE"))

# determine the average essentiality score for each gene, to be used
# when determining the color of genes in UMAP plotting
umap_colors_bone <- bias_corrected_bone |>
  mutate(avg = rowMeans(bias_corrected_bone))

print(umap_colors_bone$avg)
color_bone <- umap_colors_bone$avg

# determine the gene and index of the highest average essentiality score
most_essential_bone <- as.data.frame(color_bone) |>
  mutate(index = 1:length(color_bone)) |>
  slice_max(color_bone)
most_essential_index_bone <- most_essential_bone$index

# determine the gene and index of the lowest average essentiality score
least_essential_bone <- as.data.frame(color_bone) |>
  mutate(index = 1:length(color_bone)) |>
  slice_min(color_bone)
least_essential_index_bone <- least_essential_bone$index

# perform the UMAP analysis, reducing the genes to 2 dimensions
bone_umap <- umap(bias_corrected_bone)
head(bone_umap$layout, 3)
umap_plot_bone <- data.frame(bone_umap$layout)
head(umap_plot_bone, 3)
print(nrow(umap_plot_bone))

# find the coordinates for the most essential gene in the UMAP
umap_most_essential_bone <- umap_plot_bone |>
  slice(most_essential_index_bone)

# find the coordiantes for the least essential gene in the UMAP
umap_least_essential_bone <- umap_plot_bone |>
  slice(least_essential_index_bone)

# create the plot of the UMAP
ggplot(
  umap_plot_bone,
  aes(
    x = X1,
    y = X2,
    color = color_bone
  )
) +
  geom_point() +
  # add in the tri-color scale for essentiality scores
  scale_color_gradient2(high = "#ff0015", mid = "white", low = "#001eff") +
  # add a label for the most essential gene in the UMAP analysis
  geom_label(
    label = paste0("Gene with Highest \n Co-Essentiality: ",
                   rownames(umap_most_essential_bone)),
    x = umap_most_essential_bone$X1,
    y = (umap_most_essential_bone$X2 + 0.5),
    label.padding = unit(0.75, "lines"),
    color = "black",
    alpha = 0.3
  ) +
  # add a label for the least essential gene in the UMAP analysis
  geom_label(
    label = paste0("Gene with Lowest \n Co-Essentiality: ",
                   rownames(umap_least_essential_bone)),
    x = umap_least_essential_bone$X1,
    y = (umap_least_essential_bone$X2 + 0.5),
    label.padding = unit(0.75, "lines"),
    color = "black",
    alpha = 0.3
  ) + labs(title = "2D Co-Essentiality Map of
           \n Genes Found in Bone Cancer Cell Lines")



### plotting gene co-essentiality in lung cancer

# filter the total list of genes to just include lung cancer lines
bias_corrected_lung <- bias_corrected |>
  select(contains("LUNG"))

# determine the average essentiality score for each gene, to be used
# when determining the color of genes in UMAP plotting
umap_colors_lung <- bias_corrected_lung |>
  mutate(avg = rowMeans(bias_corrected_lung))

print(umap_colors_lung$avg)
color_lung <- umap_colors_lung$avg

# determine the gene and index of the highest average essentiality score
most_essential_lung <- as.data.frame(color_lung) |>
  mutate(index = 1:length(color_lung)) |>
  slice_max(color_lung)
most_essential_index_lung <- most_essential_lung$index

# determine the gene and index of the lowest average essentiality score
least_essential_lung <- as.data.frame(color_lung) |>
  mutate(index = 1:length(color_lung)) |>
  slice_min(color_lung)
least_essential_index_lung <- least_essential_lung$index

# perform the UMAP analysis, reducing the genes to 2 dimensions
lung_umap <- umap(bias_corrected_lung)
head(lung_umap$layout, 3)
umap_plot_lung <- data.frame(lung_umap$layout)
head(umap_plot_lung, 3)
print(nrow(umap_plot_lung))

# find the coordinates for the most essential gene in the UMAP
umap_most_essential_lung <- umap_plot_lung |>
  slice(most_essential_index_lung)

# find the coordinates for the least essential gene in the UMAP
umap_least_essential_lung <- umap_plot_lung |>
  slice(least_essential_index_lung)

length(color_lung)

# create the plot of the UMAP
ggplot(
  umap_plot_lung,
  aes(
    x = X1,
    y = X2,
    color = color_lung
  )
) +
  geom_point() +
  # add in the tri-color scale for essentiality scores
  scale_color_gradient2(high = "#ff0015", mid = "white", low = "#001eff") +
  # add a label for the most essential gene in the UMAP analysis
  geom_label(
    label = paste0("Gene with Highest \n Co-Essentiality: ",
                   rownames(umap_most_essential_lung)),
    x = umap_most_essential_lung$X1,
    y = (umap_most_essential_lung$X2 + 0.5),
    label.padding = unit(0.75, "lines"),
    color = "black",
    alpha = 0.3
  ) +
  # add a label for the least essential gene in the UMAP analysis
  geom_label(
    label = paste0("Gene with Lowest \n Co-Essentiality: ",
                   rownames(umap_least_essential_lung)),
    x = umap_least_essential_lung$X1,
    y = (umap_least_essential_lung$X2 + 0.5),
    label.padding = unit(0.75, "lines"),
    color = "black",
    alpha = 0.3
  ) + labs(title = "2D Co-Essentiality Map of
           \n Genes Found in Lung Cancer Cell Lines")


################# MY NEW FIGURE 7 #######################

# autonomic ganglia
bias_corrected_ag <- bias_corrected |>
  select(contains("AUTONOMIC_GANGLIA"))
umap_colors_ag <- bias_corrected_ag |>
  mutate(avg = rowMeans(bias_corrected_ag))
color_ag <- umap_colors_ag$avg
top_3_ag <- as.data.frame(color_ag) |>
  mutate(index = 1:length(color_ag)) |>
  slice_max(color_ag, n = 3) |>
  rownames_to_column("genes")
top_3_ag_genes <- (top_3_ag$genes)
ag_table <- select(hs, 
       keys = top_3_ag_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# bone
top_3_bone <- as.data.frame(color_bone) |>
  mutate(index = 1:length(color_bone)) |>
  slice_max(color_bone, n = 3) |>
  rownames_to_column("genes")
top_3_bone_genes <- (top_3_bone$genes)
bone_table <- select(hs, 
       keys = top_3_bone_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# breast
bias_corrected_breast <- bias_corrected |>
  select(contains("BREAST"))
umap_colors_breast <- bias_corrected_breast |>
  mutate(avg = rowMeans(bias_corrected_breast))
color_breast <- umap_colors_breast$avg
top_3_breast <- as.data.frame(color_breast) |>
  mutate(index = 1:length(color_breast)) |>
  slice_max(color_breast, n = 3) |>
  rownames_to_column("genes")
top_3_breast_genes <- (top_3_breast$genes)
breast_table <- select(hs, 
       keys = top_3_breast_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# central nervous system
bias_corrected_cns <- bias_corrected |>
  select(contains("CENTRAL_NERVOUS_SYSTEM"))
umap_colors_cns <- bias_corrected_cns |>
  mutate(avg = rowMeans(bias_corrected_cns))
color_cns <- umap_colors_cns$avg
top_3_cns <- as.data.frame(color_cns) |>
  mutate(index = 1:length(color_cns)) |>
  slice_max(color_cns, n = 3) |>
  rownames_to_column("genes")
top_3_cns_genes <- (top_3_cns$genes)
cns_table <- select(hs, 
       keys = top_3_cns_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# endometrium
bias_corrected_endo <- bias_corrected |>
  select(contains("ENDOMETRIUM"))
umap_colors_endo <- bias_corrected_endo |>
  mutate(avg = rowMeans(bias_corrected_endo))
color_endo <- umap_colors_endo$avg
top_3_endo <- as.data.frame(color_endo) |>
  mutate(index = 1:length(color_endo)) |>
  slice_max(color_endo, n = 3) |>
  rownames_to_column("genes")
top_3_endo_genes <- (top_3_endo$genes)
endo_table <- select(hs, 
       keys = top_3_endo_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# hematopoietic and lymphoid tissue
bias_corrected_hlt <- bias_corrected |>
  select(contains("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
umap_colors_hlt <- bias_corrected_hlt |>
  mutate(avg = rowMeans(bias_corrected_hlt))
color_hlt <- umap_colors_hlt$avg
top_3_hlt <- as.data.frame(color_hlt) |>
  mutate(index = 1:length(color_hlt)) |>
  slice_max(color_hlt, n = 3) |>
  rownames_to_column("genes")
top_3_hlt_genes <- (top_3_hlt$genes)
hlt_table <- select(hs, 
       keys = top_3_hlt_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# kidney
bias_corrected_kidney <- bias_corrected |>
  select(contains("KIDNEY"))
umap_colors_kidney <- bias_corrected_kidney |>
  mutate(avg = rowMeans(bias_corrected_kidney))
color_kidney <- umap_colors_kidney$avg
top_3_kidney <- as.data.frame(color_kidney) |>
  mutate(index = 1:length(color_kidney)) |>
  slice_max(color_kidney, n = 3) |>
  rownames_to_column("genes")
top_3_kidney_genes <- (top_3_kidney$genes)
kidney_table <- select(hs, 
       keys = top_3_kidney_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# large intestine
bias_corrected_largeint <- bias_corrected |>
  select(contains("LARGE_INTESTINE"))
umap_colors_largeint <- bias_corrected_largeint |>
  mutate(avg = rowMeans(bias_corrected_largeint))
color_largeint <- umap_colors_largeint$avg
top_3_largeint <- as.data.frame(color_largeint) |>
  mutate(index = 1:length(color_largeint)) |>
  slice_max(color_largeint, n = 3) |>
  rownames_to_column("genes")
top_3_largeint_genes <- (top_3_largeint$genes)
largeint_table <- select(hs, 
       keys = top_3_largeint_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# liver
bias_corrected_liver <- bias_corrected |>
  select(contains("LIVER"))
umap_colors_liver <- bias_corrected_liver |>
  mutate(avg = rowMeans(bias_corrected_liver))
color_liver <- umap_colors_liver$avg
top_3_liver <- as.data.frame(color_liver) |>
  mutate(index = 1:length(color_liver)) |>
  slice_max(color_liver, n = 3) |>
  rownames_to_column("genes")
top_3_liver_genes <- (top_3_liver$genes)
liver_table <- select(hs, 
       keys = top_3_liver_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# lung
top_3_lung <- as.data.frame(color_lung) |>
  mutate(index = 1:length(color_lung)) |>
  slice_max(color_lung, n = 3) |>
  rownames_to_column("genes")
top_3_lung_genes <- (top_3_lung$genes)
lung_table <- select(hs, 
       keys = top_3_lung_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# esophagus
bias_corrected_eso <- bias_corrected |>
  select(contains("OESOPHAGUS"))
umap_colors_eso <- bias_corrected_eso |>
  mutate(avg = rowMeans(bias_corrected_eso))
color_eso <- umap_colors_eso$avg
top_3_eso <- as.data.frame(color_eso) |>
  mutate(index = 1:length(color_eso)) |>
  slice_max(color_eso, n = 3) |>
  rownames_to_column("genes")
top_3_eso_genes <- (top_3_eso$genes)
eso_table <- select(hs, 
       keys = top_3_eso_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# ovary
bias_corrected_ovary <- bias_corrected |>
  select(contains("OVARY"))
umap_colors_ovary <- bias_corrected_ovary |>
  mutate(avg = rowMeans(bias_corrected_ovary))
color_ovary <- umap_colors_ovary$avg
top_3_ovary <- as.data.frame(color_ovary) |>
  mutate(index = 1:length(color_ovary)) |>
  slice_max(color_ovary, n = 3) |>
  rownames_to_column("genes")
top_3_ovary_genes <- (top_3_ovary$genes)
ovary_table <- select(hs, 
       keys = top_3_ovary_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# pancreas
bias_corrected_pan <- bias_corrected |>
  select(contains("PANCREAS"))
umap_colors_pan <- bias_corrected_pan |>
  mutate(avg = rowMeans(bias_corrected_pan))
color_pan <- umap_colors_pan$avg
top_3_pan <- as.data.frame(color_pan) |>
  mutate(index = 1:length(color_pan)) |>
  slice_max(color_pan, n = 3) |>
  rownames_to_column("genes")
top_3_pan_genes <- (top_3_pan$genes)
pan_table <- select(hs, 
       keys = top_3_pan_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# pleura
bias_corrected_ple <- bias_corrected |>
  select(contains("PLEURA"))
umap_colors_ple <- bias_corrected_ple |>
  mutate(avg = rowMeans(bias_corrected_ple))
color_ple <- umap_colors_ple$avg
top_3_ple <- as.data.frame(color_ple) |>
  mutate(index = 1:length(color_ple)) |>
  slice_max(color_ple, n = 3) |>
  rownames_to_column("genes")
top_3_ple_genes <- (top_3_ple$genes)
ple_table <- select(hs, 
       keys = top_3_ple_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# skin
bias_corrected_skin <- bias_corrected |>
  select(contains("SKIN"))
umap_colors_skin <- bias_corrected_skin |>
  mutate(avg = rowMeans(bias_corrected_skin))
color_skin <- umap_colors_skin$avg
top_3_skin <- as.data.frame(color_skin) |>
  mutate(index = 1:length(color_skin)) |>
  slice_max(color_skin, n = 3) |>
  rownames_to_column("genes")
top_3_skin_genes <- (top_3_skin$genes)
skin_table <- select(hs, 
       keys = top_3_skin_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# soft tissue
bias_corrected_soft <- bias_corrected |>
  select(contains("SOFT_TISSUE"))
umap_colors_soft <- bias_corrected_soft |>
  mutate(avg = rowMeans(bias_corrected_soft))
color_soft <- umap_colors_soft$avg
top_3_soft <- as.data.frame(color_soft) |>
  mutate(index = 1:length(color_soft)) |>
  slice_max(color_soft, n = 3) |>
  rownames_to_column("genes")
top_3_soft_genes <- (top_3_soft$genes)
soft_table <- select(hs, 
       keys = top_3_soft_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# stomach
bias_corrected_sto <- bias_corrected |>
  select(contains("STOMACH"))
umap_colors_sto <- bias_corrected_sto |>
  mutate(avg = rowMeans(bias_corrected_sto))
color_sto <- umap_colors_sto$avg
top_3_sto <- as.data.frame(color_sto) |>
  mutate(index = 1:length(color_sto)) |>
  slice_max(color_sto, n = 3) |>
  rownames_to_column("genes")
top_3_sto_genes <- (top_3_sto$genes)
sto_table <- select(hs, 
       keys = top_3_sto_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# thyroid
bias_corrected_thy <- bias_corrected |>
  select(contains("THYROID"))
umap_colors_thy <- bias_corrected_thy |>
  mutate(avg = rowMeans(bias_corrected_thy))
color_thy <- umap_colors_thy$avg
top_3_thy <- as.data.frame(color_thy) |>
  mutate(index = 1:length(color_thy)) |>
  slice_max(color_thy, n = 3) |>
  rownames_to_column("genes")
top_3_thy_genes <- (top_3_thy$genes)
thy_table <- select(hs, 
       keys = top_3_thy_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# upper aerodigestive tract
bias_corrected_uat <- bias_corrected |>
  select(contains("UPPER_AERODIGESTIVE_TRACT"))
umap_colors_uat <- bias_corrected_uat |>
  mutate(avg = rowMeans(bias_corrected_uat))
color_uat <- umap_colors_uat$avg
top_3_uat <- as.data.frame(color_uat) |>
  mutate(index = 1:length(color_uat)) |>
  slice_max(color_uat, n = 3) |>
  rownames_to_column("genes")
top_3_uat_genes <- (top_3_uat$genes)
uat_table <- select(hs, 
       keys = top_3_uat_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# urinary tract
bias_corrected_urine <- bias_corrected |>
  select(contains("URINARY_TRACT"))
umap_colors_urine <- bias_corrected_urine |>
  mutate(avg = rowMeans(bias_corrected_urine))
color_urine <- umap_colors_urine$avg
top_3_urine <- as.data.frame(color_urine) |>
  mutate(index = 1:length(color_urine)) |>
  slice_max(color_urine, n = 3) |>
  rownames_to_column("genes")
top_3_urine_genes <- (top_3_urine$genes)
urine_table <- select(hs, 
       keys = top_3_urine_genes,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

columns(org.Hs.eg.db)

# get the entrez ID for each of the genes
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

hs <- org.Hs.eg.db
symbols <- top_3_urine_genes

select(hs, 
       keys = symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")














############################# FIGURE 2a ###############################

data_2a <- read.csv("data/modules_d_0.5.csv")

## DoRothEA data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dorothea")
BiocManager::install("OmnipathR")

library(dorothea)
library(decoupleR)

human_regulons <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))
head(human_regulons)

colnames(human_regulons)
hr_source <- human_regulons$source
hr_confidence <- vector(, nrow(human_regulons))
for (i in 1:length(hr_confidence)) {
  if (human_regulons$confidence[i] == "A") {
    hr_confidence[i] <- 1
  } else if(human_regulons$confidence[i] == "B") {
    hr_confidence[i] <- 0.7
  } else if(human_regulons$confidence[i] == "C") {
    hr_confidence[i] <- 0.4
  } else if(human_regulons$confidence[i] == "D") {
    hr_confidence[i] <- 0.1
  }
}
hr_target <- human_regulons$target

dorothea_pairs <- bind_cols(hr_source, hr_target, hr_confidence)
colnames(dorothea_pairs) <- c("source", "target", "confidence")


## hu.MAP data
hu.MAP <- read.csv("data/hu.MAP.pairsWprob", sep = '\t')

## coxpres data
coxpres <- read.csv("data/coxpres_db.csv", header = TRUE)
coxpres_clean <- coxpres[ , -c(1)]
rownames(coxpres_clean) <- row_col_names
colnames(coxpres_clean) <- row_col_names

diag(coxpres_clean) <- NA
coxpres_clean[upper.tri(coxpres_clean, diag = TRUE)] <- NA

# the following code was provided by ChatGPT because of memory constraints on my computer
library(Matrix)
library(data.table)
library(parallel)

output_file <- "triplets.csv"
fwrite(data.table(gene1=character(), gene2=character(), score=numeric()),
       output_file)  # create header

# Step 3: Process each row one by one
for (i in seq_len(nrow(coxpres_clean))) {
  
  # Extract row
  row_vals <- coxpres_clean[i, ]
  
  # Identify columns with non-NA values
  non_na_cols <- which(!is.na(row_vals))
  
  # Only proceed if there are valid entries
  if (length(non_na_cols) > 0) {
    # Create triplet table for this row
    triplets <- data.table(
      gene1 = rownames(coxpres_clean)[i],
      gene2 = colnames(coxpres_clean)[non_na_cols],
      score = row_vals[non_na_cols]
    )
    
    # Immediately append to disk
    fwrite(triplets, output_file, append = TRUE)
  }
  
  # Free memory
  rm(row_vals, triplets)
  gc()
}
# end of code provided by ChatGPT for memory purposes
