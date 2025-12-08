# put in terminal "git pull" every time I write code BEFORE I START

# install the necessary packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("umap")
install.packages("purrr")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("umap")
library("purrr")



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


################# MY NEW FIGURE 7a/7c #######################

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
overall_genes <- as.data.frame(top_3_ag_genes) |>
  `colnames<-`(c("genes"))

bottom_3_ag <- as.data.frame(color_ag) |>
  mutate(index = 1:length(color_ag)) |>
  slice_min(color_ag, n = 3) |>
  rownames_to_column("genes")
bottom_3_ag_genes <- (bottom_3_ag$genes)
overall_genes_bottom <- as.data.frame(bottom_3_ag_genes) |>
  `colnames<-`(c("genes"))

# bone
top_3_bone <- as.data.frame(color_bone) |>
  mutate(index = 1:length(color_bone)) |>
  slice_max(color_bone, n = 3) |>
  rownames_to_column("genes")
top_3_bone_genes <- as.data.frame(top_3_bone$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_bone_genes)

bottom_3_bone <- as.data.frame(color_bone) |>
  mutate(index = 1:length(color_bone)) |>
  slice_min(color_bone, n = 3) |>
  rownames_to_column("genes")
bottom_3_bone_genes <- as.data.frame(bottom_3_bone$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_bone_genes)

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
top_3_breast_genes <- as.data.frame(top_3_breast$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_breast_genes)

bottom_3_breast <- as.data.frame(color_breast) |>
  mutate(index = 1:length(color_breast)) |>
  slice_min(color_breast, n = 3) |>
  rownames_to_column("genes")
bottom_3_breast_genes <- as.data.frame(bottom_3_breast$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_breast_genes)

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
top_3_cns_genes <- as.data.frame(top_3_cns$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_cns_genes)

bottom_3_cns <- as.data.frame(color_cns) |>
  mutate(index = 1:length(color_cns)) |>
  slice_min(color_cns, n = 3) |>
  rownames_to_column("genes")
bottom_3_cns_genes <- as.data.frame(bottom_3_cns$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_cns_genes)

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
top_3_endo_genes <- as.data.frame(top_3_endo$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_endo_genes)

bottom_3_endo <- as.data.frame(color_endo) |>
  mutate(index = 1:length(color_endo)) |>
  slice_min(color_endo, n = 3) |>
  rownames_to_column("genes")
bottom_3_endo_genes <- as.data.frame(bottom_3_endo$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_endo_genes)

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
top_3_hlt_genes <- as.data.frame(top_3_hlt$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_hlt_genes)

bottom_3_hlt <- as.data.frame(color_hlt) |>
  mutate(index = 1:length(color_hlt)) |>
  slice_min(color_hlt, n = 3) |>
  rownames_to_column("genes")
bottom_3_hlt_genes <- as.data.frame(bottom_3_hlt$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_hlt_genes)

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
top_3_kidney_genes <- as.data.frame(top_3_kidney$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_kidney_genes)

bottom_3_kidney <- as.data.frame(color_kidney) |>
  mutate(index = 1:length(color_kidney)) |>
  slice_min(color_kidney, n = 3) |>
  rownames_to_column("genes")
bottom_3_kidney_genes <- as.data.frame(bottom_3_kidney$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_kidney_genes)

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
top_3_largeint_genes <- as.data.frame(top_3_largeint$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_largeint_genes)

bottom_3_largeint <- as.data.frame(color_largeint) |>
  mutate(index = 1:length(color_largeint)) |>
  slice_min(color_largeint, n = 3) |>
  rownames_to_column("genes")
bottom_3_largeint_genes <- as.data.frame(bottom_3_largeint$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_largeint_genes)

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
top_3_liver_genes <- as.data.frame(top_3_liver$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_liver_genes)

bottom_3_liver <- as.data.frame(color_liver) |>
  mutate(index = 1:length(color_liver)) |>
  slice_min(color_liver, n = 3) |>
  rownames_to_column("genes")
bottom_3_liver_genes <- as.data.frame(bottom_3_liver$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_liver_genes)

# lung
top_3_lung <- as.data.frame(color_lung) |>
  mutate(index = 1:length(color_lung)) |>
  slice_max(color_lung, n = 3) |>
  rownames_to_column("genes")
top_3_lung_genes <- as.data.frame(top_3_lung$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_lung_genes)

bottom_3_lung <- as.data.frame(color_lung) |>
  mutate(index = 1:length(color_lung)) |>
  slice_min(color_lung, n = 3) |>
  rownames_to_column("genes")
bottom_3_lung_genes <- as.data.frame(bottom_3_lung$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_lung_genes)

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
top_3_eso_genes <- as.data.frame(top_3_eso$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_eso_genes)

bottom_3_eso <- as.data.frame(color_eso) |>
  mutate(index = 1:length(color_eso)) |>
  slice_min(color_eso, n = 3) |>
  rownames_to_column("genes")
bottom_3_eso_genes <- as.data.frame(bottom_3_eso$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_eso_genes)

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
top_3_ovary_genes <- as.data.frame(top_3_ovary$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_ovary_genes)

bottom_3_ovary <- as.data.frame(color_ovary) |>
  mutate(index = 1:length(color_ovary)) |>
  slice_min(color_ovary, n = 3) |>
  rownames_to_column("genes")
bottom_3_ovary_genes <- as.data.frame(bottom_3_ovary$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_ovary_genes)

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
top_3_pan_genes <- as.data.frame(top_3_pan$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_pan_genes)

bottom_3_pan <- as.data.frame(color_pan) |>
  mutate(index = 1:length(color_pan)) |>
  slice_min(color_pan, n = 3) |>
  rownames_to_column("genes")
bottom_3_pan_genes <- as.data.frame(bottom_3_pan$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_pan_genes)

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
top_3_ple_genes <- as.data.frame(top_3_ple$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_ple_genes)

bottom_3_ple <- as.data.frame(color_ple) |>
  mutate(index = 1:length(color_ple)) |>
  slice_min(color_ple, n = 3) |>
  rownames_to_column("genes")
bottom_3_ple_genes <- as.data.frame(bottom_3_ple$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_ple_genes)

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
top_3_skin_genes <- as.data.frame(top_3_skin$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_skin_genes)

bottom_3_skin <- as.data.frame(color_skin) |>
  mutate(index = 1:length(color_skin)) |>
  slice_min(color_skin, n = 3) |>
  rownames_to_column("genes")
bottom_3_skin_genes <- as.data.frame(bottom_3_skin$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_skin_genes)

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
top_3_soft_genes <- as.data.frame(top_3_soft$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_soft_genes)

bottom_3_soft <- as.data.frame(color_soft) |>
  mutate(index = 1:length(color_soft)) |>
  slice_min(color_soft, n = 3) |>
  rownames_to_column("genes")
bottom_3_soft_genes <- as.data.frame(bottom_3_soft$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_soft_genes)

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
top_3_sto_genes <- as.data.frame(top_3_sto$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_sto_genes)

bottom_3_sto <- as.data.frame(color_sto) |>
  mutate(index = 1:length(color_sto)) |>
  slice_min(color_sto, n = 3) |>
  rownames_to_column("genes")
bottom_3_sto_genes <- as.data.frame(bottom_3_sto$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_sto_genes)

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
top_3_thy_genes <- as.data.frame(top_3_thy$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_thy_genes)

bottom_3_thy <- as.data.frame(color_thy) |>
  mutate(index = 1:length(color_thy)) |>
  slice_min(color_thy, n = 3) |>
  rownames_to_column("genes")
bottom_3_thy_genes <- as.data.frame(bottom_3_thy$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_thy_genes)

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
top_3_uat_genes <- as.data.frame(top_3_uat$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_uat_genes)

bottom_3_uat <- as.data.frame(color_uat) |>
  mutate(index = 1:length(color_uat)) |>
  slice_min(color_uat, n = 3) |>
  rownames_to_column("genes")
bottom_3_uat_genes <- as.data.frame(bottom_3_uat$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_uat_genes)

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
top_3_urine_genes <- as.data.frame(top_3_urine$genes) |>
  `colnames<-`(c("genes"))
overall_genes <- rbind(overall_genes, top_3_urine_genes)

bottom_3_urine <- as.data.frame(color_urine) |>
  mutate(index = 1:length(color_urine)) |>
  slice_min(color_urine, n = 3) |>
  rownames_to_column("genes")
bottom_3_urine_genes <- as.data.frame(bottom_3_urine$genes) |>
  `colnames<-`(c("genes"))
overall_genes_bottom <- rbind(overall_genes_bottom, bottom_3_urine_genes)


summarized_genes_top <- overall_genes |>
  count(genes, sort = TRUE)
summarized_genes_bottom <- overall_genes_bottom |>
  count(genes, sort = TRUE)

gene_count_plot_top <- ggplot(data = summarized_genes_top,
                              aes(x = reorder(genes, n), y = n)) +
  geom_bar(stat = "identity", fill = "#ed2727") +
  labs(title = "Genes Most Highly Co-Expressed Across Cancer Types",
       x = "Genes", y = "Count") + coord_flip()

gene_count_plot_bottom <- ggplot(data = summarized_genes_bottom,
                                 aes(x = reorder(genes, n), y = n)) +
  geom_bar(stat = "identity", fill = "#4545fc") +
  labs(title = "Genes Least Co-Expressed Across Cancer Types",
       x = "Genes", y = "Count") + coord_flip()


################# MY NEW FIGURE 7b #######################

cluster1 <- data.frame(genes = c("MYC", "TAF5L"), clusters = 1)
cluster2 <- data.frame(genes = c("EP300", "TADA2B", "MED1"), clusters = 2)
cluster3 <- data.frame(genes = c("SCAP"), clusters = 3)
cluster4 <- data.frame(genes = c("ATP1A1"), clusters = 4)
cluster5 <- data.frame(genes = c("CDK6", "GRB2", "PTEN", "FERMT2", "EFR3A"),
                       clusters = 5)
cluster6 <- data.frame(genes = c("FURIN", "PSMB5", "VHL"), clusters = 6)
cluster7 <- data.frame(genes = c("TSC2", "RB1CC1", "PPP2R1A", "KNTC1", "RB1",
                                 "CDKN1A"), clusters = 7)
cluster8 <- data.frame(genes = c("RPL21"), clusters = 8)
cluster9 <- data.frame(genes = c("TYMS", "ADSL"), clusters = 9)
cluster10 <- data.frame(genes = c("CRKL", "YRDC"), clusters = 10)
cluster11 <- data.frame(genes = c("TP53"), clusters = 11)
cluster12 <- data.frame(genes = c("GPX4", "SEPHS2"), clusters = 12)
cluster13 <- data.frame(genes = c("CHMP4B", "ELMO2", "RAB10"), clusters = 13)
cluster14 <- data.frame(genes = c("NF2", "ITGAV", "TUBB"), clusters = 14)
cluster15 <- data.frame(genes = c("RPP25L"), clusters = 15)
cluster16 <- data.frame(genes = c("CDAN1"), clusters = 16)
cluster17 <- data.frame(genes = c("PDCD10"), clusters = 17)

cluster_list <- list(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6,
                     cluster7, cluster8, cluster9, cluster10, cluster11,
                     cluster12, cluster13, cluster14, cluster15,
                     cluster16, cluster17)
df_join <- purrr::reduce(cluster_list, full_join)

summarized_clusters <- df_join |>
  count(clusters, sort = FALSE)

gene_count_plot <- ggplot(data = summarized_clusters,
                          aes(x = clusters, y = n)) +
  geom_bar(stat = "identity", fill = c("#fc8c8c", "#f1b13a", "#a29f9f",
                                       "#605d5d", "#f3f36b", "#6dfb6d",
                                       "#5681f9", "#a49d9d", "#c889ef",
                                       "#fc8c8c", "#6d6b6b", "#f1b13a",
                                       "#f3f36b", "#6dfb6d", "#818080",
                                       "#484646", "#737272")) +
  labs(title = "Clusters Most Highly Co-Expressed Across Cancer Types",
       x = "Clusters", y = "Count") + coord_flip() +
  geom_text(aes(label = clusters), hjust = -1)
