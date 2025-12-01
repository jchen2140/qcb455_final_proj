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
# mini_coxpres <- as(coxpres_clean, "sparseMatrix")
mini_coxpres <- which(!is.na(coxpres_clean), arr.ind = TRUE)

library(reshape2)
subset(melt(mini_coxpres))
