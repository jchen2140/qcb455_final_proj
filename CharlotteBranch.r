print("Charlotte learning to use git!")
# put in terminal "git pull" every time I write code BEFORE I START

print("trying push again")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("umap")
library("umap")
library("dplyr")
library("ggplot2")
library("tidyverse")

gene_effect <- read.csv("data/gene_effect.csv", header = TRUE)

print(nrow(gene_effect))

line_names <- read.csv("data/DepMap-2018q3-celllines.csv",
                       header = TRUE, sep = ",", fill = TRUE)

bias_corrected <- read.csv("data/bias_corrected_data.csv") |>
  column_to_rownames("X")

bias_corrected_bone <- bias_corrected |>
  select(contains("BONE"))

umap_colors <- bias_corrected_bone |>
  mutate(avg = rowMeans(bias_corrected_bone))

print(umap_colors$avg)
color <- umap_colors$avg


gene_colors <- as.data.frame(t(bias_corrected_bone))
nrow(bias_corrected_bone_trans)

bone_umap <- umap(bias_corrected_bone)
head(bone_umap$layout, 3)
umap_plot_bone <- data.frame(bone_umap$layout)
print(nrow(umap_plot_bone))

ggplot(
  umap_plot_bone,
  aes(
    x = X1,
    y = X2,
    color = color
  )
) +
  geom_point() + scale_color_gradient(high = "blue", low = "yellow")
