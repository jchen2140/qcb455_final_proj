print("Charlotte learning to use git!")
# put in terminal "git pull" every time I write code BEFORE I START

print("trying push again")

install.packages("dplyr")
install.packages("devtools")
library("dplyr")
library("devtools")


# the following packages and instructions for downloading is
# copied directly from https://schochastics.github.io/netVizR/
install.packages(c("igraph", "graphlayouts", "ggraph", "ggforce"))
install.packages("networkdata")
install.packages("ggplot2")
devtools::install_github("schochastics/networkdata")
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(ggplot2)

gene_effect <-
  read.csv("/Users/charlottecox/Desktop/qcb455_final_proj/data/gene_effect.csv",
           header = TRUE)

print(nrow(gene_effect))

line_names <-
  read.csv("data/DepMap-2018q3-celllines.csv",
           header = TRUE, sep = ",", fill = TRUE)

bone_samples_find <- grepl("BONE", line_names$CCLE_Name, fixed = TRUE)
bone_samples <- cbind(line_names, bone <- bone_samples_find)
bone_samples <- bone_samples |>
  filter(bone == TRUE)

cell_lines <- gene_effect |>
  filter(Broad_ID %in% bone_samples$Broad_ID)


dim(cell_lines)
plot(cell_lines)
print(cell_lines[2])



install.packages("umap")
library("umap")

cell_data <- read.csv("data/modules_d_0.9.csv", header = TRUE)
head(cell_data, 3)
cell_data_labels <- cell_data[, "Members"]
cell_data_data <- cell_data[c("Cluster", "Size", "Density", "Internal.weight",
                              "External.weight", "Quality", "P.value")]
cell_umap <- umap(cell_data_data)
head(cell_umap$layout, 3)

umap_plot_df <- data.frame(cell_umap$layout)

ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2
  )
) +
  geom_point()







# this function was copied by code in the UMAP CRAN documentation
plot_umap <- function(x, labels,
                      main = "A UMAP visualization of the dataset",
                      colors = c("#ff7f00", "#e377c2", "#17becf"),
                      pad = 0.1, cex = 0.6, pch = 19,
                      add = FALSE, legend.suffix = "",
                      cex.main = 1, cex.legend = 0.85) {

  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2] - xylim[1]) * pad) * c(-0.5, 0.5)
  if (!add) {
    par(mar = c(0.2, 0.7, 1.2, 0.7), ps = 10)
    plot(xylim, xylim, type = "n", axes = FALSE, frame = FALSE)
    rect(xylim[1], xylim[1], xylim[2], xylim[2],
         border = "#aaaaaa", lwd = 0.25)
  }
  points(layout[, 1], layout[, 2], col = colors[as.integer(labels)],
         cex = cex, pch = pch)
  mtext(side = 3, main, cex = cex.main)
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend = legend.text, inset = 0.03,
         col = colors[as.integer(labels.u)],
         bty = "n", pch = pch, cex = cex.legend)
}
## <bytecode: 0x5568b6bf7680>