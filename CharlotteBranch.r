print("Charlotte learning to use git!")
# put in terminal "git pull" every time I write code BEFORE I START

print("trying push again")

install.packages("dplyr")
library("dplyr")

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


# gene_effect_table <- as.data.frame(gene_effect)
# head(gene_effect_table)