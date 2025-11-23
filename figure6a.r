library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor) 

# Read the TSV file
df <- read_tsv("cancer_type_dependencies.tsv") %>%
  clean_names()  # converts column names to lowercase and underscores

# Check column names
colnames(df)
# Should now be: "rank", "module_genes", "cancer_type", "p", "fdr"

# Prepare for plotting
plot_df <- df %>%
  rename(
    ModuleGenes = module_genes,
    CancerType = cancer_type,
    Pvalue = p,
    FDR = fdr
  ) %>%
  mutate(
    NegLogP = -log10(Pvalue),
    NegLogFDR = -log10(FDR)
  )

# Compute average FDR threshold per cancer type
fdr_df <- plot_df %>%
  group_by(CancerType) %>%
  summarise(FDR_thresh = mean(NegLogFDR, na.rm = TRUE), .groups = "drop")

# Select top 20 tissues by median -log10(p)
top20_tissues <- plot_df %>%
  group_by(CancerType) %>%
  summarise(medP = median(NegLogP, na.rm = TRUE)) %>%
  arrange(desc(medP)) %>%
  slice_head(n = 20) %>%
  pull(CancerType)

# Filter to top tissues
plot_df <- plot_df %>% filter(CancerType %in% top20_tissues)
fdr_df <- fdr_df %>% filter(CancerType %in% top20_tissues)

# Factor levels for consistent ordering
plot_df <- plot_df %>%
  mutate(CancerType = factor(CancerType, levels = sort(unique(CancerType))))
fdr_df <- fdr_df %>%
  mutate(CancerType = factor(CancerType, levels = sort(unique(CancerType))))

# Plot
p <- ggplot(plot_df, aes(x = CancerType, y = NegLogP)) +
  geom_jitter(aes(color = CancerType), width = 0.25, size = 1.3, alpha = 0.7, show.legend = FALSE) +
  geom_segment(
    data = fdr_df,
    aes(x = CancerType, xend = CancerType, y = FDR_thresh, yend = FDR_thresh),
    color = "red", size = 1.2
  ) +
  scale_y_continuous(
    breaks = seq(0, 18, by = 2),
    limits = c(0, 18),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = NULL,
    y = expression(-log[10](P)),
    title = "Differential essentiality of modules by tissue type"
  )

print(p)

# Save figure
ggsave("figure6a_alphabetized.png", plot = p, width = 10, height = 6, dpi = 300)
