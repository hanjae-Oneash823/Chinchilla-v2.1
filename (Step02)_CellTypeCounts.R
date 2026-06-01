# Read the metadata CSV file
metadata <- read.csv("seurat_metadata.csv", row.names = 1)

# Get cell type counts
cell_type_counts <- table(metadata$SingleR_label)

# Display results
print("Cell Type Counts:")
print(cell_type_counts)

# Total cells
total_cells <- nrow(metadata)
print(paste("Total cells:", total_cells))

# Optional: Create a nice summary table
summary_df <- data.frame(
  Cell_Type = names(cell_type_counts),
  Count = as.numeric(cell_type_counts),
  Percentage = round(as.numeric(cell_type_counts) / total_cells * 100, 2)
)

print(summary_df)

# Save summary to file
write.csv(summary_df, "cell_type_summary.csv", row.names = FALSE)





# =========================================
# summary



## Comprehensive review of cell type annotations across multiple samples

# Install required packages if missing
required_packages <- c("ggplot2", "dplyr", "readr", "stringr", "tidyr", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(scales)

# Locate all metadata CSV files (e.g., GSM*_seurat_metadata.csv)
metadata_files <- list.files(
  path = ".",
  pattern = "_seurat_metadata\\.csv$",
  full.names = TRUE
)

if (length(metadata_files) == 0) {
  stop("No metadata files found matching *_seurat_metadata.csv in the current directory.")
}

# Read and combine all metadata files
all_metadata <- lapply(metadata_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df$sample_id <- str_replace(basename(f), "_seurat_metadata\\.csv$", "")
  df
}) %>% bind_rows()

if (!"SingleR_label" %in% colnames(all_metadata)) {
  stop("SingleR_label column not found in metadata files.")
}

# Compute counts per sample and cell type
cell_type_counts <- all_metadata %>%
  count(sample_id, SingleR_label, name = "count") %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(count),
    percent = 100 * count / total_cells
  ) %>%
  ungroup()

# Save combined summary table
write.csv(cell_type_counts, "combined_cell_type_summary.csv", row.names = FALSE)

# Create plots directory
plots_dir <- "plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Neon palette settings (turbo-like)
neon_gradient <- c("#050505", "#00F5FF", "#7CFF6B", "#FEE440", "#FF2D95")

# Dark/minimal scientific theme (true black background)
theme_dark_minimal <- theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#000000", color = NA),
    plot.background = element_rect(fill = "#000000", color = NA),
    legend.background = element_rect(fill = "#000000", color = NA),
    legend.key = element_rect(fill = "#000000", color = NA),
    text = element_text(color = "#EDEDED"),
    axis.text = element_text(color = "#CFCFCF"),
    axis.title = element_text(color = "#EDEDED"),
    plot.title = element_text(color = "#FFFFFF", face = "bold"),
    plot.subtitle = element_text(color = "#CFCFCF"),
    panel.grid.major = element_line(color = "#1a1a1a"),
    panel.grid.minor = element_line(color = "#101010")
  )

# Plot 1: Stacked bar (counts per sample)
p1 <- ggplot(cell_type_counts, aes(x = sample_id, y = count, fill = SingleR_label)) +
  geom_col(width = 0.85) +
  labs(
    title = "Cell Type Counts by Sample",
    x = "Sample",
    y = "Cell Count",
    fill = "Cell Type"
  ) +
  scale_fill_viridis_d(option = "turbo", end = 0.95) +
  theme_dark_minimal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(plots_dir, "cell_type_counts_by_sample.png"),
  plot = p1, width = 11, height = 6, dpi = 300, bg = "#000000")

# Plot 2: Stacked bar (percent per sample)
p2 <- ggplot(cell_type_counts, aes(x = sample_id, y = percent, fill = SingleR_label)) +
  geom_col(width = 0.85) +
  labs(
    title = "Cell Type Proportions by Sample",
    x = "Sample",
    y = "Percent of Cells",
    fill = "Cell Type"
  ) +
  scale_fill_viridis_d(option = "turbo", end = 0.95) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_dark_minimal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(plots_dir, "cell_type_proportions_by_sample.png"),
  plot = p2, width = 11, height = 6, dpi = 300, bg = "#000000")

# Plot 3: Heatmap of counts (log scale)
heatmap_data <- cell_type_counts %>%
  mutate(log_count = log10(count + 1))

p3 <- ggplot(heatmap_data, aes(x = sample_id, y = SingleR_label, fill = log_count)) +
  geom_tile(color = "#000000") +
  scale_fill_gradientn(colors = neon_gradient) +
  labs(
    title = "Cell Type Abundance (log10 counts)",
    x = "Sample",
    y = "Cell Type",
    fill = "log10(count+1)"
  ) +
  theme_dark_minimal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(plots_dir, "cell_type_heatmap_logcounts.png"),
  plot = p3, width = 10, height = 7, dpi = 300, bg = "#000000")

# Plot 4: Overall composition across all samples
overall_counts <- all_metadata %>%
  count(SingleR_label, name = "count") %>%
  mutate(percent = 100 * count / sum(count))

p4 <- ggplot(overall_counts, aes(x = reorder(SingleR_label, count), y = count, fill = SingleR_label)) +
  geom_col(width = 0.85, show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Overall Cell Type Composition",
    x = "Cell Type",
    y = "Cell Count"
  ) +
  scale_fill_viridis_d(option = "turbo", end = 0.95) +
  theme_dark_minimal

ggsave(filename = file.path(plots_dir, "overall_cell_type_composition.png"),
  plot = p4, width = 9, height = 7, dpi = 300, bg = "#000000")





