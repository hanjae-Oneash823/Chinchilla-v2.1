#!/usr/bin/env Rscript
# Step 03: Feature Extraction and Visualization
# Complete working script with all visualizations and filtering

library(Matrix)
library(Seurat)
library(ggplot2)

cat("\n=== Step 03: Feature Extraction ===\n")

# Dark theme function
theme_dark_custom <- function() {
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#000000", color = NA),
    panel.background = element_rect(fill = "#000000", color = NA),
    plot.title = element_text(color = "white", size = 14, face = "bold"),
    axis.title = element_text(color = "white", size = 12),
    axis.text = element_text(color = "white", size = 10),
    legend.background = element_rect(fill = "#000000", color = "white"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )
}

main <- function() {
  tryCatch({
    # STEP 1: Manual folder selection
    cat("\n[Step 1] Folder Selection\n")
    cat("Available folders in HCA_liverdata:\n")
    
    base_path <- "./HCA_liverdata"
    available_folders <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
    available_folders <- available_folders[available_folders != ""]
    
    if (length(available_folders) == 0) {
      stop("No data folders found in HCA_liverdata!")
    }
    
    for (i in seq_along(available_folders)) {
      cat(sprintf("  %d. %s\n", i, available_folders[i]))
    }
    
    cat("\nEnter folder name or number: ")
    user_input <- readLines(n = 1)
    
    # Check if input is a number
    if (grepl("^[0-9]+$", user_input)) {
      idx <- as.numeric(user_input)
      if (idx < 1 || idx > length(available_folders)) {
        stop(sprintf("Invalid selection: %d", idx))
      }
      selected_folder <- available_folders[idx]
    } else {
      # Check if input is a valid folder name
      if (!(user_input %in% available_folders)) {
        stop(sprintf("Folder not found: %s", user_input))
      }
      selected_folder <- user_input
    }
    
    data_folder <- file.path(base_path, selected_folder)
    cat(sprintf("\n  Selected: %s\n", selected_folder))
    
    # STEP 2: Load MTX files
    cat("\n[Step 2] Loading expression matrix...\n")
    
    matrix_file <- file.path(data_folder, "matrix.mtx.gz")
    features_file <- file.path(data_folder, "features.tsv.gz")
    barcodes_file <- file.path(data_folder, "barcodes.tsv.gz")
    
    if (!file.exists(matrix_file)) stop(sprintf("matrix.mtx.gz not found!"))
    if (!file.exists(features_file)) stop(sprintf("features.tsv.gz not found!"))
    if (!file.exists(barcodes_file)) stop(sprintf("barcodes.tsv.gz not found!"))
    
    expr_matrix <- readMM(matrix_file)
    features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
    barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
    
    rownames(expr_matrix) <- features
    colnames(expr_matrix) <- barcodes
    
    cat(sprintf("  Loaded: %d genes x %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))
    
    # STEP 3: Load metadata
    cat("\n[Step 3] Loading metadata...\n")
    
    metadata_files <- list.files(data_folder, pattern = "_seurat_metadata\\.csv$", full.names = TRUE)
    if (length(metadata_files) == 0) stop("Metadata file not found!")
    
    metadata <- read.csv(metadata_files[1], row.names = 1)
    cat(sprintf("  Loaded: %d cells\n", nrow(metadata)))
    
    # STEP 4: Filter cells to metadata
    cat("\n[Step 4] Filtering cells to metadata...\n")
    
    cells_keep <- intersect(colnames(expr_matrix), rownames(metadata))
    cat(sprintf("  Before: %d cells\n", ncol(expr_matrix)))
    cat(sprintf("  After: %d cells\n", length(cells_keep)))
    
    expr_matrix <- expr_matrix[, cells_keep]
    metadata <- metadata[cells_keep, , drop = FALSE]
    
    # STEP 5: Feature counts histogram
    cat("\n[Step 5] Creating feature counts histogram...\n")
    
    gene_counts <- rowSums(expr_matrix > 0)
    
    p1 <- ggplot(data.frame(cells = gene_counts), aes(x = cells)) +
      geom_histogram(bins = 50, fill = "#FF6B6B", color = "white", alpha = 0.8) +
      geom_vline(xintercept = 10, color = "#FFD700", linetype = "dashed", size = 1.2) +
      annotate("text", x = 10, y = Inf, label = "10-cell cutoff",
               color = "#FFD700", vjust = 1.5, hjust = -0.1, size = 4) +
      labs(
        title = "Feature Counts: Cells Expressing Each Gene",
        x = "Number of Cells",
        y = "Number of Genes"
      ) +
      scale_x_log10() +
      theme_dark_custom()
    
    out1 <- file.path(data_folder, "01_feature_counts_histogram.png")
    ggsave(out1, plot = p1, width = 10, height = 6, bg = "black", dpi = 300)
    cat(sprintf("  ✓ Saved: %s\n", out1))
    
    # STEP 6: Filter genes by cutoff
    cat("\n[Step 6] Filtering genes (>10 cells expression)...\n")
    
    genes_keep <- gene_counts > 10
    cat(sprintf("  Genes before: %d\n", length(genes_keep)))
    cat(sprintf("  Genes after: %d\n", sum(genes_keep)))
    
    expr_matrix <- expr_matrix[genes_keep, ]
    
    # STEP 7: Create Seurat object
    cat("\n[Step 7] Creating Seurat object...\n")
    
    seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = metadata)
    
    # STEP 8: Normalize and find HVGs
    cat("\n[Step 8] Normalizing and finding HVGs...\n")
    
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                        nfeatures = 3000, verbose = FALSE)
    
    hvg_features <- VariableFeatures(seurat_obj)
    cat(sprintf("  Top 3000 HVGs selected: %d\n", length(hvg_features)))
    
    # STEP 9: Mean-Variance plot
    cat("\n[Step 9] Creating Mean-Variance plot...\n")
    
    # Get normalized data using layer argument for Seurat 5.0
    norm_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    
    # Convert to matrix if sparse
    if (is(norm_data, "sparseMatrix")) {
      norm_data <- as.matrix(norm_data)
    }
    
    feature_stats <- data.frame(
      feature = rownames(seurat_obj),
      mean = rowMeans(norm_data),
      variance = apply(norm_data, 1, var)
    )
    
    feature_stats$is_hvg <- feature_stats$feature %in% hvg_features
    
    p2 <- ggplot(feature_stats, aes(x = log10(mean + 1), y = log10(variance + 1))) +
      geom_point(data = subset(feature_stats, !is_hvg),
                 color = "#888888", size = 1, alpha = 0.5) +
      geom_point(data = subset(feature_stats, is_hvg),
                 color = "#FF3333", size = 2, alpha = 0.7) +
      labs(
        title = "Feature Selection: Mean-Variance Plot",
        x = expression(log[10]("Mean Expression")),
        y = expression(log[10]("Variance"))
      ) +
      theme_dark_custom() +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Red: %d HVGs", sum(feature_stats$is_hvg)),
               color = "#FF3333", vjust = 1.5, hjust = 1, size = 4)
    
    out2 <- file.path(data_folder, "02_mean_variance_plot.png")
    ggsave(out2, plot = p2, width = 10, height = 8, bg = "black", dpi = 300)
    cat(sprintf("  ✓ Saved: %s\n", out2))
    
    # STEP 10: Filter to top 3000 HVGs only
    cat("\n[Step 10] Filtering to top 3000 HVGs..\n")
    
    cat(sprintf("  Genes before: %d\n", nrow(seurat_obj)))
    seurat_obj <- seurat_obj[hvg_features, ]
    cat(sprintf("  Genes after (top 3000 HVGs): %d\n", nrow(seurat_obj)))
    cat(sprintf("  Cells kept: %d\n", ncol(seurat_obj)))
    
    # STEP 11: UMAP plots - before/after gene filtering
    cat("\n[Step 11] Computing UMAP for before/after comparison...\n")
    
    # For "before" plot - all cells, all genes after >10 filter
    cat("  Processing original data with all genes...\n")
    original_seurat <- CreateSeuratObject(counts = expr_matrix)
    original_seurat <- NormalizeData(original_seurat, verbose = FALSE)
    original_seurat <- FindVariableFeatures(original_seurat, nfeatures = 3000, verbose = FALSE)
    
    # Mark which genes will be kept (HVGs)
    hvg_genes_all <- VariableFeatures(original_seurat)
    original_seurat$hvg_status <- ifelse(rownames(original_seurat) %in% hvg_genes_all, "HVG", "Other")
    
    original_seurat <- ScaleData(original_seurat, verbose = FALSE)
    original_seurat <- RunPCA(original_seurat, npcs = 50, verbose = FALSE)
    original_seurat <- RunUMAP(original_seurat, dims = 1:50, verbose = FALSE)
    
    # For "after" plot - all cells, only top 3000 HVG genes
    cat("  Processing filtered data with top 3000 HVGs...\n")
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, verbose = FALSE)
    
    p3 <- DimPlot(original_seurat, reduction = "umap", group.by = "hvg_status",
                  cols = c("HVG" = "#FF3333", "Other" = "#444444")) +
      ggtitle("Before Filtering: HVGs (Red) vs Other Genes (Gray)") +
      theme_dark_custom() +
      theme(legend.position = "bottom", legend.text = element_text(size = 10))
    
    p4 <- DimPlot(seurat_obj, reduction = "umap") +
      ggtitle("After Filtering: Top 3000 HVGs Only") +
      theme_dark_custom() +
      theme(legend.position = "none")
    
    out3 <- file.path(data_folder, "03_umap_before_filtering.png")
    out4 <- file.path(data_folder, "04_umap_after_filtering.png")
    
    ggsave(out3, plot = p3, width = 10, height = 8, bg = "black", dpi = 300)
    cat(sprintf("  ✓ Saved: %s\n", out3))
    
    ggsave(out4, plot = p4, width = 8, height = 8, bg = "black", dpi = 300)
    cat(sprintf("  ✓ Saved: %s\n", out4))
    
    # STEP 12: Save filtered genes list
    cat("\n[Step 12] Saving filtered genes list...\n")
    
    filtered_genes <- data.frame(
      gene_id = rownames(seurat_obj),
      rank = 1:nrow(seurat_obj)
    )
    
    genes_csv_path <- file.path(data_folder, "filtered_top3000_genes.csv")
    write.csv(filtered_genes, file = genes_csv_path, row.names = FALSE)
    cat(sprintf("  ✓ Saved: %s\n", genes_csv_path))
    
    # SUMMARY
    cat("\n" %+% strrep("=", 60) %+% "\n")
    cat("ANALYSIS COMPLETE\n")
    cat(strrep("=", 60) %+% "\n")
    cat(sprintf("Folder: %s\n", selected_folder))
    cat(sprintf("Initial genes: %d\n", nrow(expr_matrix)))
    cat(sprintf("Genes after >10 cell filter: %d\n", sum(gene_counts > 10)))
    cat(sprintf("Top 3000 HVGs selected: %d\n", nrow(seurat_obj)))
    cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
    cat("\nPlots saved in data folder:\n")
    cat("  1. 01_feature_counts_histogram.png\n")
    cat("  2. 02_mean_variance_plot.png\n")
    cat("  3. 03_umap_before_filtering.png\n")
    cat("  4. 04_umap_after_filtering.png\n")
    cat(strrep("=", 60) %+% "\n\n")
    
    return(seurat_obj)
    
  }, error = function(e) {
    cat("\n❌ ERROR:\n")
    cat(sprintf("  %s\n\n", e$message))
    return(NULL)
  })
}

# Execute
result <- main()

if (is.null(result)) {
  quit(status = 1)
}
