#!/usr/bin/env Rscript
# Step 04: Filtered Genes Voting and Selection
# =============================================
# This script compares HVG lists across samples and creates consensus gene lists

library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)

cat("\n=== Step 04: Gene List Comparison and Voting ===\n\n")

# Dark theme matching previous scripts
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
    
    # STEP 1: Load all filtered gene lists
    cat("[Step 1] Loading filtered gene lists...\n")
    
    base_path <- "./FilteredFeatures"
    
    if (!dir.exists(base_path)) {
      stop(sprintf("FilteredFeatures folder not found!"))
    }
    
    # Find all filtered_top3000_genes.csv files in the folder
    csv_files <- list.files(base_path, pattern = "filtered_top3000_genes\\.csv$", full.names = TRUE)
    
    if (length(csv_files) == 0) {
      stop("No filtered_top3000_genes.csv files found in FilteredFeatures folder!")
    }
    
    gene_lists <- list()
    rank_lists <- list()
    
    for (csv_file in csv_files) {
      # Extract sample name from filename (e.g., GSM6416567_filtered_top3000_genes.csv -> GSM6416567)
      sample_name <- gsub("_filtered_top3000_genes\\.csv$", "", basename(csv_file))
      
      df <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Strip Ensembl version numbers (remove .X from ENSG00000243485.2 -> ENSG00000243485)
      df$gene_id <- gsub("\\.\\d+$", "", df$gene_id)
      
      gene_lists[[sample_name]] <- df$gene_id
      rank_lists[[sample_name]] <- setNames(df$rank, df$gene_id)
      cat(sprintf("  ✓ Loaded %s: %d genes\n", sample_name, nrow(df)))
    }
    
    n_samples <- length(gene_lists)
    cat(sprintf("\nTotal samples loaded: %d\n\n", n_samples))
    
    if (n_samples < 2) {
      stop("Need at least 2 samples to compare!")
    }
    
    # STEP 2: UpSet plot of original samples
    cat("[Step 2] Creating UpSet plot for original samples...\n")
    
    # Create binary matrix for UpSet
    all_genes <- unique(unlist(gene_lists))
    upset_matrix <- matrix(0, nrow = length(all_genes), ncol = n_samples)
    rownames(upset_matrix) <- all_genes
    colnames(upset_matrix) <- names(gene_lists)
    
    for (i in seq_along(gene_lists)) {
      upset_matrix[gene_lists[[i]], i] <- 1
    }
    
    # Convert to combination matrix for ComplexHeatmap
    m <- make_comb_mat(upset_matrix)
    
    # Create UpSet plot with clean black and white theme
    png(file.path(base_path, "01_upset_original_samples.png"), 
        width = 14, height = 6.5, units = "in", res = 300, bg = "white")
    
    ht <- UpSet(m,
                set_order = names(gene_lists),
                comb_col = viridis(1, option = "magma", end = 0.7),
                comb_order = order(colSums(m), decreasing = TRUE),
                bg_col = "#FFFFFF",
                bg_pt_col = "#DDDDDD",
                pt_size = unit(5, "mm"),
                lwd = 2.5,
                gap = unit(0.3, "cm"),
                top_annotation = upset_top_annotation(m, ylim = c(0, max(comb_size(m)) * 1.15), gp = gpar(fill = viridis(1, option = "magma", end = 0.7), col = NA), annotation_name_gp = gpar(col = "black", fontsize = 12), axis_param = list(gp = gpar(col = "black", fontsize = 11))),
                right_annotation = upset_right_annotation(m, ylim = c(0, max(set_size(m)) * 1.15), gp = gpar(fill = "#333333", col = NA), annotation_name_gp = gpar(col = "black", fontsize = 12), axis_param = list(gp = gpar(col = "black", fontsize = 11))))
    
    draw(ht)
    dev.off()
    cat("  ✓ Saved: 01_upset_original_samples.png\n\n")
    
    # STEP 3: RRHO comparisons between samples (skipped - optional visualization)
    cat("[Step 3] Skipping RRHO sample comparison\n\n")
    
    # STEP 4: Create combined gene lists using 3 methods
    cat("[Step 4] Creating combined gene lists using 3 methods...\n\n")
    
    # Method 1: Intersection (genes in ALL samples)
    cat("  Method 1: Intersection (all samples)...\n")
    intersection_genes <- Reduce(intersect, gene_lists)
    cat(sprintf("    Genes in all %d samples: %d\n", n_samples, length(intersection_genes)))
    
    # Method 2: N-sample voting (genes in at least 3 samples)
    cat("  Method 2: N-sample voting (≥3 samples)...\n")
    gene_counts <- table(unlist(gene_lists))
    voting_genes <- names(gene_counts[gene_counts >= min(3, n_samples)])
    cat(sprintf("    Genes in ≥3 samples: %d\n", length(voting_genes)))
    
    # Method 3: Rank aggregation (average rank, top 3000)
    cat("  Method 3: Rank aggregation (average rank)...\n")
    
    # Calculate average rank for each gene
    all_unique_genes <- unique(unlist(gene_lists))
    avg_ranks <- sapply(all_unique_genes, function(gene) {
      ranks <- sapply(rank_lists, function(rl) {
        if (gene %in% names(rl)) {
          return(rl[gene])
        } else {
          return(3001)  # Penalize genes not in list
        }
      })
      return(mean(ranks))
    })
    
    # Sort by average rank and take top 3000
    sorted_genes <- names(sort(avg_ranks))
    rank_agg_genes <- head(sorted_genes, 3000)
    cat(sprintf("    Top 3000 by average rank selected\n"))
    
    cat("\n")
    
    # STEP 5: Save combined gene lists
    cat("[Step 5] Saving combined gene lists...\n")
    
    # Save intersection list
    if (length(intersection_genes) > 0) {
      write.csv(data.frame(gene_id = intersection_genes),
                file.path(base_path, "combined_intersection.csv"),
                row.names = FALSE)
      cat("  ✓ Saved: combined_intersection.csv\n")
    }
    
    # Save voting list
    write.csv(data.frame(gene_id = voting_genes),
              file.path(base_path, "combined_voting.csv"),
              row.names = FALSE)
    cat("  ✓ Saved: combined_voting.csv\n")
    
    # Save rank aggregation list
    write.csv(data.frame(gene_id = rank_agg_genes,
                         avg_rank = avg_ranks[rank_agg_genes]),
              file.path(base_path, "combined_rank_aggregation.csv"),
              row.names = FALSE)
    cat("  ✓ Saved: combined_rank_aggregation.csv\n\n")
    
    # STEP 6: Compare the 3 methods using UpSet
    cat("[Step 6] Creating UpSet plot to compare methods...\n")
    
    method_lists <- list(
      Intersection = intersection_genes,
      Voting = voting_genes,
      RankAgg = rank_agg_genes
    )
    
    # Create binary matrix
    all_method_genes <- unique(c(intersection_genes, voting_genes, rank_agg_genes))
    method_matrix <- matrix(0, nrow = length(all_method_genes), ncol = 3)
    rownames(method_matrix) <- all_method_genes
    colnames(method_matrix) <- c("Intersection", "Voting (≥3)", "Rank Aggregation")
    
    method_matrix[intersection_genes, 1] <- 1
    method_matrix[voting_genes, 2] <- 1
    method_matrix[rank_agg_genes, 3] <- 1
    
    m2 <- make_comb_mat(method_matrix)
    
    png(file.path(base_path, "03_upset_method_comparison.png"),
        width = 12, height = 5.5, units = "in", res = 300, bg = "white")
    
    ht2 <- UpSet(m2,
                 comb_col = viridis(1, option = "magma", end = 0.7),
                 comb_order = order(colSums(m2), decreasing = TRUE),
                 bg_col = "#FFFFFF",
                 bg_pt_col = "#DDDDDD",
                 pt_size = unit(5, "mm"),
                 lwd = 2.5,
                 gap = unit(0.3, "cm"),
                 top_annotation = upset_top_annotation(m2,ylim = c(0, max(comb_size(m2)) * 1.15), gp = gpar(fill = viridis(1, option = "magma", end = 0.7), col = NA), annotation_name_gp = gpar(col = "black", fontsize = 12), axis_param = list(gp = gpar(col = "black", fontsize = 11))),
                 right_annotation = upset_right_annotation(m2, ylim = c(0, max(set_size(m2)) * 1.15), gp = gpar(fill = "#333333", col = NA), annotation_name_gp = gpar(col = "black", fontsize = 12), axis_param = list(gp = gpar(col = "black", fontsize = 11))))
    
    draw(ht2)
    dev.off()
    cat("  ✓ Saved: 03_upset_method_comparison.png\n\n")
    
    # STEP 7: Intersection Heatmap (Consensus vs. Originals)
    cat("[Step 7] Creating intersection heatmap (Consensus vs. Originals)...\n")
    
    # Create 3x5 matrix: methods x samples
    # Each cell = percentage of original sample genes captured by each method
    intersection_matrix <- matrix(0, nrow = 3, ncol = n_samples)
    rownames(intersection_matrix) <- c("Rank Aggregation", "Voting", "Intersection")
    colnames(intersection_matrix) <- names(gene_lists)
    
    # Calculate percentage overlap
    intersection_matrix[1, ] <- sapply(gene_lists, function(sample_genes) {
      length(intersect(sample_genes, rank_agg_genes)) / length(sample_genes) * 100
    })
    intersection_matrix[2, ] <- sapply(gene_lists, function(sample_genes) {
      length(intersect(sample_genes, voting_genes)) / length(sample_genes) * 100
    })
    intersection_matrix[3, ] <- sapply(gene_lists, function(sample_genes) {
      length(intersect(sample_genes, intersection_genes)) / length(sample_genes) * 100
    })
    
    # Create heatmap
    png(file.path(base_path, "04_intersection_heatmap.png"),
        width = 8, height = 5, units = "in", res = 300, bg = "white")
    
    ht_intersection <- Heatmap(intersection_matrix,
                               name = "% Overlap",
                               col = viridis(256, option = "viridis"),
                               rect_gp = gpar(col = "white", lwd = 1),
                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 grid.text(sprintf("%.1f%%", intersection_matrix[i, j]),
                                           x, y, gp = gpar(col = "white", fontsize = 10))
                               },
                               column_names_gp = gpar(col = "black", fontsize = 11),
                               row_names_gp = gpar(col = "black", fontsize = 11),
                               heatmap_legend_param = list(title_gp = gpar(col = "black", fontsize = 11),
                                                            labels_gp = gpar(col = "black", fontsize = 10)))
    
    draw(ht_intersection)
    dev.off()
    cat("  ✓ Saved: 04_intersection_heatmap.png\n\n")
    
    # STEP 8: Multi-Set UpSet Plot (All 8 lists)
    cat("[Step 8] Creating multi-set UpSet plot (3 methods + 5 samples)...\n")
    
    # Create matrix with all 8 lists
    all_8_genes <- unique(unlist(c(gene_lists, method_lists)))
    upset_matrix_all <- matrix(0, nrow = length(all_8_genes), ncol = 8)
    rownames(upset_matrix_all) <- all_8_genes
    colnames(upset_matrix_all) <- c(names(gene_lists), "Intersection", "Voting", "Rank Agg")
    
    # First 5 columns: original samples
    for (i in 1:n_samples) {
      sample_genes <- gene_lists[[i]]
      matching_rows <- rownames(upset_matrix_all) %in% sample_genes
      upset_matrix_all[matching_rows, i] <- 1
    }
    
    # Last 3 columns: consensus methods
    upset_matrix_all[rownames(upset_matrix_all) %in% intersection_genes, 6] <- 1
    upset_matrix_all[rownames(upset_matrix_all) %in% voting_genes, 7] <- 1
    upset_matrix_all[rownames(upset_matrix_all) %in% rank_agg_genes, 8] <- 1
    
    m_all <- make_comb_mat(upset_matrix_all)
    
    png(file.path(base_path, "05_upset_all_sets.png"),
        width = 16, height = 8, units = "in", res = 300, bg = "white")
    
    ht_all <- UpSet(m_all,
                    comb_col = viridis(1, option = "magma", end = 0.7),
                    comb_order = order(colSums(m_all), decreasing = TRUE),
                    bg_col = "#FFFFFF",
                    bg_pt_col = "#DDDDDD",
                    pt_size = unit(4, "mm"),
                    lwd = 2,
                    gap = unit(0.3, "cm"),
                    top_annotation = upset_top_annotation(m_all,
                                                           ylim = c(0, max(comb_size(m_all)) * 1.15),
                                                           gp = gpar(fill = viridis(1, option = "magma", end = 0.7), col = NA),
                                                           annotation_name_gp = gpar(col = "black", fontsize = 12),
                                                           axis_param = list(gp = gpar(col = "black", fontsize = 10))),
                    right_annotation = upset_right_annotation(m_all,
                                                               ylim = c(0, max(set_size(m_all)) * 1.15),
                                                               gp = gpar(fill = "#333333", col = NA),
                                                               annotation_name_gp = gpar(col = "black", fontsize = 11),
                                                               axis_param = list(gp = gpar(col = "black", fontsize = 10))))
    
    draw(ht_all)
    dev.off()
    cat("  ✓ Saved: 05_upset_all_sets.png\n\n")
    
    # STEP 9: Summary statistics
    cat("[Step 9] Generating summary statistics...\n\n")
    cat("ANALYSIS COMPLETE\n")
    cat(paste(strrep("=", 60), "\n", sep=""))
    cat(sprintf("Samples analyzed: %d\n", n_samples))
    cat(sprintf("Total unique genes across all samples: %d\n", length(all_genes)))
    cat("\nCombined gene list sizes:\n")
    cat(sprintf("  Intersection (all samples): %d genes\n", length(intersection_genes)))
    cat(sprintf("  Voting (≥3 samples): %d genes\n", length(voting_genes)))
    cat(sprintf("  Rank Aggregation (top 3000): %d genes\n", length(rank_agg_genes)))
    cat("\nOverlap between methods:\n")
    cat(sprintf("  Intersection ∩ Voting: %d genes\n", 
                length(intersect(intersection_genes, voting_genes))))
    cat(sprintf("  Intersection ∩ Rank Agg: %d genes\n",
                length(intersect(intersection_genes, rank_agg_genes))))
    cat(sprintf("  Voting ∩ Rank Agg: %d genes\n",
                length(intersect(voting_genes, rank_agg_genes))))
    cat(sprintf("  All 3 methods: %d genes\n",
                length(Reduce(intersect, method_lists))))
    
    cat("\nPlots saved:\n")
    cat("  1. 01_upset_original_samples.png\n")
    cat("  2. 03_upset_method_comparison.png\n")
    cat("  3. 04_intersection_heatmap.png\n")
    cat("  4. 05_upset_all_sets.png\n")
    cat(paste(strrep("=", 60), "\n\n", sep=""))
    
    return(method_lists)
    
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


