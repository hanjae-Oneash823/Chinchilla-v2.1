
# ============================================================================
# Training Data Generation from Filtered Seurat Object
# GISTomics-Chinchilla Project
# ============================================================================
# This script:
# 1. Allows user to select a subfolder from HCA_liverdata
# 2. Loads matrix/features/barcodes as a Seurat object
# 3. Filters for specified genes and cells
# 4. Generates training data CSV files (gene counts + cell type labels)
# ============================================================================

library(Seurat)
library(Matrix)

# Function to list available subfolders
list_available_subfolders <- function(base_path) {
  subfolders <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
  return(sort(subfolders))
}

# Function to check if required files exist (handles .gz compressed files)
check_required_files <- function(sample_dir) {
  # Look for either compressed (.gz) or uncompressed versions
  matrix_file <- NULL
  features_file <- NULL
  barcodes_file <- NULL
  
  # Check for matrix files
  if (file.exists(file.path(sample_dir, "matrix.mtx.gz"))) {
    matrix_file <- file.path(sample_dir, "matrix.mtx.gz")
  } else if (file.exists(file.path(sample_dir, "matrix.mtx"))) {
    matrix_file <- file.path(sample_dir, "matrix.mtx")
  }
  
  # Check for features files
  if (file.exists(file.path(sample_dir, "features.tsv.gz"))) {
    features_file <- file.path(sample_dir, "features.tsv.gz")
  } else if (file.exists(file.path(sample_dir, "features.tsv"))) {
    features_file <- file.path(sample_dir, "features.tsv")
  }
  
  # Check for barcodes files
  if (file.exists(file.path(sample_dir, "barcodes.tsv.gz"))) {
    barcodes_file <- file.path(sample_dir, "barcodes.tsv.gz")
  } else if (file.exists(file.path(sample_dir, "barcodes.tsv"))) {
    barcodes_file <- file.path(sample_dir, "barcodes.tsv")
  }
  
  return(list(
    matrix = matrix_file,
    features = features_file,
    barcodes = barcodes_file
  ))
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

generateTrainingData <- function() {
  
  # Display Chinchilla Logo
  cat("\n")
  cat("    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\n")
  cat("    ▓                                               ▓\n")
  cat("    ▓   🐭 GISTomics-CHINCHILLA 🐭                 ▓\n")
  cat("    ▓                                               ▓\n")
  cat("    ▓   ➤ (´◔ω◔`) TRAINING DATA GENERATOR ◀        ▓\n")
  cat("    ▓                                               ▓\n")
  cat("    ▓   >> UNLEASHING RNA-SEQ FURY <<              ▓\n")
  cat("    ▓                                               ▓\n")
  cat("    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\n")
  cat("\n")
  
  cat("\n🐭 Chinchilla is ready to process your data! 🐭\n\n")
  cat("========================================\n")
  cat("Training Data Generation Script\n")
  cat("========================================\n\n")

  # Set base directories
  base_dir <- getwd()
  hca_data_dir <- file.path(base_dir, "HCA_liverdata")
  training_data_dir <- file.path(base_dir, "Training Data")

# Create Training Data directory if it doesn't exist
if (!dir.exists(training_data_dir)) {
  dir.create(training_data_dir, showWarnings = FALSE)
  cat("Created 'Training Data' directory\n\n")
}

# List available subfolders
available_subfolders <- list_available_subfolders(hca_data_dir)

if (length(available_subfolders) == 0) {
  stop("No subfolders found in HCA_liverdata directory!")
}

cat("Available samples:\n")
for (i in seq_along(available_subfolders)) {
  cat(sprintf("%d. %s\n", i, available_subfolders[i]))
}

cat("\n")
choice <- as.integer(readline(prompt = "Select a sample (enter number): "))

if (is.na(choice) || choice < 1 || choice > length(available_subfolders)) {
  stop("Invalid selection!")
}

selected_sample <- available_subfolders[choice]
sample_dir <- file.path(hca_data_dir, selected_sample)

cat(sprintf("\nSelected: %s\n", selected_sample))

# Check for required files
file_check <- check_required_files(sample_dir)

if (is.null(file_check$matrix) || is.null(file_check$features) || is.null(file_check$barcodes)) {
  stop("Required files (matrix.mtx, features.tsv, barcodes.tsv) not found in ", sample_dir)
}

cat("Loading matrix data...\n")

# Load the matrix, features, and barcodes
# readMM handles .gz files automatically
matrix_data <- readMM(file_check$matrix)
features <- read.delim(file_check$features, header = FALSE, row.names = 1)
barcodes <- read.delim(file_check$barcodes, header = FALSE)

colnames(matrix_data) <- barcodes[, 1]
rownames(matrix_data) <- rownames(features)

cat(sprintf("Matrix dimensions: %d genes x %d cells\n", nrow(matrix_data), ncol(matrix_data)))
cat("🐭 Chinchilla approves of this data! 🐭\n\n")

# Create Seurat object
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = matrix_data)

# Strip Ensembl gene version numbers (e.g., ENSG00000000003.14 -> ENSG00000000003)
rownames(seurat_obj) <- gsub("\\..*$", "", rownames(seurat_obj))
cat("Stripped Ensembl gene version numbers\n")

# Load filtered genes list
filtered_genes_path <- file.path(hca_data_dir, "filteredfeatures.csv")

if (!file.exists(filtered_genes_path)) {
  stop("filteredfeatures.csv not found in HCA_liverdata directory!")
}

cat("Loading filtered genes list...\n")
filtered_genes <- read.csv(filtered_genes_path, stringsAsFactors = FALSE)$gene_id

cat(sprintf("Total filtered genes on list: %d\n", length(filtered_genes)))
cat("🐭 Chinchilla's fluffy brain is processing 3000 genes! 🐭\n\n")

# Find matching genes between filtered list and seurat object
genes_in_seurat <- intersect(filtered_genes, rownames(seurat_obj))
cat(sprintf("Genes found in seurat object: %d (out of %d filtered genes)\n", 
            length(genes_in_seurat), length(filtered_genes)))

if (length(genes_in_seurat) == 0) {
  stop("No matching genes found between matrix and filtered gene list!")
}

# Keep all filtered genes (will fill missing ones with 0 in training data)
genes_to_keep <- filtered_genes

# Load cell metadata (filtered cells with cell type annotations)
# The metadata file has the GSM ID prefix (e.g., GSM6416567_seurat_metadata.csv)
gsm_id <- basename(sample_dir)
metadata_path <- file.path(sample_dir, paste0(gsm_id, "_seurat_metadata.csv"))

if (!file.exists(metadata_path)) {
  stop("seurat_metadata.csv not found in ", sample_dir, ". Looking for: ", paste0(gsm_id, "_seurat_metadata.csv"))
}

cat("Loading cell metadata...\n")
cell_metadata <- read.csv(metadata_path, row.names = 1)

# Check for cell type column - look for columns like "cell_type", "celltype", "type", "ident", etc.
celltype_col <- NULL
possible_names <- c("SingleR_label", "cell_type", "celltype", "type", "ident", "seurat_clusters")

for (col_name in possible_names) {
  if (col_name %in% colnames(cell_metadata)) {
    celltype_col <- col_name
    break
  }
}

if (is.null(celltype_col)) {
  # If no common name found, display available columns and ask user
  cat("Available columns in metadata:\n")
  for (i in seq_along(colnames(cell_metadata))) {
    cat(sprintf("%d. %s\n", i, colnames(cell_metadata)[i]))
  }
  col_choice <- as.integer(readline(prompt = "Which column contains the cell type? (enter number): "))
  celltype_col <- colnames(cell_metadata)[col_choice]
}

cat(sprintf("Using column '%s' for cell types\n", celltype_col))

# Filter cells that are in the metadata
cells_to_keep <- intersect(colnames(seurat_obj), rownames(cell_metadata))
cat(sprintf("Found %d cells with annotations\n", length(cells_to_keep)))

if (length(cells_to_keep) == 0) {
  stop("No matching cells found between matrix and metadata!")
}

# Subset Seurat object (only for cells, keep all genes for now)
seurat_filtered <- subset(seurat_obj, cells = cells_to_keep)
cat(sprintf("Filtered Seurat object: %d genes x %d cells\n", 
            nrow(seurat_filtered), ncol(seurat_filtered)))

# Get the gene expression matrix (counts)
# Use LayerData() for compatibility with Assay5 objects
expr_matrix <- as.matrix(LayerData(seurat_filtered, assay = "RNA", layer = "counts"))

# Create a new matrix with all 3000 filtered genes in the correct order
# Initialize with zeros
full_expr_matrix <- matrix(0, nrow = length(genes_to_keep), ncol = ncol(expr_matrix),
                            dimnames = list(genes_to_keep, colnames(expr_matrix)))

# Fill in the counts for genes that exist in the seurat object
genes_available <- intersect(genes_to_keep, rownames(expr_matrix))
full_expr_matrix[genes_available, ] <- expr_matrix[genes_available, ]

expr_matrix <- full_expr_matrix

cat(sprintf("Full training matrix dimensions: %d genes x %d cells\n",
            nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("Matrix includes %d genes (with 0 counts for missing genes)\n", nrow(expr_matrix)))
cat("🐭 Chinchilla has prepared the perfect training matrix! 🐭\n\n")


# Get cell type information
cell_types <- cell_metadata[colnames(expr_matrix), ]
cell_type_col_name <- colnames(cell_types)[1]

# Option to filter for major cell types
major_cell_types <- c(
  "Hepatocytes", "Endothelial_cells", "Macrophage", "Smooth_muscle_cells",
  "Tissue_stem_cells", "T_cells", "Monocyte", "Fibroblasts", "NK_cell",
  "Neutrophils", "B_cell"
)
cat("\nWould you like to filter to only the 11 major cell types?\n")
cat(paste(major_cell_types, collapse = ", "), "\n")
filter_major <- readline(prompt = "Filter to major cell types only? (y/n): ")
if (tolower(filter_major) == "y") {
  keep_cells <- rownames(cell_types)[cell_types[["SingleR_label"]] %in% major_cell_types]
  cell_types <- cell_types[keep_cells, , drop = FALSE]
  expr_matrix <- expr_matrix[, keep_cells, drop = FALSE]
  cat(sprintf("Filtered to %d cells of major types.\n", ncol(expr_matrix)))
}

# Display cell type statistics
cat("\n📊 Cell Type Distribution:\n")
cell_type_counts <- table(cell_types[["SingleR_label"]])
for (ct in names(cell_type_counts)) {
  cat(sprintf("  - %s: %d cells\n", ct, cell_type_counts[ct]))
}

# Ask user for sampling method
cat("\n")
cat("Select sampling method:\n")
cat("1. Random sampling (enter total number or 'all')\n")
cat("2. Balanced sampling (equal samples per cell type)\n")
sampling_method <- readline(prompt = "Choose method (1 or 2): ")

if (sampling_method == "2") {
  # Balanced sampling - equal number per cell type
  cat("\nHow many samples per cell type? (enter number or 'all'): ")
  samples_per_type_input <- readline(prompt = "")
  
  if (tolower(samples_per_type_input) == "all") {
    samples_per_type <- 200
  } else {
    samples_per_type <- as.integer(samples_per_type_input)
  }
  
   selected_cells <- c()
   cell_type_actual <- c()
   per_type_counts <- c()
   
   for (ct in unique(cell_types[["SingleR_label"]])) {
     cells_of_type <- colnames(expr_matrix)[cell_types[["SingleR_label"]] == ct]
     available_count <- length(cells_of_type)
     
     if (is.infinite(samples_per_type)) {
       samples_to_take <- available_count
     } else {
       samples_to_take <- min(samples_per_type, available_count)
     }
     
     if (samples_to_take < samples_per_type && !is.infinite(samples_per_type)) {
       cat(sprintf("⚠️  %s: only %d cells available (requested %d)\n", ct, available_count, samples_per_type))
     }
     
     if (samples_to_take > 0) {
       selected_indices <- sample(1:length(cells_of_type), samples_to_take, replace = FALSE)
       selected_cells <- c(selected_cells, cells_of_type[selected_indices])
       cell_type_actual <- c(cell_type_actual, rep(ct, samples_to_take))
       per_type_counts <- c(per_type_counts, samples_to_take)
     }
   }
   
   cat(sprintf("Balanced sampling: requested up to %d per cell type. Actual per type: %s\n", samples_per_type, paste(per_type_counts, collapse=", ")))
   num_samples <- sum(per_type_counts)
} else {
  # Random sampling - original method
  cat("\nHow many training samples? (enter number or 'all'): ")
  num_samples_input <- readline(prompt = "")
  
  if (tolower(num_samples_input) == "all") {
    num_samples <- ncol(expr_matrix)
    selected_cells <- colnames(expr_matrix)
  } else {
    num_samples <- as.integer(num_samples_input)
    
    if (is.na(num_samples) || num_samples < 1) {
      stop("Invalid number of samples!")
    }
    
    if (num_samples > ncol(expr_matrix)) {
      cat(sprintf("Warning: Requested %d samples but only %d cells available. Using all.\n", 
                  num_samples, ncol(expr_matrix)))
      num_samples <- ncol(expr_matrix)
    }
    
    # Randomly select cells
    selected_indices <- sample(1:ncol(expr_matrix), num_samples, replace = FALSE)
    selected_cells <- colnames(expr_matrix)[selected_indices]
  }
}

# Ask user for output format
cat("\n")
cat("Select output format:\n")
cat("1. Separate files (cellbarcode.csv + cellbarcode_type.csv)\n")
cat("2. Combined files (all_counts.csv + all_types.csv)\n")
cat("3. Single file (training_data.csv with counts and type in last column)\n")
output_format <- readline(prompt = "Choose format (1, 2, or 3): ")

cat(sprintf("\nGenerating training data for %d cells...\n", length(selected_cells)))

# Initialize data structures for combined/single file formats
all_counts_matrix <- NULL
all_types_vector <- c()

# Generate training data files
num_generated <- 0

for (cell_barcode in selected_cells) {
  # Get gene counts for this cell
  gene_counts <- expr_matrix[, cell_barcode]
  
  # Get cell type from the correct column
  cell_type <- as.character(cell_metadata[cell_barcode, celltype_col])
  
  if (output_format == "1") {
    # Format 1: Separate files per cell
    counts_filename <- file.path(training_data_dir, paste0(cell_barcode, ".csv"))
    write.table(as.matrix(gene_counts), file = counts_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    type_filename <- file.path(training_data_dir, paste0(cell_barcode, "_type.csv"))
    write.table(cell_type, file = type_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  } else if (output_format == "2" || output_format == "3") {
    # Format 2 & 3: Build matrices with cells as columns
    all_counts_matrix <- cbind(all_counts_matrix, gene_counts)
    all_types_vector <- c(all_types_vector, cell_type)
  }
  
  num_generated <- num_generated + 1
  
  if (num_generated %% 100 == 0) {
    cat(sprintf("Generated training data for %d cells...\n", num_generated))
  }
}

# Save combined/single file formats
if (output_format == "2") {
  # Format 2: Separate combined files
  # all_counts.csv: genes (rows) x cells (columns)
  all_counts_filename <- file.path(training_data_dir, "all_counts.csv")
  write.table(all_counts_matrix, file = all_counts_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # all_types.csv: one row with all cell types
  all_types_filename <- file.path(training_data_dir, "all_types.csv")
  write.table(all_types_vector, file = all_types_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE, nrow = 1)
  
  cat(sprintf("Saved combined files:\n  - %s (genes x cells)\n  - %s (cell types row)\n", all_counts_filename, all_types_filename))
  
} else if (output_format == "3") {
  # Format 3: Single file with counts and types
  # Combine counts matrix with types row at the bottom
  combined_data <- rbind(all_counts_matrix, all_types_vector)
  
  single_filename <- file.path(training_data_dir, "training_data.csv")
  write.table(combined_data, file = single_filename, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat(sprintf("Saved single file: %s (genes x cells, with cell types row at bottom)\n", single_filename))
}

cat(sprintf("\n========================================\n"))
cat(sprintf("Successfully generated %d training samples!\n", num_generated))
cat(sprintf("Saved to: %s\n", training_data_dir))
cat(sprintf("========================================\n"))
cat("🐭 Chinchilla's work is complete! *happy squeaks* 🐭\n\n")

# Summary statistics
cat("Summary Statistics:\n")
cat(sprintf("Sample: %s\n", selected_sample))
cat(sprintf("Total cells processed: %d\n", num_generated))
cat(sprintf("Genes per cell: %d\n", length(genes_to_keep)))
cat(sprintf("Cell types: %s\n", paste(unique(cell_types[, 1]), collapse = ", ")))

}

# ============================================================================
# To run the script, call:
generateTrainingData()
# ============================================================================
