# ===================================================
# ===================================================
# scRNAseq_pipeline.R
# Install required packages if not already installed
cran_packages <- c("Seurat", "Matrix")
bioc_packages <- c("SingleR", "celldex", "SingleCellExperiment")
# Install CRAN packages
for (pkg in cran_packages) {
	if (!requireNamespace(pkg, quietly = TRUE)) {
		install.packages(pkg, repos = "https://cloud.r-project.org/")
	}
}
# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
	install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}
# Install Bioconductor packages
for (pkg in bioc_packages) {
	if (!requireNamespace(pkg, quietly = TRUE)) { 
		BiocManager::install(pkg, ask = FALSE)
	}
}
# Load required libraries
library(Seurat)
library(SingleR)
library(celldex)
library(Matrix)
# ===================================================
# ===================================================

# Set data directory for 10X Genomics files
data_dir <- "HCA_liverdata/GSM6416567/"  # Path to folder containing barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz

# Read 10X Genomics data
seurat_obj <- CreateSeuratObject(counts = Read10X(data.dir = data_dir))

# Extract project ID from data directory
project_id <- sub(".*/", "", sub("/$", "", data_dir))

# Export sparse matrix data for feature filtering
# 1. Save sparse counts matrix
counts_matrix <- GetAssayData(seurat_obj, layer = "counts")
saveRDS(counts_matrix, file = paste0(project_id, "_counts_matrix.rds"))
# 2. Export cell metadata (barcodes + project)
cell_metadata <- data.frame(
  barcode = colnames(seurat_obj),
  project = project_id,
  stringsAsFactors = FALSE
)
write.csv(cell_metadata, file = paste0(project_id, "_cell_metadata.csv"), row.names = FALSE)
# 3. Export features (gene names)
features_df <- data.frame(feature = rownames(seurat_obj))
write.csv(features_df, file = paste0(project_id, "_features.csv"), row.names = FALSE)

# QC: Calculate percent mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics (optional)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Dimensionality reduction
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Cell type annotation using SingleR and celldex reference
ref <- celldex::HumanPrimaryCellAtlasData()
seurat_obj_sce <- as.SingleCellExperiment(seurat_obj)
annots <- SingleR(test = seurat_obj_sce, ref = ref, labels = ref$label.main)
seurat_obj$SingleR_label <- annots$labels

# Save results
saveRDS(seurat_obj, file = "seurat_annotated.rds")
write.csv(seurat_obj@meta.data, file = "seurat_metadata.csv")

# Optional: UMAP plot colored by annotation
DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE)
