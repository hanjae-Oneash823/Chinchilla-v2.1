# Training Data Documentation

## Overview
This document describes the training datasets used for hepatic cell type annotation in the scRNA-seq analysis pipeline.

## Primary Training Dataset

**Dataset ID:** GSE1854577

**Title:** Single-Cell, Single-Nucleus and Spatial RNA Sequencing of the Human Liver Identifies Hepatic Stellate Cell and Cholangiocyte Heterogeneity

**Description:**
- Comprehensive single-cell RNA sequencing data from human liver tissue
- Multiple donors (MacParland et al., 2,3,4 tissues) with varying tissue types
- Combined single-cell (scRNA-seq) and single-nucleus (snRNA-seq) sequencing approaches
- Spatial transcriptomics data integration (Visium technology) for spatial context

**Cell Types Identified:**
- ~73,295 total cells/nuclei sequenced
- Major cell types include:
  - **Hepatocytes** - primary liver parenchymal cells
  - **Cholangiocytes** - biliary epithelial cells
  - **Mesenchymal cells** - stromal cell populations
  - Additional immune and endothelial populations

**Data Characteristics:**
- High-quality coverage across major liver cell types
- Comprehensive hepatocyte and cholangiocyte heterogeneity
- Mesenchymal cell diversity well-represented
- Suitable as a reference for liver-specific cell type annotation

## Secondary Training Dataset

**Dataset Source:** Human Cell Atlas (HCA)

**Sample IDs:** 
- GSM6416567
- GSM6416569
- GSM8493744
- GSM8493745
- GSM8493746

**Purpose:**
- Complementary liver samples for cross-validation
- Broader representation of human liver cell diversity
- Integration with HCA standardized annotations

## Data Format

**File Format:** CSV (Comma-Separated Values)

**Structure:**
- Each file represents a single cell
- Columns: Gene features (expression counts)
- Rows: Individual genes/features
- Values: Normalized or raw expression counts

**Specifications:**
- Format: CSV
- Features: Multiple (n features per cell)
- Values: Expression counts or normalized values
- One sample per file for modular analysis

## Usage in Pipeline

1. **Reference Construction:** Primary and secondary datasets serve as reference for cell type annotation
2. **Annotation Method:** SingleR (Single-cell Recognition) uses these datasets for probabilistic cell type assignment
3. **Cell Type Resolution:** 14 major liver cell types annotated:
   - Hepatocytes
   - Cholangiocytes
   - LSECs (Liver Sinusoidal Endothelial Cells)
   - Macro. Endothelial cells
   - Stellate Cells
   - Fibroblasts
   - Kupffer Cells
   - Macrophages
   - Dendritic Cells
   - T Cells
   - NK/NKT Cells
   - B Cells
   - Plasma Cells
   - Erythroid Cells

## References

- GSE1854577: MacParland et al. - Human liver single-cell atlas
- HCA (Human Cell Atlas): Standardized reference for human cell types
- SingleR: Aran et al. - Automated cell type annotation reference-based method
