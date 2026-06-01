import pandas as pd
import scipy.io
import scipy.sparse
import numpy as np
import scanpy as sc
import celltypist

# Load barcodes, features, and matrix
data_dir = 'data/'
barcodes = pd.read_csv(f'{data_dir}barcodes.tsv.gz', header=None)[0].tolist()
features = pd.read_csv(f'{data_dir}features.tsv.gz', sep='\t', header=None)
matrix = scipy.io.mmread(f'{data_dir}matrix.mtx.gz').tocsc()

# Load data into AnnData
adata = sc.read_10x_mtx(data_dir, var_names='gene_symbols', cache=True)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells and genes
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]
adata = adata[:, adata.var.n_cells_by_counts > 3]

# Normalize, log-transform, and find highly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]
sc.pp.scale(adata, max_value=10)

# --- scANVI cell type annotation (unsupervised demo) ---
import scvi
import anndata

# Prepare AnnData for scvi-tools
print("\nRunning scANVI (scvi-tools) for cell type annotation...")

# Deep cleaning for NaN/inf before scANVI

# --- CellTypist annotation (following official tutorial) ---
import celltypist

# Deep cleaning for NaN/inf before CellTypist
adata_ct = adata.copy()
print('\nBefore cleaning for CellTypist:')
print(f'  Cells: {adata_ct.shape[0]}, Genes: {adata_ct.shape[1]}')
print(f'  Any NaN in X: {np.isnan(adata_ct.X).any()}')
print(f'  Any inf in X: {np.isinf(adata_ct.X).any()}')

# Remove cells with any NaN or inf values
cell_mask = ~np.isnan(adata_ct.X).any(axis=1) & ~np.isinf(adata_ct.X).any(axis=1)
adata_ct = adata_ct[cell_mask, :]
# Remove genes with any NaN or inf values
gene_mask = ~np.isnan(adata_ct.X).any(axis=0) & ~np.isinf(adata_ct.X).any(axis=0)
adata_ct = adata_ct[:, gene_mask]

print('After removing NaN/inf cells/genes for CellTypist:')
print(f'  Cells: {adata_ct.shape[0]}, Genes: {adata_ct.shape[1]}')
print(f'  Any NaN in X: {np.isnan(adata_ct.X).any()}')
print(f'  Any inf in X: {np.isinf(adata_ct.X).any()}')

# Remove any remaining NaN/inf (rare, but deep check)
if np.isnan(adata_ct.X).any() or np.isinf(adata_ct.X).any():
	print('Deep cleaning for CellTypist: removing all rows/columns with any NaN/inf (may reduce data size)')
	mask_rows = ~(np.isnan(adata_ct.X).any(axis=1) | np.isinf(adata_ct.X).any(axis=1))
	adata_ct = adata_ct[mask_rows, :]
	mask_cols = ~(np.isnan(adata_ct.X).any(axis=0) | np.isinf(adata_ct.X).any(axis=0))
	adata_ct = adata_ct[:, mask_cols]
	print(f'  After deep cleaning: Cells: {adata_ct.shape[0]}, Genes: {adata_ct.shape[1]}')
	print(f'  Any NaN in X: {np.isnan(adata_ct.X).any()}')
	print(f'  Any inf in X: {np.isinf(adata_ct.X).any()}')

# Final check
if np.isnan(adata_ct.X).any() or np.isinf(adata_ct.X).any():
	raise ValueError('Data still contains NaN or inf after all cleaning steps for CellTypist. Please check your input data.')

# Ensure adata_ct is a copy, not a view
adata_ct = adata_ct.copy()

# --- CellTypist expects log1p normalized expression to 10,000 counts per cell ---
import scanpy as sc
print('Normalizing total counts per cell to 10,000 and log1p transforming for CellTypist...')
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)

print('\nRunning CellTypist for cell type annotation...')
# Download the latest model (or use a specific model if desired)
model = celltypist.models.download_models('Immune_All_Low.pkl')

# Run CellTypist prediction
pred = celltypist.annotate(adata_ct, model=model, majority_voting=True)

# Add predictions to AnnData object
adata_ct = pred.to_adata()

# Save predictions to CSV
celltypist_df = adata_ct.obs[['cell_type', 'majority_voting']]
celltypist_df.to_csv('celltypist_annotations.csv')
print('CellTypist annotation complete. Results saved to celltypist_annotations.csv')

