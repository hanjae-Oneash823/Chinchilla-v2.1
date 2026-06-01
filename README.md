# Chinchilla-v2.1
Dual-stream VAE for liver scRNA-seq cell type annotation — trained on Human Cell Atlas data, classifies 14 hepatic cell types.

A PyTorch model that jointly learns a compressed latent representation of gene expression profiles (VAE) and predicts cell types (MLP classifier) in a single training loop. Trained on Human Cell Atlas liver data (5 samples, ~73k cells) preprocessed through a 5-step R pipeline (Seurat → SingleR → feature filtering → voting-based gene selection). Classifies 14 hepatic cell types including Hepatocytes, Cholangiocytes, LSECs, Kupffer Cells, and Stellate Cells.

<img width="3300" height="1800" alt="cell_type_counts_by_sample" src="https://github.com/user-attachments/assets/614c4e3b-2ea4-4d6d-9896-002775df30fd" />
