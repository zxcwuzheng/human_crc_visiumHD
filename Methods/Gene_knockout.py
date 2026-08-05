import os
import scanpy as sc
from GenKI import GenKI


# using scanpy to load ST dataset before
adata_merged_T = ...
# Save the raw counts in layers
adata_merged_T.layers['counts'] = adata_merged_T.X
# Filter cells with at least 50 counts
sc.pp.filter_cells(adata_merged_T, min_counts=50)
# Filter genes detected in at least 2% of cells
sc.pp.filter_genes(adata_merged_T, min_cells=adata_merged_T.shape[0]*0.02) # adata_merged_T.shape[0]*0.05
# Mark mitochondrial genes
adata_merged_T.var['mt'] = adata_merged_T.var_names.str.startswith(('MT-', 'mt-'))
# Remove mitochondrial genes
adata_merged_T = adata_merged_T[:, ~adata_merged_T.var['mt']].copy()
# Normalize total counts to 1e6
sc.pp.normalize_total(adata_merged_T, target_sum=1e6)
# Log transform the data
sc.pp.log1p(adata_merged_T)
# Identify highly variable genes, 5000 genes selected
sc.pp.highly_variable_genes(adata_merged_T,n_top_genes=5000,flavor="seurat",batch_key="batch")

# Check if target gene CD44 is among highly variable genes
'CD44' in adata_merged_T.var.index[adata_merged_T.var['highly_variable']]

# Restore raw counts to X for GenKI analysis
adata_merged_T.X = adata_merged_T.layers['counts']


# Initialize GenKI for gene knockout analysis, target is CD44
gk = GenKI.from_adata(adata_merged_T[:,adata_merged_T.var['highly_variable']].copy(), target_gene=["CD44"], preprocess=True)
# Run analysis with 100 permutations, using KL divergence for ranking, random seed 2026
ranked = gk.run(seed=2026,n_permutations=100,by="KL")
# Print the ranked results
print(ranked)