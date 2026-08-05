library(Seurat)
library(tidyverse)
library(WGCNA)
library(hdWGCNA)
# Enable multi-threading for WGCNA to speed up computations
enableWGCNAThreads(nThreads = 12)
# Set random seed for reproducibility of results
set.seed(12345)


# Combine spatial coordinates from multiple tissue slices
coords <- rbind(GetTissueCoordinates(combined_obj,'slice1'),
                GetTissueCoordinates(combined_obj,'slice1.2'),
                GetTissueCoordinates(combined_obj,'slice1.3'))
colnames(coords) <- c('spatial_1', 'spatial_2')
coords <- as.matrix(coords)
combined_obj@reductions$spatial <- CreateDimReducObject(
  embeddings = coords,
  assay = 'RNA',
  key = 'spatial_'  
)

# Normalize data before WGCNA analysis
combined_obj <- NormalizeData(combined_obj)

# Initialize hdWGCNA object with specified features (bulk-derived biomarkers)
combined_obj <- SetupForWGCNA(
  combined_obj,
  features = marker, # Bulk biomarkers
  wgcna_name = "All",
)

# Create metacells by aggregating spots 
# Group by both slice (sample) and cell type
combined_obj <- MetacellsByGroups(
  seurat_obj = combined_obj,
  group.by = c("slice","DeconvolutionLabel1"), # sample, cell type
  ident.group = 'DeconvolutionLabel1',
  reduction = 'spatial', # Use combined spatial coordinates for clustering
  k = 50, min_cells = 50,max_iter=2000, 
  max_shared = 15 # Maximum overlap allowed between metacells
)

# Prepare expression matrix for WGCNA 
combined_obj  <- SetDatExpr(combined_obj)

# Test different soft-thresholding powers to select optimal power
combined_obj <- TestSoftPowers(combined_obj)

# Construct the gene co-expression network and identify modules
combined_obj <- ConstructNetwork(
  combined_obj,
  overwrite_tom = TRUE, 
  tom_name='All_marker', 
  minModuleSize = 10, # Minimum number of genes required to form a module
  detectCutHeight = 0.9995 # Cut height for module detection
)

# Calculate module eigengenes
combined_obj <- ModuleEigengenes(combined_obj)

# Calculate connectivity (kME)
combined_obj <- ModuleConnectivity(
  combined_obj
)

# Plot the hierarchical clustering dendrogram
PlotDendrogram(combined_obj, main='Visium HD hdWGCNA dendrogram')


# Integration with protein-protein interaction (PPI) network is documented in Fig4.ipynb