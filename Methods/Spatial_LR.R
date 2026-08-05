# Load required libraries
library(SpaGene)
library(Seurat)

# Load ligand-receptor pair database for human
load("LRpair_human.rds")

# Extract count matrix from the spatial assay of the Seurat object
count <- GetAssayData(obj, assay = "spatial", slot = "counts")

# Get spatial coordinates of each spot and convert to matrix
loc <- as.matrix(GetTissueCoordinates(obj))

# Run SpaGene to identify spatially ligand-receptor pairs
lr <- SpaGene_LR(count, loc, LRpair=LRpair)  

# Filter results to keep only significant interactions 
lr <- lr[lr$adj<0.05,]

# Sort by significance and take top 50 interactions
lr <- head(lr[order(lr$adj),], 50)