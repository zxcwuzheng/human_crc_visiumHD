library(Seurat)
library(spacexr)

#----------------------------
# 1. Deconvolution using RCTD 
#----------------------------
# Load reference single-cell RNA-seq count matrix 
ref <- Read10X_h5('HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5')
# Load metadata containing cell type annotations 
MetaData<-read.csv('HumanColonCancer_VisiumHD-main/MetaData/SingleCell_MetaData.csv.gz')

# Filter reference to keep only cell types with more than 25 cells
KpIdents<-names(which(table(MetaData$Level1)>25))
MetaData<-MetaData[MetaData$Level1%in%KpIdents,]
ref<-ref[,MetaData$Barcode]

# Fix cell type labels: replace '/' with '_' to avoid path issues
CTRef<-MetaData$Level1
CTRef<-gsub("/","_",CTRef)
CTRef<-as.factor(CTRef)
# Name the vector with cell barcodes
names(CTRef)<-MetaData$Barcode

# Build RCTD reference object from single-cell data
reference <- Reference(ref[,names(CTRef)], CTRef , colSums(ref),n_max_cells = 30000)

# Prepare Spatial transcriptomics  query object
counts_hd <- object[["Spatial.Polygons"]]$counts
cortex_cells_hd <- colnames(object[["Spatial.Polygons"]])
coords <- GetTissueCoordinates(object)[cortex_cells_hd, 1:2]

# Create RCTD SpatialRNA object for query
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))


# Run RCTD deconvolution with doublet mode 
RCTD <- create.RCTD(query, reference, max_cores = 30)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# Add deconvolution results to Seurat object metadata
object <- AddMetaData(object, metadata = RCTD@results$results_df)


# Filter low-quality deconvolution results
object$first_type <- as.character(object$first_type)
object$first_type[is.na(object$first_type)] <- "Unknown"
# Keep only singlet predictions
filter_obj <- object[,object$spot_class=='singlet']
# Remove QC filtered spots
filter_obj <- filter_obj[,filter_obj$first_type !='QC_Filtered']


#---------------------------------
# Identification of spatial regions/niche
#--------------------------------
# Build niche assay using Seurat's BuildNicheAssay function
# group.by: use deconvolved cell type labels
# niches.k: number of spatial niches to identify
# neighbors.k: number of neighbors for constructing spatial graph
filter_obj2 <- BuildNicheAssay(object = filter_obj, fov = "slice1.polygons", group.by = "first_type",
    niches.k = 5, neighbors.k = 30)