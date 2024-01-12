### Load in the data ###
library(Seurat)
library(ggplot2)
library(patchwork)

# Load in the RNA UMI matrix
cbmc.rna <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))

cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file="GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))


all.equal(colnames(cbmc.rna), colnames(cbmc.adt))


### Setup a Seurat object, add the RNA and protein data ### 
# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

Assays(cbmc)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)

# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])

# List the current default assay
DefaultAssay(cbmc)

# Switch the default to ADT
DefaultAssay(cbmc) <- "ADT"
DefaultAssay(cbmc)


### Cluster cells on the basis on their scRNA-seq profiles ###
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)

# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(cbmc, label = TRUE)


### Visualize multiple modalities side-by-side ###
# Normalize ADT data,
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
# cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2

# Alternative method
# Identify the key for the RNA and protein assays
Key(cbmc[["RNA"]])
Key(cbmc[["ADT"]])

# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2


### Identify cell surface markers for scRNA-seq clusters ###
# CD19 is a B cell marker, identify cluster 6 as expressing CD19 on the surface
VlnPlot(cbmc, "adt_CD19")

# also identify alternative protein and RNA markers for this cluster through
# differential expression
adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")

head(adt_markers)
head(rna_markers)


### Additional visualizations of multimodal data ###
# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")

# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")

FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")

# Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# particularly in comparison to RNA values. This is due to the significantly higher protein
# copy number in cells, which significantly reduces 'drop-out' in ADT data
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")