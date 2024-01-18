integrative_analysis
================
Matthew Esqueda
2024-01-18

## Introduction

Integration of single-cell datasets, for example across experimental
batches, donors, or conditions, is often an important step in scRNA-seq
workflows. Itegrative analysis can help to match shared cell types and
states across datasets, which can boost statistical power, and most
importantly, facilitate accurate comparative analysis across datasets.
Here, a [dataset of human PBMC profiled with seven different
technologies](https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data),
profiled as part of a systematic comparative analysis (`pbmcsa`). The
data is available as part of the
[SeuratData](https://github.com/satijalab/seurat-data) package.

## Layers in the Seurat v5 object

Seurat v5 assays store data in layers. These layers can store raw,
un-normalized (`layer='counts'`), normalized data (`layer='data'`), or
z-scored/variance-stabilized data (`layer='scale.data'`). Load the data,
remove low-quality cells, and obtain predicted cell annotations with the
[Azimuth
pipeline](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html).

``` r
# load in pbmc systematic comparative analysis dataset
obj <- LoadData("pbmcsca")
```

    ## Validating object structure

    ## Updating object slots

    ## Ensuring keys are in the proper structure

    ## Warning: Assay RNA changing from Assay to Assay

    ## Ensuring keys are in the proper structure

    ## Ensuring feature names don't have underscores or pipes

    ## Updating slots in RNA

    ## Validating object structure for Assay 'RNA'

    ## Object representation is consistent with the most current Seurat version

    ## Warning: Assay RNA changing from Assay to Assay5

``` r
obj <- subset(obj, nFeature_RNA > 1000)
obj <- RunAzimuth(obj, reference = "pbmcref")
```

    ## Warning: Overwriting miscellanous data for model

    ## Warning: Adding a dimensional reduction (refUMAP) without the associated assay
    ## being present

    ## Warning: Adding a dimensional reduction (refUMAP) without the associated assay
    ## being present

    ## detected inputs from HUMAN with id type Gene.name

    ## reference rownames detected HUMAN with id type Gene.name

    ## Normalizing query using reference SCT model

    ## Warning: 735 features of the features specified were not present in both the reference query assays. 
    ## Continuing with remaining 4265 features.

    ## Projecting cell embeddings

    ## Finding query neighbors

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 13407 anchors

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Predicting cell labels
    ## Predicting cell labels

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Predicting cell labels

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## 
    ## Integrating dataset 2 with reference dataset

    ## Finding integration vectors

    ## Integrating data

    ## Warning: Keys should be one or more alphanumeric characters followed by an
    ## underscore, setting key from integrated_dr_ to integrateddr_

    ## Computing nearest neighbors

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## Running UMAP projection

    ## Warning in RunUMAP.default(object = neighborlist, reduction.model =
    ## reduction.model, : Number of neighbors between query and reference is not equal
    ## to the number of neighbors within reference

    ## 00:17:58 Read 10434 rows

    ## 00:17:58 Processing block 1 of 1

    ## 00:17:58 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 20
    ## 00:17:58 Initializing by weighted average of neighbor coordinates using 1 thread
    ## 00:17:58 Commencing optimization for 67 epochs, with 208680 positive edges
    ## 00:18:01 Finished

    ## Warning: No assay specified, setting assay as RNA by default.

    ## Projecting reference PCA onto query
    ## Finding integration vector weights
    ## Projecting back the query cells into original PCA space
    ## Finding integration vector weights
    ## Computing scores:
    ##     Finding neighbors of original query cells
    ##     Finding neighbors of transformed query cells
    ##     Computing query SNN
    ##     Determining bandwidth and computing transition probabilities
    ## Total elapsed time: 7.91480302810669

``` r
# currently, the object has two layers in the RNA assay: counts and data
obj
```

    ## An object of class Seurat 
    ## 33789 features across 10434 samples within 4 assays 
    ## Active assay: RNA (33694 features, 0 variable features)
    ##  2 layers present: counts, data
    ##  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
    ##  2 dimensional reductions calculated: integrated_dr, ref.umap

The object contains data from nine different batches (stores in the
`Method` column in the object metadata), representing seven different
technologies. In Seurat v5, integrate the different batches together and
keep all the data in one object, but split into layers. After splitting,
there are now 18 layers (a `counts` and `data` layer for each batch). We
can also run a standard scRNA-seq analysis (i.e. without integration).
Note that since the data is split into layers, normalization and
variable feature identification is performed for each batch
independently (a consensus set of variable features is automatically
identified).

``` r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method) 
obj
```

    ## An object of class Seurat 
    ## 33789 features across 10434 samples within 4 assays 
    ## Active assay: RNA (33694 features, 0 variable features)
    ##  18 layers present: counts.Smart-seq2, counts.CEL-Seq2, counts.10x_Chromium_v2_A, counts.10x_Chromium_v2_B, counts.10x_Chromium_v3, counts.Drop-seq, counts.Seq-Well, counts.inDrops, counts.10x_Chromium_v2, data.Smart-seq2, data.CEL-Seq2, data.10x_Chromium_v2_A, data.10x_Chromium_v2_B, data.10x_Chromium_v3, data.Drop-seq, data.Seq-Well, data.inDrops, data.10x_Chromium_v2
    ##  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
    ##  2 dimensional reductions calculated: integrated_dr, ref.umap

``` r
obj <- NormalizeData(obj)
```

    ## Normalizing layer: counts.Smart-seq2

    ## Normalizing layer: counts.CEL-Seq2

    ## Normalizing layer: counts.10x_Chromium_v2_A

    ## Normalizing layer: counts.10x_Chromium_v2_B

    ## Normalizing layer: counts.10x_Chromium_v3

    ## Normalizing layer: counts.Drop-seq

    ## Normalizing layer: counts.Seq-Well

    ## Normalizing layer: counts.inDrops

    ## Normalizing layer: counts.10x_Chromium_v2

``` r
obj <- FindVariableFeatures(obj)
```

    ## Finding variable features for layer counts.Smart-seq2

    ## Finding variable features for layer counts.CEL-Seq2

    ## Finding variable features for layer counts.10x_Chromium_v2_A

    ## Finding variable features for layer counts.10x_Chromium_v2_B

    ## Finding variable features for layer counts.10x_Chromium_v3

    ## Finding variable features for layer counts.Drop-seq

    ## Finding variable features for layer counts.Seq-Well

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at -2.5648

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 0.49783

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 2.2129e-14

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 0.090619

    ## Finding variable features for layer counts.inDrops

    ## Finding variable features for layer counts.10x_Chromium_v2

``` r
obj <- ScaleData(obj)
```

    ## Centering and scaling data matrix

``` r
obj <- RunPCA(obj)
```

    ## PC_ 1 
    ## Positive:  LYZ, FCN1, CST3, S100A9, SERPINA1, CTSS, VCAN, CD68, AIF1, PSAP 
    ##     S100A8, LST1, GRN, SPI1, CLEC7A, FCER1G, NCF2, CSTA, TYROBP, TYMP 
    ##     CD14, FGL2, CSF3R, MPEG1, CPVL, COTL1, MS4A6A, CFD, HCK, MAFB 
    ## Negative:  IL32, MALAT1, TRBC2, CD3D, TRAC, CCL5, TRBC1, IL7R, RPS18, CTSW 
    ##     CST7, CD247, CD3G, CD2, LDHB, GZMA, SYNE2, RORA, LTB, NKG7 
    ##     CD69, SPOCK2, CD7, ITM2A, LYAR, RPS8, GIMAP7, PRF1, CD8A, CXCR4 
    ## PC_ 2 
    ## Positive:  NKG7, S100A4, SRGN, IL32, CST7, CCL5, ANXA1, CTSW, GZMA, GZMH 
    ##     FGFBP2, PRF1, ITGB2, GNLY, CD247, ID2, CD3D, KLRD1, ARL4C, GZMB 
    ##     GAPDH, PFN1, CD7, S100A6, CD8A, LYAR, CD3G, RORA, TRBC1, HOPX 
    ## Negative:  MS4A1, CD79A, BANK1, IGHM, HLA-DQA1, CD79B, HLA-DQB1, IGKC, LINC00926, HLA-DRA 
    ##     CD74, CD22, IGHD, VPREB3, TNFRSF13C, MEF2C, IGHG1, RALGPS2, HLA-DMB, HLA-DPA1 
    ##     BLK, HLA-DPB1, JCHAIN, ADAM28, HLA-DRB1, FAM129C, P2RX5, FCRL1, CD24, IGLC2 
    ## PC_ 3 
    ## Positive:  IL7R, RPLP1, RPS8, RPS18, LTB, MAL, LEF1, RCAN3, LDHB, EEF1A1 
    ##     CCR7, RPS2, JUNB, NOSIP, CAMK4, PIK3IP1, TCF7, RPLP0, VIM, TRAC 
    ##     TRABD2A, NELL2, TRAT1, RGCC, AQP3, ZFP36L2, INPP4B, RGS10, CD40LG, CD27 
    ## Negative:  PRF1, NKG7, GZMB, FGFBP2, KLRD1, GNLY, FCGR3A, CST7, GZMH, KLRF1 
    ##     SPON2, GZMA, ADGRG1, CTSW, MYOM2, CCL4, S1PR5, CLIC3, CCL5, HOPX 
    ##     TTC38, TRDC, IL2RB, PRSS23, CX3CR1, FCRL6, SH2D1B, IFITM2, ITGB2, PTGDR 
    ## PC_ 4 
    ## Positive:  MTRNR2L12, MTRNR2L1, GPM6A, KCNJ3, AP000769.1, XIST, MGST2, SLFN5, OAS2, GIMAP7 
    ##     GIMAP4, FYB, LINC00861, IL7R, IKZF2, OAS3, SAMHD1, CX3CR1, LBH, SPN 
    ##     GBP5, MYOM2, MTRNR2L8, CASP8AP2, SELL, NCL, FOXP1, TCF7, SYNE2, MLLT6 
    ## Negative:  RPS2, RPS4Y1, RPLP1, FTH1, RPS18, TMSB10, LY6E, MALAT1, FTL, SUB1 
    ##     GAPDH, LGALS1, HERPUD1, JUNB, SNX3, SEC61B, S100A4, S100A6, C12orf75, HLA-DPB1 
    ##     RPS8, H1FX, PPP1R14B, DDIT4, HLA-DRB1, HLA-DPA1, MAP3K8, S100A10, MT2A, PRELID1 
    ## PC_ 5 
    ## Positive:  NFKBIA, CXCL8, TNFAIP3, THBS1, EGR1, NR4A2, SLC2A3, G0S2, IER3, IL1B 
    ##     VCAN, DUSP1, MTRNR2L8, CXCR4, SGK1, PTGS2, RPS4Y1, FOSB, BHLHE40, CD83 
    ##     DDX3Y, MTRNR2L12, JUN, CYP1B1, KLF6, HIF1A, CD14, KLF10, ZFP36, MT-CO1 
    ## Negative:  TMSB4X, CDKN1C, CSF1R, FCGR3A, LINC01272, MS4A7, CKB, HES4, CSTB, RHOC 
    ##     AIF1, LST1, TMSB10, LRRC25, FAM110A, XIST, GIMAP7, GIMAP4, PYCARD, RPS2 
    ##     ABI3, BATF3, RPS18, TCF7L2, RNASET2, CORO1A, SIGLEC10, ARPC1B, LYPD2, ZNF703

Visualize the results of a standard analysis without integration. Cells
are grouping both by cell type and by underlying method. A UMAP analysis
is a visualization of this, clustering the dataset would return
predominantly batch-specific clusters. Especially if previous cell-type
annotations were not available, this would make downstream analysis
extremely challenging.

``` r
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10434
    ## Number of edges: 412660
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8981
    ## Number of communities: 48
    ## Elapsed time: 0 seconds

``` r
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
```

    ## 00:18:48 UMAP embedding parameters a = 0.9922 b = 1.112

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:18:48 Read 10434 rows and found 30 numeric columns

    ## 00:18:48 Using Annoy for neighbor search, n_neighbors = 30

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:18:48 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 00:18:49 Writing NN index file to temp file C:\Users\matta\AppData\Local\Temp\Rtmp8cVbQ1\filea76c2d157ba5
    ## 00:18:49 Searching Annoy index using 1 thread, search_k = 3000
    ## 00:18:51 Annoy recall = 100%
    ## 00:18:52 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 00:18:55 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 00:18:55 Commencing optimization for 200 epochs, with 428980 positive edges
    ## 00:19:04 Optimization finished

``` r
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))
```

![](seurat5_integration_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Perform streamlined (one-line) integrative analysis

Seurat v5 enables streamlined integrative analysis using the
`IntegrateLaters` function. The method currently supports five
integration methods. Each of these methods performs integration in
low-dimensional space, and returns a dimensional reduction
(i.e. `integrated.rcpa`) that aims to co-embed shared cell types across
batches: - Anchor-based CCA integration (method=CCAIntegration) -
Anchor-based RPCA integration (method=RPCAIntegration) - Harmony
(method=HarmonyIntegration) - FastMNN (method= FastMNNIntegration) -
scVI (method=scVIIntegration) The anchor-based RPCA integration
represents a faster and more conservative (less correction) method for
integration.

Each of the following lines perform a new integration using a single
line of code

``` r
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
```

``` r
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
```

``` r
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
```

    ## Warning: HarmonyMatrix is deprecated and will be removed in the future from the
    ## API in the future

    ## Warning: Warning: The parameters do_pca and npcs are deprecated. They will be ignored for this function call and please remove parameters do_pca and npcs and pass to harmony cell_embeddings directly.
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter tau is deprecated. It will be ignored for this function call and please remove parameter tau in future function calls. Advanced users can set value of parameter tau by using parameter .options and function harmony_options().
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter block.size is deprecated. It will be ignored for this function call and please remove parameter block.size in future function calls. Advanced users can set value of parameter block.size by using parameter .options and function harmony_options().
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter max.iter.harmony is replaced with parameter max_iter. It will be ignored for this function call and please use parameter max_iter in future function calls.
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter max.iter.cluster is deprecated. It will be ignored for this function call and please remove parameter max.iter.cluster in future function calls. Advanced users can set value of parameter max.iter.cluster by using parameter .options and function harmony_options().
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter epsilon.cluster is deprecated. It will be ignored for this function call and please remove parameter epsilon.cluster in future function calls. Advanced users can set value of parameter epsilon.cluster by using parameter .options and function harmony_options().
    ## This warning is displayed once per session.

    ## Warning: Warning: The parameter epsilon.harmony is deprecated. It will be ignored for this function call and please remove parameter epsilon.harmony in future function calls. If users want to control if harmony would stop early or not, use parameter early_stop. Advanced users can set value of parameter epsilon.harmony by using parameter .options and function harmony_options().
    ## This warning is displayed once per session.

Visualize and cluster the datasets for CCI integration and RPCA
integration.

``` r
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10434
    ## Number of edges: 614214
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8056
    ## Number of communities: 25
    ## Elapsed time: 1 seconds

``` r
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
```

    ## 00:23:24 UMAP embedding parameters a = 0.9922 b = 1.112

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:23:24 Read 10434 rows and found 30 numeric columns

    ## 00:23:24 Using Annoy for neighbor search, n_neighbors = 30

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:23:24 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 00:23:25 Writing NN index file to temp file C:\Users\matta\AppData\Local\Temp\Rtmp8cVbQ1\filea76c15337b73
    ## 00:23:25 Searching Annoy index using 1 thread, search_k = 3000
    ## 00:23:27 Annoy recall = 100%
    ## 00:23:29 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 00:23:32 Found 2 connected components, falling back to 'spca' initialization with init_sdev = 1
    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'
    ## Also defined by 'BiocGenerics'
    ## 00:23:32 Using 'irlba' for PCA
    ## 00:23:32 PCA: 2 components explained 58.74% variance
    ## 00:23:32 Scaling init to sdev = 1
    ## 00:23:32 Commencing optimization for 200 epochs, with 492740 positive edges
    ## 00:23:41 Optimization finished

``` r
p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
  combine = FALSE, label.size = 2
)
```

``` r
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
obj <- FindClusters(obj, resolution = 2, cluster.name = "rpca_clusters")
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10434
    ## Number of edges: 487050
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7953
    ## Number of communities: 26
    ## Elapsed time: 1 seconds

``` r
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
```

    ## 00:23:45 UMAP embedding parameters a = 0.9922 b = 1.112

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:23:45 Read 10434 rows and found 30 numeric columns

    ## 00:23:45 Using Annoy for neighbor search, n_neighbors = 30

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:23:45 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 00:23:46 Writing NN index file to temp file C:\Users\matta\AppData\Local\Temp\Rtmp8cVbQ1\filea76c3a71647c
    ## 00:23:46 Searching Annoy index using 1 thread, search_k = 3000
    ## 00:23:48 Annoy recall = 100%
    ## 00:23:50 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 00:23:52 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 00:23:53 Commencing optimization for 200 epochs, with 476490 positive edges
    ## 00:24:02 Optimization finished

``` r
p2 <- DimPlot(
  obj,
  reduction = "umap.rpca",
  group.by = c("Method", "predicted.celltype.l2", "rpca_clusters"),
  combine = FALSE, label.size = 2
)
```

``` r
wrap_plots(c(p1, p2), ncol = 2, byrow = F)
```

![](seurat5_integration_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Visualize and cluster the datasets for harmony integration.

``` r
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10434
    ## Number of edges: 460218
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7903
    ## Number of communities: 24
    ## Elapsed time: 1 seconds

``` r
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
```

    ## 00:24:11 UMAP embedding parameters a = 0.9922 b = 1.112

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:24:11 Read 10434 rows and found 30 numeric columns

    ## 00:24:11 Using Annoy for neighbor search, n_neighbors = 30

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:24:11 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 00:24:11 Writing NN index file to temp file C:\Users\matta\AppData\Local\Temp\Rtmp8cVbQ1\filea76c452936eb
    ## 00:24:12 Searching Annoy index using 1 thread, search_k = 3000
    ## 00:24:14 Annoy recall = 100%
    ## 00:24:15 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 00:24:18 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 00:24:19 Commencing optimization for 200 epochs, with 466698 positive edges
    ## 00:24:28 Optimization finished

``` r
p3 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("Method", "predicted.celltype.l2", "harmony_clusters"),
  combine = FALSE, label.size = 2
)
```

Compare the expression of biological markers based on different
clustering solutions, or visualize the method’s clustering solution on
different UMAP visualizations

``` r
p1 <- VlnPlot(
  obj,
  features = "rna_CD8A", group.by = "unintegrated_clusters"
) + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")

p2 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "cca_clusters"
) + NoLegend() + ggtitle("CD8A - CCA Clusters")

p3 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "rpca_clusters"
) + NoLegend() + ggtitle("CD8A - scVI Clusters")

p1 | p2 | p3
```

![](seurat5_integration_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
```

    ## 00:24:32 UMAP embedding parameters a = 0.9922 b = 1.112

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:24:32 Read 10434 rows and found 30 numeric columns

    ## 00:24:32 Using Annoy for neighbor search, n_neighbors = 30

    ## Found more than one class "dist" in cache; using the first, from namespace 'spam'

    ## Also defined by 'BiocGenerics'

    ## 00:24:32 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 00:24:33 Writing NN index file to temp file C:\Users\matta\AppData\Local\Temp\Rtmp8cVbQ1\filea76c787c2bad
    ## 00:24:33 Searching Annoy index using 1 thread, search_k = 3000
    ## 00:24:35 Annoy recall = 100%
    ## 00:24:36 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 00:24:39 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 00:24:40 Commencing optimization for 200 epochs, with 476490 positive edges
    ## 00:24:49 Optimization finished

``` r
p4 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("cca_clusters"))
p5 <- DimPlot(obj, reduction = "umap.rpca", group.by = c("cca_clusters"))
p6 <- DimPlot(obj, reduction = "umap.harmony", group.by = c("cca_clusters"))

p4 | p5 | p6
```

![](seurat5_integration_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Once integrative analysis is complete you can rejoin the layers - which
collapses the individual datasets together and recreates the original
`counts` and `data` layers. This needs to be done before performing any
differential expression analysis. You can always resplit the layers in
case you would like to perform integrative analysis again.

``` r
obj <- JoinLayers(obj)
obj
```

    ## An object of class Seurat 
    ## 33789 features across 10434 samples within 4 assays 
    ## Active assay: RNA (33694 features, 2000 variable features)
    ##  3 layers present: data, counts, scale.data
    ##  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
    ##  10 dimensional reductions calculated: integrated_dr, ref.umap, pca, umap.unintegrated, integrated.cca, integrated.rpca, harmony, umap.cca, umap.rpca, umap.harmony

Integration can also be performed using sctransform-normalized data by
first running SCTransform normalization, and then setting the
`normalization.method` argument in `IntegrateLayers`.
