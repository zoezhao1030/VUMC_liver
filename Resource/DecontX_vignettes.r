---
title: "Decontamination of ambient RNA in single-cell genomic data with DecontX" 
author:
- name: Shiyi (Iris) Yang
  affiliation: &id Boston University School of Medicine
- name: Zhe Wang
  affiliation: *id
- name: Yuan Yin
  affiliation: *id
- name: Joshua Campbell
  affiliation: *id
  email: camp@bu.edu
date: "`r Sys.Date()`"
output: 
  BiocStyle::html_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{Estimate and remove cross-contamination from ambient RNA in single-cell data with DecontX}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, dev = "png")
```

# Introduction
Droplet-based microfluidic devices have become widely used to perform single-cell RNA sequencing (scRNA-seq). However, ambient RNA present in the cell suspension can be aberrantly counted along with a cell’s native mRNA and result in cross-contamination of transcripts between different cell populations. DecontX is a Bayesian method to estimate and remove contamination in individual cells. DecontX assumes the observed expression of a cell is a mixture of counts from two multinomial distributions: (1) a distribution of native transcript counts from the cell’s actual population and (2) a distribution of contaminating transcript counts from all other cell populations captured in the assay. Overall, computational decontamination of single cell counts can aid in downstream clustering and visualization. 

The package can be loaded using the `library` command.

```{r load, eval=TRUE, message=FALSE}
library(celda)
```

# Importing data

DecontX can take either a [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object or a counts matrix as input. `decontX` will attempt to convert any input matrix to class `dgCMatrix` from package [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html) before starting the analysis.

To import datasets directly into an SCE object, the [singleCellTK](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html) package has several importing functions for different preprocessing tools including CellRanger, STARsolo, BUStools, Optimus, DropEST, SEQC, and Alevin/Salmon. For example, the following code can be used as a template to read in the filtered and raw matrices for multiple samples processed with CellRanger:

```{r sce_import, eval = FALSE}
library(singleCellTK)
sce <- importCellRanger(sampleDirs = c("path/to/sample1/", "path/to/sample2/"))
```

Within each sample directory, there should be subfolders called `"outs/filtered_feature_bc_matrix/"` or `"outs/raw_feature_bc_matrix/"` with files called `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`. If these files are in different subdirectories, the `importCellRangerV3Sample` function can be used to import data from a different directory instead. 

Optionally, the "raw" or "droplet" matrix can also be easily imported by setting the `dataType` argument to "raw":

```{r sce_import_raw, eval = FALSE}
sce.raw <- importCellRanger(sampleDirs = c("path/to/sample1/", "path/to/sample2/"), dataType = "raw")
```

The raw matrix can be passed to the `background` parameter in `decontX` as described below. If using Seurat, go to the [Working with Seurat](#seurat) section for details on how to convert between SCE and Seurat objects.

# Load PBMC4k data from 10X

We will utilize the 10X PBMC 4K dataset as an example in this vignette. This data can be easily retrieved from the package [TENxPBMCData](http://bioconductor.org/packages/release/data/experiment/html/TENxPBMCData.html). Make sure the the column names are set before running decontX.

```{r load_10X, eval=TRUE, message=FALSE}
# Load PBMC data
library(TENxPBMCData)
sce <- TENxPBMCData("pbmc4k")
colnames(sce) <- paste(sce$Sample, sce$Barcode, sep = "_")
rownames(sce) <- rowData(sce)$Symbol_TENx
counts(sce) <- as(counts(sce), "dgCMatrix")
```


# Running decontX

A SingleCellExperiment (SCE) object or a sparse matrix containing the counts for filtered cells can be passed to decontX via the `x` parameter. The matrix to use in an SCE object can be specified with the `assayName` parameter, which is set to `"counts"` by default. There are two major ways to run decontX: with and without the raw/droplet matrix containing empty droplets. Here is an example of running decontX without supplying the background:

```{r decontX, eval=TRUE, message=FALSE}
sce <- decontX(sce)
```

In this scenario, `decontX` will estimate the contamination distribution for each cell cluster based on the profiles of the other cell clusters in the filtered dataset. The estimated contamination results can be found in the `colData(sce)$decontX_contamination` and the decontaminated counts can be accessed with `decontXcounts(sce)`. `decontX` will perform heuristic clustering to quickly define major cell clusters. However if you have your own cell cluster labels, they can be specified with the `z` parameter. These results will be used throughout the rest of the vignette. 

The raw/droplet matrix can be used to empirically estimate the distribution of ambient RNA, which is especially useful when cells that contributed to the ambient RNA are not accurately represented in the filtered count matrix containing the cells. For example, cells that were removed via flow cytometry or that were more sensitive to lysis during dissociation may have contributed to the ambient RNA but were not measured in the filtered/cell matrix. The raw/droplet matrix can be input as an SCE object or a sparse matrix using the `background` parameter:

```{r decontX_background, eval=FALSE, message=FALSE}
sce <- decontX(sce, background = sce.raw)
```

Only empty droplets in the background matrix should be used to estimate the ambient RNA. If any cell ids (i.e. `colnames`) in the raw/droplet matrix supplied to the `background` parameter are also found in the filtered counts matrix (`x`), decontX will automatically remove them from the raw matrix. However, if the cell ids are not available for the input matrices, decontX will treat the entire `background` input as empty droplets. All of the outputs are the same as when running decontX without setting the `background` parameter.

> Note: If the input object is just a matrix and not an SCE object, make sure to save the output into a variable with a different name (e.g. `result <- decontX(mat)`). The result object will be a list with contamination in `result$contamination` and the decontaminated counts in `result$decontXcounts`. 

# Plotting DecontX results

## Cluster labels on UMAP
DecontX creates a UMAP which we can use to plot the cluster labels automatically identified in the analysis. Note that the clustering approach used here is designed to find "broad" cell types rather than individual cell subpopulations within a cell type. 

```{r UMAP_Clusters}
umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
```

## Contamination on UMAP
The percentage of contamination in each cell can be plotting on the UMAP to visualize what what clusters may have higher levels of ambient RNA.

```{r plot_decon}
plotDecontXContamination(sce)
```

## Expression of markers on UMAP
Known marker genes can also be plotted on the UMAP to identify the cell types for each cluster. We will use CD3D and CD3E for T-cells, LYZ, S100A8, and S100A9 for monocytes, CD79A, CD79B, and MS4A1 for B-cells, GNLY for NK-cells, and PPBP for megakaryocytes.

```{r plot_feature, message=FALSE}
library(scater)
sce <- logNormCounts(sce)
plotDimReduceFeature(as.matrix(logcounts(sce)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = c("CD3D", "CD3E", "GNLY",
        "LYZ", "S100A8", "S100A9",
        "CD79A", "CD79B", "MS4A1"),
    exactMatch = TRUE)
```

## Barplot of markers detected in cell clusters
The percetage of cells within a cluster that have detectable expression of marker genes can be displayed in a barplot. Markers for cell types need to be supplied in a named list. First, the detection of marker genes in the original `counts` assay is shown:

```{r barplotCounts}
markers <- list(Tcell_Markers = c("CD3E", "CD3D"),
    Bcell_Markers = c("CD79A", "CD79B", "MS4A1"),
    Monocyte_Markers = c("S100A8", "S100A9", "LYZ"),
    NKcell_Markers = "GNLY")
cellTypeMappings <- list(Tcells = 2, Bcells = 5, Monocytes = 1, NKcells = 6)
plotDecontXMarkerPercentage(sce,
    markers = markers,
    groupClusters = cellTypeMappings,
    assayName = "counts")
```

We can then look to see how much decontX removed aberrant expression of marker genes in each cell type by changing the `assayName` to `decontXcounts`:

```{r barplotDecontCounts}
plotDecontXMarkerPercentage(sce,
    markers = markers,
    groupClusters = cellTypeMappings,
    assayName = "decontXcounts")
```

Percentages of marker genes detected in other cell types were reduced or completely removed. For example, the percentage of cells that expressed Monocyte marker genes was greatly reduced in T-cells, B-cells, and NK-cells.
The original counts and decontamined counts can be plotted side-by-side by listing multiple assays in the `assayName` parameter. This option is only available if the data is stored in `SingleCellExperiment` object.

```{r barplotBoth}
plotDecontXMarkerPercentage(sce,
    markers = markers,
    groupClusters = cellTypeMappings,
    assayName = c("counts", "decontXcounts"))
```

Some helpful hints when using `plotDecontXMarkerPercentage`:

1. Cell clusters can be renamed and re-grouped using the `groupCluster` parameter, which also needs to be a named list. If `groupCluster` is used, cell clusters not included in the list will be excluded in the barplot. For example, if we wanted to group T-cells and NK-cells together, we could set `cellTypeMappings <- list(NK_Tcells = c(2,6), Bcells = 5, Monocytes = 1)`
2. The level a gene that needs to be expressed to be considered detected in a cell can be adjusted using the `threshold` parameter.
3. If you are not using a `SingleCellExperiment`, then you will need to supply the original counts matrix or the decontaminated counts matrix as the first argument to generate the barplots. 

## Violin plot to compare the distributions of original and decontaminated counts
Another useful way to assess the amount of decontamination is to view the expression of marker genes before and after `decontX` across cell types. Here we view the monocyte markers in each cell type. The violin plot shows that the markers have been removed from T-cells, B-cells, and NK-cells, but are largely unaffected in monocytes.

```{r plotDecontXMarkerExpression}
plotDecontXMarkerExpression(sce,
    markers = markers[["Monocyte_Markers"]],
    groupClusters = cellTypeMappings,
    ncol = 3)
```

Some helpful hints when using `plotDecontXMarkerExpression`:

1. `groupClusters` works the same way as in `plotDecontXMarkerPercentage`.
2. This function will plot each pair of markers and clusters (or cell type specified by `groupClusters`). Therefore, you may want to keep the number of markers small in each plot and call the function multiple times for different sets of marker genes. 
3. You can also plot the individual points by setting `plotDots = TRUE` and/or log transform the points on the fly by setting `log1p = TRUE`. 
4. This function can plot any assay in a `SingleCellExperiment`. Therefore you could also examine normalized expression of the original and decontaminated counts. For example:


```{r plot_norm_counts, eval = TRUE}
library(scater)
sce <- logNormCounts(sce,
    exprs_values = "decontXcounts",
    name = "decontXlogcounts")

plotDecontXMarkerExpression(sce,
    markers = markers[["Monocyte_Markers"]],
    groupClusters = cellTypeMappings,
    ncol = 3,
    assayName = c("logcounts", "decontXlogcounts"))
```


# Other important notes

## Choosing appropriate cell clusters
The ability of DecontX to accurately identify contamination is dependent on the cell cluster labels. DecontX assumes that contamination for a cell cluster comes from combination of counts from all other clusters. The default clustering approach used by DecontX tends to select fewer clusters that represent broader cell types. For example, all T-cells tend to be clustered together rather than splitting naive and cytotoxic T-cells into separate clusters. Custom cell type labels can be suppled via the `z` parameter if some cells are not being clustered appropriately by the default method.

## Adjusting the priors to influence contamination estimates
There are ways to force `decontX` to estimate more or less contamination across a dataset by manipulating the priors. The `delta` parameter is a numeric vector of length two. It is the concentration parameter for the Dirichlet distribution which serves as the prior for the proportions of native and contamination counts in each cell. The first element is the prior for the proportion of native counts while the second element is the prior for the proportion of contamination counts. These essentially act as pseudocounts for the native and contamination in each cell. If `estimateDelta = TRUE`, `delta` is only used to produce a random sample of proportions for an initial value of contamination in each cell. Then `delta` is updated in each iteration. If `estimateDelta = FALSE`, then `delta` is fixed with these values for the entire inference procedure. Fixing `delta` and setting a high number in the second element will force `decontX` to be more aggressive and estimate higher levels of contamination in each cell at the expense of potentially removing native expression. For example, in the previous PBMC example, we can see what the estimated `delta` was by looking in the estimates:

```{r findDelta}
metadata(sce)$decontX$estimates$all_cells$delta
```

Setting a higher value in the second element of delta and `estimateDelta = FALSE` will force `decontX` to estimate higher levels of contamination per cell:

```{r newDecontX, eval=TRUE, message=FALSE}
sce.delta <- decontX(sce, delta = c(9, 20), estimateDelta = FALSE)

plot(sce$decontX_contamination, sce.delta$decontX_contamination,
     xlab = "DecontX estimated priors",
     ylab = "Setting priors to estimate higher contamination")
abline(0, 1, col = "red", lwd = 2)
```

## Working with Seurat {#seurat}

If you are using the [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) package for downstream analysis, the following code can be used to read in a matrix and convert between Seurat and SCE objects:

```{r seurat_create, eval = FALSE}
# Read counts from CellRanger output
library(Seurat)
counts <- Read10X("sample/outs/filtered_feature_bc_matrix/")

# Create a SingleCellExperiment object and run decontX
sce <- SingleCellExperiment(list(counts = counts))
sce <- decontX(sce)

# Create a Seurat object from a SCE with decontX results
seuratObject <- CreateSeuratObject(round(decontXcounts(sce)))
```

Optionally, the "raw" matrix can be also be imported and used as the background:

```{r seurat_raw, eval = FALSE}
counts.raw <- Read10X("sample/outs/raw_feature_bc_matrix/")
sce.raw <- SingleCellExperiment(list(counts = counts.raw))
sce <- decontX(sce, background = sce.raw)
```

Note that the decontaminated matrix of decontX consists of floating point numbers and must be rounded to integers before adding it to a Seurat object. If you already have a Seurat object containing the counts matrix and would like to run decontX, you can retrieve the count matrix, create a SCE object, and run decontX, and then add it back to the Seurat object:

```{r seurat_create2, eval = FALSE}
counts <- GetAssayData(object = seuratObject, slot = "counts")
sce <- SingleCellExperiment(list(counts = counts))
sce <- decontX(sce)
seuratObj[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))
```


# Session Information

```{r}
sessionInfo()
```
