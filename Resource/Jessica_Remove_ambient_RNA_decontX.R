library(Seurat)
library(SingleCellExperiment)
library(celda)
source("~/Documents/Source/scFunctions.R")

# Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6
# Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/celda/inst/doc/decontX.pdf


#### Run with just filtered cells ####
seuratObject <- readRDS("~/Documents/WT_5XFAD/Final_HP_WT_5XFAD_Renamed.rds")
seuratObject <- subset(seuratObject, subset=orig.ident!="TGW-3-HP")
cds <- SingleCellExperiment(assays=list(counts=seuratObject@assays$RNA@counts,
                                        norm=seuratObject@assays$RNA@data),
                            colData = seuratObject@meta.data)
cds <- decontX(x = cds)
saveRDS(cds, "~/Documents/WT_5XFAD/decontX_HP_WT_5XFAD.rds")

seuratObject@assays$decontX <- CreateAssayObject(counts = decontXcounts(cds))
DefaultAssay(seuratObject) <- "decontX"
FeaturePlot(seuratObject,"Tyrobp", reduction = "tsne")

umap <- reducedDim(cds, "decontX_UMAP")
plotDimReduceCluster(dim1 = umap[, 1], dim2 = umap[, 2], 
                     cluster = cds$decontX_clusters)
plotDimReduceCluster(dim1 = umap[, 1], dim2 = umap[, 2], 
                     cluster = cds$condition)
plotDimReduceCluster(dim1 = umap[, 1], dim2 = umap[, 2], 
                     cluster = cds$Cell_typev2)
plotDecontXContamination(cds)

# visualize on Seurat tsne
decontam <- decontXcounts(cds)
seuratObject$decontX_contamination <- cds$decontX_contamination

pdf("~/Documents/WT_5XFAD/decontX/decontX_contamination.pdf")
FeaturePlot(seuratObject,"decontX_contamination", reduction = "tsne") + 
  scale_colour_gradientn(
    colors = c("blue","green","yellow","orange","red"),
    #name = "Contamination",
    limits = c(0, 1)
  ) +
  theme_cleanDimPlotLegend()
dev.off()

seuratObject$decontX_contamination[seuratObject$decontX_contamination>0.5] <- 0.5
pdf("~/Documents/WT_5XFAD/decontX/decontX_contamination_max50.pdf")
FeaturePlot(seuratObject,"decontX_contamination", reduction = "tsne") + 
  scale_colour_gradientn(
    colors = c("blue","green","yellow","orange","red"),
    #name = "Contamination",
    limits = c(0, 0.5)
  ) +
  theme_cleanDimPlotLegend()
dev.off()

seuratObject$decontX_contamination[seuratObject$decontX_contamination>0.1] <- 0.1
pdf("~/Documents/WT_5XFAD/decontX/decontX_contamination_max10.pdf")
FeaturePlot(seuratObject,"decontX_contamination", reduction = "tsne") + 
  theme_cleanDimPlotLegend()
dev.off()

# Normalize corrected data
seuratObject@assays$RNA@decontX_data <- decontXcounts(cds)
seuratObject@assays$decontX <- CreateAssayObject(counts = decontXcounts(cds))

DefaultAssay(seuratObject) <- "decontX"

seuratObject <- NormalizeData(seuratObject)
seuratObject <- FindVariableFeatures(seuratObject)
length(intersect(seuratObject@assays$RNA@var.features, seuratObject@assays$decontX@var.features)) # 1972
setdiff(seuratObject@assays$RNA@var.features, seuratObject@assays$decontX@var.features)
seuratObject <- ScaleData(seuratObject)
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject, 
                                                                 assay = "decontX"))
nDims = 30
resUse = seq(0.5,4,by=0.5)
seuratObject <- FindNeighbors(seuratObject, dims = 1:nDims, reduction="pca",k.param = 25)
seuratObject <- FindClusters(object = seuratObject, resolution = resUse, verbose = T, reduction = "pca")
seuratObject <- RunUMAP(seuratObject, dims = 1:nDims, reduction = "pca")
seuratObject <- RunTSNE(object = seuratObject, reduction = "pca", dims = 1:nDims)

DimPlot(seuratObject, reduction="tsne", group.by="condition")
saveRDS(seuratObject, "~/Documents/WT_5XFAD/decontX/Rerun_clustering_HP_WT_5XFAD_decontX.rds")

# Rerun DEGs for astrocytes to see if microglia marker genes are still strong DEGs
Idents(seuratObject) <- "Cell_type_condition"
degs <- FindMarkers(seuratObject, ident.1="Astrocyte_TGW", ident.2 = "Astrocyte_WTW")
DefaultAssay(seuratObject) <- "RNA"
degs_RNA <- FindMarkers(seuratObject, ident.1="Astrocyte_TGW", ident.2 = "Astrocyte_WTW")


# using the background did not work well for me
#### Run with raw matrix UMI > 50 #### 
#374708 background droplets

# clear space
rm(cds)
# get cells to remove from raw matrix to get just background

# get unfiltered data from 10x (outs/raw_feature_bc_matrix)

# worked in hoffman2
setwd("/u/flashscratch/j/jading/raw_counts_scAD/")
HP <- readRDS("/u/flashscratch/j/jading/10X_Single_Cell/HP_HYP/Final_HP_WT_5XFAD_Renamed.rds")
filtered_cells <- colnames(HP)
sample_dirs = c("TGW-3-HP","TGW-5-HP","TGW-6-HP","WTW-1-HP","WTW-2-HP","WTW-3-HP")

for(sample in sample_dirs){
  
  count.data <- Read10X(data.dir = sample)
  count.data@Dimnames[[2]] <- gsub("-1","",count.data@Dimnames[[2]])
  count.data@Dimnames[[2]] = paste0(sample,"_",count.data@Dimnames[[2]])
  cat(sum(count.data@Dimnames[[2]] %in% filtered_cells),"\n")
  count.data <- count.data[,!(count.data@Dimnames[[2]] %in% filtered_cells)]
  sample.seurat <- CreateSeuratObject(counts = count.data, project = sample)
  rm(count.data)
  sample.seurat <- subset(sample.seurat, subset=nCount_RNA>50)
  
  # meta data
  descriptions = unlist(strsplit(sample,split = "-"))
  meta_data = c("condition","replicate","tissue")
  for(iter in 1:length(meta_data)){
    sample.seurat@meta.data[,meta_data[iter]] <- descriptions[iter]
  }
  
  # merge
  if(!exists("seuratObject")){
    seuratObject <- sample.seurat
    firstSampleName = sample
    firstSample = TRUE
  } 
  else{
    if(firstSample==TRUE){
      seuratObject <- merge(x = seuratObject, y = sample.seurat, project = "AD_Metabolic")
      firstSample = FALSE
    } 
    else{
      seuratObject <- merge(x = seuratObject, y = sample.seurat, project = "AD_Metabolic")
    }
  }
  
  #clean up environment
  rm(sample.seurat)
  
}
saveRDS(seuratObject, 
        "/u/flashscratch/j/jading/10X_Single_Cell/AD_HP/HP_WT_5XFAD_over50UMI.rds")


cds <- SingleCellExperiment(assays=list(counts=HP@assays$RNA@counts),
                            colData = HP@meta.data)
rm(HP)
cds_background <- SingleCellExperiment(assays=list(counts=seuratObject@assays$RNA@counts),
                                       colData = seuratObject@meta.data)
rm(seuratObject)
cds <- decontX(x = cds, background = cds_background)
saveRDS(cds, "/u/flashscratch/j/jading/10X_Single_Cell/AD_HP/HP_WT_5XFAD_decontX_over50UMIbkgd.rds")

cds_UMIgr8r50 <- readRDS("~/Documents/WT_5XFAD/decontX/HP_WT_5XFAD_decontX_over50UMIbkgd.rds")

#### Run with background cells from filtered matrix ####
# did the same workflow but for "background cells" of the filtered matrix
filtered_removed_cells <- readRDS("~/Documents/WT_5XFAD/decontX/HP_WT_5XFAD_filtered_removed_cells.rds")
dim(filtered_removed_cells)
# 13208 background droplets

cds <- readRDS("~/Documents/WT_5XFAD/decontX_HP_WT_5XFAD.rds")
cds_UMIgr8r50 <- readRDS("~/Documents/WT_5XFAD/decontX/HP_WT_5XFAD_decontX_over50UMIbkgd.rds")
cds_filtered_bkgd <- readRDS("~/Documents/WT_5XFAD/decontX/HP_WT_5XFAD_decontX_filtered_bkgd.rds")

seuratObject <- readRDS("~/Documents/WT_5XFAD/Final_HP_WT_5XFAD_Renamed.rds")

seuratObject$decontX_contamination <- cds$decontX_contamination
seuratObject$decontX_contamination_umi50 <- cds_UMIgr8r50$decontX_contamination
seuratObject$decontX_contamination_filtered <- cds_filtered_bkgd$decontX_contamination

FeaturePlot(seuratObject, c("decontX_contamination","decontX_contamination_umi50",
                            "decontX_contamination_filtered"), reduction="tsne")

seuratObject$decontX_contamination[seuratObject$decontX_contamination>0.5] <- 0.5
seuratObject$decontX_contamination_umi50[seuratObject$decontX_contamination_umi50>0.1] <- 0.1
seuratObject$decontX_contamination_filtered[seuratObject$decontX_contamination_filtered>0.1] <- 0.1

plots <- list()
plots[[1]] <- FeaturePlot(seuratObject,"decontX_contamination", reduction = "tsne") + 
  scale_colour_gradientn(
    colors = c("blue","green","yellow","orange","red"),
    #name = "Contamination",
    limits = c(0, 0.5)
  ) +
  theme_cleanDimPlotLegend() +
  ggtitle("No background")
plots[[2]] <- FeaturePlot(seuratObject,"decontX_contamination_umi50", reduction = "tsne") + 
  scale_colour_gradientn(
    colors = c("blue","green","yellow","orange","red"),
    #name = "Contamination",
    limits = c(0, max(seuratObject$decontX_contamination_umi50))
  ) +
  theme_cleanDimPlotLegend() +
  ggtitle("With raw droplets >50 UMI\nas background")
plots[[3]] <- FeaturePlot(seuratObject,"decontX_contamination_filtered", reduction = "tsne") + 
  scale_colour_gradientn(
    colors = c("blue","green","yellow","orange","red"),
    #name = "Contamination",
    limits = c(0, max(seuratObject$decontX_contamination_filtered))
  ) +
  theme_cleanDimPlotLegend()+
  ggtitle("With filtered (not included) droplets\nas background")

pdf("~/Documents/WT_5XFAD/decontX/contamination_results_maxcutoff.pdf", width = 21)
ggarrange(plotlist = plots, nrow = 1, ncol = 3)
dev.off()

cds_redo <- SingleCellExperiment(assays=list(counts=seuratObject@assays$RNA@counts,
                                        norm=seuratObject@assays$RNA@data),
                                 colData = seuratObject@meta.data)
cds_background <- SingleCellExperiment(assays=list(counts=filtered_removed_cells@assays$RNA@counts),
                                       colData = filtered_removed_cells@meta.data)
# doesn't work on my local computer?
cds_redo <- decontX(x = cds_redo, background = cds_background)


seuratObject <- readRDS("~/Documents/WT_5XFAD/decontX/Rerun_clustering_HP_WT_5XFAD_decontX.rds")

pdf("~/Documents/WT_5XFAD/decontX/Rerun_clustering.pdf", width = 14)
ggarrange(DimPlot(seuratObject, reduction="tsne", group.by="condition") + theme_cleanDimPlotLegend(),
          DimPlot(seuratObject, reduction="tsne", group.by="Cell_typev2", label = TRUE, label.size = 5) + theme_cleanDimPlot(),
          ncol = 2)
dev.off()

DefaultAssay(seuratObject) <- "decontX"
FeaturePlot(seuratObject, "Tyrobp", reduction = "tsne")
DefaultAssay(seuratObject) <- "RNA"
FeaturePlot(seuratObject, "Tyrobp", reduction = "tsne")

plots <- list()
DefaultAssay(seuratObject) <- "decontX"
plots[[1]] <- FeaturePlot(seuratObject, "Tyrobp", reduction = "tsne", pt.size = 1) + 
  theme_cleanDimPlot() + ggtitle("decontX assay")
DefaultAssay(seuratObject) <- "RNA"
plots[[2]] <- FeaturePlot(seuratObject, "Tyrobp", reduction = "tsne", pt.size = 1) + 
  theme_cleanDimPlot() + ggtitle("RNA/unadjusted assay")

pdf("~/Documents/WT_5XFAD/decontX/Tyrobp_res.pdf", width = 14, height = 8)
plt <- ggarrange(plotlist = plots,
                ncol = 2)
annotate_figure(p = plt, top = text_grob("Tyrobp expression", size = 30))
dev.off()


Idents(seuratObject) <- "Cell_type_condition"
DefaultAssay(seuratObject) <- "decontX"
degs <- FindMarkers(seuratObject, ident.1="Oligodendrocyte_TGW", ident.2 = "Oligodendrocyte_WTW", logfc.threshold = 0)
write.table(degs, "Oligo_AD_DEGs_decontX_assay.txt", row.names = FALSE, sep = "\t", quote = FALSE)
DefaultAssay(seuratObject) <- "RNA"
degs_RNA <- FindMarkers(seuratObject, ident.1="Oligodendrocyte_TGW", ident.2 = "Oligodendrocyte_WTW")

degs$Cell_type <- "Oligodendrocyte"
degs$GENE <- rownames(degs)
library(ggrepel)
All_HP_DEGs <- read.delim("~/Documents/AD_Metabolic/DEGs/HP_DEGs_2019_2020.txt")
v1 <- volcano_plot(DEG_df = degs, 
                   cell_type = "Oligodendrocyte", 
                   lfc = .5, 
                   #padjusted = 1e-25) + 
  padjusted = 1e-100) + 
  ggtitle("2019 decontX") +
  theme(plot.title = element_text(hjust = 0.5))
v2 <- volcano_plot(DEG_df = All_HP_DEGs[All_HP_DEGs$Comparison=="TGW2020.v.WTW2019",], 
                   cell_type = "Oligodendrocyte", 
                   lfc = .5, 
                   #padjusted = 1e-25) + 
  padjusted = 1e-100) + 
  ggtitle("2020") +
  theme(plot.title = element_text(hjust = 0.5))
pdf("~/Documents/WT_5XFAD/decontX/Oligo_DEGs.pdf", width = 14)
fig <- ggarrange(v1,v2, nrow = 1)
print(annotate_figure(fig,
                      top = text_grob("Oligodendrocyte, AD effect", 
                                      face = "bold", size = 18)))
dev.off()
