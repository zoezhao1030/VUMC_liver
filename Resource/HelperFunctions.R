########################################################
#Load required libraries for Seurat Analysis
loadLibraries <- function(){
  library(Seurat) # v 3.0.0.9000
  library(Matrix)
  library(dplyr)
  library(xlsx)
  library(varhandle)
  library(scMCA)
  library(DropSeq.util)
  library(DESeq2)
  library(MAST)
  library(metap)
  library(DoubletFinder)
  library(plyr)
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  library(ggpubr)
  library(sctransform)
  library(foreach)
  library(doParallel)
  library(doSNOW)
  library(reshape2)
  library(ggrepel)
  library(caret)
  source("./FindConservedMarkers.R")
}
########################################################
#Set threshold for filtering cells
setThresholds <- function(MinCells=3,MinGenes=400,MinUMIs=700,MinPercentMT=-Inf,MaxGenes=3000,
                          MaxPercentMT=0.1,MaxUMIs=5000,MaxPercentRibo=0.1,MaxPercentPred=0.05){
  assign("mincells", MinCells, envir = .GlobalEnv)               #min number of cells that a gene has to be detected in to be kept
  assign("mingenes", MinGenes, envir = .GlobalEnv)               #min number of genes to keep the cell
  assign("minUMIs", MinUMIs, envir = .GlobalEnv)                 #min number of UMIs to keep the cell
  assign("minPercentMT", MinPercentMT, envir = .GlobalEnv)       #min mitochondrial percentage to keep a cell
  assign("maxgenes", MaxGenes, envir = .GlobalEnv)               #max number of genes to keep a cell
  assign("maxPercentMT", MaxPercentMT, envir = .GlobalEnv)       #max mitochondrial percentage to keep a cell
  assign("maxUMIs", MaxUMIs, envir = .GlobalEnv)                 #max UMIs to keep a cell
  assign("maxPercentRibo", MaxPercentRibo, envir = .GlobalEnv)   #max ribosomal percentage to keep a cell
  assign("maxPercentPred", MaxPercentPred, envir = .GlobalEnv)   #max predicted genes percentage to keep a cell
}
########################################################
#Read in the dropEST data
readDropESTdata = function(tissue=NULL,directoryPath="./",timePointINT=NULL,
                           conditionINT=5,sampleNameINT=6,projectName=NULL){
  #Get dropEST matrices directories
  directories = list.dirs(path = directoryPath, full.names = TRUE, recursive = TRUE)
  #Get Tissue Samples
  directories = directories[grep(tissue,directories)]
  #iterate through the samples
  for(sample in directories){
    #get meta data
    splitSample = unlist(strsplit(sample,"/"))
    if(!is.null(timePointINT)){
      timePoint = splitSample[timePointINT]
    }
    condition = splitSample[conditionINT]
    sampleName = splitSample[sampleNameINT]
    print(sampleName)
    animal = unlist(strsplit(sampleName,tissue))[1]
    #load the data
    count.data <- Read10X(data.dir = sample)
    count.data@Dimnames[[2]] = paste0(sampleName,"_",count.data@Dimnames[[2]])
    sample.seurat <- CreateSeuratObject(counts = count.data, project = sampleName)
    rm(count.data)
    #define meta data
    if(!is.null(timePointINT)){
      sample.seurat@meta.data$timpoint = timePoint
    }
    sample.seurat@meta.data$condition = condition
    sample.seurat@meta.data$sampleName = sampleName
    sample.seurat@meta.data$animal = animal
    sample.seurat@meta.data$tissue = tissue
    if(!is.null(timePointINT)){
      sample.seurat@meta.data$timpoint.condition = paste0(timePoint,condition)
    }
    #merge seurat objects
    if(!exists("seuratObject")){
      seuratObject <- sample.seurat
      firstSampleName = sampleName
      firstSample = TRUE
    } else{
      if(firstSample==TRUE){
        seuratObject <- merge(x = seuratObject, y = sample.seurat, project = projectName)
        firstSample = FALSE
      } else{
        seuratObject <- merge(x = seuratObject, y = sample.seurat, project = projectName)
      }
    }
    #clean up environment
    rm(sample.seurat)
  }
  return(seuratObject)
}
########################################################
#Get mito, ribo, and predicted % and add to seurat object meta data
getCellQuality <- function(seuratObject=dropEST.combined,mtPrefix="^mt-",riboPrefix1="^Rps",riboPrefix2="^Rpl",
                           predictedFeaturePrefix=c("^Gm1","^Gm2","^Gm3","^Gm4","^Gm5","^Gm6","^Gm7","^Gm8","^Gm9")){
  # get mito %
  mito.features <- grep(pattern = mtPrefix, x = rownames(x = seuratObject), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts"))
  # get ribosome %
  ribo.features <- grep(pattern = riboPrefix1, x = rownames(x = seuratObject), value = TRUE)
  ribo.features <- c(ribo.features, grep(pattern = riboPrefix2, x = rownames(x = seuratObject), value = TRUE))
  percent.ribo <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts")[ribo.features,])/Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts"))
  
  # get predicted genes %
  pred.features <- grep(pattern = predictedFeaturePrefix[1], x = rownames(x = seuratObject), value = TRUE)
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[2], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[3], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[4], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[5], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[6], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[7], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[8], x = rownames(x = seuratObject), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = predictedFeaturePrefix[9], x = rownames(x = seuratObject), value = TRUE))
  percent.pred <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts")[pred.features,])/Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts"))
  
  # Add % mito, ribo and pred to the meta data
  seuratObject[["percent.mito"]] <- percent.mito
  seuratObject[["percent.ribo"]] <- percent.ribo
  seuratObject[["percent.pred"]] <- percent.pred
  
  return(seuratObject)
}

########################################################
cellQualityPlot <- function(seuratObject=dropEST.combined,fileName=NULL,H=9,W=40,
                            featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","percent.pred"),
                            identPlot="orig.ident",pointSize=0){
  #Make directory
  filePath = unlist(strsplit(x=fileName,split="/"))
  filePath = paste0(filePath[1:(length(filePath)-1)],collapse = "/")
  ifelse(!dir.exists(file.path(filePath)), dir.create(file.path(filePath),recursive = T), FALSE)
  #Generate QCPlots
  seuratObject = SetIdent(seuratObject, value = identPlot)
  pdf(file=fileName,height=H,width=W)
  print(VlnPlot(object = seuratObject, features = featuresPlot ,pt.size = pointSize))
  dev.off()
}

########################################################
normalizeAndScale <- function(seuratObject=dropEST.combined.filtered,varsRegress = c("nCount_RNA","percent.mito","percent.ribo")){
  seuratObject <- NormalizeData(object = seuratObject, normalization.method = "LogNormalize",scale.factor = 10000)
  seuratObject <- FindVariableFeatures(object = seuratObject)
  #Regress out number of UMI and percent mito and percent ribo
  seuratObject <- ScaleData(object = seuratObject, features = rownames(x = seuratObject),
                            vars.to.regress = varsRegress)
  
  #Regress out number of UMI and percent mito and percent ribo
  # seuratObject <- SCTransform(object = seuratObject, vars.to.regress = varsRegress, verbose = TRUE)
  
  return(seuratObject)
}

########################################################
runPCAandICA <- function(seuratObject=dropEST.combined.filtered,do.plot=T,tissue=NULL,identPlot="orig.ident"){
  #Perform PCA on the scaled data (uses the highly var genes)
  seuratObject <- RunPCA(object = seuratObject, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)
  #Perform ICA on the scaled data (uses the highly var genes)
  seuratObject <- RunICA(object = seuratObject, verbose = T, nics = 75, ndims.print = 1:5, nfeatures.print = 10)
  #Plot PCA Plots
  ifelse(!dir.exists(file.path("./Plots/QCPlots/PCA")), dir.create(file.path("./Plots/QCPlots/PCA"),recursive = T), FALSE)
  seuratObject = SetIdent(seuratObject, value = identPlot)
  pdf(file=paste0("./Plots/QCPlots/PCA/",tissue,".pdf"))
  print(DimPlot(object = seuratObject, dims = c(1,2), reduction = "pca"))
  dev.off()
  #Plot Heatmaps of PCs
  ifelse(!dir.exists(file.path("./Plots/QCPlots/PCHeatmap")), dir.create(file.path("./Plots/QCPlots/PCHeatmap"),recursive = T), FALSE)
  pdf(file=paste0("./Plots/QCPlots/PCHeatmap/",tissue,".pdf"),height=20,width=20)
  print(DimHeatmap(object = seuratObject, dims = 1:12, balanced = TRUE, cells = 100, reduction = "pca"))
  dev.off()
  #Plot ICA Plots
  ifelse(!dir.exists(file.path("./Plots/QCPlots/ICA")), dir.create(file.path("./Plots/QCPlots/ICA"),recursive = T), FALSE)
  seuratObject = SetIdent(seuratObject, value = identPlot)
  pdf(file=paste0("./Plots/QCPlots/ICA/",tissue,".pdf"))
  print(DimPlot(object = seuratObject, dims = c(1,2), reduction = "ica"))
  dev.off()
  #Plot Heatmaps of ICs
  ifelse(!dir.exists(file.path("./Plots/QCPlots/ICHeatmap")), dir.create(file.path("./Plots/QCPlots/ICHeatmap"),recursive = T), FALSE)
  pdf(file=paste0("./Plots/QCPlots/ICHeatmap/",tissue,".pdf"),height=20,width=20)
  print(DimHeatmap(object = seuratObject, dims = 1:12, balanced = TRUE, cells = 100, reduction = "ica"))
  dev.off()
  return(seuratObject)
}

########################################################
runSampleCCA <- function(seuratObject=dropEST.combined.filtered,combineLevel="timpoint.condition",numFeatures=2000,numDims=30,
                         tissue=NULL,kParam=25,resUse=seq(0.5,4,by=0.5),kWeight = 100, kAnchor = 5, kFilter = 200){
  seuratObject = SetIdent(seuratObject,value=combineLevel)
  seurat.list = list()
  for(i in 1:length(levels(seuratObject@active.ident))){
    seurat.list[[i]] = subset(seuratObject,idents=levels(seuratObject@active.ident)[i])
  }
  for (i in 1:length(x = seurat.list)) {
    # seurat.list[[i]] <- SCTransform(object = seurat.list[[i]] , verbose = TRUE)
    seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]], verbose = FALSE)
    seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", nfeatures = numFeatures, verbose = FALSE)
  }
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:numDims, k.anchor = kAnchor, k.filter = kFilter)
  seuratIntegrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:numDims, k.weight = kWeight)
  library(ggplot2)
  library(cowplot)
  DefaultAssay(object = seuratIntegrated) <- "integrated"
  # Run the standard workflow for visualization and clustering
  seuratIntegrated <- ScaleData(object = seuratIntegrated, verbose = FALSE)
  seuratIntegrated <- RunPCA(object = seuratIntegrated, npcs = numDims, verbose = FALSE)
  seuratIntegrated <- RunUMAP(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- RunTSNE(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- FindNeighbors(object = seuratIntegrated, reduction = "pca", dims = 1:numDims, k.param = kParam)
  seuratIntegrated <- FindClusters(object = seuratIntegrated, resolution = resUse, verbose = T, reduction = "pca")
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEcondition"))), dir.create(file.path(paste0("./Plots/ccaTSNEcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPcondition"))), dir.create(file.path(paste0("./Plots/ccaUMAPcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEresolution/",tissue))), dir.create(file.path(paste0("./Plots/ccaTSNEresolution/",tissue)),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPresolution/",tissue))), dir.create(file.path(paste0("./Plots/ccaUMAPresolution/",tissue)),recursive = T), FALSE)
  pdf(file=paste0("./Plots/ccaUMAPcondition/",tissue,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = combineLevel))
  dev.off()
  pdf(file=paste0("./Plots/ccaTSNEcondition/",tissue,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = combineLevel))
  dev.off()
  resolutions <- as.character(resUse)
  for(reso in resolutions){
    pdf(file=paste0("./Plots/ccaTSNEresolution/",tissue,"/PCA_tSNEres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "tsne"))
    dev.off()
    pdf(file=paste0("./Plots/ccaUMAPresolution/",tissue,"/PCA_UMAPres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "umap"))
    dev.off()
  }
  return(seuratIntegrated)
}

########################################################
runSampleDropVizCCA <- function(seuratObject=dropEST.combined.filtered,combineLevel="timpoint.condition",numFeatures=2000,numDims=30,
                                tissue=NULL,kParam=25,resUse=seq(0.5,4,by=0.5),dropVIZdata=NULL,dropVIZclusters=NULL,
                                dropVIZsubclusters=NULL,oldNames=NULL,newNames=NULL,subsetSize=5000,fileName=NULL,subsetClusters=FALSE){
  
  ########## Set up DropVIZ
  #Read in dropVIZ data to predict the cell types
  dropVizDGE <- loadSparseDge(dropVIZdata) 
  #Read in the DropViz Cluster Assignments
  clusterAssign = readRDS(file = dropVIZclusters)
  subclusterAssign = readRDS(file = dropVIZsubclusters)
  #Subset to only cells which have a subcluster assignment
  clusterAssign = clusterAssign[which(names(clusterAssign) %in% names(subclusterAssign))]
  clusterAssign = clusterAssign[names(subclusterAssign)]
  dropVizMetaData = data.frame(cluster = clusterAssign, subcluster = subclusterAssign)
  #change from cluster numbers to cell type names
  dropVizMetaData$cluster = mapvalues(dropVizMetaData$cluster, from = oldNames, to = newNames)
  # Set up DropViz object
  dropVizSeurat = CreateSeuratObject(counts = dropVizDGE, meta.data = dropVizMetaData)
  dropVizSeurat = SetIdent(dropVizSeurat,value = "cluster")
  # Keep only a subset of the cell populations?
  if(!isFALSE(subsetClusters)){
    dropVizSeurat = subset(dropVizSeurat,idents=subsetClusters)
  }
  # Subset to make it run (runs out of memory if we don't subset)
  dropVizSeurat = subset(dropVizSeurat,max.cells.per.ident = subsetSize)
  # Clean up environment
  rm(dropVizDGE,dropVizMetaData,clusterAssign,subclusterAssign)
  ########## Predict Cell Type
  DefaultAssay(object = seuratObject) <- "RNA"
  # #Set up list of datasets
  seuratObject = SetIdent(seuratObject,value=combineLevel)
  seurat.list = list()
  for(i in 1:length(levels(seuratObject@active.ident))){
    seurat.list[[i]] = subset(seuratObject,idents=levels(seuratObject@active.ident)[i])
  }
  #Add dropViz to the list
  seurat.list[[(length(seurat.list)+1)]] = dropVizSeurat
  rm(dropVizSeurat)
  for (i in 1:length(x = seurat.list)) {
    seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]], verbose = FALSE)
    seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", nfeatures = numFeatures, verbose = FALSE)
    # seurat.list[[i]] <- SCTransform(object = seurat.list[[i]], verbose = TRUE)
  }
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:numDims)
  seuratIntegrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:numDims)
  rm(seurat.list,seurat.anchors)
  library(ggplot2)
  library(cowplot)
  DefaultAssay(object = seuratIntegrated) <- "integrated"
  # After we have the integrated data, remove the dropVIZ data before we do clustering and dimension reduction
  seuratIntegrated = SetIdent(seuratIntegrated,value="tissue")
  seuratIntegrated = subset(seuratIntegrated,idents=tissue)
  # Run the standard workflow for visualization and clustering
  seuratIntegrated <- ScaleData(object = seuratIntegrated, verbose = TRUE)
  seuratIntegrated <- RunPCA(object = seuratIntegrated, npcs = numDims, verbose = FALSE)
  seuratIntegrated <- RunUMAP(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- RunTSNE(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- FindNeighbors(object = seuratIntegrated, reduction = "pca", dims = 1:numDims, k.param = kParam)
  seuratIntegrated <- FindClusters(object = seuratIntegrated, resolution = resUse, verbose = T, reduction = "pca")
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEcondition"))), dir.create(file.path(paste0("./Plots/ccaTSNEcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPcondition"))), dir.create(file.path(paste0("./Plots/ccaUMAPcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEresolution/",fileName))), dir.create(file.path(paste0("./Plots/ccaTSNEresolution/",fileName)),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPresolution/",fileName))), dir.create(file.path(paste0("./Plots/ccaUMAPresolution/",fileName)),recursive = T), FALSE)
  pdf(file=paste0("./Plots/ccaUMAPcondition/",fileName,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = combineLevel))
  dev.off()
  pdf(file=paste0("./Plots/ccaTSNEcondition/",fileName,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = combineLevel))
  dev.off()
  resolutions <- as.character(resUse)
  for(reso in resolutions){
    pdf(file=paste0("./Plots/ccaTSNEresolution/",fileName,"/PCA_tSNEres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "tsne"))
    dev.off()
    pdf(file=paste0("./Plots/ccaUMAPresolution/",fileName,"/PCA_UMAPres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "umap"))
    dev.off()
  }
  return(seuratIntegrated)
}

########################################################
# NewNames = c("Interneuron_Gad2","Neuron_Subiculum_Slc17a6","Neuron_Subiculum_Entorhinal_Nxph3","Neuron_Dentate_C1ql2",
#              "Neuron_CA1_Subiculum_Postsubiculum_Entorhinal_Fibcd1","Neuron_CA2CA3_Pvrl3","Astrocyte_Gja1",
#              "Oligodendrocyte_Tfr","Polydendrocyte_Tnr","Microglia_Macrophage_C1qb","Ependyma","Choroid_Plexus_Ttr",
#              "Neurogenesis_Sox4","Neuron_CajalRetzius_Lhx1","Endothelial_Flt1","Mural_Rgs5Acta2","Fibroblast-Like_Dcn")
# clusterSubset = c("Interneuron_Gad2","Neuron_Subiculum_Slc17a6","Neuron_Subiculum_Entorhinal_Nxph3","Neuron_Dentate_C1ql2",
#                   "Neuron_CA1_Subiculum_Postsubiculum_Entorhinal_Fibcd1","Neuron_CA2CA3_Pvrl3")
# seuratObject=neuron.subset
# dropVIZdata="../Seurat/DropVizData/Hippocampus/F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"
# dropVIZclusters="../Seurat/DropVizData/Hippocampus/F_GRCm38.81.P60Hippocampus.cluster.assign.RDS"
# dropVIZsubclusters="../Seurat/DropVizData/Hippocampus/F_GRCm38.81.P60Hippocampus.subcluster.assign.RDS"
# oldNames=seq(1,17)
# newNames=NewNames
# subsetSize=5000
# numFeatures=2000
# numDims=20
# clusterConfidence=0.5
# subclusterConfidence=0.3
# tissue="Hip.Neurons"
# defaultAssay="RNA"
# subsetClusters=clusterSubset

PredictCellType <- function(seuratObject=dropEST.combined.filtered,dropVIZdata=NULL,
                            dropVIZclusters=NULL,dropVIZsubclusters=NULL,oldNames=NULL,
                            newNames=NULL,subsetSize=5000,numFeatures=5000,numDims=20,
                            clusterConfidence=0.5,subclusterConfidence=0.3,tissue=NULL,
                            defaultAssay="RNA",subsetClusters=FALSE){
  ########## Set up DropVIZ
  #Read in dropVIZ data to predict the cell types
  dropVizDGE <- loadSparseDge(dropVIZdata) 
  #Read in the DropViz Cluster Assignments
  clusterAssign = readRDS(file = dropVIZclusters)
  subclusterAssign = readRDS(file = dropVIZsubclusters)
  #Subset to only cells which have a subcluster assignment
  clusterAssign = clusterAssign[which(names(clusterAssign) %in% names(subclusterAssign))]
  clusterAssign = clusterAssign[names(subclusterAssign)]
  dropVizMetaData = data.frame(cluster = clusterAssign, subcluster = subclusterAssign)
  #change from cluster numbers to cell type names
  dropVizMetaData$cluster = mapvalues(dropVizMetaData$cluster, from = oldNames, to = newNames)
  # Set up DropViz object
  dropVizSeurat = CreateSeuratObject(counts = dropVizDGE, meta.data = dropVizMetaData)
  dropVizSeurat = SetIdent(dropVizSeurat,value = "cluster")
  # Keep only a subset of the cell populations?
  if(!isFALSE(subsetClusters)){
    dropVizSeurat = subset(dropVizSeurat,idents=subsetClusters)
  }
  # Subset to make it run (runs out of memory if we don't subset)
  dropVizSeurat = subset(dropVizSeurat,max.cells.per.ident = subsetSize)
  # Clean up environment
  rm(dropVizDGE,dropVizMetaData,clusterAssign,subclusterAssign)
  ########## Predict Cell Type
  DefaultAssay(object = seuratObject) <- defaultAssay
  #Set up list of datasets
  data.list = list()
  data.list[["dropVIZ"]] = dropVizSeurat
  data.list[["dropEST"]] = seuratObject
  #Find variable genes
  for (i in 1:length(x = data.list)) {
    # data.list[[i]] <- SCTransform(object = data.list[[i]] , verbose = TRUE)
    data.list[[i]] <- NormalizeData(object = data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(object = data.list[[i]], selection.method = "vst", nfeatures = numFeatures, verbose = FALSE)
  }
  
  
  #MODIFY???? -- SET GENES TO 4K?
  # #set the overlapping ones as the variable genes
  # overlappingGenes = intersect(data.list[['dropEST']]@assays$RNA@var.features,data.list[['dropVIZ']]@assays$RNA@var.features)
  # length(overlappingGenes)
  # # #remove a particular gene?
  # overlappingGenes = overlappingGenes[-which(overlappingGenes == "Ttr")]
  # length(overlappingGenes) #1823
  # #set the variable genes to the overlapping genes
  # data.list[['dropEST']]@assays$RNA@var.features = overlappingGenes
  # data.list[['dropVIZ']]@assays$RNA@var.features = overlappingGenes
  
  
  
  data.query <- data.list[["dropEST"]]
  #dims 1:20, reduction = cca from biorxiv paper
  data.anchors <- FindTransferAnchors(reference = data.list[["dropVIZ"]], query = data.query, dims = 1:numDims, npcs = numDims)
  #dims 1:20, reduction = cca from biorxiv paper
  predictions1 <- TransferData(anchorset = data.anchors, refdata = data.list[["dropVIZ"]]$cluster, dims = 1:numDims)
  predictions2 <- TransferData(anchorset = data.anchors, refdata = data.list[["dropVIZ"]]$subcluster, dims = 1:numDims)
  
  predictions1.confident = predictions1[,1]
  predictions1.confident[which(predictions1$prediction.score.max<clusterConfidence)] = "unclassified"
  
  predictions2.confident = predictions2[,1]
  predictions2.confident[which(predictions2$prediction.score.max<subclusterConfidence)] = "unclassified"
  
  seuratObject <- AddMetaData(object = seuratObject, metadata = predictions1.confident, col.name = "Cluster.Prediction")
  seuratObject <- AddMetaData(object = seuratObject, metadata = predictions2.confident, col.name = "Subcluster.Prediction")
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  ifelse(!dir.exists(file.path("./Plots","./EvaluateDropVIZ")), dir.create(file.path("./Plots","./EvaluateDropVIZ")), FALSE)
  
  #check tSNE plot
  seuratObject = SetIdent(seuratObject, value = "Cluster.Prediction")
  p1 = DimPlot(seuratObject,label = T,pt.size = 0.5,reduction="umap",order = c("unclassified"),
               cols = rev(c("grey",gg_color_hue(length(levels(seuratObject@active.ident))-1))))
  pdf(file = paste0("./Plots/EvaluateDropVIZ/",tissue,"AllClusterClassification.pdf"),width = 14,height = 8)
  print(p1)
  dev.off()
  
  return(seuratObject)
}

########################################################
DetectDoublets <- function(seuratObject=dropEST.combined.filtered,numDims=10,doubletRate=0.075,cellLabels="CellType",
                           tissue=NULL){
  ## pK Identification ---------------------------------------------------------------------------------------------------------
  DefaultAssay(object = seuratObject) <- "RNA"
  # seuratObject <- SCTransform(object = seuratObject, verbose = TRUE, vars.to.regress = "nCount_RNA")
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- ScaleData(object = seuratObject, features = rownames(x = seuratObject), vars.to.regress = "nCount_RNA")
  seuratObject = FindVariableFeatures(seuratObject)
  seuratObject <- RunPCA(seuratObject)
  seuratObject <- RunTSNE(seuratObject, dims.use = 1:numDims, verbose=TRUE)
  sweep.res.list_seurat <- paramSweep_v3(seuratObject)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = FALSE)
  bcmvn_seurat <- find.pK(sweep.stats_seurat) #optimal pK=0.005
  optimal.pK = unfactor(bcmvn_seurat$pK[which(bcmvn_seurat$BCmetric==max(bcmvn_seurat$BCmetric))])
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(seuratObject@meta.data[,cellLabels])
  nExp_poi <- round(doubletRate*length(rownames(seuratObject@meta.data))) #assume doublet rate of 7.5%
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seuratObject <- doubletFinder_v3(seuratObject, pN = 0.25, pK = optimal.pK, nExp = nExp_poi, reuse.pANN = FALSE)
  pANN1 = colnames(seuratObject@meta.data)[grep("pANN",colnames(seuratObject@meta.data))][1]
  seuratObject <- doubletFinder_v3(seuratObject, pN = 0.25, pK = optimal.pK, nExp = nExp_poi.adj, reuse.pANN = pANN1)
  ## Plot results --------------------------------------------------------------------------------------------------------------
  DF.class1 = colnames(seuratObject@meta.data)[grep("DF.classifications",colnames(seuratObject@meta.data))][1]
  DF.class2 = colnames(seuratObject@meta.data)[grep("DF.classifications",colnames(seuratObject@meta.data))][2]
  seuratObject@meta.data[,"DF_hi.lo"] <- seuratObject@meta.data[,DF.class1]
  seuratObject@meta.data$DF_hi.lo[which(seuratObject@meta.data$DF_hi.lo == "Doublet" & seuratObject@meta.data[,DF.class2] == "Singlet")] <- "Doublet_lo"
  seuratObject@meta.data$DF_hi.lo[which(seuratObject@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPdoublets"))), dir.create(file.path(paste0("./Plots/ccaUMAPdoublets")),recursive = T), FALSE)
  pdf(file=paste0("./Plots/ccaUMAPdoublets/",tissue,".pdf"))
  print(DimPlot(seuratObject, group.by="DF_hi.lo", order=c("Doublet_hi","Doublet_lo","Singlet"), cols=c("black","gold","red")))
  dev.off()
  #remove doublets
  seuratObject = SetIdent(seuratObject,value="DF_hi.lo")
  seuratObject = subset(seuratObject, idents = "Singlet")
  return(seuratObject)
}

########################################################
runFeaturePlot <- function(seuratObject=dropEST.combined.filtered,markers=cellTypeMarkers,resUse="integrated_snn_res.1",
                           tissue=NULL,dataTypesPlot=c("rna","integrated","sct")){
  seuratObject = SetIdent(seuratObject,value=resUse)
  for(dataType in dataTypesPlot){
    ifelse(!dir.exists(file.path(paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType))), dir.create(file.path(paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType)),recursive = T), FALSE)
    for(i in 1:length(markers[[1]])){
      if(dataType == "integrated"){
        if(markers[[1]][i] %in% rownames(seuratObject@assays$integrated@data)){
          pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/FeaturePlot_",markers[[2]][i],"MarkerGene.pdf"))
          print(FeaturePlot(seuratObject,paste0(dataType,"_",markers[[1]][i]),cols = c("lightgrey","blue")))
          dev.off()
        }
      } else{
        if(dataType == "sct"){
          if(markers[[1]][i] %in% rownames(seuratObject@assays$SCT@data)){
            pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/FeaturePlot_",markers[[2]][i],"MarkerGene.pdf"))
            print(FeaturePlot(seuratObject,paste0(dataType,"_",markers[[1]][i]),cols = c("lightgrey","blue")))
            dev.off()
          }
        } else{
          pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/FeaturePlot_",markers[[2]][i],"MarkerGene.pdf"))
          print(FeaturePlot(seuratObject,paste0(dataType,"_",markers[[1]][i]),cols = c("lightgrey","blue")))
          dev.off()
        }
      }
    }
    for(i in 1:length(markers[[1]])){
      if(dataType == "integrated"){
        if(markers[[1]][i] %in% rownames(seuratObject@assays$integrated@data)){
          pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/VlnPlot_",markers[[2]][i],"MarkerGene.pdf"),height=8,width=16)
          print(VlnPlot(seuratObject,features = paste0(dataType,"_",markers[[1]][i]),sort=T))
          dev.off()
        }
      } else{
        if(dataType == "sct"){
          if(markers[[1]][i] %in% rownames(seuratObject@assays$SCT@data)){
            pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/VlnPlot_",markers[[2]][i],"MarkerGene.pdf"),height=8,width=16)
            print(VlnPlot(seuratObject,features = paste0(dataType,"_",markers[[1]][i]),sort=T))
            dev.off()
          }
        } else{
          pdf(file=paste0("./IdentifyClusters/MarkerGenes/CCA/",tissue,"/",dataType,"/VlnPlot_",markers[[2]][i],"MarkerGene.pdf"),height=8,width=16)
          print(VlnPlot(seuratObject,features = paste0(dataType,"_",markers[[1]][i]),sort=T))
          dev.off()
        }
      }
    }
  }
}

########################################################
runRBCfiltering <- function(seuratObject=dropEST.combined.filtered,tissue=NULL,rbcThresh=0.15){
  rbc.features <- grep(pattern = "Hba-", x = rownames(x = seuratObject), value = TRUE)
  rbc.features <- c(rbc.features, grep(pattern = "Hbb-", x = rownames(x = seuratObject), value = TRUE))
  rbc.features <- c(rbc.features, grep(pattern = "Hbq", x = rownames(x = seuratObject), value = TRUE))
  percent.rbc <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts")[rbc.features,])/Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts"))
  seuratObject$percent.rbc = percent.rbc
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAP.RBCfilter"))), dir.create(file.path(paste0("./Plots/ccaUMAP.RBCfilter")),recursive = T), FALSE)
  pdf(file=paste0("./Plots/ccaUMAP.RBCfilter/",tissue,".pdf"))
  print(FeaturePlot(seuratObject,"percent.rbc"))
  dev.off()
  seuratObject = subset(x = seuratObject, subset = percent.rbc < rbcThresh)
  return(seuratObject)
}

########################################################
runJackStraw <- function(seuratObject=dropEST.combined.filtered,numReplicates=100,numDims=50,tissue=NULL){
  seuratObject <- JackStraw(object = seuratObject, reduction = "pca", num.replicate = numReplicates, verbose = TRUE, dims = numDims)
  seuratObject = ScoreJackStraw(object = seuratObject, dims = 1:numDims, reduction = "pca")
  ifelse(!dir.exists(file.path("./Plots/QCPlots/JackStraw")), dir.create(file.path("./Plots/QCPlots/JackStraw"),recursive = T), FALSE)
  pdf(file=paste0("./Plots/QCPlots/JackStraw/",tissue,".pdf"),height=10,width=15)
  print(JackStrawPlot(object = seuratObject, dims = 1:50))
  dev.off()
  return(seuratObject)
}

########################################################
runClusteringAndDimReduction <- function(seuratObject=dropEST.combined.filtered,dimRed1="ica",numDims=NULL,kParam=25,
                                         resUse=seq(0.5,4,by=0.5),tissue=NULL,sampleMetaLabel="orig.ident",
                                         treatmentMetaLabel="timpoint.condition",dataLevel="SCT"){
  
  #Find clusters using graph-based method (check with many different resolutions which we will plot)
  seuratObject <- FindNeighbors(object = seuratObject, reduction = dimRed1, dims = 1:numDims, k.param = kParam)
  seuratObject <- FindClusters(object = seuratObject, resolution = resUse, verbose = T, reduction = dimRed1)

  #Use tSNE and UMAP to visualize in 2 dimensions
  ifelse(!dir.exists(file.path(paste0("./Plots/FirstTSNE/",tissue))), dir.create(file.path(paste0("./Plots/FirstTSNE/",tissue)),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/FirstUMAP/",tissue))), dir.create(file.path(paste0("./Plots/FirstUMAP/",tissue)),recursive = T), FALSE)
  ifelse(!dir.exists(file.path("./Plots/TreatmentUMAP")), dir.create(file.path("./Plots/TreatmentUMAP"),recursive = T), FALSE)
  ifelse(!dir.exists(file.path("./Plots/SampleUMAP")), dir.create(file.path("./Plots/SampleUMAP"),recursive = T), FALSE)
  ifelse(!dir.exists(file.path("./Plots/TreatmentTSNE")), dir.create(file.path("./Plots/TreatmentTSNE"),recursive = T), FALSE)
  ifelse(!dir.exists(file.path("./Plots/SampleTSNE")), dir.create(file.path("./Plots/SampleTSNE"),recursive = T), FALSE)

  seuratObject <- RunUMAP(object = seuratObject, reduction = dimRed1, dims = 1:numDims)
  seuratObject <- RunTSNE(object = seuratObject, reduction = dimRed1, dims = 1:numDims, check_duplicates = FALSE)
  resolutions <- as.character(resUse)
  for(reso in resolutions){
    if(dimRed1=="ica"){
      pdf(file=paste0("./Plots/FirstTSNE/",tissue,"/ICA_tSNEres",reso,".pdf"))
    } else{
      pdf(file=paste0("./Plots/FirstTSNE/",tissue,"/PCA_tSNEres",reso,".pdf"))
    }
    seuratObject <- SetIdent(seuratObject, value = paste0(dataLevel,"_snn_res.",reso))
    print(DimPlot(seuratObject, label = T, reduction = "tsne"))
    dev.off()
    if(dimRed1=="ica"){
      pdf(file=paste0("./Plots/FirstUMAP/",tissue,"/ICA_UMAPres",reso,".pdf"))
    } else{
      pdf(file=paste0("./Plots/FirstUMAP/",tissue,"/PCA_UMAPres",reso,".pdf"))
    }
    seuratObject <- SetIdent(seuratObject, value = paste0(dataLevel,"_snn_res.",reso))
    print(DimPlot(seuratObject, label = T, reduction = "umap"))
    dev.off()
  }

  seuratObject <- SetIdent(seuratObject, value = treatmentMetaLabel)
  if(dimRed1=="ica"){
    pdf(file=paste0("./Plots/TreatmentUMAP/ICA_",tissue,"_TreatmentUMAP.pdf"))
  } else{
    pdf(file=paste0("./Plots/TreatmentUMAP/PCA_",tissue,"_TreatmentUMAP.pdf"))
  }
  print(DimPlot(seuratObject,do.label = F,pt.size = 0.5, reduction = "umap"))
  dev.off()
  if(dimRed1=="ica"){
    pdf(file=paste0("./Plots/TreatmentTSNE/ICA_",tissue,"_TreatmentTSNE.pdf"))
  } else{
    pdf(file=paste0("./Plots/TreatmentTSNE/PCA_",tissue,"_TreatmentTSNE.pdf"))
  }
  print(DimPlot(seuratObject,do.label = F,pt.size = 0.5, reduction = "tsne"))
  dev.off()

  seuratObject <- SetIdent(seuratObject, value = sampleMetaLabel)
  if(dimRed1=="ica"){
    pdf(file=paste0("./Plots/SampleUMAP/ICA_",tissue,"_SampleUMAP.pdf"))
  } else{
    pdf(file=paste0("./Plots/SampleUMAP/PCA_",tissue,"_SampleUMAP.pdf"))
  }
  print(DimPlot(seuratObject,do.label = F,pt.size = 0.5, reduction = "umap"))
  dev.off()
  if(dimRed1=="ica"){
    pdf(file=paste0("./Plots/SampleTSNE/ICA_",tissue,"_SampleTSNE.pdf"))
  } else{
    pdf(file=paste0("./Plots/SampleTSNE/PCA_",tissue,"_SampleTSNE.pdf"))
  }
  print(DimPlot(seuratObject,do.label = F,pt.size = 0.5, reduction = "tsne"))
  dev.off()
  return(seuratObject)
}

########################################################
plotTSNEandDendrogram <- function(seuratObject=dropEST.combined.filtered,tissue=NULL,cellMetaDataLabel="CellType",lineageGroups=NULL,
                                  lineageNumberSize=6,lineageTextSize=6,lineageFill=NULL,reductionType="tsne",clusterColors=NULL,
                                  clusterAnnotationSize=6,cellNumSize=6,legendTextSize=14,clusterPointSize=2,legendPtSize=10,
                                  legendSpacing=1.5,outputFile=NULL,fileHeight,fileWidth){
  seuratObject = SetIdent(seuratObject,value = cellMetaDataLabel)
  #Get the mean value of all highly variable genes for each cluster to compute euclidean distance
  allHVGmeans = NULL
  #Iterate through each cell type
  for(cellType in levels(seuratObject@meta.data[,cellMetaDataLabel])){
    #Get all single cells for each cell type
    cellTypeCells = WhichCells(seuratObject,idents = cellType)
    #Extract the data matrix of highly variable genes for the single cells from the cell type
    cellTypeHVG = data.matrix(seuratObject@assays$RNA@data[seuratObject@assays$RNA@var.features,cellTypeCells])
    HVGmeans = apply(cellTypeHVG, 1, FUN = mean)
    if(is.null(allHVGmeans)){
      allHVGmeans = HVGmeans
    } else{
      allHVGmeans = rbind(allHVGmeans,HVGmeans)
    }
  }
  rownames(allHVGmeans) = levels(seuratObject@active.ident)
  dd <- dist(allHVGmeans, method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  dend <- as.dendrogram(hc)
  dend_data <- dendro_data(dend, type = "rectangle")
  #Add cluster counts to dendrogram (will use as labels)
  ClusterCounts = table(seuratObject@meta.data[,cellMetaDataLabel])
  labelsReplace = data.frame(ClusterCounts[unfactor(dend_data$labels$label)])
  dend_data$labels$counts = labelsReplace$Freq
  
  # Draw colored rectangles on the dendrogram
  # Labels and boxes on the dendrogram
  labelCells = lineageGroups
  #Get the min and max x-axis values of the cell types in each group
  rectangleCoords = NULL
  for(lineage in names(labelCells)){
    cellTypes = labelCells[[lineage]]
    #Get the x and y limits for the rectangle for the lineages
    xaxisRange = c(min(dend_data$labels$x[which(dend_data$labels$label %in% cellTypes)]),max(dend_data$labels$x[which(dend_data$labels$label %in% cellTypes)]))
    yaxisRange = c(min(dend_data$segments[which(dend_data$segments$x==xaxisRange[1] | dend_data$segments$xend==xaxisRange[1]),c("y","yend")],
                       dend_data$segments[which(dend_data$segments$x==xaxisRange[2] | dend_data$segments$xend==xaxisRange[2]),c("y","yend")]),
                   max(dend_data$segments[which(dend_data$segments$x==xaxisRange[1] | dend_data$segments$xend==xaxisRange[1]),c("y","yend")],
                       dend_data$segments[which(dend_data$segments$x==xaxisRange[2] | dend_data$segments$xend==xaxisRange[2]),c("y","yend")]))
    #Need to be careful with lineage text when it is only one cell type
    if(length(cellTypes)==1){
      #Get the location of where the lineage text should be written
      lineageSegementY = (yaxisRange[2] - yaxisRange[1])/3
      lineageSegementX = xaxisRange[1]
      if(is.null(rectangleCoords)){
        rectangleCoords = data.frame(x1=xaxisRange[1],x2=xaxisRange[2],y1=yaxisRange[1],y2=yaxisRange[2],lineage=lineage,Lx=lineageSegementX,Ly=lineageSegementY)
      }
      else{
        rectangleCoords = rbind(rectangleCoords,data.frame(x1=xaxisRange[1],x2=xaxisRange[2],y1=yaxisRange[1],y2=yaxisRange[2],lineage=lineage,Lx=lineageSegementX,Ly=lineageSegementY))
      }
    } else{
      #Get the location of where the lineage text should be written
      lineageSegement = dend_data$segments[which(dend_data$segments$y == yaxisRange[2] | dend_data$segments$yend == yaxisRange[2]),]
      maxVal = max(max(lineageSegement$y),max(lineageSegement$yend))
      lineageSegement = lineageSegement[which(lineageSegement$y==maxVal | lineageSegement$yend==maxVal),]
      if(is.null(rectangleCoords)){
        rectangleCoords = data.frame(x1=xaxisRange[1],x2=xaxisRange[2],y1=yaxisRange[1],y2=yaxisRange[2],lineage=lineage,Lx=lineageSegement$x,Ly=min(lineageSegement$y,lineageSegement$yend))
      }
      else{
        rectangleCoords = rbind(rectangleCoords,data.frame(x1=xaxisRange[1],x2=xaxisRange[2],y1=yaxisRange[1],y2=yaxisRange[2],lineage=lineage,Lx=lineageSegement$x,Ly=min(lineageSegement$y,lineageSegement$yend)))
      }
    }
  }
  #Make the boxes a little outside of the lines
  rectangleCoords$x1 = rectangleCoords$x1 - (max(c(dend_data$segments$x,dend_data$segments$xend))-min(c(dend_data$segments$x,dend_data$segments$xend)))/100
  rectangleCoords$x2 = rectangleCoords$x2 + (max(c(dend_data$segments$x,dend_data$segments$xend))-min(c(dend_data$segments$x,dend_data$segments$xend)))/100
  rectangleCoords$y1 = rectangleCoords$y1 - (max(c(dend_data$segments$y,dend_data$segments$yend))-min(c(dend_data$segments$y,dend_data$segments$yend)))/100
  rectangleCoords$y2 = rectangleCoords$y2 + (max(c(dend_data$segments$y,dend_data$segments$yend))-min(c(dend_data$segments$y,dend_data$segments$yend)))/100
  
  labelAdjust = (max(c(dend_data$segments$y,dend_data$segments$yend))-min(c(dend_data$segments$y,dend_data$segments$yend)))/50
  lineageLabelAdjustX = (max(c(dend_data$segments$x,dend_data$segments$xend))-min(c(dend_data$segments$x,dend_data$segments$xend)))/25
  lineageLabelAdjustY = (max(c(dend_data$segments$y,dend_data$segments$yend))-min(c(dend_data$segments$y,dend_data$segments$yend)))/50
  plotHeight = max(c(dend_data$segments$y,dend_data$segments$yend))+2
  
  
  # Plot line segments and add labels
  p1 <- ggplot(dend_data$segments) + 
    geom_rect(data=rectangleCoords, mapping = aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=lineage)) + #Draw a rectangle
    scale_fill_manual(values=lineageFill) + #set the colors for the rectangles
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+ #Draw the dendrogram
    geom_text(data = dend_data$labels, aes(x, y-labelAdjust, label = counts), #Label the clusters with the number of cells in the cluster
              hjust = 1, angle = 0, size = lineageNumberSize)+
    geom_text(data = rectangleCoords, aes(Lx + lineageLabelAdjustX, Ly + lineageLabelAdjustY, label = lineage), #Label the lineage boxes
              hjust = 0, angle = 0, size = lineageTextSize)+
    ylim(-4, plotHeight) + #Set limits so we can see the cluster labels
    coord_flip() + #make horizontal
    theme_bw() + #remove grey background
    theme(panel.grid.major = element_blank(), #remove gridlines 
          panel.grid.minor = element_blank(), #remove gridlines 
          axis.title.x=element_blank(), #remove x axis title
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), #remove y axis title
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none", #remove the legend
          panel.border = element_blank()) #remove the border from the panel

  #Reverse order so we can plot correctly
  dend_data$labels$label = rev(as.character(dend_data$labels$label))
  dend_data$labels$label = as.factor(dend_data$labels$label)
  dend_data$labels$counts = rev(dend_data$labels$counts)
  
  #Combine with tSNE plot or UMAP
  
  #Use tSNE or UMAP and extract cell embeddings
  if(reductionType=="umap"){
    reductionPlot = data.frame(seuratObject@reductions$umap@cell.embeddings)
  } else{
    reductionPlot = data.frame(seuratObject@reductions$tsne@cell.embeddings)
  }
  #Add cell type labels
  reductionPlot$group = seuratObject@meta.data[rownames(reductionPlot),cellMetaDataLabel]
  #Add in numbers for the clusters
  ClusterNum = seuratObject@meta.data[,cellMetaDataLabel]
  ClusterNum = unfactor(ClusterNum)
  names(ClusterNum) = rownames(seuratObject@meta.data)
  for(i in 1:length(unfactor(dend_data$labels$label))){
    ClusterNum[which(ClusterNum == dend_data$labels$label[i])] = i
  }
  seuratObject$ClusterNum = ClusterNum
  reductionPlot$cluster = seuratObject$ClusterNum[rownames(reductionPlot)]
  
  #Reorder the cell types to correspond to the dendrogram
  for(lineage in names(lineageGroups)){
    for(cellType in lineageGroups[[lineage]])
      reductionPlot$group <- relevel(reductionPlot$group, cellType)
  }
  reductionPlot$group <- factor(reductionPlot$group, levels=rev(levels(reductionPlot$group)))
  #Put cluster number in the correct order
  clustOrder = as.character(unique(as.numeric(reductionPlot$cluster))[order(unique(as.numeric(reductionPlot$cluster)))])
  reductionPlot$cluster <- factor(reductionPlot$cluster, levels=clustOrder)
  
  #Find where to plot the cluster numbers
  allClusterMedians = NULL
  for(cluster in unique(reductionPlot$cluster)){
    if(reductionType=="umap"){
      dim1median = median(reductionPlot$UMAP_1[which(reductionPlot$cluster==cluster)])
      dim2median = median(reductionPlot$UMAP_2[which(reductionPlot$cluster==cluster)])
    } else{
      dim1median = median(reductionPlot$tSNE_1[which(reductionPlot$cluster==cluster)])
      dim2median = median(reductionPlot$tSNE_2[which(reductionPlot$cluster==cluster)])
    }
    if(is.null(allClusterMedians)){
      allClusterMedians = data.frame(dim1median,dim2median,cluster)
    } else{
      allClusterMedians = rbind(allClusterMedians,data.frame(dim1median,dim2median,cluster))
    }
  }
  
  #Number of cells
  annotations  <- data.frame(xpos = c(Inf),
                             ypos = c(-Inf),
                             annotateText = c(paste0("Cell Number:\n",nrow(seuratObject@meta.data)," ")),
                             hjustvar = c(1.1),
                             vjustvar = c(-0.25))
  if(reductionType=="umap"){
    p2 <- ggplot(reductionPlot, aes(x=UMAP_1, y=UMAP_2)) 
  } else{
    p2 <- ggplot(reductionPlot, aes(x=tSNE_1, y=tSNE_2)) 
  }
  p2 <- p2 +
    geom_point(aes(color=cluster),size=clusterPointSize) +
    scale_color_manual(values=clusterColors, labels = c(paste(dend_data$labels$x,dend_data$labels$label,sep = "   "))) +
    geom_text(data = allClusterMedians, aes(dim1median, dim2median, label = cluster), #Label the clusters the cluster number
              hjust = 1, angle = 0, size = clusterAnnotationSize) +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = cellNumSize) +
    theme_bw() + #remove grey background
    theme(panel.grid.major = element_blank(), #remove gridlines 
          panel.grid.minor = element_blank(), #remove gridlines 
          legend.title=element_blank(), #remove legend title
          legend.text=element_text(size=legendTextSize), #increaese legend font size
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.border = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=legendPtSize))) + #Change size of points in legend
    theme(legend.key = element_rect(size = 6),
          legend.key.height = unit(legendSpacing, "cm"), #Change the vertical spacing between elements in the legend
          legend.key.width = unit(1, "cm"))
  
  #arrange plots for output
  p3 <- ggarrange(p2, p1, 
                  labels = c("A", ""),
                  ncol = 2, nrow = 1,
                  widths = c(2, 1),
                  heights = c(1))
  
  ifelse(!dir.exists(file.path("./Plots/Figure1")), dir.create(file.path("./Plots/Figure1"),recursive = T), FALSE)
  pdf(file=paste0("./Plots/Figure1/",tissue,".",reductionType,".pdf"),height=fileHeight,width=fileWidth)
  print(p3)
  dev.off()
  
  return(seuratObject)
}

########################################################
# Main DEG function -- can wrap it to include additional variable levels (e.g. timepoint and condition)
runCellTypeDEGsMain <- function(seuratObject=metaObject1Subset,metaData1=NULL,cellMetaDataLabel="CellTypeRefined",metaObject1=NULL,
                                groupsCompare=NULL,sampleNameMetaEntry=NULL,compareMetaDataLabel=NULL,minCounts=5,
                                cellTypeDEGsNovelMethod=NULL,cellTypeDEGsStandard=NULL){
  
  metaObject1Subset = seuratObject
  
  #iterate through each cell type
  for(cellType in levels(metaObject1Subset@active.ident)){
    print(cellType)
    #generate a cell type subset
    cellTypeSubset = subset(metaObject1Subset,idents = cellType)
    #get Group1 sample names
    group1Samples = unique(cellTypeSubset@meta.data[which(cellTypeSubset@meta.data[,compareMetaDataLabel]==groupsCompare[1]),sampleNameMetaEntry])
    #get Group2 sample names
    group2Samples = unique(cellTypeSubset@meta.data[which(cellTypeSubset@meta.data[,compareMetaDataLabel]==groupsCompare[2]),sampleNameMetaEntry])
    
    #set identity to sample level
    cellTypeSubset = SetIdent(cellTypeSubset, value = sampleNameMetaEntry)
    
    #remove samples which have less than 5 cells
    sampleCounts = table(cellTypeSubset@active.ident)
    samplesRemove = names(sampleCounts)[which(sampleCounts<minCounts)]
    
    group1Overlap = which(group1Samples %in% samplesRemove)
    if(length(group1Overlap)>0){
      group1Samples = group1Samples[-group1Overlap]
    }
    
    group2Overlap = which(group2Samples %in% samplesRemove)
    if(length(group2Overlap)>0){
      group2Samples = group2Samples[-group2Overlap]
    }
    
    if(length(group1Samples) > 0 & length(group2Samples) > 0){
      
      #list for overlapping genes
      allGenes = list()
      
      #get DEGs between each group1 sample and all group2 samples using Wilcoxon rank sum test
      group1DEGs = list()
      for(i in 1:length(group1Samples)){
        group1DEGs[[group1Samples[i]]] = FindMarkers(cellTypeSubset, ident.1 = group1Samples[i], ident.2 = group2Samples)
        colnames(group1DEGs[[group1Samples[i]]]) = paste0(group1Samples[i],colnames(group1DEGs[[group1Samples[i]]]))
        allGenes[[group1Samples[i]]] = rownames(group1DEGs[[group1Samples[i]]])
      }
      
      #get DEGs between each group2 sample and all group1 samples using Wilcoxon rank sum test
      group2DEGs = list()
      for(i in 1:length(group2Samples)){
        group2DEGs[[group2Samples[i]]] = FindMarkers(cellTypeSubset, ident.1 = group2Samples[i], ident.2 = group1Samples)
        colnames(group2DEGs[[group2Samples[i]]]) = paste0(group2Samples[i],colnames(group2DEGs[[group2Samples[i]]]))
        allGenes[[group2Samples[i]]] = rownames(group2DEGs[[group2Samples[i]]])
      }
      
      #get overlapping genes from all the statistical tests at the sample level
      overlappingGenes = Reduce(intersect, allGenes)
      
      #generate a matrix with the pvals, FC, and adj p vals
      for(i in 1:length(group1Samples)){
        if(i==1){
          combinedMatrix = group1DEGs[[group1Samples[i]]][overlappingGenes,c(1,2,5)]
        } else{
          combinedMatrix = cbind(combinedMatrix,group1DEGs[[group1Samples[i]]][overlappingGenes,c(1,2,5)])
        }
      }
      for(i in 1:length(group2Samples)){
        combinedMatrix = cbind(combinedMatrix,group2DEGs[[group2Samples[i]]][overlappingGenes,c(1,2,5)])
      }
      
      colsUseGroup1 = seq(from = 1, to = length(group1Samples)*3, by = 3)
      colsUseGroup2 = seq(from = length(group1Samples)*3+1, to = length(group1Samples)*3+1+(length(group2Samples)-1)*3, by = 3)
      
      minP = apply(combinedMatrix, 1, function(x) minimump(x[c(colsUseGroup1,colsUseGroup2)])$p)
      maxP = apply(combinedMatrix, 1, function(x) maximump(x[c(colsUseGroup1,colsUseGroup2)])$p)
      
      minPadj = p.adjust(minP,method = 'bonferroni',n = nrow(cellTypeSubset@assays$RNA@data))
      maxPadj = p.adjust(maxP,method = 'bonferroni',n = nrow(cellTypeSubset@assays$RNA@data))
      #add to the combined matrix
      combinedMatrix$minP = minP
      combinedMatrix$maxP = maxP
      combinedMatrix$minPadj = minPadj
      combinedMatrix$maxPadj = maxPadj
      
      #get the conversative p-val (max of minP and maxP)
      ConsP = apply(combinedMatrix, 1, function(x) max(x[c('minP','maxP')]))
      combinedMatrix$ConsP = ConsP
      
      #only keep genes which show the FC in the same direction
      colsUseGroup1 = seq(from = 2, to = length(group1Samples)*3, by = 3)
      colsUseGroup2 = seq(from = length(group1Samples)*3+2, to = length(group1Samples)*3+2+(length(group2Samples)-1)*3, by = 3)
      
      if(length(colsUseGroup1)==1){
        goodGroup1 = rownames(combinedMatrix)
      } else{
        goodGroup1 = names(which(rowSums(combinedMatrix[,colsUseGroup1] > 0)==length(group1Samples) | rowSums(combinedMatrix[,colsUseGroup1] < 0)==length(group1Samples)))
      }
      if(length(colsUseGroup2)==1){
        goodGroup2 = rownames(combinedMatrix)
      } else{
        goodGroup2 = names(which(rowSums(combinedMatrix[,colsUseGroup2] > 0)==length(group2Samples) | rowSums(combinedMatrix[,colsUseGroup2] < 0)==length(group2Samples)))
      }
      
      correctDirection = intersect(goodGroup1,goodGroup2)
      combinedMatrix = combinedMatrix[correctDirection,]
      #sort the matrix by the max pvalue
      combinedMatrix = combinedMatrix[order(combinedMatrix$ConsP),]
      
      if(!(is.null(metaData1))){
        cellTypeDEGsNovelMethod[[metaObject1]][[cellType]] = combinedMatrix
      } else{
        cellTypeDEGsNovelMethod[[cellType]] = combinedMatrix
      }
    }
    #do traditional Seurat DEG identification
    cellTypeSubset = SetIdent(cellTypeSubset, value = compareMetaDataLabel)
    #only proceed if both have at least 3 cells
    statusCounts = table(cellTypeSubset@active.ident)
    if(length(which(statusCounts>=3))==2){
      StandardDEGs = FindMarkers(cellTypeSubset, ident.1 = groupsCompare[1], ident.2 = groupsCompare[2])
      #add cell type DEGs to matrices
      if(!(is.null(metaData1))){
        cellTypeDEGsStandard[[metaObject1]][[cellType]] = StandardDEGs
      } else{
        cellTypeDEGsStandard[[cellType]] = StandardDEGs
      }
    }
  }
  return(list(cellTypeDEGsNovelMethod,cellTypeDEGsStandard))
}

########################################################
#Wrapper for the main DEG function -- allows us to add an additional meta data level (e.g. timepoint + condition rather than just condition)
runCellTypeDEGsNovelMethod <- function(seuratObject=dropEST.combined.filtered,metaData1=NULL,cellMetaDataLabel="CellTypeRefined",
                                       groupsCompare=NULL,sampleNameMetaEntry=NULL,compareMetaDataLabel=NULL,minCounts=5,tissue=NULL){
  
  #initialize a list to hold marker genes for all cell types
  cellTypeDEGsStandard = list()
  cellTypeDEGsNovelMethod = list()
  
  #iterate through the first supplied level of meta data
  if(!(is.null(metaData1))){
    seuratObject = SetIdent(seuratObject, value = metaData1)
    for(metaObject1 in levels(seuratObject@active.ident)){
      print(metaObject1)
      metaObject1Subset = subset(seuratObject,idents = metaObject1)
      metaObject1Subset = SetIdent(metaObject1Subset, value = cellMetaDataLabel)
      outputList = runCellTypeDEGsMain(seuratObject=metaObject1Subset,metaData1=metaData1,cellMetaDataLabel=cellMetaDataLabel,metaObject1=metaObject1,
                                       groupsCompare=groupsCompare,sampleNameMetaEntry=sampleNameMetaEntry,compareMetaDataLabel=compareMetaDataLabel,
                                       minCounts=minCounts,cellTypeDEGsNovelMethod=cellTypeDEGsNovelMethod,cellTypeDEGsStandard=cellTypeDEGsStandard)
      cellTypeDEGsNovelMethod = outputList[[1]]
      cellTypeDEGsStandard = outputList[[2]]
    }
  } else{
    metaObject1Subset = SetIdent(seuratObject, value = cellMetaDataLabel)
    outputList = runCellTypeDEGsMain(seuratObject=metaObject1Subset,metaData1=metaData1,cellMetaDataLabel=cellMetaDataLabel,metaObject1=NULL,
                                     groupsCompare=groupsCompare,sampleNameMetaEntry=sampleNameMetaEntry,compareMetaDataLabel=compareMetaDataLabel,
                                     minCounts=minCounts,cellTypeDEGsNovelMethod=cellTypeDEGsNovelMethod,cellTypeDEGsStandard=cellTypeDEGsStandard)
    cellTypeDEGsNovelMethod = outputList[[1]]
    cellTypeDEGsStandard = outputList[[2]]
  }
  #save DEGs
  saveRDS(cellTypeDEGsNovelMethod, file = paste0("./Saves/",tissue,".DEGsNovelMethod.rds"))
  saveRDS(cellTypeDEGsStandard, file = paste0("./Saves/",tissue,".DEGsStandardMethod.rds"))
}

########################################################
iterateDEGsWrapper <- function(degList = NULL, degMatrix = NULL, novel = TRUE, metaData1 = NULL, tissue = NULL, group = NULL, meta1 = NULL){
  for(cellType in names(degList)){
    cellTypeDEGs = degList[[cellType]]
    if(isTRUE(novel)){
      cellTypeDEGs = cellTypeDEGs[,c(2,(ncol(cellTypeDEGs)-4):ncol(cellTypeDEGs))]
      colnames(cellTypeDEGs)[1] = paste0(tissue,".",group[1],"_logFC")
    }
    cellTypeDEGs$CellType = cellType
    if(!(is.null(metaData1))){
      cellTypeDEGs[,metaData1] = meta1
    }
    cellTypeDEGs$Gene = rownames(cellTypeDEGs)
    rownames(cellTypeDEGs) = seq(1,nrow(cellTypeDEGs))
    if(is.null(degMatrix)){
      degMatrix = cellTypeDEGs
    } else{
      degMatrix = rbind(degMatrix,cellTypeDEGs)
    }
  }
  return(degMatrix)
}

#Write DEGs to File
writeDEGs <- function(cellTypeDEGsNovelMethod=NULL,cellTypeDEGsStandard=NULL,metaData1=NULL,tissue=NULL,groupsCompare=NULL){
  cellTypeDEGsNovelMethod = readRDS(file = cellTypeDEGsNovelMethod)
  cellTypeDEGsStandard = readRDS(file = cellTypeDEGsStandard)
  
  novelMethodDEGmatrix = NULL
  standardMethodDEGmatrix = NULL
  
  if(!is.null(metaData1)){
    for(meta1 in names(cellTypeDEGsNovelMethod)){
      novelMethodDEGmatrix = iterateDEGsWrapper(degList = cellTypeDEGsNovelMethod[[meta1]], degMatrix = novelMethodDEGmatrix, novel = TRUE, metaData1 = metaData1,
                                                tissue = tissue, group = groupsCompare[1], meta1 = meta1)
      standardMethodDEGmatrix = iterateDEGsWrapper(degList = cellTypeDEGsStandard[[meta1]], degMatrix = standardMethodDEGmatrix, novel = FALSE, metaData1 = metaData1,
                                                   tissue = tissue, group = groupsCompare[1], meta1 = meta1)
    }
  } else{
    novelMethodDEGmatrix = iterateDEGsWrapper(degList = cellTypeDEGsNovelMethod, degMatrix = novelMethodDEGmatrix, novel = TRUE, metaData1 = NULL,
                                              tissue = tissue, group = groupsCompare[1], meta1 = NULL)
    standardMethodDEGmatrix = iterateDEGsWrapper(degList = cellTypeDEGsStandard, degMatrix = standardMethodDEGmatrix, novel = FALSE, metaData1 = NULL,
                                                 tissue = tissue, group = groupsCompare[1], meta1 = NULL)
  }

  ifelse(!dir.exists(file.path("./CellTypeDEGs")), dir.create(file.path("./CellTypeDEGs"),recursive = T), FALSE)
  write.table(novelMethodDEGmatrix,file = paste0("./CellTypeDEGs/",tissue,".CellTypeDEGsNovelMethod.txt"),sep = "\t",quote = F,row.names = F)
  write.table(standardMethodDEGmatrix,file = paste0("./CellTypeDEGs/",tissue,".CellTypeDEGsStandardMethod.txt"),sep = "\t",quote = F,row.names = F)
}

########################################################
#Euclidean Distance

runEuclideanDistance <- function(seuratObject=dropEST.combined.filtered,metaData1=NULL,metaDataSubset=NULL,cellTypeMetaData="NeuronCellType",
                                       nPerm=1000,do.par=TRUE,num.cores="auto",comparisonMetaData="condition",group1="TBI",group2="Sham",
                                       distanceScorePlot=NULL,fcPlot=NULL,distanceSave=NULL,topGenesTrim=20,bottomFractionTrim=0.75,zScore=F){
  
  #Subset to metadata slot if specified
  if(!is.null(metaData1)){
    seuratObject = SetIdent(seuratObject, value = metaData1)
    seuratObject = subset(seuratObject, idents = metaDataSubset)
  }

  #Set ident to meta data level which holds the cell identities
  seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)
  
  #Number of permutations
  numPerm = nPerm
  
  #set up parallel processing
  if(do.par==TRUE){
    #auto detect number of available cores
    if(num.cores == "auto"){
      cores=detectCores()
      #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
      cores = max(1,cores[1]-1)
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    } else{
      #if a number of cores for parallel processing is specified, use this number of cores
      if(is.numeric(num.cores)){
        cores = num.cores
        #setup and start parallel processing cluster
        cl <- makeCluster(cores)
        registerDoSNOW(cl)
      }
    }
  }
  #if we are not using parallel processing, just use one core
  if(do.par==FALSE){
    cores = 1
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  
  #output list and vector for distances and permutation distances --> using euclidian distance
  cellTypePermutationDistances = list()
  cellTypeDistances = c()
  
  #iterate through each cell type
  for(cellType in levels(seuratObject@active.ident)){
    print(cellType)
    #subset to cell type of interest
    cellTypeSubset = subset(seuratObject,idents = cellType)
    cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
    
    #Get the top most highly expressed genes
    avgGeneExpression = rowMeans(data.matrix(cellTypeSubset@assays$RNA@data))
    avgGeneExpression = avgGeneExpression[order(avgGeneExpression,decreasing = T)]
    #Remove the lowly expressed genes (just add noise, mostly zero, and speeds up computation)
    if(bottomFractionTrim != 0){
      # trimValue = quantile(avgGeneExpression, bottomFractionTrim)
      # avgGeneExpression = avgGeneExpression[which(avgGeneExpression>trimValue)]
      #try num genes to keep
      avgGeneExpression = avgGeneExpression[1:bottomFractionTrim]
    }
    #Trim top genes (the top few genes are extremely highly expressed and can skew the results)
    if(topGenesTrim!=0){
      avgGeneExpression = avgGeneExpression[-seq(1,topGenesTrim)]
    }
    #Update our matrix
    cellTypeSubset@assays$RNA@data = cellTypeSubset@assays$RNA@data[names(avgGeneExpression),]
    
    #See if we should use Z-scores
    if(zScore==T){
      cellTypeSubset = ScaleData(cellTypeSubset)
    }
    
    group1Length = length(which(cellTypeSubset@active.ident==group1))
    group2length = length(which(cellTypeSubset@active.ident==group2))
    
    #Only continue if we have at least 3 cells of the cell type in each group
    if(group1Length > 2 & group2length > 2){
      
      #Label group1 and group2 cells
      group1Cells = WhichCells(cellTypeSubset,idents = group1)
      group2Cells = WhichCells(cellTypeSubset,idents = group2)
      #Get the average expression of group1 and group2 cells
      #See if we are using zscore
      if(zScore==T){
        avgExp = suppressWarnings(AverageExpression(cellTypeSubset,use.scale = T,assay = "RNA"))
        avgExp = avgExp$RNA
      } else{
        #compute in non-log space
        avgExp = suppressWarnings(AverageExpression(cellTypeSubset,assay = "RNA"))
        avgExp = avgExp$RNA
      }
      
      #calculate the distance between the group1 and group2 cells
      cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
      
      cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
      
      #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
      allCells = c(group1Cells,group2Cells)
      
      #set up a progress bar
      pb <- txtProgressBar(max = numPerm, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
        set.seed(j)
        permSubset = cellTypeSubset
        sampleCells = sample(allCells,length(allCells),replace = FALSE)
        #place the cells into either group1 or group2 (same size as original group1 and group2)
        group1Synth = sampleCells[1:length(group1Cells)]
        group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
        
        permStatus = permSubset@active.ident
        permStatus[group1Synth] = group1
        permStatus[group2Synth] = group2
        permSubset$permStatus = permStatus
        permSubset = SetIdent(permSubset, value = "permStatus")
        
        #See if we are using zscore
        if(zScore==T){
          permAvgExp = suppressWarnings(AverageExpression(permSubset,use.scale = T,assay = "RNA"))
          permAvgExp = permAvgExp$RNA
        } else{
          #compute in non-log space
          permAvgExp = suppressWarnings(AverageExpression(permSubset,assay = "RNA"))
          permAvgExp = permAvgExp$RNA
        }
        
        # #compute in non-log space
        # permAvgExp = suppressWarnings(AverageExpression(permSubset))
        # permAvgExp = permAvgExp$RNA
        
        permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
        permCellTypeDistance
      }
      cellTypePermutationDistances[[cellType]] = permutationDistances
      close(pb)
    }
  }
  
  #stop parallel processing
  stopCluster(cl)
  
  names(cellTypeDistances) = levels(seuratObject@active.ident)
  cellTypeDistancePvals = c()
  for(cellType in levels(seuratObject@active.ident)){
    numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
    cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
  }
  names(cellTypeDistancePvals) = levels(seuratObject@active.ident)
  
  melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
  melted.cellTypePvals = rep(cellTypeDistancePvals,each=numPerm)
  melted.cellTypeDistances = rep(cellTypeDistances,each=numPerm)
  melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)
  
  #########
  #Bonferroni correction
  allDistancePvals = c(cellTypeDistancePvals)
  # all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
  allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm
  
  #########
  #Combine neuronal subtype and cell types to plot
  colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
  melted.combinedPermutationDistances = melted.cellTypePermutationDistances
  melted.combinedPermutationDistances$pvals = allDistancePvals
  melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))
  
  #########
  #Plot
  ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
  p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
    geom_density(adjust = 1) +
    facet_wrap( ~ L1, ncol=5,scales="free") +
    xlab("Euclidian Distance") +
    ylab("Density") +
    # Add line
    geom_vline(aes(xintercept=distances),
               color="blue", linetype="dashed", size=1) +
    annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)
  
  pdf(file = paste0("./Plots/EuclideanDistance/",distanceScorePlot),height = 16,width = 16)
  print(p1)
  dev.off()
  
  #########
  #get the log FC of the euclidean distance versus the median of the background distribution
  
  medianBackground = c()
  for(cellType in levels(seuratObject@active.ident)){
    medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
  }
  names(medianBackground) = levels(seuratObject@active.ident)
  
  logFC = log10(cellTypeDistances) - log10(medianBackground)
  pval = -log10(allDistancePvals)
  toPlot = as.data.frame(cbind(logFC,pval))
  
  p2 <- ggplot(toPlot, aes(pval, logFC)) +
    geom_point(color = 'black') +
    theme_classic(base_size = 10) +
    xlab("-log10 pval") + 
    geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 
  
  pdf(file = paste0("./Plots/EuclideanDistance/",fcPlot))
  print(p2)
  dev.off()
  
  #Save incase we need to alter the figures, don't want to have to rerun everything
  save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Saves/",distanceSave))
}

########################################################
#Wrapper for the main DEG function -- allows us to add an additional meta data level (e.g. timepoint + condition rather than just condition)
# seuratObject=dropEST.combined.filtered
# metaData1="timpoint"
# cellMetaDataLabel="NeuronCellType"
# metaDataClassify="condition"
# numGenesUse=1000
# numBootStrap = 5
# do.par=TRUE
# num.cores="auto"
# percentTrain=0.7
# tissue="Hip"
# kNumber=10
# kRepeats=3

runClassifier <- function(seuratObject=dropEST.combined.filtered,metaData1=NULL,cellMetaDataLabel="CellTypeRefined",
                          metaDataClassify=NULL,numGenesUse=1000,numBootStrap=1000,do.par=TRUE,num.cores="auto",
                          percentTrain=0.7,kNumber=10,kRepeats=3,tissue=NULL){
  
  #initialize a list to hold the accuracy and variable importance for all cell types
  allCellTypeAccuracy = list()
  # cellTypeVariableImportance = list()
  
  #iterate through the first supplied level of meta data
  if(!(is.null(metaData1))){
    seuratObject = SetIdent(seuratObject, value = metaData1)
    for(metaObject1 in levels(seuratObject@active.ident)){
      print(metaObject1)
      metaObject1Subset = subset(seuratObject,idents = metaObject1)
      metaObject1Subset = SetIdent(metaObject1Subset, value = cellMetaDataLabel)
      
      cellTypeAccuracy = runClassifierMain(seuratObject=metaObject1Subset,cellMetaDataLabel=cellMetaDataLabel,metaDataClassify=metaDataClassify,
                                     numGenesUse=numGenesUse,numBootStrap=numBootStrap,do.par=do.par,num.cores=num.cores,percentTrain=percentTrain,
                                     kNumber=kNumber,kRepeats=kRepeats)
      
      allCellTypeAccuracy[[metaObject1]] = cellTypeAccuracy
    }
  } else{
    metaObject1Subset = SetIdent(seuratObject, value = cellMetaDataLabel)
    cellTypeAccuracy = runClassifierMain(seuratObject=metaObject1Subset,cellMetaDataLabel=cellMetaDataLabel,metaDataClassify=metaDataClassify,
                                         numGenesUse=numGenesUse,numBootStrap=numBootStrap,do.par=do.par,num.cores=num.cores,percentTrain=percentTrain,
                                         kNumber=kNumber,kRepeats=kRepeats)
    allCellTypeAccuracy = cellTypeAccuracy
  }
  #save classifier accuracy
  saveRDS(allCellTypeAccuracy, file = paste0("./Saves/",tissue,".ClassifierAccuracy.rds"))
}

########################################################
runClassifierMain <- function(seuratObject=metaObject1Subset,cellMetaDataLabel="CellTypeRefined",metaDataClassify=NULL,numGenesUse=1000,
                              numBootStrap=1000,do.par=TRUE,num.cores="auto",percentTrain=0.7,kNumber=10,kRepeats=3){
  
  metaObject1Subset = seuratObject
  
  #list to hold accuracy for each cell type
  cellTypeAccuracy = list()
  
  #Initialize list that has classifier accuracy of each cell type
  
  #set up parallel processing
  if(do.par==TRUE){
    #auto detect number of available cores
    if(num.cores == "auto"){
      cores=detectCores()
      #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
      cores = max(1,cores[1]-1)
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    } else{
      #if a number of cores for parallel processing is specified, use this number of cores
      if(is.numeric(num.cores)){
        cores = num.cores
        #setup and start parallel processing cluster
        cl <- makeCluster(cores)
        registerDoSNOW(cl)
      }
    }
  }
  #if we are not using parallel processing, just use one core
  if(do.par==FALSE){
    cores = 1
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  
  #iterate through each cell type
  for(cellType in levels(metaObject1Subset@active.ident)){
    print(cellType)
    #generate a cell type subset
    cellTypeSubset = subset(metaObject1Subset,idents = cellType)
    
    #get the data labels which will be used for the classifier
    classify = cellTypeSubset@meta.data[,metaDataClassify]
    
    #set up a gene expression matrix with the identity labels for use in svm classifier
    cellTypeSubsetData = data.frame(t(data.matrix(cellTypeSubset@assays$RNA@data)))
    #subset to top N most highly expressed genes in that cell type
    avgGeneExpression = colMeans(cellTypeSubsetData)
    avgGeneExpression = sort(avgGeneExpression,decreasing = T)
    #Keep the top N most highly expressed genes
    cellTypeSubsetData = cellTypeSubsetData[,names(avgGeneExpression)[1:numGenesUse]]
    cellTypeSubsetData$classify = classify
    
    #set up a progress bar
    pb <- txtProgressBar(max = numBootStrap, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    classifierAccuracy <- foreach(j=1:numBootStrap, .combine=c, .options.snow=opts, .packages='caret') %dopar% {
      set.seed(j)
      
      #Partition data into training and testing
      intrain <- createDataPartition(y = cellTypeSubsetData$classify, p= percentTrain, list = FALSE)
      training <- cellTypeSubsetData[intrain,]
      testing <- cellTypeSubsetData[-intrain,]
      training$classify = factor(training$classify)
      testing$classify = factor(testing$classify)
      
      #K-fold cross-validation (10-fold with 3 repeats)
      trctrl <- trainControl(method = "repeatedcv", number = kNumber, repeats = kRepeats)

      svm_Linear <- train(classify ~., data = training, method = "svmLinear",
                          trControl=trctrl,
                          preProcess = c("center", "scale"),
                          tuneLength = 10)
      
      test_pred <- predict(svm_Linear, newdata = testing)
      output = confusionMatrix(test_pred, testing$classify)
      output$overall["Accuracy"]
    }
    cellTypeAccuracy[[cellType]] = classifierAccuracy
    close(pb)
    
  }
  #stop parallel processing
  stopCluster(cl)
  
  return(cellTypeAccuracy)
}
