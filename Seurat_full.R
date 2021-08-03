library(Seurat)
library(dplyr)
library(cowplot)
library(scExtras)
library(readxl)

#Set working directory to project directory
setwd("~/Desktop/")

####################### Set parameters ############################
# All results will be saved in this directory including plots and final seurat object
outdir<-'2975371125_8347_4_WT/Seurat'

# specify project name. This will also be the RDS file name
projectname<-'WT'

# dir of the 10x output files, genes.tsv,barcodes.tsv and matrix.mtx. If there are reps, comma seperate the files
files <- c('2975371125_8347_4_WT/filtered_feature_bc_matrix/')

#provide Organism (only works for mouse and human)
org<-'mouse'

#Remind me to give you this file. Changes mouse gene names to human
mouseorthologfile <- '~/Desktop/Apoorva/NGSshare/homologs/mouse_human.csv'

#How many inital Principle Component dimensions to compute.
npcs<-40

#This for nearest neighbors, 30 is default
k=30

#Set cut-offs for filteration
LowerFeatureCutoff=200
UpperFeatureCutoff="MAD"
UpperMitoCutoff=5

#Specify if you want to run doublet detection and scale cell cycle genes (T/F)
doubletdetection = T
ccscale = T
######################## QC DATA########################################
## Preprocess the data and create a seurat object

#Create the output directory
dir.create(outdir,recursive = T)

#If only one sample, read file and change column names to project name
inputdata <- Read10X(data.dir =files[1])

#Use this command if your input file is a h5 file instead of the matrix, barcode and feature files
#inputdata <- Read10X_h5(filename =files[1])
colnames(inputdata) <- paste0(colnames(inputdata), '-',projectname)

# Initialize the Seurat object with the raw (non-normalized data).
#Remove genes expressed in less than 200 cells and remove cells expressing less than 10 genes
scrna <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200,project = projectname)

####################################if you have reps####################################
#If you have reps, initialize the first object with the raw (non-normalized data) and add rest of the data
inputdata <- Read10X(data.dir =files[1])
#inputdata <- Read10X_h5(filename =files[1])
colnames(inputdata) <- paste0(colnames(inputdata), '-',name, '-rep1')
scrna <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200, project = name)
   for(i in 2:length(files)){
      tmp.data <- Read10X(data.dir =files[1])
      colnames(tmp.data) <- paste0(colnames(tmp.data), '-',name, '-rep',i)
      tmp.object <- CreateSeuratObject(counts= tmp.data, min.cells = 10, min.features = 200, project = name)
      scrna <- merge(scrna, tmp.object, do.normalize = FALSE, min.cells = 0, min.features = 0)
    }
################################ DOUBLET DETECTION ###################################################
  #Doublet detection (this section is optional to run)
if(doubletdetection ==T){
  sce= as.SingleCellExperiment(scrna)
  sce = cxds(sce)
  sce = bcds(sce)
  sce = cxds_bcds_hybrid(sce)
  scrna=as.Seurat(sce)
}
######################################## FILTER DATA ##################################################
  #Filter mitochondrial genes
  if(org=='mouse'){
    scrna[["percent.mito"]] <- PercentageFeatureSet(scrna, pattern = "^mt-")
  }else{
    scrna[["percent.mito"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
  }
  write.csv(object@meta.data, file=paste(dir,"/percentmito.csv",sep=""))

  #plot qc violin plot
  png(paste(dir,"/QC_Vlnplot.png",sep=""), width=10, height=6, units="in", res=300)
  print({VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
                 ncol = 3)})
  dev.off()

# save filter stats to seurat object
  scrna@misc[["filterstats"]] <- list()
  scrna@misc[["filterstats"]][['TotalCellsbeforefilteration']] <- dim(scrna)[2]

  if(UpperFeatureCutoff=="MAD"){
      UpperFeatureCutoff <- median(scrna$nFeature_RNA) + 3*mad(scrna$nFeature_RNA)
    }

    scrna@misc[["filterstats"]][['TotalSamples']] <- dim(scrna[[]][1])[1]

    cells.use <- colnames(scrna)[which(scrna[[]]['percent.mito'] < UpperMitoCutoff)]
    scrna@misc[["filterstats"]][['Mitofilter']] <- dim(scrna[[]][1])[1] - length(cells.use)
    scrna <- subset(scrna, cells = cells.use)

    cells.use <- colnames(scrna)[which(scrna[[]]['nFeature_RNA'] > LowerFeatureCutoff)]
    scrna@misc[["filterstats"]][['LowFeatureFilter']] <- dim(scrna[[]][1])[1] - length(cells.use)
    scrna <- subset(scrna, cells = cells.use)

    cells.use <- colnames(scrna)[which(scrna[[]]['nFeature_RNA'] < UpperFeatureCutoff)]
    scrna@misc[["filterstats"]][['HighFeatureFilter']] <- dim(scrna[[]][1])[1] - length(cells.use)
    scrna <- subset(scrna, cells = cells.use)

################################### PREPROCESS DATA #######################################################
return_var_genes = F
if(ccscale==T){
    scrna <- NormalizeData(scrna = scrna)
    #detection of variable genes
    #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
    scrna <-FindVariableFeatures(scrna = scrna,selection.method = "vst", nfeatures = 3000, verbose = FALSE)

    if(org=='human'){
      #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in scrna@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
      scrna <- CellCycleScoring(scrna = scrna, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }else{
      m2h <- readr::read_csv(mouseorthologfile)
      cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
      cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
      scrna <- CellCycleScoring(scrna = scrna, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }
      scrna <- SCTransform(scrna, verbose = FALSE, vars.to.regress=c("nCount_RNA","percent.mito","S.Score", "G2M.Score"),return.only.var.genes = return_var_genes)
    }else{
      scrna <- SCTransform(scrna, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = return_var_genes, verbose = FALSE)
    }
scrna <- LogSeuratCommand(scrna = scrna)

################################## RUN PCA ########################################################

PCATools <- function(object,npcs=npcs,jackstraw=T,plotdir=plotdir){
  object <- RunPCA(object = object, npcs = npcs, verbose = FALSE)
  if(jackstraw==TRUE){
    object <- JackStraw(object = object, num.replicate = 100,dims = npcs)
    object <- ScoreJackStraw(object = object,dims=1:npcs)
  }
  object <- LogSeuratCommand(object = object)
  object
}
scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = outdir)

png(file=paste0(plotdir,"/elbowplot.png",sep=""), height = 15, width = 15, units = "in", res=500)
ElbowPlot(scrna, ndims = npcs)
dev.off()

png(file=paste0(outdir,"/jackstraw.png",sep=""), height = 15, width = 15, units = "in", res=500)
JackStrawPlot(scrna, dims = 1:npcs)
dev.off()

########################### TSNE, UMAP AND OTHER DIMENSION REDUCTION #####################################
#Based on above plots, select dimensions
npcs=30
ClusterDR <-function(object,dims,n.neighbors=30,k=20,algorithm=1,DM=F,UMAP=T,TSNE=T,findallmarkers=T,resolution=0.5,n.components=2,min.dist=0.3){

  if(TSNE==TRUE){
    object <- RunTSNE(object = object, reduction = "pca",dims = dims)
  }
  if(UMAP==TRUE){
    object <- RunUMAP(object = object, reduction = "pca", n.neighbors = n.neighbors,n.components = n.components,dims = dims,min.dist=min.dist)
  }
  if(DM==TRUE){
    object <- RunDiffusion(object = object,dims=dims)
  }
  object <- FindNeighbors(object = object,dims=dims,k.param = k)
  object <- FindClusters(object = object,resolution=resolution)
  object$var_cluster <- object@active.ident
  if(findallmarkers==TRUE){
    object@misc[["findallmarkers"]] <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }

  object
}
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =k)
################################# SAVE DATASET #########################################################
saveRDS(scrna,file=paste0(outdir,"/",projectname,'.RDS'))

