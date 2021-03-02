library(Seurat)
library(dplyr)
library(cowplot)
library(scExtras)
library(readxl)

setwd("~/Desktop/sengupta/")
cpallette =cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7",
                       "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                       "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                       "#8A7C64", "#599861")
####################### Global Vars ############################
outdir<-'2975371125_8347_4_WT/Seurat' # All results will be saved here plots, RData file
projectname<-'WT' # specify project name,this will also be the Rdata file name
input10x <- c('2975371125_8347_4_WT/filtered_feature_bc_matrix/') # dir(s) of the 10x output files, genes.tsv,barcodes.tsv
org<-'mouse'

mouseorthologfile <- '~/Desktop/Apoorva/NGSshare/homologs/mouse_human.csv'
npcs<-30 #How many inital PC dimensions to compute.
k=30 #This for nearest neighbors, 30 is default

################################################################
## Preprocess the data and create a seurat object
dir.create(outdir,recursive = T)
scrna= RunQC(dir=outdir,org=org,name=projectname,files=input10x ,filter=T, doubletdetection = T)
scrna = processExper(scrna ,ccscale = T, sc.transform = T)

scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = outdir)

png(file=paste0(plotdir,"/elbowplot.png",sep=""), height = 15, width = 15, units = "in", res=500)
ElbowPlot(scrna, ndims = npcs)
dev.off()

png(file=paste0(outdir,"/jackstraw.png",sep=""), height = 15, width = 15, units = "in", res=500)
JackStrawPlot(scrna, dims = 1:npcs)
dev.off()

npcs =30
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =k)
scrna = RunDoubletfinder(scrna, dims=npcs,group="var_cluster",doublet.formrate =0.75, sct=T)

scrna <- RunLigRec(scrna,org=org, group.by = "var_cluster")

saveRDS(scrna,file=paste0(outdir,"/",projectname,'.RDS'))

