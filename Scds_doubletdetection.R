library(Seurat)
library(scds)

scrna= readRDS("~/Desktop/AGIQ424_RUL_D_CD45neg_20190922.RDS")
sce= as.SingleCellExperiment(scrna)
sce = cxds(sce)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)
seurat=as.Seurat(sce)
scrna=seurat
saveRDS(scrna, file="~/Desktop/Maria_scds_test.RDS")

process<- function(scrna){
  sce= as.SingleCellExperiment(scrna)
  sce = cxds(sce)
  sce = bcds(sce,retRes = TRUE,verb=TRUE)
  sce = cxds_bcds_hybrid(sce)
  seurat=as.Seurat(sce)
  meta= seurat@meta.data
  return(meta)
  }

setwd("~/Desktop/labproject/projects/Morrisey/Maria/lungMAP/")
scrna=readRDS(file = 'R710000511_COPD_D_CD45neg_20181109/Seurat/R710000511_COPD_D_CD45neg_20181109.RDS')
sce= as.SingleCellExperiment(scrna)
sce = cxds(sce)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)
seurat=as.Seurat(sce)
R710000511_meta= seurat@meta.data

scrna=readRDS(file = 'R710000546_COPD_D_CD45neg_20190517/Seurat/R710000546_COPD_D_CD45neg_20190517.RDS')
R710000546_meta= process(scrna)


scrna=readRDS(file = 'R710000548_COPD_D_CD45neg_20190528/Seurat/R710000548_COPD_D_CD45neg_20190528.RDS')
R710000548_meta= process(scrna)

scrna=readRDS(file = 'AEL4146_RUL_D_CD45neg_20180101/Seurat/AEL4146_RUL_D_CD45neg_20180101.RDS')
AEL4146_meta= process(scrna)

scrna= readRDS(file = 'AFAE328_RML_D_CD45neg_20180117/Seurat/AFAE328_RML_D_CD45neg_20180117.RDS')
AFAE328_meta= process(scrna)

scrna= readRDS(file = 'AFEW459_LUL_D_CD45neg_20180525/Seurat/AFEW459_LUL_D_CD45neg_20180525.RDS')
AFEW459_meta= process(scrna)

scrna= readRDS(file = 'AFJZ055_RML_D_CD45neg_20181029/Seurat/AFJZ055_RML_D_CD45neg_20181029.RDS')
AFJZ055_meta= process(scrna)

scrna <- readRDS(file = 'AGIQ424_RUL_D_CD45neg_20190922/scRNA/Seurat/AGIQ424_RUL_D_CD45neg_20190922.RDS')
AGIQ424_meta = process(scrna)

meta.data=rbind(R710000511_meta,R710000546_meta,R710000548_meta,AEL4146_meta,AFAE328_meta,AFEW459_meta,AFJZ055_meta,AGIQ424_meta)
