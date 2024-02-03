
library(dplyr)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(stats)
library(mgcv)
library(tidyr)

library(Seurat)
library(scuttle)
library(umap)
library(harmony)
library(DoubletFinder)

##############Read_cell_matrix###############

MscT_matrix <- Read10X(data.dir = './MscT_matrix/', gene.column = 1)
MscT_primary <- CreateSeuratObject(MscT_matrix, assay = 'RNA', project = 'CreateSeuratObject')

################filtering##################

qc.genes  <- isOutlier(MscT_primary$nFeature_RNA , nmads = 3, log=TRUE, type="both")
qc.counts <- isOutlier(MscT_primary$nCount_RNA   , nmads = 3, log=TRUE, type="both")

attr(qc.genes,  "thresholds")
attr(qc.counts, "thresholds")

MscT_primary@meta.data %>% 
  ggplot(aes( x=nFeature_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")+
  geom_vline(xintercept = c(39,973))

MscT_primary@meta.data %>% 
  ggplot(aes( x=nCount_RNA)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")+
  geom_vline(xintercept = c(367,61098))

MscT_primary@meta.data$QC_1 <- 'L'
MscT_primary@meta.data$QC_1[which(MscT_primary$nCount_RNA>367&
                            MscT_primary$nCount_RNA<61098&
                            MscT_primary$nFeature_RNA>39&
                            MscT_primary$nFeature_RNA<973)] <- 'H'

MscT_primary_sub <- subset(MscT_primary , nCount_RNA   > 367 & nCount_RNA  < 61098 &
                              nFeature_RNA > 39 & nFeature_RNA < 973 )

saveRDS(MscT_primary_sub, './MscT_primary_sub.rds')
MscT_primary <- readRDS('./MscT_primary_sub.rds')

####################Remove_Doublets################

set.seed(1234)
MscT_primary <- NormalizeData(MscT_primary,normalization.method = "LogNormalize", scale.factor = 10000)
MscT_primary <- FindVariableFeatures(object = MscT_primary,selection.method = "vst",nfeatures = 2000)
MscT_primary <- ScaleData(object = MscT_primary)

MscT_primary <- RunPCA(object = MscT_primary)
ElbowPlot(MscT_primary, ndims=50, reduction="pca")
MscT_primary <- RunHarmony(MscT_primary , group.by.vars="orig.ident" , plot_convergence = FALSE)
MscT_primary <- RunUMAP(MscT_primary, reduction = "harmony",dims = 1:25)
MscT_primary <- FindNeighbors(MscT_primary, reduction = "harmony",dims = 1:25)
MscT_primary <- FindClusters(MscT_primary, resolution = 0.2)

saveRDS(MscT_primary,'./MscT_primary_sub_normalize.rds')
MscT_primary <- readRDS('./MscT_primary_sub_normalize.rds')

set.seed(1234)
sweep.res.list <- paramSweep_v3(MscT_primary, PCs = 1:25, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
DoubletRate = 0.117

MscT_primary@meta.data$seurat_clusters <- as.character(MscT_primary@meta.data$seurat_clusters)
annotations <- MscT_primary@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(DoubletRate*length(MscT_primary$seurat_clusters))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gc()
MscT_primary <- doubletFinder_v3(MscT_primary, PCs = 1:25, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
gc()

MscT_primary <- subset(MscT_primary , DF.classifications_0.25_0.3_20082%in%'Singlet' )
saveRDS(MscT_primary,'./MscT_primary_sub_normal_doublet.rds')


