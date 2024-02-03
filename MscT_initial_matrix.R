
library(rstatix)
library(ggsci)
library(dplyr)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(stats)
library(mgcv)
library(tidyr)
library(readxl)
library(paletteer)

library(Seurat)
library(scuttle)
library(patchwork)
library(umap)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(ggprism)
library(ggpubr)
library(pheatmap)
library(vegan)
library(DoubletFinder)
library(KernSmooth)

library(igraph)
library(Hmisc)
library(ggsignif)
library(clustree)

library(viridis)
library(SCopeLoomR)
library(Matrix)

##################color########################

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', '#D6E7A3')

fashion_color_1 <- c('#6699CC','#663366','#CCCC99','#990033','#006633',
                     '#CCCC00','#CC9900','#FF9933','#666699','#660033',
                     '#99CC99','#993366','#666633','#FF0033','#FFFF66')

fashion_color_2 <- c('#CC9999','#009966','#666699','#CC6600','#999999',
                     '#CCCC33','#66CCCC','#FFCC99','#CC0066','#0099CC',
                     '#009999','#CCCC99','#FF6666','#006699','#CC6633')

fashion_color_3 <- paletteer_d("palettesForR::Tango")
fashion_color_4 <- paletteer_d("ggthemes::Classic_20",n=20)

color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
           '#f032e6', '#bcf60c', '#e6beff', 
           '#9a6324', '#fabebe', '#008080', '#808000', '#fffac8', '#800000', '#aaffc3',
           '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

########################MscT_benchmarking######################

MscT_initial <- readRDS('./MscT_primary.rds')
MscT_initial_test <- MscT_initial

set.seed(1234)
MscT_initial_test <- NormalizeData(MscT_initial_test,normalization.method = "LogNormalize", scale.factor = 10000)
MscT_initial_test <- FindVariableFeatures(object = MscT_initial_test,selection.method = "vst",nfeatures = 2000)
MscT_initial_test <- ScaleData(object = MscT_initial_test)

MscT_initial_test <- RunPCA(object = MscT_initial_test)
ElbowPlot(MscT_initial_test, ndims=50, reduction="pca")
MscT_initial_test <- RunHarmony(MscT_initial_test , group.by.vars="orig.ident" , plot_convergence = FALSE)
MscT_initial_test <- RunUMAP(MscT_initial_test, reduction = "harmony",dims = 1:15)
MscT_initial_test <- FindNeighbors(MscT_initial_test, reduction = "harmony",dims = 1:15)

seq <- seq(0.1, 0.9, by=0.1)
for (i in seq) {
  MscT_initial_test <- FindClusters(MscT_initial_test, resolution = i)
}

clustree(MscT_initial_test, prefix = 'RNA_snn_res.') 
saveRDS(MscT_initial_test, './MscT_initial_test.rds')

##############MscT_Optimal_resolution_value_0.3#############

MscT_initial <- readRDS('./MscT_initial.rds')

set.seed(1234)
MscT_initial <- NormalizeData(MscT_initial,normalization.method = "LogNormalize", scale.factor = 10000)
MscT_initial <- FindVariableFeatures(object = MscT_initial,selection.method = "vst",nfeatures = 2000)
MscT_initial <- ScaleData(object = MscT_initial)

MscT_initial <- RunPCA(object = MscT_initial)
ElbowPlot(MscT_initial, ndims=50, reduction="pca")
MscT_initial <- RunHarmony(MscT_initial , group.by.vars="orig.ident" , plot_convergence = FALSE)
MscT_initial <- RunUMAP(MscT_initial, reduction = "harmony",dims = 1:15)
MscT_initial <- FindNeighbors(MscT_initial, reduction = "harmony",dims = 1:15)
MscT_initial <- FindClusters(MscT_initial, resolution = 0.3)

###################MscT_Findmarkers###############

Ann_cog <- read.table('./list_cog.tsv', header = T , sep = "\t", quote = '')
Ann_go <- read.table('./list_go.tsv', header = T , sep = "\t", quote = '')
Ann_kegg <- read.table('./list_kegg.tsv', header = T , sep = "\t", quote = '')
Ann_cazy <- read.table('./list_cazy.tsv', header = T , sep = "\t", quote = '')

Ann_MscT <- merge(Ann_cog, Ann_go, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_kegg, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_cazy, by='cowSGB_Genes', all=T)

Ann_MscT$cowSGB_Genes <- gsub(pattern = '_', replacement = '-', Ann_MscT$cowSGB_Genes)
Ann_MscT$cowSGB_Genes <- gsub('.{2}$', '', Ann_MscT$cowSGB_Genes)

Idents(MscT_initial) <- 'seurat_clusters'
MscT_initial_markers <- FindAllMarkers(MscT_initial, only.pos = TRUE, logfc.threshold = 0.25 , min.pct = 0.1)
MscT_initial_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
MscT_initial_markers <- MscT_initial_markers[,c(6,7,2,3,4)]
colnames(MscT_initial_markers)[2] <- 'cowSGB_Genes'

MscT_initial_markers_ann <- merge(Msc_hope_markers, Ann_MscT, by='cowSGB_Genes', all = F)
write.csv(MscT_initial_markers, './MscT_initial_markers.csv')
write.csv(MscT_initial_markers_ann, './MscT_initial_markers_ann.csv')

########################MscT_identify_cluster####################

MscT_initial@meta.data$clusters_0.3 <- as.character(MscT_initial@meta.data$seurat_clusters)
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('0','22','23','24')] <- 'HSP90+ HMACs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('1','5','7','10','13','18','20','21','25')] <- 'Integrase/amylase+ MCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('2')] <- 'Formate/nitrite TCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('3','15','16')] <- 'RPA+ lipid MCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('4')] <- 'GTPase+ proliferating cells'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('6')] <- 'Pyruvate/glucan MCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('8')] <- 'Signal communication cells'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('9')] <- 'Lipid/ion MCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('11')] <- 'Unknown cells'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('12','14')] <- 'ARH+ cells'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('17')] <- 'Sulfur MCs'
MscT_initial@meta.data$clusters_0.3[MscT_initial@meta.data$clusters_0.3%in%c('19')] <- 'Transposase+ cells'

DimPlot(MscT_initial, reduction = "umap", group.by = "clusters_0.3" , cols = my36colors , pt.size=0.5, label = F , label.size = 4 , raster=FALSE)+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))
saveRDS(MscT_initial, './MscT_initial.rds')

