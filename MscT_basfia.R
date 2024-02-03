
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
library(KernSmooth)

library(igraph)
library(Hmisc)
library(ggsignif)
library(clustree)

##################color########################

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175', '#D6E7A3')

fashion_color_1 <- c('#6699CC','#663366','#CCCC99','#990033','#006633',
                     '#CCCC00','#CC9900','#FF0033','#666699','#660033',
                     '#99CC99','#993366','#666633','#FF9933','#FFFF66')

fashion_color_2 <- c('#999999','#CC9999','#009966','#CC6600','#666699','#CCCC33',
                     '#CC0066','#66CCCC','#FFCC99','#0099CC',
                     '#009999','#CCCC99','#FF6666','#006699','#CC6633')


fashion_color_3 <- paletteer_d("palettesForR::Tango")
fashion_color_4 <- paletteer_d("ggthemes::Classic_20",n=20)

################MscT_subset_basfia######################

MscT_initial <- readRDS('./MscT_initial.rds')
MscT_basfia_temp <- subset(MscT_initial, Species %in% 's__Basfia succiniciproducens')
MscT_basfia <- CreateSeuratObject(MscT_basfia_temp@assays$RNA@counts, object='CreateSeuratObject', min.cells = 1)

set.seed(1234)
MscT_basfia <- NormalizeData(MscT_basfia,normalization.method = "LogNormalize", scale.factor = 10000)
MscT_basfia <- FindVariableFeatures(object = MscT_basfia,selection.method = "vst",nfeatures = 2000)
MscT_basfia <- ScaleData(object = MscT_basfia)

MscT_basfia <- RunPCA(object = MscT_basfia)
ElbowPlot(MscT_basfia, ndims=50, reduction="pca")
MscT_basfia <- RunHarmony(MscT_basfia , group.by.vars="orig.ident" , plot_convergence = FALSE)
MscT_basfia <- RunUMAP(MscT_basfia, reduction = "harmony",dims = 1:20)
MscT_basfia <- FindNeighbors(MscT_basfia, reduction = "harmony",dims = 1:20)

seq <- seq(0.1, 1, by=0.1)
for (i in seq) {
  MscT_basfia_test <- FindClusters(MscT_basfia_test, resolution = i)
}

colnames(MscT_basfia_test@meta.data)
clustree(MscT_basfia_test, prefix = 'RNA_snn_res.') 
saveRDS(MscT_basfia_test, './Finaldata_Pre_20231219/pbmc_Msc_hope/RDS/MscT_basfia_test(0.1_1).rds')

#Optimal resolution value: 0.8

Msc_hope_basfia <- FindClusters(Msc_hope_basfia, resolution = 0.8)
DimPlot(Msc_hope_basfia, reduction = "umap", group.by = "seurat_clusters" , cols = fashion_color_2 , pt.size=0.5, label = T , label.size = 4 , raster=FALSE)+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))

Msc_hope_basfia <- readRDS('./Finaldata_Pre_20231219/pbmc_Msc_hope/RDS/Msc_hope_basfia.rds')

######################MscT_basfia_Findmarkers##################

Ann_cog <- read.table('./list_cog.tsv', header = T , sep = "\t", quote = '')
Ann_go <- read.table('./list_go.tsv', header = T , sep = "\t", quote = '')
Ann_kegg <- read.table('./list_kegg.tsv', header = T , sep = "\t", quote = '')
Ann_cazy <- read.table('./list_cazy.tsv', header = T , sep = "\t", quote = '')

Ann_MscT <- merge(Ann_cog, Ann_go, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_kegg, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_cazy, by='cowSGB_Genes', all=T)

Ann_MscT$cowSGB_Genes <- gsub(pattern = '_', replacement = '-', Ann_MscT$cowSGB_Genes)
Ann_MscT$cowSGB_Genes <- gsub('.{2}$', '', Ann_MscT$cowSGB_Genes)

Idents(MscT_basfia) <- 'seurat_clusters'
MscT_basfia_markers <- FindAllMarkers(MscT_basfia, only.pos = TRUE, logfc.threshold = 0.25 , min.pct = 0.1)
MscT_basfia_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
MscT_basfia_markers <- MscT_basfia_markers[,c(6,7,2,3,4)]
colnames(MscT_basfia_markers)[2] <- 'cowSGB_Genes'

MscT_basfia_markers_ann <- merge(MscT_basfia_markers, Ann_Msc_hope, by='cowSGB_Genes', all = F)
table(MscT_basfia_markers_ann$cluster)
write.csv(MscT_basfia_markers, './Finaldata/MscT_basfia_markers.csv')
write.csv(MscT_basfia_markers_ann, './Finaldata/MscT_basfia_markers_ann.csv')
saveRDS(MscT_basfia, './Finaldata/MscT_basfia.rds')

##################MscT_basfia_identify_clusters##################

MscT_basfia@meta.data$clusters_0.8 <- as.character(MscT_basfia@meta.data$seurat_clusters)
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('0')] <- 'Multi signal cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('1')] <- 'Integrase+ cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('2','3','7')] <- 'Domain protein+ cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('5')] <- 'Signal transduction cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('6')] <- 'Stress response cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('8','9')] <- 'Transposase+ Formate/nitrite TCs'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('10')] <- 'DNA binding protein+ cells'
MscT_basfia@meta.data$clusters_0.8[MscT_basfia@meta.data$clusters_0.8%in%c('4','11','12')] <- 'Unknown cells'

DimPlot(MscT_basfia, reduction = "umap", group.by = "clusters_0.8" , cols =fashion_color_2 , pt.size=0.5, label = F , label.size = 4 , raster=FALSE)+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))

saveRDS(MscT_basfia, './Finaldata_Pre_20231219/pbmc_Msc_hope/RDS/MscT_basfia.rds')

