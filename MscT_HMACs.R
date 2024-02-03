
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
library(DoubletFinder)
library(KernSmooth)

library(igraph)
library(Hmisc)
library(ggsignif)
library(clustree)
library(cowplot)
library(paletteer)

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

fashion_color_2 <- c('#CC9999','#009966','#666699','#CC6600','#999999',
                     '#CCCC33','#66CCCC','#FFCC99','#CC0066','#0099CC',
                     '#009999','#CCCC99','#FF6666','#006699','#CC6633')

fashion_color_3 <- paletteer_d("palettesForR::Tango")
fashion_color_4 <- paletteer_d("ggthemes::Classic_20",n=20)

cb_palette <- c("#246b93", "#2baeb5", "#ff66fc", "#e07233", "#d561dd", "#c93f00", "#ddd53e",
                "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#8249aa", "#cc8e12", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                "#1a918f", "#ed1299", "#2927c4", "#7149af" ,"#57e559", "#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                "#dd27ce", "#07a301", "#167275", "#391c82", "#09f9f5", "#925bea", "#63ff4f")

########################MscT_subset_HMACs##########################

MscT_HMACs_temp <- subset(Msc_hope, clusters_0.3=='HSP90+ HMACs')
MscT_HMACs <- CreateSeuratObject(MscT_HMACs_temp@assays$RNA@counts,
                                    object = 'CreateSeuratObject', min.cells = 1)

set.seed(1234)
MscT_HMACs <- NormalizeData(MscT_HMACs,normalization.method = "LogNormalize", scale.factor = 10000)
MscT_HMACs <- FindVariableFeatures(object = MscT_HMACs,selection.method = "vst",nfeatures = 2000)
MscT_HMACs <- ScaleData(object = MscT_HMACs)

MscT_HMACs <- RunPCA(object = MscT_HMACs)
ElbowPlot(MscT_HMACs, ndims=50, reduction="pca")
MscT_HMACs <- RunHarmony(MscT_HMACs , group.by.vars="orig.ident" , plot_convergence = FALSE)
MscT_HMACs <- RunUMAP(MscT_HMACs, reduction = "harmony",dims = 1:20)
MscT_HMACs <- FindNeighbors(MscT_HMACs, reduction = "harmony",dims = 1:20)

MscT_HMACs_test <- MscT_HMACs
seq <- seq(0.1, 0.3, by=0.02)
for (i in seq) {
  MscT_HMACs_test <- FindClusters(MscT_HMACs_test, resolution = i)
}

clustree(MscT_HMACs_test, prefix = 'RNA_snn_res.') 
saveRDS(MscT_HMACs_test, './MscT_HMACs_test.rds')

#Optimal resolution value: 0.14

MscT_HMACs <- FindClusters(MscT_HMACs, resolution = 0.14)
colnames(MscT_HMACs@meta.data)

DimPlot(MscT_HMACs, reduction = "umap", group.by = "seurat_clusters" , cols = fashion_color_1 , pt.size=0.5, label = T , label.size = 4 , raster=FALSE)+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))

###################MscT_HMACs_Findmarkers###############

Ann_cog <- read.table('./list_cog.tsv', header = T , sep = "\t", quote = '')
Ann_go <- read.table('./list_go.tsv', header = T , sep = "\t", quote = '')
Ann_kegg <- read.table('./list_kegg.tsv', header = T , sep = "\t", quote = '')
Ann_cazy <- read.table('./list_cazy.tsv', header = T , sep = "\t", quote = '')

Ann_MscT <- merge(Ann_cog, Ann_go, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_kegg, by='cowSGB_Genes', all=T)
Ann_MscT <- merge(Ann_MscT, Ann_cazy, by='cowSGB_Genes', all=T)

Ann_MscT$cowSGB_Genes <- gsub(pattern = '_', replacement = '-', Ann_MscT$cowSGB_Genes)
Ann_MscT$cowSGB_Genes <- gsub('.{2}$', '', Ann_MscT$cowSGB_Genes)

Idents(MscT_HMACs) <- 'seurat_clusters'
MscT_HMACs_markers <- FindAllMarkers(MscT_HMACs, only.pos = TRUE, logfc.threshold = 0.25 , min.pct = 0.1)
MscT_HMACs_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
MscT_HMACs_markers <- MscT_HMACs_markers[,c(6,7,2,3,4)]
colnames(MscT_HMACs_markers)[2] <- 'cowSGB_Genes'

MscT_HMACs_markers_ann <- merge(MscT_HMACs_markers, Ann_Msc_hope, by='cowSGB_Genes', all = F)
write.csv(MscT_HMACs_markers, './Finaldata/MscT_HMACs_markers.csv')
write.csv(MscT_HMACs_markers_ann, './Finaldata/MscT_HMACs_markers_ann.csv')
saveRDS(MscT_HMACs, './Finaldata/MscT_HMACs.rds')

################MscT_HMACs_identify_clusters#################

MscT_HMACs@meta.data$clusters_0.14 <- as.character(MscT_HMACs@meta.data$seurat_clusters)
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('0','4','5')] <- 'HTH+ HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('1')] <- 'Integrase+ HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('2','7','8')] <- 'Peptide transporter HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('3','11')] <- 'Motility HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('6')] <- 'Membrane protein+ HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('9')] <- 'Lipid metabolism HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('10')] <- 'Secretion HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('12')] <- 'Multi-signal HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('13')] <- 'TonB-linked protein+ HMACs'
MscT_HMACs@meta.data$clusters_0.14[MscT_HMACs@meta.data$clusters_0.14%in%c('14')] <- 'His Kinase A+ HMACs'

DimPlot(MscT_HMACs, reduction = "umap", group.by = "clusters_0.14" , cols = my36colors , pt.size=0.5, label = F , label.size = 4 , raster=FALSE)+
  theme(panel.border = element_rect(color = "black",fill = NA,size = 1))
saveRDS(MscT_HMACs, './MscT_HMACs.rds')

