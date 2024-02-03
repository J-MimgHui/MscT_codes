
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
library(patchwork)
library(harmony)
library(umap)
library(Biobase)
library(monocle)

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

################Create_cds_file######################

pbmc <- CreateSeuratObject(MscT_basfia@assays$RNA@counts, subject = 'CreateSeuratObject', min.cells = 3, min.features = 150)
pbmc$clusters_0.8 <- MscT_basfia$clusters_0.8[which(  rownames(MscT_basfia@meta.data)%in%
                                                      rownames(pbmc@meta.data)   )]

data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily = negbinomial.size())

###################run_monocle###################

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

saveRDS(monocle_cds,'./monocle_cds_bas.rds')
monocle_cds <- readRDS('./monocle_cds_bas.rds')

monocle_cds <- detectGenes(monocle_cds , min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds)))
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~clusters_0.8")
deg <- subset(diff_test_res , qval < 0.01)
ordering_genes <- row.names(deg)[order(deg$qval)]
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds)
write.csv(diff_test_res, './diff_test_res.csv')

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, reduction_method = 'DDRTree', residualModelFormulaStr = "~orig.ident")
monocle_cds <- orderCells(monocle_cds)

##############plot_cell_trajectories################

plot_cell_trajectory(monocle_cds, color_by = "clusters_0.8", cell_size = 1.5)
plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1.5)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")

time_diff <- differentialGeneTest(monocle_cds[ordering_genes,],cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
time_genes <- time_diff %>% pull(gene_short_name) %>% as.character()
time_genes <- as.data.frame(time_genes)

write.csv(time_genes, './time_genes.csv')
time_genes <- read.csv('./time_genes.csv')

###############plot_differential_genes###################

P_1 <- plot_pseudotime_heatmap(Msc_hope_monocle[time_genes$time_genes,],
                               num_clusters = 3,
                               cores = 1,
                               show_rownames = FALSE,
                               return_heatmap = T)

clusters <- cutree(P_1$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Gene <- rownames(clustering)
clustering$Gene_Clusters <- paste('state', clustering$Gene_Clusters)
write.csv(clustering , './Pseudotime_Genes.csv')

