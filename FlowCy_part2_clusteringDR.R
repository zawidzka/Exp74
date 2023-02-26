rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)
library(tidyverse)
library(data.table)
library(ggpubr)
library(flowAI)
library(PeacoQC)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# set workingDir
workingDir <- "230223_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
setwd(workingDirPath)


sce <- readRDS("SCE_3timepts_BMTUM_3.rds")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

MDSplot <- CATALYST::pbMDS(sce, features = type_markers(sce), fun = "median", label_by = "condition")
MDSplot

MDSplot + scale_color_manual(values = c("blue", "blue", "blue", "green", "green" ,"green")) + geom_point() + stat_ellipse()

# ggsave("pbMDS_3timepoints.pdf", height = 7, width = 7, plot = last_plot())


pdf("MDS_3timepoints.pdf", width = 10, height = 7)
CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median", label_by = "condition")
dev.off()

# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 5000
n_events <- min(n_cells(sce))
ifelse(n_cells > n_events, n_cells <- n_events, n_cells <- n_cells)
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
# sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
#              distMethod = "euclidean",
#              PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type", assay = "exprs")
# sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")
sce2 <- runDR(sce, dr = "PCA", cells = n_cells, features = "type", assay = "exprs")

saveRDS(sce, file = "SCE_3timepts_BMTUM_UMAP_3_DR.rds")

display.brewer.all(colorblindFriendly = TRUE)
delta_area(sce)
cdx2 <- type_markers(sce)
cdx <- state_markers(sce)
Heatmap_meta6 <- plotMultiHeatmap(sce, k = "meta6",
                 hm1 = cdx2, hm2 = "abundances", 
                 bars = FALSE, perc = TRUE, row_anno = TRUE, scale = "last")


heatmap_meta5 <- 
  plotMultiHeatmap(sce, k = "meta6",
                 hm1 = cdx2, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = TRUE, scale = "last")



plotExprHeatmap(sce, features = type_markers(sce), k = "meta6", by = "cluster_id",  fun = "mean", scale = "first")


UMAP1 <- plotDR(sce, dr = "UMAP", color_by = "condition") + geom_point(size = 2)
UMAP2 <- plotDR(sce, dr = "UMAP", color_by = "condition") + 
  geom_point(size = 2) + scale_color_manual(values = c("salmon", "salmon", "salmon", "darkturquoise", "darkturquoise" ,"darkturquoise")) 
PCA <- plotDR(sce, dr = "PCA", color_by = "condition") + geom_point()
plotDR(sce, dr = "UMAP", color_by = "condition")
PCA + stat_ellipse()



pdf("MDS.pdf", height = 8, width = 12, useDingbats = FALSE)
MDSplot
dev.off()

tiff("MDS.tiff", height = 8, width = 12, units = "in", res = 600)
MDSplot
dev.off()

pdf("UMAP_colorbycondition.pdf", height = 6, width = 7, useDingbats = FALSE)
UMAP1
dev.off()

tiff("UMAP_colorbycondition.tiff", height = 6, width = 7, units = "in", res = 600)
UMAP1
dev.off()

pdf("HeatmapMeta6_2.pdf", height = 8, width = 14, useDingbats = FALSE)
Heatmap_meta6
dev.off()

tiff("HeatmapMeta6_2.tiff", height = 8, width = 14, units = "in", res = 600)
Heatmap_meta6
dev.off()


Abundances <- plotAbundances(sce, k = "meta6", by = "cluster_id", group_by = "condition")
ggsave("abundances.pdf", plot = last_plot(), height = 6, width = 10)

pdf("Clusters_meta6.pdf", height = 8, width = 12, useDingbats = FALSE)
Abundances + 
  labs(title = "Frequency of Cluster Expression by Condition") +
  theme(plot.title = element_text(face = "bold", size = "14"))
dev.off()

tiff("Clusters_meta6.tiff", height = 8, width = 12, units = "in", res = 600)
Abundances +
labs(title = "Frequency of Cluster Expression by Condition") +
  theme(plot.title = element_text(face = "bold", size = "14"))
dev.off()




# plot expression of type markers by cluster
plotClusterHeatmap(sce,
                   hm2 = NULL, k = "meta10", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
plotClusterExprs(sce, k = "meta10", features = "type")

CATALYST::plotDR(sce2, dr = "UMAP", color_by = "condition") + geom_point(size = 1)
                                                                        



CATALYST::plotDR(sce, dr = "UMAP", color_by = "meta5", facet_by = "condition")
CATALYST::plotDR(sce, dr = "UMAP", color_by = c("PD1", "TIM3", "TOX"))


UMAP_singlemarkers <- CATALYST::plotDR(sce, dr = "UMAP", 
                 color_by = c("Ly108", "CX3CR1", "CD62L", "TCF1", "CD44", "CD69", "CD127", "KLRG1", 
                              "Ki67", "BCL2", "PD1", "TIM3", "TOX", "CD101"))#, 
                 facet_by = "condition") +
                 geom_point(size = 1)
plotDR(sce, dr = "UMAP", color_by = "TCF1", facet_by = "condition") +  
  geom_density2d(binwidth = 0.006, colour = "grey")
plotDR(sce2, dr = "UMAP", color_by = "meta9", facet_by = "condition", ncol = NULL) +  
  geom_density2d(binwidth = 0.016, colour = "grey") + geom_point(size = 1)



pdf("UMAP_singlemarkers_2.pdf", height = 8, width = 12, useDingbats = FALSE)
UMAP_singlemarkers
dev.off()

tiff("UMAP_singlemarkers_2.tiff", height = 8, width = 12, units = "in", res = 600)
UMAP_singlemarkers
dev.off()
                                                                                        
