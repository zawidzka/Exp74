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
workingDir <- "220915_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
setwd(workingDirPath)

sce <- readRDS("SCE_4timepts_BMTUM_DR.rds")
CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", 
                by = "cluster_id", scale = "first", bars = TRUE, perc = TRUE)


# annotate clusters

annotation_table <- as.data.frame(cbind(c(1:8), c(1:8)))

colnames(annotation_table) <- c("meta8", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta8", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
# filtered_sce <- filterSCE(sce, cluster_id %in% c(paste0("C", c(1:12))), k = "cluster_annotation")

# store original_sce
# old_sce <- sce
# sce <- filtered_sce

FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")

# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrastBMvBM <- createContrast(c(0, 1, 0))
contrastBMvTUM <-createContrast(c(0, 0, 1))

nrow(contrastBMvBM) == ncol(design)
nrow(contrastBMvTUM) == ncol(design)


# prelim plots

pdf("Multiheatmap_4tps_TEFF_BMTUM.pdf", width = 14, height = 6)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()

pdf("abundances_4tps_TEFF_BMTUM.pdf", width = 14, height = 6)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
dev.off()

# contrast <- createContrast(c(0,1))

# nrow(contrast) == ncol(design)

out_DA1 <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrastBMvBM,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE, min_samples = 1)

da1 <- rowData(out_DA1$res)
plotDiffHeatmap(sce, da1, top_n = 9, all = TRUE, fdr = FDR_cutoff)

out_DA2 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBMvTUM,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da2 <- rowData(out_DA2$res)
plotDiffHeatmap(sce, da2, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# create a different design matrix for contrast of other groups, reference level switch
tissue <- factor(c("BM_D7", "BM_D28", "BM_D7", "BM_D28", 
                   "BM_D7", "BM_D28", "BM_D7", "BM_D28", 
                   "BM_D28","BM_D28", "BM_D28",
                   "TUM_D28", "TUM_D28", "TUM_D28", 
                   "TUM_D28", "TUM_D28","TUM_D28", "TUM_D28"))
model.matrix(~tissue)
tissue <- relevel(tissue, ref = "BM_D28")
design2<-model.matrix(~tissue)
design2

# create another contrast
contrastBMlatevTUM <- createContrast(c(0, 0, 1))


out_DA3 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastBMlatevTUM,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da3 <- rowData(out_DA3$res)
plotDiffHeatmap(sce, da3, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# design2 <-createDesignMatrix(ei, cols_design = c("condition", "patient_id"))
# design2



save(list = ls(), file = "workspaceDA.rds")









