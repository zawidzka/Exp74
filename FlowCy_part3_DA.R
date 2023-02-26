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

sce <- readRDS("SCE_3timepts_BMTUM_UMAP_3_DR.rds")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta6", 
                by = "cluster", scale = "first", bars = TRUE, perc = TRUE)



# annotate clusters

annotation_table <- as.data.frame(cbind(c(1:6), c(1:6)))

colnames(annotation_table) <- c("meta6", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta6", 
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
# design <- createDesignMatrix(ei,
#                             cols_design = c("condition"))

sample_condition <- factor(c("BM_D14", "TUM_D14","TUM_D14","BM_D14", "BM_D14","BM_D28",
                             "TUM_D28","TUM_D28","TUM_D28", "BM_D28", "BM_D28","BM_D21",
                             "TUM_D21","TUM_D21","TUM_D21","BM_D21","BM_D21"))

model.matrix(~sample_condition)
design<-model.matrix(~sample_condition)
design

contrastBM_D14vBM_D21 <- createContrast(c(0, 1, 0, 0, 0, 0))
contrastBM_D14vBM_D28 <- createContrast(c(0, 0, 1, 0, 0, 0))
contrastBM_D14vTUM_D14 <- createContrast(c(0, 0, 0, 1, 0, 0))
contrastBM_D14vTUM_D21 <- createContrast(c(0, 0, 0, 0, 1, 0))
contrastBM_D14vTUM_D28 <- createContrast(c(0, 0, 0, 0, 0, 1))

nrow(contrastBM_D14vBM_D21) == ncol(design)
nrow(contrastBM_D14vBM_D28) == ncol(design)
nrow(contrastBM_D14vTUM_D14) == ncol(design)
nrow(contrastBM_D14vTUM_D21) == ncol(design)
nrow(contrastBM_D14vTUM_D28) == ncol(design)
# prelim plots

pdf("Multiheatmap_4tps_TEFF_BMTUM.pdf", width = 14, height = 6)
plotMultiHeatmap(sce, k = "meta6", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()

pdf("abundances_4tps_TEFF_BMTUM.pdf", width = 14, height = 6)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
dev.off()


out_DA1 <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrastBM_D14vBM_D21,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE, min_samples = 1)

da1 <- rowData(out_DA1$res)
plotDiffHeatmap(sce, da1, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA2 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBM_D14vBM_D28,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da2 <- rowData(out_DA2$res)
plotDiffHeatmap(sce, da2, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA3 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBM_D14vTUM_D14,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da3 <- rowData(out_DA3$res)
plotDiffHeatmap(sce, da3, top_n = 6, all = TRUE, fdr = FDR_cutoff)


out_DA4 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBM_D14vTUM_D21,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da4 <- rowData(out_DA4$res)
plotDiffHeatmap(sce, da4, top_n = 6, all = TRUE, fdr = FDR_cutoff)


out_DA5 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBM_D14vTUM_D28,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da5 <- rowData(out_DA5$res)
plotDiffHeatmap(sce, da5, top_n = 6, all = TRUE, fdr = FDR_cutoff)


### create a different design matrix for contrast of other groups, reference level switch
# sample_condition <- factor(c("BM_D14", "TUM_D14","TUM_D14","BM_D14", "BM_D14","BM_D28",
#                              "TUM_D28","TUM_D28","TUM_D28", "BM_D28", "BM_D28","BM_D21",
#                              "TUM_D21","TUM_D21","TUM_D21","BM_D21","BM_D21"))
#                              
#   
#   
#   factor(c("BM_D7", "BM_D28", "BM_D7", "BM_D28", 
#                    "BM_D7", "BM_D28", "BM_D7", "BM_D28", 
#                    "BM_D28","BM_D28", "BM_D28",
#                    "TUM_D28", "TUM_D28", "TUM_D28", 
#                    "TUM_D28", "TUM_D28","TUM_D28", "TUM_D28"))
sample_condition <- relevel(sample_condition, ref = "BM_D21")
design2<-model.matrix(~sample_condition)
design2

# create another contrast set
contrastBM_D21vBM_D28 <- createContrast(c(0, 0, 1, 0, 0, 0))
contrastBM_D21vTUM_D14 <- createContrast(c(0, 0, 0, 1, 0, 0))
contrastBM_D21vTUM_D21 <- createContrast(c(0, 0, 0, 0, 1, 0))
contrastBM_D21vTUM_D28 <- createContrast(c(0, 0, 0, 0, 0, 1))

nrow(contrastBM_D21vBM_D28) == ncol(design2)
nrow(contrastBM_D21vTUM_D14) == ncol(design2)
nrow(contrastBM_D21vTUM_D21) == ncol(design2)
nrow(contrastBM_D21vTUM_D28) == ncol(design2)


out_DA6 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastBM_D21vBM_D28,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da6 <- rowData(out_DA6$res)
plotDiffHeatmap(sce, da6, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA7 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastBM_D21vTUM_D14,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da7 <- rowData(out_DA7$res)
plotDiffHeatmap(sce, da7, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA8 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastBM_D21vTUM_D21,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da8 <- rowData(out_DA8$res)
plotDiffHeatmap(sce, da8, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA9 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastBM_D21vTUM_D28,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da9 <- rowData(out_DA9$res)
plotDiffHeatmap(sce, da9, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# design2 <-createDesignMatrix(ei, cols_design = c("condition", "patient_id"))
# design2


### re-level design matrix again for contrast of other groups, reference level switch
sample_condition <- relevel(sample_condition, ref = "BM_D28")
design3<-model.matrix(~sample_condition)
design3

# create another contrast set
contrastBM_D28vTUM_D14 <- createContrast(c(0, 0, 0, 1, 0, 0))
contrastBM_D28vTUM_D21 <- createContrast(c(0, 0, 0, 0, 1, 0))
contrastBM_D28vTUM_D28 <- createContrast(c(0, 0, 0, 0, 0, 1))

nrow(contrastBM_D28vTUM_D14) == ncol(design3)
nrow(contrastBM_D28vTUM_D21) == ncol(design3)
nrow(contrastBM_D28vTUM_D28) == ncol(design3)

out_DA10 <- diffcyt(sce,
                   experiment_info = ei, design = design3, contrast = contrastBM_D28vTUM_D14,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da10 <- rowData(out_DA10$res)
plotDiffHeatmap(sce, da10, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA11 <- diffcyt(sce,
                    experiment_info = ei, design = design3, contrast = contrastBM_D28vTUM_D21,
                    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                    clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                    transform = FALSE, normalize = FALSE, min_samples = 1)

da11 <- rowData(out_DA11$res)
plotDiffHeatmap(sce, da11, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA12 <- diffcyt(sce,
                    experiment_info = ei, design = design3, contrast = contrastBM_D28vTUM_D28,
                    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                    clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                    transform = FALSE, normalize = FALSE, min_samples = 1)

da12 <- rowData(out_DA12$res)
plotDiffHeatmap(sce, da12, top_n = 6, all = TRUE, fdr = FDR_cutoff)


### relevel next to last time
sample_condition <- relevel(sample_condition, ref = "TUM_D14")
design4<-model.matrix(~sample_condition)
design4

# create another contrast set
contrastTUM_D14vTUM_D21 <- createContrast(c(0, 0, 0, 0, 1, 0))
contrastTUM_D14vTUM_D28 <- createContrast(c(0, 0, 0, 0, 0, 1))

nrow(contrastTUM_D14vTUM_D21) == ncol(design4)
nrow(contrastTUM_D14vTUM_D28) == ncol(design4)

out_DA13 <- diffcyt(sce,
                    experiment_info = ei, design = design4, contrast = contrastTUM_D14vTUM_D21,
                    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                    clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                    transform = FALSE, normalize = FALSE, min_samples = 1)

da13 <- rowData(out_DA13$res)
plotDiffHeatmap(sce, da13, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA14 <- diffcyt(sce,
                    experiment_info = ei, design = design4, contrast = contrastTUM_D14vTUM_D28,
                    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                    clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                    transform = FALSE, normalize = FALSE, min_samples = 1)

da14 <- rowData(out_DA14$res)
plotDiffHeatmap(sce, da14, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# final relevel
sample_condition <- relevel(sample_condition, ref = "TUM_D21")
design5<-model.matrix(~sample_condition)
design5

# create another contrast set
contrastTUM_D21vTUM_D28 <- createContrast(c(0, 0, 0, 0, 0, 1))

nrow(contrastTUM_D21vTUM_D28) == ncol(design5)

out_DA15 <- diffcyt(sce,
                    experiment_info = ei, design = design5, contrast = contrastTUM_D21vTUM_D28,
                    analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                    clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                    transform = FALSE, normalize = FALSE, min_samples = 1)

da15 <- rowData(out_DA15$res)
plotDiffHeatmap(sce, da15, top_n = 6, all = TRUE, fdr = FDR_cutoff)

save(list = ls(), file = "workspaceDA.rds")










