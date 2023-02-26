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

load(file = "workspaceDA.rds")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "first", bars = TRUE, perc = TRUE)

# set outputPath
plotsPath <- paste(getwd(), "plots", sep = "/")
dir.create(plotsPath)

display.brewer.all(colorblindFriendly = FALSE)

# BM D14 v D21
plotDiffHeatmap(sce, da1, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test1 <- as_tibble(da1)
p.adj.signif_1 <- c("ns","ns","ns","*","ns","ns")
group1 <- (rep("BM_D14", nrow(stat_test1)))
group2 <- (rep("BM_D21", nrow(stat_test1)))
y.position <- c(10,10,10,10,10,10)
stat.test1 <- cbind(stat_test1, group1, group2, p.adj.signif_1, y.position)
stat.test1 <- as_tibble(stat.test1)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp1 <- bxp + stat_pvalue_manual(stat.test1, label = "p.adj.signif_1", tip.length = 0, size = 2.5)
bxp1

# BM 14 v D28

plotDiffHeatmap(sce, da2, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test2 <- as_tibble(da2)
p.adj.signif_2 <- c(rep("ns", 6))
group1 <- (rep("BM_D14", nrow(stat_test2)))
group2 <- (rep("BM_D28", nrow(stat_test2)))
y.position <- c(15,15,15,15,15,15)
stat.test2 <- cbind(stat_test2, group1, group2, p.adj.signif_2, y.position)
stat.test2 <- as_tibble(stat.test2)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp2 <- bxp1 + stat_pvalue_manual(stat.test2, label = "p.adj.signif_2", tip.length = 0, size = 2.5)
bxp2
# 
# # BM D14 v TUM D14
plotDiffHeatmap(sce, da3, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test3 <- as_tibble(da3)
p.adj.signif_3 <- c("****", "ns", "*", "****", "****", "****")
group1 <- (rep("BM_D14", nrow(stat_test3)))
group2 <- (rep("TUM_D14", nrow(stat_test3)))
#y.position <- c(100, 30, 15, 10, 10, 10)
y.position <- c(20,20,20,20,20,20)
stat.test3 <- cbind(stat_test3, group1, group2, p.adj.signif_3, y.position)
stat.test3 <- as_tibble(stat.test3)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp3 <- bxp2 + stat_pvalue_manual(stat.test3, label = "p.adj.signif_3", tip.length = 0,
                                    size = 2.5)
bxp3


# # BM D14 v TUM D21
plotDiffHeatmap(sce, da4, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test4 <- as_tibble(da4)
p.adj.signif_4 <- c("****", "ns", "****", "****", "****", "****")
group1 <- (rep("BM_D14", nrow(stat_test4)))
group2 <- (rep("TUM_D21", nrow(stat_test4)))
y.position <- c(25,25,25,25,25,25)
stat.test4 <- cbind(stat_test4, group1, group2, p.adj.signif_4, y.position)
stat.test4 <- as_tibble(stat.test4)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp4 <- bxp3 + stat_pvalue_manual(stat.test4, label = "p.adj.signif_4", tip.length = 0,
                                  size = 2.5)
bxp4

# BM D14 v TUM D28
plotDiffHeatmap(sce, da5, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test5 <- as_tibble(da5)
p.adj.signif_5 <- rep("****", 6)
group1 <- (rep("BM_D14", nrow(stat_test5)))
group2 <- (rep("TUM_D28", nrow(stat_test5)))
y.position <- c(30,30,30,30,30,30)
stat.test5 <- cbind(stat_test5, group1, group2, p.adj.signif_5, y.position)
stat.test5 <- as_tibble(stat.test5)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp5 <- bxp4 + stat_pvalue_manual(stat.test5, label = "p.adj.signif_5", tip.length = 0,
                                  size = 2.5)
bxp5


# BM D21 v D28
plotDiffHeatmap(sce, da6, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test6 <- as_tibble(da6)
p.adj.signif_6 <- c(rep("ns", 6))
group1 <- (rep("BM_D21", nrow(stat_test6)))
group2 <- (rep("BM_D28", nrow(stat_test6)))
y.position <- c(35,35,35,35,35,35)
stat.test6 <- cbind(stat_test6, group1, group2, p.adj.signif_6, y.position)
stat.test6 <- as_tibble(stat.test6)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp6 <- bxp5 + stat_pvalue_manual(stat.test6, label = "p.adj.signif_6", tip.length = 0,
                                  size = 2.5)
bxp6


# BM D21 v TUM D14
plotDiffHeatmap(sce, da7, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test7 <- as_tibble(da7)
p.adj.signif_7 <- c("****", "*", "**", "****", "****", "****")
group1 <- (rep("BM_D21", nrow(stat_test7)))
group2 <- (rep("TUM_D14", nrow(stat_test7)))
y.position <- c(40,40,40,40,40,40)
stat.test7 <- cbind(stat_test7, group1, group2, p.adj.signif_7, y.position)
stat.test7 <- as_tibble(stat.test7)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp7 <- bxp6 + stat_pvalue_manual(stat.test7, label = "p.adj.signif_7", tip.length = 0,
                                  size = 2.5)
bxp7


# BM D21 v TUM D21
plotDiffHeatmap(sce, da8, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test8 <- as_tibble(da8)
p.adj.signif_8 <- c(rep("****", 5), "*")
group1 <- (rep("BM_D21", nrow(stat_test8)))
group2 <- (rep("TUM_D21", nrow(stat_test8)))
y.position <- c(45,45,45,45,45,45)
stat.test8 <- cbind(stat_test8, group1, group2, p.adj.signif_8, y.position)
stat.test8 <- as_tibble(stat.test8)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp8 <- bxp7 + stat_pvalue_manual(stat.test8, label = "p.adj.signif_8", tip.length = 0,
                                  size = 2.5)
bxp8


# BM D21 v TUM D28
plotDiffHeatmap(sce, da9, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test9 <- as_tibble(da9)
p.adj.signif_9 <- c(rep("****", 6))
group1 <- (rep("BM_D21", nrow(stat_test9)))
group2 <- (rep("TUM_D28", nrow(stat_test9)))
y.position <- c(50,50,50,50,50,50)
stat.test9 <- cbind(stat_test9, group1, group2, p.adj.signif_9, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp9 <- bxp8 + stat_pvalue_manual(stat.test9, label = "p.adj.signif_9", tip.length = 0,
                                  size = 2.5)
bxp9

# BM D28 v TUM D14
plotDiffHeatmap(sce, da10, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test10 <- da10 %>% as_tibble()
p.adj.signif_10 <- c("****", "**", "**", "****", "****", '****')
group1 <- (rep("BM_D28", nrow(stat_test10)))
group2 <- (rep("TUM_D14", nrow(stat_test10)))
y.position <- c(55,55,55,55,55,55)
stat.test10 <- cbind(stat_test10, group1, group2, p.adj.signif_10, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp10 <- bxp9 + stat_pvalue_manual(stat.test10, label = "p.adj.signif_10", tip.length = 0,
                                  size = 2.5)
bxp10


# BM D28 v TUM D21
plotDiffHeatmap(sce, da11, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test11 <- da11 %>% as_tibble()
p.adj.signif_11 <- c("****", "**", "****", "****", "****", '****')
group1 <- (rep("BM_D28", nrow(stat_test11)))
group2 <- (rep("TUM_D21", nrow(stat_test11)))
y.position <- c(60,60,60,60,60,60)
stat.test11 <- cbind(stat_test11, group1, group2, p.adj.signif_11, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp11 <- bxp10 + stat_pvalue_manual(stat.test11, label = "p.adj.signif_11", tip.length = 0,
                                   size = 2.5)
bxp11

# BM D28 v TUM D28
plotDiffHeatmap(sce, da12, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test12 <- da12 %>% as_tibble()
p.adj.signif_12 <- c(rep("****", 6))
group1 <- (rep("BM_D28", nrow(stat_test12)))
group2 <- (rep("TUM_D28", nrow(stat_test12)))
y.position <- c(65,65,65,65,65,65)
stat.test12 <- cbind(stat_test12, group1, group2, p.adj.signif_12, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp12 <- bxp11 + stat_pvalue_manual(stat.test12, label = "p.adj.signif_12", tip.length = 0,
                                    size = 2.5)
bxp12


# TUM D14 v TUM D21
plotDiffHeatmap(sce, da13, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test13 <- da13 %>% as_tibble()
p.adj.signif_13 <- c(rep("ns", 6))
group1 <- (rep("TUM_D14", nrow(stat_test13)))
group2 <- (rep("TUM_D21", nrow(stat_test13)))
y.position <- c(70,70,70,70,70,70)
stat.test13 <- cbind(stat_test13, group1, group2, p.adj.signif_13, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp13 <- bxp12 + stat_pvalue_manual(stat.test13, label = "p.adj.signif_13", tip.length = 0,
                                    size = 2.5)
bxp13


# TUM D14 v TUM D28
plotDiffHeatmap(sce, da14, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test14 <- da14 %>% as_tibble()
p.adj.signif_14 <- c("**", "ns", "ns", "****", "**", "*")
group1 <- (rep("TUM_D14", nrow(stat_test14)))
group2 <- (rep("TUM_D28", nrow(stat_test14)))
y.position <- c(75,75,75,75,75,75)
stat.test14 <- cbind(stat_test14, group1, group2, p.adj.signif_14, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp14 <- bxp13 + stat_pvalue_manual(stat.test14, label = "p.adj.signif_14", tip.length = 0,
                                    size = 2.5)
bxp14


# TUM D21 v TUM D28
plotDiffHeatmap(sce, da15, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test15 <- da15 %>% as_tibble()
p.adj.signif_15 <- c("ns", "*", "ns", "****", "ns", "****")
group1 <- (rep("TUM_D21", nrow(stat_test15)))
group2 <- (rep("TUM_D28", nrow(stat_test15)))
y.position <- c(80,80,80,80,80,80)
stat.test15 <- cbind(stat_test15, group1, group2, p.adj.signif_15, y.position) %>% as_tibble()
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp15 <- bxp14 + stat_pvalue_manual(stat.test15, label = "p.adj.signif_15", tip.length = 0,
                                    size = 2.5)
bxp15
# view by sample to see outliers
# bxp2 <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "sample_id")
# bxp2 <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)
# bxp2

# ggsave("abundances_stat.svg", plot = last_plot(), dpi = 300)
# ggsave("abundances_stat.png", plot = last_plot(), dpi = 300, width = 6, height = 5)
# 

pdf("abundances_boxplot_withstats.pdf", width = 20, height = 10)
bxp15 + ylim(0,120)
dev.off()


svg("ExprHeatmap.svg", width = 10, height = 10)
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "median",
                scale = "last", bars = FALSE, perc = FALSE)
dev.off()

pdf("ExprHeatmap.pdf", width = 5, height = 6)
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "median",
                scale = "first", bars = FALSE, perc = TRUE)
dev.off()


plotDR(sce2, dr = "UMAP", color_by = "condition") + geom_point(size = 1.5)
# + geom_density2d(binwidth = 0.016, colour = "black")  # + geom_point(aes(size = 0.5))
ggsave("UMAP_wContours.svg", plot = last_plot())
ggsave("UMAP.pdf", plot = last_plot())
ggsave("UMAP.pdf", height = 5, width = 6, plot = last_plot())



CATALYST::plotDR(sce2, dr = "UMAP", color_by = c("Ly108", "CX3CR1", "CD62L", "CD44", "CD127", "KLRG1", "BCL2", "Ki67", "CXCR3", "CXCR5", "SCA1", "TIM3", "PD1", "CD101"), facet_by = "condition")

plotDR(sce, dr = "UMAP", color_by = c("CTV", "CD107a", "IL2", "IFNg", "TNFa","OVATet", "Ly108", "TIM3", "PD1", "CD101", "CX3CR1"), 
                              facet_by = "condition") # +  geom_density2d(binwidth = 0.016, colour = "grey")

ggsave("UMAP_wContours_bymarker.svg", height = 4, width = 20, plot = last_plot())
ggsave("UMAP_bymarker.pdf", height = 8, width = 20, plot = last_plot())
ggsave("UMAP_bymarker.png", height = 8, width = 18, dpi = 300, plot = last_plot())





CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")
ggsave("pbMDS.svg", plot = last_plot())
ggsave("pbMDS.pdf", plot = last_plot())


svg("Multiheatmap.svg", width = 12, height = 10)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "last")
dev.off()


pdf("Multiheatmap_4tps_BMTUM.pdf", width = 14, height = 6)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()
ggsave("Multiheatmap_4tps_BMTUM.png", width = 14, height = 6, dpi = 300, plot = last_plot())

pdf("Multiheatmap2.pdf", width = 8, height = 10)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "state", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()




plotClusterHeatmap(sce,
                   hm2 = NULL, k = "meta6", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
plotClusterExprs(sce, k = "meta6", features = "state")
ggsave("clusterexp_CTV.pdf", plot = last_plot())

plotClusterExprs(sce, k = "meta6", features = "state")
ggsave("clusterexp_state.pdf", plot = last_plot())
