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
library(renv)

renv::init()

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create folders//directories

dirFlowFiles <- "flow_files"
dirPath <- paste(PrimaryDirectory, dirFlowFiles, sep = "/")
dir.create(dirPath)

dirCSVfiles <- "csv_files2"
CSVdirPath <- paste(dirPath, dirCSVfiles, sep = "/")
dir.create(CSVdirPath)

dirFCSfiles <- "fcs_files2"
FCSdirPath <- paste(dirPath, dirFCSfiles, sep = "/")
dir.create(FCSdirPath)

workingDir <- "230223_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
dir.create(workingDirPath)

# List files
CSVfiles <- list.files(CSVdirPath, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

# create sample_metadata.csv
fwrite(list(fileName), "sample_metadata2.csv")


# err_vals <- c(15, 44)
  
convertCSVtoFCS <- function(CSVfiles, csvDirectory, csv2fcsDirectory, fileName){
  for(i in c(1:length(CSVfiles))){
    data <- fread(paste(csvDirectory, CSVfiles[i], sep = "/"))
    print(CSVfiles[i])
    print(fileName[i])
    cytofCore.write.FCS(as.matrix(data), 
                        filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                        what = "numeric")
  }
}

convertCSVtoFCS(CSVfiles = CSVfiles[-2], csvDirectory = CSVdirPath, csv2fcsDirectory = FCSdirPath, fileName = fileName[-2])

# Create flowSet from FCSfiles
# filter selected FCS files out
FCSfiles <- list.files(FCSdirPath, pattern = ".fcs$", full = FALSE)
FCSfiles <- FCSfiles[!grepl('BM', FCSfiles)]
FCSfiles <- FCSfiles[!grepl('LN', FCSfiles)]
FCSfiles <- FCSfiles[!grepl('SP', FCSfiles)]
FCSfiles <- FCSfiles[!grepl('TEFF', FCSfiles)]
FCSfiles <- FCSfiles[!grepl('TPEX', FCSfiles)]

FCSfiles <- FCSfiles[grepl('*_[23456]', FCSfiles)]

FCSfiles <- FCSfiles


FCSfiles <- FCSfiles[!grepl('BM_[1234567]_TTERM', FCSfiles)]


FCSfiles

flowSet <- read.flowSet(files = FCSfiles, path = FCSdirPath, truncate_max_range = FALSE)
colnames(flowSet)

# prep sample_md
sample_md <- fread(file = paste(PrimaryDirectory, "sample_metadata2.csv", sep = "/"), header = TRUE)
sample_md

# remove selected samples

sample_md <- sample_md[-err_vals]

sample_md <- sample_md[!grepl('BM', file_name)]
sample_md <- sample_md[!grepl('LN', sample_id)]
sample_md <- sample_md[!grepl('SP', sample_id)]
sample_md <- sample_md[!grepl('TEFF', sample_id)]
sample_md <- sample_md[!grepl('TPEX', sample_id)]

sample_md <- sample_md[grepl('*_[23456]', sample_id)]

sample_md <- sample_md[-2]
sample_md


is.data.table(sample_md)
sample_md[, sample_id := factor(sample_id)]
                     
sample_md[ , condition := factor(condition, levels = c("BM_D14", "BM_D21", "BM_D28", "TUM_D14",
                                                        "TUM_D21","TUM_D28"))]





levels(sample_md$condition)
colnames(sample_md)
# fwrite(sample_md, "sample_metadata.csv")


# create panel_md.csv
#fwrite(list(colnames(flowSet)), "panel_md.csv")

panel_md <- fread(file = paste(PrimaryDirectory, "panel_md.csv", sep = "/"), header = TRUE)
all(colnames(flowSet) == panel_md$fcs_colname)

# preprocessing - QC

setwd(workingDirPath)

QC_dir <- "QC"
if(!dir.exists(QC_dir)){
  dir.create(QC_dir)
  dir.create(file.path(QC_dir, "flowAI"))
  dir.create(file.path(QC_dir, "PeacoQC"))
}
fs_AI <- flowCore::fsApply(flowSet, function(ff){
  resQC_AI <- flow_auto_qc(fcsfiles = ff,
                           folder_results = file.path(QC_dir, paste("flowAI", gsub(".fcs", " Folder", ff@description$GUID), sep = "/")),
                           output = 1)
  return(resQC_AI)
})

fs_PeacoQC <- flowCore::fsApply(flowSet, function(ff){
  resQC <- PeacoQC(ff = ff,
                   determine_good_cells = "all",
                   channels = c(10:28),
                   plot = TRUE,
                   output_folder = file.path(QC_dir, "PeacoQC"))
  ff <- ff[resQC$GoodCells, ]
  return(ff)
})

flowSet[[1]]@description$`$CYT`
flowSet[[1]]@description$`$CYT` <- "FACS"

chs_of_interest <- colnames(flowSet)[10:28]
plot_aggregate(flowSet, channels = chs_of_interest, output_image = "FCSpreNorm2.png")
plot_aggregate(fs_AI, channels = chs_of_interest, output_image = "FCSpreNormpostAI2.png")
plot_aggregate(fs_PeacoQC, channels = chs_of_interest, output_image = "FCSpreNormpostPeacoQC.png")
normFlowSet <- warpSet(flowSet, c("Comp-Alexa Fluor 488-A :: TCF1", "Comp-Alexa Fluor 700-A :: CD101", 
                                  "Comp-BV570-A :: CD44", "Comp-PE-Cy7-A :: BCL2", "Comp-BV510-A :: CD62L"))
plot_aggregate(normFlowSet, channels = chs_of_interest, output_image = "FCSpostNorm2.png")

# try fs_AI if no norm needed or apply norm
fs_AI <- flowSet
fs_AI <- normFlowSet

fs_AI[[1]]@description$`$CYT`
fs_AI[[1]]@description$`$CYT` <- "FACS"

# prepData - BM vs PB!
colnames(sample_md)
all(colnames(fs_AI) == panel_md$fcs_colname)
flowSet[[1]]@description$`$CYT`

sce <- CATALYST::prepData(fs_AI, panel = panel_md, md = sample_md, transform = FALSE)
assay(sce, "exprs") <- assay(sce, "counts")
assays(sce)

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 9
p

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")

saveRDS(sce, file = "SCE_3timepts_BMTUM_3.rds")
# file set 3 is post OTI cleaning

# ggplot(mtcars, aes(x = disp, y = hp)) + 
#   geom_point() +
#   theme_bw() +
#   geom_smooth()

renv::snapshot()


