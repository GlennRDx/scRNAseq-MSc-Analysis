BiocManager::install("SingleCellExperiment")
library("SingleCellExperiment")
BiocManager::install("Seurat")
library("Seurat")
BiocManager::install("scater")
library("scater")
BiocManager::install("SCINA")
library("SCINA")
BiocManager::install("devtools")
library("devtools")
BiocManager::install("dplyr")
library("dplyr")
BiocManager::install("scmap")
library("scmap")
BiocManager::install("celldex")
library("celldex")
BiocManager::install("SingleR")
library("SingleR")
BiocManager::install("ggplot2")
library("ggplot2")
BiocManager::install("Harmony")
library("Harmony")
BiocManager::install("cerebroApp")
library("cerebroApp")
BiocManager::install("msigdb")
library("msigdb")


atlas_metadata = read.table('/home/glennrd/Downloads/Atlas/atlas_metadata_regional.txt') # 	Study metadata file for all cells.
colnames(atlas_metadata) <- atlas_metadata[1, ]
atlas_metadata <- subset(atlas_metadata, atlas_metadata$Region == "Ileum")
UMIs_metadata = as.character(atlas_metadata$NAME)

atlas_UMIcounts = read.table('/home/glennrd/Downloads/Atlas/GSE92332_Regional_UMIcounts.txt')
UMIs = as.character(colnames(atlas_UMIcounts))

matching_columns <- UMIs %in% UMIs_metadata
atlas_UMIcounts_subset <- atlas_UMIcounts[, matching_columns]

# Change column names to only cell labels
new_column_names <- gsub(".*?_([^_]*_){2}", "", colnames(atlas_UMIcounts_subset))
unique(new_column_names)

# Convert the data into a SingleCellExperiment object
ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=Matrix::Matrix(assays(ref)$logcounts)), 
                                                      colData=colData(ref), rowData=rowData(ref))
