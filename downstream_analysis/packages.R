# Install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

# bioconductor packages
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DOSE")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("gage")
BiocManager::install("gageData")
BiocManager::install("KEGGREST")
BiocManager::install("AnnotationDbi")

# regular packages
pak('ggplot2')
pak('stringr')
pak('tools')
pak('ggridges')
pak('viridis')
pak('tidyverse')
pak('ggtext')
pak('dplyr')

# Load packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(stringr)
library(gage)
library(gageData)
library(tools)
library(ggridges)
library(viridis)
library(tidyverse)
library(KEGGREST)
library(ggtext)

