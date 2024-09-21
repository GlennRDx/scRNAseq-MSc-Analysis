# Install and load packages
# BiocManager::valid()
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   BiocManager::install(version = "3.14")
# 
# BiocManager::install("ggtree")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("pathview")
# 
# options(repos = c(
#   PPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest", 
#   CRAN = "https://cloud.r-project.org"
# ))
# pak::pkg_install("clusterProfiler")
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)

# Load data
df = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_Enterocyte.csv')

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = df$index, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# Get significant genes
significant_genes <- df$index[df$P.Value < 0.05]
significant_entrez_ids <- entrez_ids[significant_genes]
significant_entrez_ids <- na.omit(significant_entrez_ids)
x <- paste0("mmu:",significant_entrez_ids[[1]])

# KEGG over-representation analysis
kegg_ora <- enrichKEGG(gene = x, organism = 'mmu', pvalueCutoff = 0.05)

# View and visualize results
head(kegg_ora)
dotplot(kegg_ora)