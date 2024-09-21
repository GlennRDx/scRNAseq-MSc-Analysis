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

df = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_ISC.csv')

entrez_ids <- mapIds(org.Mm.eg.db, keys = df$index, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# Create a named vector of logFC values
geneList <- df$logFC
names(geneList) <- entrez_ids[df$index]
geneList <- na.omit(geneList)

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Perform KEGG GSEA
kegg_gsea_result <- gseKEGG(geneList = geneList, organism = 'mmu', pvalueCutoff = 0.05, eps = 0)

# Append KEGG pathway IDs to the descriptions
kegg_gsea_result@result$Description <- paste0(kegg_gsea_result@result$Description, " (", kegg_gsea_result@result$ID, ")")

# Remove unnecessary text from the descriptions
kegg_gsea_result@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_gsea_result@result$Description)

# Assign dotplot to an object
p <- enrichplot::dotplot(kegg_gsea_result, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)  +
  theme(axis.text.y = element_text(size = 8), 
        plot.title = element_text(size = 14), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))

print(p)
