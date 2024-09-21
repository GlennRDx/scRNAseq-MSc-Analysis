# Load required libraries
library(Seurat)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Define the gene lists
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

################################################################################

# Basic function to convert human to mouse gene names
mapfun <- function(humangenes) {
  # Map human gene symbols to ENTREZ IDs
  hgns <- mapIds(org.Hs.eg.db, humangenes, "ENTREZID", "SYMBOL")
  
  # Map human ENTREZ IDs to mouse ENTREZ IDs
  mapped <- select(Orthology.eg.db, hgns, "Mus_musculus", "Homo_sapiens")
  
  # Remove rows with no mouse counterpart
  naind <- is.na(mapped$Mus_musculus)
  
  # Map mouse ENTREZ IDs to mouse gene symbols
  msymb <- mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus[!naind]), "SYMBOL", "ENTREZID")
  
  # Prepare output
  out <- data.frame(Human_symbol = humangenes, mapped, Mouse_symbol = NA)
  out$Mouse_symbol[!naind] <- msymb
  out
}

# Convert s.genes and g2m.genes from human to mouse
s.genes.mouse <- mapfun(s.genes)
g2m.genes.mouse <- mapfun(g2m.genes)

s.genes.mouse = s.genes.mouse$Mouse_symbol
g2m.genes.mouse = g2m.genes.mouse$Mouse_symbol

s.genes.mouse <- na.omit(s.genes.mouse)
g2m.genes.mouse <- na.omit(g2m.genes.mouse)

# Write s.genes.mouse to a text file
writeLines(s.genes.mouse, "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/s_genes_mouse.txt")

# Write g2m.genes.mouse to a text file
writeLines(g2m.genes.mouse, "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/g2m_genes_mouse.txt")

