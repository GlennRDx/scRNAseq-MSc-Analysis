############################# Download Packages ################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.19")

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('edgeR')

# pak::pkg_install("dplyr")

############################## Import Packages #################################

library(biomaRt)
library(topGO)
library(dplyr)
library(limma)
library(KEGGREST)
library(pathview)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

############################### Import Data ####################################

df_ent = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_Enterocyte.csv')
df_isc = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_ISC.csv')

.ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") # Initialise the Ensembl Mart for the mouse dataset

# Function to convert gene symbols in a data frame to Ensembl gene IDs
convert_to_ensembl <- function(df, ensembl = .ensembl, sort_col='P.Value') {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }
  # Retrieve the mapping from gene symbols to Ensembl Gene IDs
  genes_ent <- getBM(
    attributes = c("mgi_symbol", "ensembl_gene_id"),
    filters = "mgi_symbol",
    values = df[['index']],
    mart = ensembl
  )
  
  # Create a named vector of Ensembl Gene IDs
  genemap_ent <- setNames(genes_ent$ensembl_gene_id, genes_ent$mgi_symbol)
  
  # Map the gene identifiers in your data frame to Ensembl Gene IDs with a new column
  df$gene_id <- genemap_ent[df[['index']]]
  
  # Sort the data frame by P.Value
  df <- df %>% arrange(sort_col)
  
  return(df)
}
df_ent <- convert_to_ensembl(df_ent)
df_isc <- convert_to_ensembl(df_isc)

length(which(df_ent$adj.P.Val < 0.05)) # number of DE genes
head(df_ent)

########################### GO and KEGG analysis ###############################

create_topGOdata <- function(df, ontology = "ALL", mapping = "org.Mm.eg.db") {
  # Create the gene list
  geneList <- df$P.Value
  
  # Get the list of ENSEMBL to Entrez Gene mappings
  ensembl_to_entrez <- as.list(org.Mm.egENSEMBL2EG)
  
  # Extract the first part of the gene_id and use it to name geneList
  names(geneList) <- ensembl_to_entrez[sapply(strsplit(df$gene_id, split = "\\."), "[[", 1L)]
  
  # Create the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                geneSelectionFun = function(x) x,
                annot = annFUN.org,
                mapping = mapping)
  
  return(GOdata)
}
# GOdata_ent <- create_topGOdata(df_ent)
GOdata_isc <- create_topGOdata(df_isc)


# The topGOdata object is then used as input for enrichment testing (Kolmogorov-Smirnov)
resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
head(tab, 15)

# plot the GO graph for the 3 most significant terms and their parents
par(cex = 0.3)
showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 3, useInfo = "def")

#topGO Example Using Fisher’s Exact Test
GOdata <- new("topGOdata", # Create topGOData object
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x) (x < 0.05),
              annot = annFUN.org , mapping = "org.Mm.eg.db")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                numChar = 120)
head(tab)

# KEGG Pathway Enrichment Testing (KEGGREST)
# Pull all pathways for mmu
pathways.list <- keggList("pathway", "mmu")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list))
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[paste0("path:",outdat$pathway.code)]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)

# Plotting Pathways
foldChangeList <- df$logFC
xx <- as.list(org.Mm.egENSEMBL2EG)
names(foldChangeList) <- xx[sapply(strsplit(df$gene_id,split="\\."),"[[", 1L)]
head(foldChangeList)

########################## KEGG Pathway Analysis ###############################

# View KEGG pathway (Wnt signaling pathway)
mmu04930 <- pathview(gene.data  = FILL_IN_HERE,
                     pathway.id = "mmu04930",
                     species    = "mmu",
                     limit      = list(gene=max(abs(FILL_IN_HERE)), cpd=1))

# View KEGG pathway (PI3K / AKT disease)
entercode <- pathview(gene.data  = foldChangeList,
                     pathway.id = "mmu04151",
                     species    = "mmu",
                     limit      = list(gene=max(abs(foldChangeList)), cpd=1))

# View KEGG pathway (pathway name)
entercode <- pathview(gene.data  = foldChangeList,
                     pathway.id = "entercode",
                     species    = "mmu",
                     limit      = list(gene=max(abs(foldChangeList)), cpd=1))
















################################################################################
# # To illustrate the KS test, we plot probability distributions of p-values that are and that are not annotated with the term GO:0046661 “male sex differentiation”
# rna.pp.terms <- genesInTerm(GOdata)[["GO:0043491"]] # get genes associated with term
# p.values.in <- geneList[names(geneList) %in% rna.pp.terms]
# p.values.out <- geneList[!(names(geneList) %in% rna.pp.terms)]
# plot.ecdf(p.values.in, verticals = T, do.points = F, col = "red", lwd = 2, xlim = c(0,1),
#           main = "Empirical Distribution of DE P-Values by Annotation with 'PI3K/AKT signal transduction'",
#           cex.main = 0.9, xlab = "p", ylab = "Probabilty(P-Value < p)")
# ecdf.out <- ecdf(p.values.out)
# xx <- unique(sort(c(seq(0, 1, length = 201), knots(ecdf.out))))
# lines(xx, ecdf.out(xx), col = "black", lwd = 2)
# legend("bottomright", legend = c("Genes Annotated with 'PI3K/AKT signal transduction'", "Genes Not Annotated with 'PI3K/AKT signal transduction'"), lwd = 2, col = 2:1, cex = 0.9)



### REMOVED CODE ###

########################## Sanbiomics Tutorial #################################

df_ent_sig = df_ent$gene_id[which(df_ent$adj.P.Val < 0.05)]
GO_results_ent <- enrichGO(gene = df_ent_sig, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results_ent)

barplot(GO_results_ent, showCategory = 50)
dotplot(GO_results_ent, showCategory = 50)

########################## Gencore Tutorial #################################

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# we want the log2 fold change 
original_gene_list <- df_ent$abs.log2FC

# name the vector
names(original_gene_list) <- df_ent$gene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

gse <- pairwise_termsim(gse)
emapplot(gse, showCategory = 10)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

