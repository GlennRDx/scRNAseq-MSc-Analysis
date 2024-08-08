library(aPEAR)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(stringr)

perform_GO_enrichment <- function(df, ONTOLOGY = NULL, orgdb = org.Mm.eg.db, showCategory = 50, nPermSimple = 10000) {
  set.seed(1)
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(orgdb, keys = df$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Remove duplicates
  df <- df[!duplicated(entrez_ids[df$X]), ]
  entrez_ids <- entrez_ids[!duplicated(entrez_ids)]
  
  # Create a named vector of logFC values
  geneList <- df$logFC
  names(geneList) <- entrez_ids[df$X]
  geneList <- na.omit(geneList)
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform GSEA
  set.seed(1)
  enrich <- gseGO(geneList = geneList, OrgDb = orgdb, ont = ONTOLOGY, eps = 0, nPermSimple = nPermSimple)
  
  # Extract the top 50 GO terms by p-value
  top50_go_terms <- head(enrich@result[order(enrich@result$pvalue), ], 50)
  
  # Print the top 50 GO terms
  print(top50_go_terms)
  
  # Assign dotplot to an object
  print('Creating dotplot...')
  enrich_dotplot <- enrichplot::dotplot(enrich, showCategory = showCategory, split = ".sign") + 
    facet_grid(. ~ .sign) +
    theme(axis.text.y = element_text(size = 8), 
          plot.title = element_text(size = 14), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))
  
  print(enrich_dotplot)
  
  print('Creating enrichment map')
  set.seed(1)
  settings = aPEAR.methods
  settings$minClusterSize = 3
  enrichnet = enrichmentNetwork(enrich@result, 
                    drawEllipses = TRUE,
                    repelLabels = T, 
                    outerCutoff = 0.6, # 0.5 Decreasing this value results in greater connectivity between the nodes in different clusters
                    innerCutoff = 0.005, # 0.1 Decreasing this value results in greater connectivity within the nodes in the same cluster.
                    fontSize = 4,
                    
                    )
  print(enrichnet)
  return(list(enrich = enrich, enrichnet = enrichnet))
}

enrichnet = perform_GO_enrichment(df = df_enp, ONTOLOGY = 'BP')
