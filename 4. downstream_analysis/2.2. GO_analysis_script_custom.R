# Load necessary libraries
library(org.Mm.eg.db)
library(GO.db)
library(gage)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(stringdist)
library(dendextend)

# Load necessary data
data(go.sets.mm)
data(go.subs.mm)

# Function to retrieve gene list for GO or KEGG terms
get_gene_list <- function(identifier, type) {
  go_id <- strsplit(identifier, " ")[[1]][1]
  
  if (type == "GO") {
    genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "GOALL", 
                                   keys = go_id, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else if (type == "KEGG") {
    kegg_genes <- gsub("mmu:", "", keggLink("mmu", go_id))
    genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "ENTREZID", 
                                   keys = kegg_genes, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else {
    stop("Invalid type. Use 'GO' or 'KEGG'.")
  }
  return(genes)
}

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(df, ont) {
  gene_symbols <- df$X
  gene_entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  gene_entrez_ids <- as.character(gene_entrez_ids)
  
  foldchanges <- df$abs.log2FC
  names(foldchanges) <- gene_entrez_ids
  foldchanges <- na.omit(foldchanges)
  
  if (ont == "BP") {
    gosets <- go.sets.mm[go.subs.mm$BP]
  } else if (ont == "CC") {
    gosets <- go.sets.mm[go.subs.mm$CC]
  } else if (ont == "MF") {
    gosets <- go.sets.mm[go.subs.mm$MF]
  } else {
    stop("Invalid ontology. Use 'BP', 'CC', or 'MF'.")
  }
  
  gores <- gage(exprs = foldchanges, gsets = gosets, same.dir = TRUE)
  return(gores)
}

# Function to extract top GO terms
extract_top_go_terms <- function(gores, n) {
  up <- data.frame(id = rownames(gores$greater), gores$greater) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(p.val) %>%
    dplyr::slice_head(n = n)
  
  down <- data.frame(id = rownames(gores$less), gores$less) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(p.val) %>%
    dplyr::slice_head(n = n)
  
  up$direction <- "Up"
  down$direction <- "Down"
  
  return(rbind(up, down))
}

# Function to calculate gene ratio for GO term
calculate_gene_ratio <- function(go_id, df) {
  go_id <- strsplit(go_id, " ")[[1]][1]
  
  if (!(go_id %in% keys(GOTERM))) {
    message("Invalid GO term: ", go_id)
    return(NA)
  }
  
  significant_genes <- df$X[df$adj.P.Val < 0.05]
  
  go_genes <- tryCatch({
    get_gene_list(go_id, "GO")
  }, error = function(e) {
    message("Error retrieving genes for GO term: ", go_id)
    return(NULL)
  })
  
  if (is.null(go_genes)) return(NA)
  
  go_gene_symbols <- unique(go_genes$SYMBOL)
  intersection <- intersect(significant_genes, go_gene_symbols)
  gene_ratio <- length(intersection) / length(go_gene_symbols)
  
  return(gene_ratio)
}

# Function to create a combined dataframe with GO terms and gene ratios
create_combined_df <- function(df, ont_list, n) {
  result_df <- data.frame()
  
  for (ont in ont_list) {
    gores <- perform_go_enrichment(df, ont)
    top_terms <- extract_top_go_terms(gores, n)
    
    top_terms$gene_ratio <- sapply(top_terms$id, calculate_gene_ratio, df = df)
    top_terms$ont <- ont
    
    top_terms$Term <- sapply(top_terms$id, function(id) {
      tryCatch(Term(GOTERM[[id]]), error = function(e) id)
    })
    
    # Add gene_list column
    top_terms$gene_list <- sapply(top_terms$id, function(id) {
      go_genes <- get_gene_list(id, "GO")
      significant_genes <- df$X[df$adj.P.Val < 0.05]
      intersection <- intersect(significant_genes, unique(go_genes$SYMBOL))
      paste(intersection, collapse = ", ")
    })
    
    result_df <- rbind(result_df, top_terms)
  }
  
  # Remove rows where gene_ratio is NA
  result_df <- result_df[!is.na(result_df$gene_ratio), ]
  
  return(result_df)
}

# Function to cluster GO terms
cluster_go_terms <- function(result_df, n_clusters = 5) {
  # Create a matrix of gene lists
  gene_lists <- strsplit(result_df$gene_list, ", ")
  names(gene_lists) <- result_df$id
  
  # Calculate Jaccard similarity matrix
  jaccard_matrix <- matrix(0, nrow = length(gene_lists), ncol = length(gene_lists))
  rownames(jaccard_matrix) <- names(gene_lists)
  colnames(jaccard_matrix) <- names(gene_lists)
  
  for (i in 1:length(gene_lists)) {
    for (j in 1:length(gene_lists)) {
      jaccard_matrix[i, j] <- stringdist::stringsim(gene_lists[[i]], gene_lists[[j]], method = "jaccard")
    }
  }
  
  # Convert similarity matrix to distance matrix
  dist_matrix <- as.dist(1 - jaccard_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "complete")
  
  # Cut the dendrogram to get clusters
  clusters <- cutree(hc, k = n_clusters)
  
  # Add cluster information to the result dataframe
  result_df$cluster <- clusters[match(result_df$id, names(clusters))]
  
  return(result_df)
}

# Function to plot GO enrichment results
plot_go_enrichment <- function(df) {
  ggplot(df, aes(x = gene_ratio, y = reorder(Term, gene_ratio))) +
    geom_point(aes(color = factor(cluster), size = -log10(p.val))) +
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_discrete(name = "Cluster") +
    scale_size_continuous(range = c(2, 8)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = "Gene Ratio", y = "GO Term", size = "-log10(p-value)")
}

# Main analysis function
perform_enrichment_analysis <- function(df, ont_list = c("BP", "CC", "MF"), n = 5, n_clusters = 5) {
  # Perform the combined analysis
  result_df <- create_combined_df(df, ont_list, n)
  
  # Cluster GO terms
  result_df <- cluster_go_terms(result_df, n_clusters)
  
  # Create the plot
  go_plot <- plot_go_enrichment(result_df)
  
  return(list(data = result_df, plot = go_plot))
}

# Run the analysis with a specific dataframe (e.g., spy_isc) and customizable top n terms
# Assuming spy_isc is your dataframe with differential expression results
# spy_isc <- read.csv("path/to/your/differential_expression_results.csv")

results <- perform_enrichment_analysis(spy_isc, n = 10, n_clusters = 5)  # Modify n and n_clusters as needed

# Print the plot
print(results$plot)

# The clustered data is available in results$data
# You can write this to a file for further analysis if needed:
# write.csv(results$data, "clustered_go_terms.csv", row.names = FALSE)