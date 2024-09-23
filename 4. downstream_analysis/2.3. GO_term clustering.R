# Function to cluster GO terms based on their gene overlap
cluster_go_terms <- function(go_data) {
  
  # 'go_data' is a dataframe with columns 'id' (GO term) and 'SYMBOL' (list of genes).
  
  # Step 1: Extract GO terms and their associated gene lists
  go_terms <- go_data$id
  genes_list <- lapply(go_data$genes, function(gene_data) gene_data)
  
  # Step 2: Filter out GO terms with no genes (already done in script 1, but extra safety check)
  valid_indices <- sapply(genes_list, function(gene_data) !is.null(gene_data) && length(gene_data) > 0)
  valid_go_terms <- go_terms[valid_indices]
  valid_genes_list <- genes_list[valid_indices]
  
  # If no valid GO terms, return an empty dataframe
  if (length(valid_go_terms) == 0) {
    return(data.frame(GO_id = character(0), cluster = integer(0)))
  }
  
  # Step 3: Define Jaccard similarity function
  jaccard_similarity <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    return(intersection / union)
  }
  
  # Step 4: Initialize similarity matrix
  n <- length(valid_go_terms)
  similarity_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(valid_go_terms, valid_go_terms))
  
  # Step 5: Calculate pairwise Jaccard similarities
  for (i in 1:n) {
    for (j in i:n) {
      sim <- jaccard_similarity(valid_genes_list[[i]], valid_genes_list[[j]])
      similarity_matrix[i, j] <- sim
      similarity_matrix[j, i] <- sim  # Symmetric matrix
    }
  }
  
  # Step 6: Convert similarity matrix to distance matrix (1 - similarity)
  distance_matrix <- 1 - similarity_matrix
  
  # Step 7: Perform hierarchical clustering
  hc <- hclust(as.dist(distance_matrix), method = "average")
  
  # Step 8: Cut the dendrogram into clusters (adjust `h` for fewer/bigger clusters)
  clusters <- cutree(hc, h = 0.9)  # You can tweak `h` for more or less clusters
  
  # Step 9: Return dataframe with valid GO terms and their clusters
  result_df <- data.frame(GO_id = valid_go_terms, cluster = clusters)
  
  return(result_df)
}

# Example of how to use the function:
# Assuming 'results' is the output from the first script, with GO terms and gene lists
cluster_result <- cluster_go_terms(results$data)
print(cluster_result)
