# Clustering function
cluster_go_terms <- function(go_data, h = 0.9) {
  # Extract GO terms and their associated gene lists
  go_terms <- go_data$id
  genes_list <- go_data$genes
  
  # Filter out GO terms with no genes
  valid_indices <- sapply(genes_list, function(gene_data) !is.null(gene_data) && length(gene_data) > 0)
  valid_go_terms <- go_terms[valid_indices]
  valid_genes_list <- genes_list[valid_indices]
  
  # If no valid GO terms, return an empty dataframe
  if (length(valid_go_terms) == 0) {
    return(data.frame(GO_id = character(0), cluster = integer(0)))
  }
  
  # Define Jaccard similarity function
  jaccard_similarity <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    return(intersection / union)
  }
  
  # Initialize similarity matrix
  n <- length(valid_go_terms)
  similarity_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(valid_go_terms, valid_go_terms))
  
  # Calculate pairwise Jaccard similarities
  for (i in 1:n) {
    for (j in i:n) {
      sim <- jaccard_similarity(valid_genes_list[[i]], valid_genes_list[[j]])
      similarity_matrix[i, j] <- sim
      similarity_matrix[j, i] <- sim  # Symmetric matrix
    }
  }
  
  # Convert similarity matrix to distance matrix (1 - similarity)
  distance_matrix <- 1 - similarity_matrix
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(distance_matrix), method = "average")
  
  # Cut the dendrogram into clusters
  clusters <- cutree(hc, h = h)
  
  # Return dataframe with valid GO terms and their clusters
  result_df <- data.frame(GO_id = valid_go_terms, cluster = clusters)
  
  return(result_df)
}

# Plotting function
plot_go_enrichment <- function(df, font_size = 8, legend_size = 12, cluster_dot_size = 8) {
  # Generate a color palette based on the number of clusters
  n_clusters <- length(unique(df$cluster))
  color_palette <- scales::hue_pal()(n_clusters)
  
  # Create labels with p-value formatting and asterisks
  df$label <- sapply(1:nrow(df), function(i) {
    pval <- df$p.val[i]
    pval_formatted <- formatC(pval, format = "e", digits = 2)
    stars <- ""
    if (pval < 0.001) {
      stars <- '<span style="font-family: monospace;">***</span>'
    } else if (pval < 0.01) {
      stars <- '<span style="font-family: monospace;">** </span>'
    } else if (pval < 0.05) {
      stars <- '<span style="font-family: monospace;"> * </span>'
    }
    paste0(df$name[i], " - ", df$id[i], " - p-val: ", pval_formatted, " ", stars)
  })
  
  # Rank by cluster, then by NES
  df <- df %>%
    arrange(cluster, desc(abs(NES))) %>%
    mutate(rank = row_number())
  
  ggplot(df, aes(x = NES, y = reorder(label, rank))) +
    geom_point(aes(color = factor(cluster), size = setSize), shape = 21, stroke = 0.9, fill = NA) +  # Outline
    geom_point(aes(color = factor(cluster), size = setSize), alpha = 0.7) +  # Filled dot
    geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold", show.legend = FALSE) +  # Cluster number
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a vertical line at x=0
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(
      axis.text.y = ggtext::element_markdown(size = font_size),
      legend.title = element_text(size = legend_size),
      legend.text = element_text(size = legend_size),
      legend.key.size = unit(cluster_dot_size, "mm")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = cluster_dot_size))
    ) +
    labs(x = "Normalized Enrichment Score (NES)", y = "GO Term", color = "Cluster", size = "Gene Set Size")
}

# Or load a specific result file later:
enrichment_results_isc <- load_enrichment_results(file.path(save_dir, "spy_ent_enrichment_results.rds"))

# Perform clustering on a specific result:
clustered_results <- cluster_go_terms(enrichment_results, h = 0.95)

# Merge clustering results with enrichment results:
final_results <- merge(enrichment_results, clustered_results, by.x = "ID", by.y = "GO_id", all.x = TRUE)

# Create and display the plot:
go_plot_isc <- plot_go_enrichment(final_results, font_size = 10)
print(go_plot_isc)








# enrichment_results <- create_combined_df(spy_isc, ont_list = c("BP", "CC", "MF"), n = 20)
# clustered_results <- cluster_go_terms(enrichment_results, h = 0.95)
# final_results <- merge(enrichment_results, clustered_results, by = "ID", all.x = TRUE)
# go_plot <- plot_go_enrichment(final_results, font_size = 10)
# print(go_plot)