# Function to perform GO enrichment analysis using gseGO
perform_go_enrichment <- function(df, ont) {
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, keys = df$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Remove duplicates
  df <- df[!duplicated(entrez_ids[df$X]), ]
  entrez_ids <- entrez_ids[!duplicated(entrez_ids)]
  
  # Create a named vector of logFC values
  geneList <- df$abs.log2FC
  names(geneList) <- entrez_ids[df$X]
  geneList <- na.omit(geneList)
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform GSEA
  set.seed(1)
  enrich <- gseGO(geneList = geneList, 
                  OrgDb = org.Mm.eg.db, 
                  ont = ont, 
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  verbose = FALSE)
  
  return(enrich)
}

# Function to extract top GO terms from the results
extract_top_go_terms <- function(enrich, n) {
  result <- as.data.frame(enrich)
  result <- result %>%
    arrange(pvalue) %>%
    slice_head(n = n)
  
  result$direction <- ifelse(result$NES > 0, "Up", "Down")
  
  return(result)
}

# Function to create a combined dataframe with GO terms, gene ratios, and gene lists
create_combined_df <- function(df, ont_list, n) {
  result_df <- data.frame()
  
  for (ont in ont_list) {
    enrich <- perform_go_enrichment(df, ont)
    top_terms <- extract_top_go_terms(enrich, n)
    
    top_terms$ont <- ont
    top_terms$gene_ratio <- top_terms$setSize / nrow(df)
    top_terms$gene_count <- top_terms$setSize
    
    result_df <- rbind(result_df, top_terms)
  }
  
  return(result_df)
}

# Clustering function
cluster_go_terms <- function(go_data, h = 0.9) {
  print("Starting clustering...")
  print(paste("Number of GO terms:", nrow(go_data)))
  
  # Extract GO terms and their associated gene lists
  go_terms <- go_data$ID
  genes_list <- strsplit(go_data$core_enrichment, "/")
  
  print(paste("Number of gene lists:", length(genes_list)))
  
  # Filter out GO terms with no genes
  valid_indices <- sapply(genes_list, function(gene_data) !is.null(gene_data) && length(gene_data) > 0)
  valid_go_terms <- go_terms[valid_indices]
  valid_genes_list <- genes_list[valid_indices]
  
  print(paste("Number of valid GO terms:", length(valid_go_terms)))
  
  # If no valid GO terms, return an empty dataframe
  if (length(valid_go_terms) == 0) {
    print("No valid GO terms found. Returning empty dataframe.")
    return(data.frame(ID = character(0), cluster = integer(0)))
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
  result_df <- data.frame(ID = valid_go_terms, cluster = clusters)
  
  print(paste("Number of clusters:", length(unique(clusters))))
  print("Clustering completed.")
  
  return(result_df)
}

# Plotting function
plot_go_enrichment <- function(df, font_size = 8, legend_size = 12, cluster_dot_size = 8, title = NULL) {
  # Check for required columns
  required_cols <- c("ID", "Description", "pvalue", "ont", "direction", "gene_ratio", "gene_count", "cluster")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure all required columns are present and non-null
  df <- df %>%
    filter(!is.na(ID) & !is.na(Description) & !is.na(pvalue) & 
             !is.na(ont) & !is.na(direction) & !is.na(gene_ratio) & 
             !is.na(gene_count) & !is.na(cluster))
  
  # Generate a color palette based on the number of clusters
  n_clusters <- length(unique(df$cluster))
  color_palette <- scales::hue_pal()(n_clusters)
  
  # Create labels with p-value formatting and asterisks
  df$label <- sapply(1:nrow(df), function(i) {
    pval <- df$pvalue[i]
    pval_formatted <- formatC(pval, format = "e", digits = 2)
    stars <- ""
    if (pval < 0.001) {
      stars <- '<span style="font-family: monospace;">***</span>'
    } else if (pval < 0.01) {
      stars <- '<span style="font-family: monospace;">** </span>'
    } else if (pval < 0.05) {
      stars <- '<span style="font-family: monospace;"> * </span>'
    }
    paste0(df$Description[i], " - ", df$ID[i], " - p-val: ", pval_formatted, " ", stars)
  })
  
  # Rank by cluster, then by descending NES value
  df <- df %>%
    arrange(cluster, abs(NES)) %>%
    mutate(rank = row_number())
  
  # Create the plot
  p <- ggplot(df, aes(x = abs(NES), y = reorder(label, rank))) +
    geom_point(aes(color = factor(cluster), size = gene_count), shape = 21, stroke = 0.9, fill = NA) +  # Outline
    geom_point(aes(color = factor(cluster), size = gene_count), alpha = 0.7) +  # Filled dot
    geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold", show.legend = FALSE) +  # Cluster number
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(range = c(2, 10)) +
    scale_x_continuous(limits = c(0, NA)) +  # Set x-axis to start from 0
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
    labs(x = "Absolute NES", y = "GO Term", color = "Cluster", size = "Gene Count", 
         title = title)  # Add the user-provided title
  
  return(p)
}


# # Function to save the plot to a specific directory
# save_plot <- function(plot, filename, output_dir, width = 2400, height = 1000, res = 200) {
#   file_path <- file.path(output_dir, paste0(filename, ".png"))
#   ggsave(file_path, plot = plot, width = width / res, height = height / res, dpi = res, units = "in")
# }

# Main workflow function (prints plots to the plotting environment)
main_workflow <- function(df, ont_list = c("BP", "CC", "MF"), n = 20, h = 0.95, title = NULL) {
  # Perform enrichment analysis
  enrichment_results <- create_combined_df(df, ont_list, n)
  
  # Perform clustering
  clustered_results <- cluster_go_terms(enrichment_results, h)
  
  # Merge enrichment and clustering results
  final_results <- merge(enrichment_results, clustered_results, by = "ID", all.x = TRUE)
  
  # Create the plot
  go_plot <- plot_go_enrichment(final_results, font_size = 10, title = title)
  
  # Print the plot to the R plotting environment
  print(go_plot)
  
  return(list(results = final_results, plot = go_plot))
}

main_workflow(spy_eec, n = 50, title = 'EEC')

df_list = dsq_list

# Loop over dataset list and apply the workflow
for (df_name in names(df_list)) {
  # Extract the dataset and set the title based on the last 3 characters of dataset, in uppercase
  dataset <- df_list[[df_name]]
  title_suffix <- toupper(substr(df_name, nchar(df_name) - 2, nchar(df_name)))  # Extract and capitalize last 3 characters
  title <- paste("GO Term Enrichment - Dot Plots -", title_suffix)  # Construct the title
  
  # Run the workflow and print the plot to the plotting environment
  main_workflow(df = dataset, title = title)
}
