pathway_heatmap <- function(df_list, 
                            pid = NULL, 
                            scale_to_one = FALSE, 
                            remove_na_rows = FALSE, 
                            order_by_sum = TRUE, 
                            custom_gene_list = NULL,
                            output_dir = ".",
                            file_name = "heatmap.jpg") {
  
  pathway_name = get_kegg_pathway_name(pid)
  
  # Step 1: Determine the gene list
  if (!is.null(custom_gene_list)) {
    gene_list <- unique(custom_gene_list)
  } else if (!is.null(pid)) {
    gene_list <- unique(get_gene_list(pid, 'KEGG')$SYMBOL)
  } else {
    stop("Either 'pid' or 'custom_gene_list' must be provided.")
  }
  
  # Initialize empty tables with appropriate dimensions
  p_val <- data.frame(matrix(ncol = length(df_list), nrow = length(gene_list)))
  logFC <- data.frame(matrix(ncol = length(df_list), nrow = length(gene_list)))
  
  # Set row names and column names
  rownames(p_val) <- gene_list
  rownames(logFC) <- gene_list
  
  # Extract the last 3 characters of each dataframe name and convert to uppercase
  colnames(p_val) <- toupper(substr(names(df_list), nchar(names(df_list)) - 2, nchar(names(df_list))))
  colnames(logFC) <- colnames(p_val)
  
  n_genes = length(gene_list)
  
  # Step 2: Populate the tables using df$P.value and df$abs.log2FC
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    
    # Ensure the df has the required columns
    if (!all(c("X", "adj.P.Val", "abs.log2FC") %in% colnames(df))) {
      stop("Dataframe missing required columns.")
    }
    
    for (gene in gene_list) {
      # Check if the gene exists in the current dataframe
      if (gene %in% df$X) {
        # Get the index of the gene
        gene_idx <- which(df$X == gene)
        
        # Populate p-value and logFC
        p_val[gene, i] <- df$adj.P.Val[gene_idx]
        
        # Set logFC to NaN if p-value is not significant
        if (df$adj.P.Val[gene_idx] < 0.05) {
          logFC[gene, i] <- df$abs.log2FC[gene_idx]
        } else {
          logFC[gene, i] <- NaN
        }
      }
    }
  }
  
  # Conditionally remove rows that are entirely NAs
  if (remove_na_rows) {
    non_na_rows <- apply(logFC, 1, function(row) !all(is.na(row)))
    p_val <- p_val[non_na_rows, ]
    logFC <- logFC[non_na_rows, ]
    n_omitted = sum(!non_na_rows)
    # Calculate the percentage of the pathway that has been omitted
    percentage_omitted = round((n_omitted / n_genes) * 100, 3)
  } else {
    # Move rows in p_val that are entirely NAs to the end
    na_rows <- apply(p_val, 1, function(row) all(is.na(row)))
    p_val <- rbind(p_val[!na_rows, ], p_val[na_rows, ])
    logFC <- logFC[rownames(p_val), ]  # Ensure the order of logFC matches p_val
    n_omitted = sum(na_rows)  # Omitted rows are those moved to the end
    percentage_omitted = 0
  }
  
  p_val_clean <- p_val
  p_val_clean[is.na(p_val_clean)] <- 1  # Assuming 1 is a non-significant p-value
  
  logFC_clean <- logFC
  logFC_clean[is.na(logFC_clean)] <- 0  # Assuming 0 is a reasonable substitute for missing logFC
  
  # Optionally sort the rows by the sum of significant logFC values
  if (order_by_sum) {
    significance_threshold <- 0.05  # Adjust this threshold if necessary
    significant_logFC <- logFC_clean
    significant_logFC[p_val_clean >= significance_threshold] <- 0  # Zero out non-significant logFC values
    logFC_sums <- rowSums(significant_logFC, na.rm = TRUE)
    
    # Order the rows by the sum of significant logFC values
    logFC_clean <- logFC_clean[order(logFC_sums, decreasing = TRUE), ]
    p_val_clean <- p_val_clean[rownames(logFC_clean), ]  # Ensure p_val matches the new order
  }
  
  # Determine the color scale range and create custom breaks
  if (scale_to_one) {
    max_abs_logFC <- max(abs(logFC_clean), na.rm = TRUE)
    color_limits <- c(-1, 1)
    if (max_abs_logFC > 1) {
      color_limits <- c(-max_abs_logFC, max_abs_logFC)
    }
  } else {
    color_limits <- c(-max(abs(logFC_clean), na.rm = TRUE), max(abs(logFC_clean), na.rm = TRUE))
  }
  
  breaks <- seq(color_limits[1], color_limits[2], length.out = 101)
  color_palette <- colorRampPalette(c("red", "white", "green"))(100)
  
  # Function to convert p-values to significance symbols
  pval_to_significance <- function(p) {
    if (p < 0.0001) {
      return("****")
    } else if (p < 0.001) {
      return("***")
    } else if (p < 0.01) {
      return("**")
    } else if (p < 0.05) {
      return("*")
    } else {
      return("")
    }
  }
  
  # Apply the function to the p-value table
  significance_symbols <- apply(p_val_clean, c(1, 2), pval_to_significance)
  
  # Generate the heatmap with significance annotations
  title_text <- if (remove_na_rows) {
    paste0(pathway_name, "\nHeatmap of logFC ", if (!is.null(pid)) pid else "Custom Genes", ' - ', n_omitted, " genes omitted (", percentage_omitted, "%)")
  } else {
    paste0(pathway_name, "\nHeatmap of logFC ", if (!is.null(pid)) pid else "Custom Genes")
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define the file path
  file_path <- file.path(output_dir, file_name)
  
  # Open a JPEG device
  png(filename = file_path, width = 1500, height = 2200, res = 200)
  
  # Plot the heatmap
  pheatmap(logFC_clean, 
           main = title_text, 
           cluster_rows = !order_by_sum,  # Cluster rows if not ordering by sum
           cluster_cols = FALSE, 
           display_numbers = significance_symbols, 
           color = color_palette, 
           breaks = breaks,
           border_color = NA,
           fontsize = 10,
           labels_col = colnames(logFC_clean))  # Use the modified column names
  
  # Close the JPEG device
  dev.off()
  
  message("Heatmap saved to ", file_path)
}


# # # Example usage with a custom gene list
# genes <- unique(unlist(lapply(df_list, function(df) df$X)))
# hsp_genes <- genes[grepl("^hsp", genes, ignore.case = TRUE)]

# Example usage with KEGG pathway
pathway_heatmap(df_list, 
                pid = 'mmu04972', 
                scale_to_one = TRUE, 
                remove_na_rows = TRUE, 
                order_by_sum = TRUE,
                output_dir = "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/downstream_analysis/KEGG_Results/crypt/Individual_Pathway_Analysis/",
                file_name = "heatmap.png")


# pathway_heatmap(df_list, custom_gene_list = hsp_genes, scale_to_one = T, remove_na_rows = T, order_by_sum = T)
