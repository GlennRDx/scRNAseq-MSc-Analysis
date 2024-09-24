# Load gene set data (make sure these datasets are available in your workspace)
data(go.sets.mm)
data(go.subs.mm)

# Function to check if the GOTERM object is available
check_goterm <- function() {
  if (!exists("GOTERM")) {
    if (!requireNamespace("GO.db", quietly = TRUE)) {
      stop("GO.db package is required but not installed. Please install it using BiocManager::install('GO.db')")
    }
    GOTERM <<- GO.db::GOTERM
  }
}

# Function to retrieve gene lists for GO terms
get_gene_list <- function(identifier, type = "GO") {
  if (type == "GO") {
    print(identifier)
    genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "GOALL", 
                                   keys = identifier, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else {
    stop("Invalid type. Currently only 'GO' is supported.")
  }
  return(genes)
}

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(df, ont) {
  # Subset only the significantly DEGs
  df <- df[df$adj.P.Val <= 0.05, ]
  gene_symbols <- df$X
  gene_entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  gene_entrez_ids <- as.character(gene_entrez_ids)
  
  foldchanges <- df$abs.log2FC
  names(foldchanges) <- gene_entrez_ids
  foldchanges <- na.omit(foldchanges)
  
  # Select GO sets based on ontology
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

# Function to extract top GO terms from the results
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

# Function to calculate gene ratio for a given GO term
calculate_gene_ratio <- function(go_id, df) {
  go_id <- strsplit(go_id, " ")[[1]][1]
  
  if (!(go_id %in% AnnotationDbi::keys(GOTERM))) {
    message("Invalid GO term: ", go_id)
    return(NA)
  }
  
  significant_genes <- df$X[df$adj.P.Val < 0.05]
  
  go_genes <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db, 
                          keys = go_id, 
                          columns = c("SYMBOL"), 
                          keytype = "GOALL")
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

# Function to count genes in each GO term
count_genes <- function(gene_list) {
  return(length(unique(unlist(gene_list))))
}

# Function to create a combined dataframe with GO terms, gene ratios, and gene lists
create_combined_df <- function(df, ont_list, n) {
  result_df <- data.frame()
  
  for (ont in ont_list) {
    gores <- perform_go_enrichment(df, ont)
    top_terms <- extract_top_go_terms(gores, n)
    
    # Calculate gene ratio and retrieve associated gene lists for each GO term
    top_terms$gene_ratio <- sapply(top_terms$id, calculate_gene_ratio, df = df)
    top_terms$ont <- ont
    
    # Split the 'id' column into 'id' and 'name'
    top_terms$name <- sub("^\\S+\\s", "", top_terms$id)
    top_terms$id <- substr(top_terms$id, 1, 10)
    
    # Retrieve GO term names from GOTERM
    top_terms$Term <- sapply(top_terms$id, function(id) {
      tryCatch(AnnotationDbi::Term(GOTERM[[id]]), error = function(e) id)
    })
    
    # Retrieve and store gene lists for each GO term
    top_terms$genes <- lapply(top_terms$id, function(go_id) {
      gene_data <- tryCatch(get_gene_list(go_id, "GO"), error = function(e) NULL)
      if (!is.null(gene_data)) return(unique(gene_data$SYMBOL)) else return(NA)
    })
    
    # Count the number of genes in each GO term
    top_terms$gene_count <- sapply(top_terms$genes, count_genes)
    
    result_df <- rbind(result_df, top_terms)
  }
  
  # Remove rows where gene_ratio or genes are NA
  result_df <- result_df[!is.na(result_df$gene_ratio) & !sapply(result_df$genes, is.null), ]
  
  return(result_df)
}

# Function to save enrichment results
save_enrichment_results <- function(results, filename) {
  saveRDS(results, file = filename)
}

# Function to load enrichment results
load_enrichment_results <- function(filename) {
  readRDS(file = filename)
}

# Main function for performing GO enrichment analysis with gene list storage and file saving
perform_and_save_enrichment_analysis <- function(df, name, save_dir, ont_list = c("BP", "CC", "MF"), n = 5) {
  # Check if GOTERM is available
  check_goterm()
  
  # Perform the combined analysis
  result_df <- create_combined_df(df, ont_list, n)
  
  # Create filename
  filename <- file.path(save_dir, paste0(name, "_enrichment_results.rds"))
  
  # Save results
  save_enrichment_results(result_df, filename)
  
  return(result_df)
}

# Function to process all datasets in spy_list
process_all_datasets <- function(spy_list, save_dir, ont_list = c("BP", "CC", "MF"), n = 5) {
  # Create save directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Process each dataset
  results_list <- lapply(names(spy_list), function(name) {
    cat("Processing", name, "...\n")
    df <- spy_list[[name]]
    perform_and_save_enrichment_analysis(df, name, save_dir, ont_list, n)
  })
  
  names(results_list) <- names(spy_list)
  return(results_list)
}


# Example usage:
# Define the save directory
save_dir <- "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/4. downstream_analysis/GO_results/crypt/GO_enrichment_files"

# Process all datasets
all_results <- process_all_datasets(spy_list, save_dir, n = 10)
