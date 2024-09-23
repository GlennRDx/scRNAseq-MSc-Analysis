get_kegg_pathway_name <- function(pid) {
  pathway_info <- keggGet(pid)
  pathway_name <- pathway_info[[1]]$NAME
  pathway_name <- gsub(" - Mus musculus \\(house mouse\\)", "", pathway_name)
  return(pathway_name)
}

# Read and parse the KEGG category file
parse_kegg_categories <- function(filepath) {
  lines <- readLines(filepath)
  categories <- data.frame(Major = character(), Sub = character(), Pathway = character(), stringsAsFactors = FALSE)
  
  current_major <- NULL
  current_sub <- NULL
  
  for (line in lines) {
    if (startsWith(line, "MAJOR_")) {
      current_major <- substr(line, 7, nchar(line))  # Remove "MAJOR_" prefix
      current_sub <- NULL  # Reset subcategory
    } else if (grepl("^[0-9]{5}", line)) {
      pathway_id <- paste0("mmu", substr(line, 1, 5))
      categories <- rbind(categories, data.frame(Major = current_major, Sub = ifelse(is.null(current_sub), "No_Subcategory", current_sub), Pathway = pathway_id))
    } else {
      current_sub <- line
    }
  }
  
  return(categories)
}
cats = parse_kegg_categories('/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/4. downstream_analysis/kegg_categories')

# Find all genes associated with a KEGG pathway
get_gene_list <- function(identifier, type) {
  library(org.Mm.eg.db)
  library(KEGGREST)
  library(AnnotationDbi)
  if (type == "GO") {
    genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "GOALL", 
                                   keys = identifier, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else if (type == "KEGG") {
    kegg_genes <- gsub("mmu:", "", keggLink("mmu", identifier))
    genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "ENTREZID", 
                                   keys = kegg_genes, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else {
    stop("Invalid type. Use 'GO' or 'KEGG'.")
  }
  return(genes)
}

# Define the function to process a list of data frames
process_df_list <- function(df_list, kegg_pathway, pval_threshold = 0.05) {
  # Initialize empty data frames for p-values and logFC values
  pval_df <- data.frame()
  logFC_df <- data.frame()
  
  for (i in seq_along(df_list)) {
    df_name <- names(df_list)[i]
    data <- df_list[[i]]
    
    # Generate the pathway report for the current data frame
    intersecting_genes <- pathway_report(data, kegg_pathway, pval_threshold)
    
    # Create temporary data frames for p-values and logFC values
    temp_pval_df <- data.frame(gene = intersecting_genes$symbol, p.val = intersecting_genes$p.val)
    temp_logFC_df <- data.frame(gene = intersecting_genes$symbol, logFC = intersecting_genes$logFC)
    
    # Rename the columns to reflect the current data frame
    colnames(temp_pval_df)[2] <- df_name
    colnames(temp_logFC_df)[2] <- df_name
    
    # Merge the temporary data frames with the main data frames
    if (nrow(pval_df) == 0) {
      pval_df <- temp_pval_df
      logFC_df <- temp_logFC_df
    } else {
      pval_df <- merge(pval_df, temp_pval_df, by = "gene", all = TRUE)
      logFC_df <- merge(logFC_df, temp_logFC_df, by = "gene", all = TRUE)
    }
  }
  
  return(list(pval_df = pval_df, logFC_df = logFC_df))
}


# Define the function to count number of KEGG pathways associated with mmu genes
count_kegg_pathways <- function(gene_symbol) {
  # Convert the gene symbol to an Entrez gene ID
  entrez_id <- mapIds(org.Mm.eg.db, keys = gene_symbol, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  if (is.na(entrez_id)) {
    stop("Gene symbol not found in Entrez database.")
  }
  
  # Create the KEGG gene identifier
  kegg_gene_id <- paste0("mmu:", entrez_id)
  
  # Retrieve the KEGG pathways for the gene
  pathways <- keggLink("pathway", kegg_gene_id)
  
  # Count the number of unique pathways
  num_pathways <- length(unique(pathways))
  
  return(num_pathways)
}

process_files <- function(input_directory, output_directory, p_val = 0.05, lfc = 0.1, export_pathway_files = TRUE, cats) {
  # List all CSV files in the input directory
  file_list <- list.files(path = input_directory, pattern = "\\.csv$", full.names = TRUE)
  
  # Loop over each file in the file list
  for (file_path in file_list) {
    # Extract the file name and cell type
    file_name <- basename(file_path)
    parts <- strsplit(file_name, "_")[[1]]
    cell_type <- paste(parts[6:length(parts)], collapse = "_")
    cell_type <- sub("\\.csv$", "", cell_type) # Remove file extension '.csv'
    cell_type <- gsub("_", " ", cell_type) # Replace underscores with spaces
    
    print(file_name)
    print(cell_type)
    
    # Read the CSV file into a data frame
    df <- read.csv(file_path)
    
    # Run the analyze_pathways function
    setwd(output_directory)
    analyze_pathways(df, cell_type = cell_type, p_val = p_val, lfc = lfc, export_pathway_files = export_pathway_files, cats = cats)
  }
}


########################## Pathway report function #############################

# Reports pathway metrics from a given deg file
pathway_report <- function(data, kegg_pathway, pval_threshold = 0.05) {
  # Load necessary libraries
  library(dplyr)
  library(ggplot2)  # Library for plotting
  
  # Ensure that the required columns are present in the data
  if (!all(c("X", "adj.P.Val", "abs.log2FC") %in% colnames(data))) {
    stop("The data must contain 'X', 'adj.P.Val', and 'abs.log2FC' columns.")
  }
  
  # Get the list of genes for the given KEGG pathway
  gene_list <- get_gene_list(kegg_pathway, "KEGG")
  
  # Ensure that the SYMBOL column is present in the gene_list
  if (!"SYMBOL" %in% colnames(gene_list)) {
    stop("The gene list must contain a 'SYMBOL' column.")
  }
  
  # Extract the gene symbols from the gene_list
  kegg_genes <- gene_list$SYMBOL
  
  # Find the intersecting genes
  intersecting_genes <- data %>% 
    filter(X %in% kegg_genes) %>% 
    select(symbol = X, p.val = adj.P.Val, logFC = abs.log2FC) %>% 
    arrange(p.val)
  
  # Calculate the number of significant genes
  significant_genes <- intersecting_genes %>% 
    filter(p.val < pval_threshold)
  
  # Calculate the gene ratio
  total_intersecting <- nrow(intersecting_genes)
  total_significant <- nrow(significant_genes)
  
  gene_ratio_decimal <- if (total_intersecting > 0) {
    total_significant / total_intersecting
  } else {
    NA  # Avoid division by zero
  }
  
  gene_ratio_fraction <- if (!is.na(gene_ratio_decimal) && total_intersecting > 0) {
    paste(total_significant, total_intersecting, sep = "/")
  } else {
    NA  # Avoid division by zero
  }
  
  # Calculate the average logFC of the significant genes
  avg_logFC <- if (total_significant > 0) {
    mean(significant_genes$logFC, na.rm = TRUE)
  } else {
    NA  # If no significant genes, average is NA
  }
  
  # Create a distribution plot of the logFC for significant genes
  if (total_significant > 0) {
    p <- ggplot(significant_genes, aes(x = logFC)) +
      geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of logFC for Significant Genes",
           x = "logFC",
           y = "Count") +
      theme_minimal()
    
    print(p)  # Print the plot to the RStudio Plots pane
  }
  
  # Return results as a list
  return(list(
    intersecting_genes = intersecting_genes, 
    gene_ratio_decimal = gene_ratio_decimal,
    gene_ratio_fraction = gene_ratio_fraction,
    avg_logFC = avg_logFC
  ))
}

# Look at KEGG pathway graph of specific pathway
specific_pathway_analysis = function(df, pid, output_directory = '/home/', p_val = 0.05){
  owd = getwd()
  setwd(output_directory)
  df <- df[df$adj.P.Val <= p_val, ]
  
  # Extract gene symbols
  gene_symbols <- df$X
  
  # Map gene symbols to Entrez IDs
  gene_entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  gene_entrez_ids <- as.character(gene_entrez_ids)
  
  # Prepare fold changes data
  foldchanges <- df$abs.log2FC
  names(foldchanges) <- gene_entrez_ids
  foldchanges <- na.omit(foldchanges)
  
  pathview(gene.data = foldchanges, 
           pathway.id = pid, 
           species = 'mmu', 
           expand.node = T,
           kegg.native = T,
           low = list(gene = "red"), 
           mid =list(gene = "gray"), 
           high = list(gene = "green"))
  setwd(owd)
}