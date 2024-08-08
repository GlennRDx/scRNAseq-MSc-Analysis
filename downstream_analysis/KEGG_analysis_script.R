setwd('/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/r_studio/')
################################### Load packages ##############################
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(stringr)
library(gage)
library(gageData)
library(tools)
library(ggridges)
library(viridis)
library(tidyverse)
library(KEGGREST)
library(ggtext)

########################## Define secondary functions ##########################

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
cats = parse_kegg_categories('/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/r_studio/kegg_categories')

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

###################### Define the 'analyse pathways' function ####################

analyze_pathways <- function(df, cell_type = 'Cell_Type', p_val = 0.05, lfc = 0.1, export_pathway_files = TRUE, cats) {
  original_wd <- getwd()
  
  # Subset only the significantly DEGs
  df <- df[df$adj.P.Val <= p_val, ]
  df <- df[df$abs.log2FC >= lfc | df$abs.log2FC <= -lfc, ]
  
  # Extract gene symbols
  gene_symbols <- df$X
  
  # Map gene symbols to Entrez IDs
  gene_entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  gene_entrez_ids <- as.character(gene_entrez_ids)
  
  # Prepare fold changes data
  foldchanges <- df$abs.log2FC
  names(foldchanges) <- gene_entrez_ids
  foldchanges <- na.omit(foldchanges)
  
  # Perform GO enrichment analysis
  data("go.sets.mm")
  data("go.subs.mm")
  gobpsets <- go.sets.mm[go.subs.mm$BP]
  gobpres <- gage(exprs = foldchanges, gsets = gobpsets, same.dir = TRUE)
  
  # Perform KEGG pathway analysis
  data("kegg.sets.mm")
  data("sigmet.idx.mm")
  kegg.sets.mm <- kegg.sets.mm[sigmet.idx.mm]
  keggres <- gage(exprs = foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)
  
  # Extract and process upregulated KEGG pathways
  upregulated_pathways <- data.frame(id = rownames(keggres$greater), keggres$greater) %>%
    tibble::as_tibble() %>%
    filter(p.val < p_val) %>%
    .[['id']] %>%
    as.character()
  upregulated_pathway_ids <- substr(upregulated_pathways, start = 1, stop = 8)
  upregulated_pvals <- keggres$greater[, "p.val"]
  
  # Extract and process downregulated KEGG pathways
  downregulated_pathways <- data.frame(id = rownames(keggres$less), keggres$less) %>%
    tibble::as_tibble() %>%
    filter(p.val < p_val) %>%
    .[['id']] %>%
    as.character()
  downregulated_pathway_ids <- substr(downregulated_pathways, start = 1, stop = 8)
  downregulated_pvals <- keggres$less[, "p.val"]
  
  # Check if there are no significant pathways
  if (length(upregulated_pathway_ids) == 0 && length(downregulated_pathway_ids) == 0) {
    cat(paste("No significant pathways found for", cell_type, "\n"))
    return(NULL)
  }
  
  # Find intersection of upregulated and downregulated pathways
  intersected_pathways <- intersect(upregulated_pathway_ids, downregulated_pathway_ids)
  if (length(intersected_pathways) > 0) {
    cat("Intersection of upregulated and downregulated pathways:\n")
    cat(paste(intersected_pathways, collapse = ", "), "\n")
  } else {
    cat("No intersection between upregulated and downregulated pathways.\n")
  }
  
  if (export_pathway_files == TRUE) {
    # Create timestamp for result folder
    timestamp <- format(Sys.time(), "%H%M%S_%d%m%Y")
    result_folder <- paste0("KEGG_Results_", cell_type, '_', timestamp)
    dir.create(result_folder)  # Create main result folder
  }
  
  # Function to export pathway files individually
  export_pathways <- function(pathway_ids, fullnames, pathway_direction) {
    pathway_folder <- file.path(result_folder, paste0("KEGG_", cell_type, '_', pathway_direction)) # creates main folder
    dir.create(pathway_folder)  # Create subfolder (increased OR decreased)
    setwd(pathway_folder)
    # Perform pathview and export each pathway
    for (pid in pathway_ids) {
      full_pathway_name <- fullnames[which(substr(fullnames, 1, 8) == pid)]
      print(full_pathway_name)
      tryCatch({
        dir.create(full_pathway_name)
        setwd(full_pathway_name)
        pathview(gene.data = foldchanges, 
                 pathway.id = pid, 
                 species = 'mmu', 
                 expand.node = F,
                 kegg.native = T,
                 low = list(gene = "red"), 
                 mid =list(gene = "gray"), 
                 high = list(gene = "green"))
      }, 
      error = function(e) {
        cat(paste("Error for pathway", full_pathway_name, ":", e$message, "\n"))
      })
      setwd("..")
    }
  }
  
  if (export_pathway_files == TRUE) {
    # Export downregulated pathways
    setwd(original_wd)
    print(downregulated_pathway_ids)
    export_pathways(downregulated_pathway_ids, downregulated_pathways, "Decreased")
    
    # Export upregulated pathways
    setwd(original_wd)
    print(upregulated_pathway_ids)
    export_pathways(upregulated_pathway_ids, upregulated_pathways, "Increased")
    setwd(original_wd)
  }
  
  # Modify the ridge_data function to use different formatting for p-values and asterisks
  ridge_data <- function(pathway_ids, pathway_pvals, pathway_direction, cats) {
    ridge_df <- data.frame()
    for (i in seq_along(pathway_ids)) {
      pid <- pathway_ids[i]
      pval <- pathway_pvals[i]
      genes_in_pathway <- get_gene_list(pid, 'KEGG')
      log2fc_values <- foldchanges[names(foldchanges) %in% genes_in_pathway$ENTREZID]
      gene_ratio <- length(log2fc_values) / nrow(genes_in_pathway)  # Calculate gene ratio
      if (length(log2fc_values) > 0) {  # Check if log2fc_values is not empty
        pathway_category <- unique(cats[cats$Pathway == pid, "Major"])
        pathway_subcategory <- unique(cats[cats$Pathway == pid, "Sub"])
        if (length(pathway_category) == 0) {
          pathway_category <- "Unknown"
        }
        if (length(pathway_subcategory) == 0) {
          pathway_subcategory <- "Unknown"
        }
        pathway_name <- get_kegg_pathway_name(pid)
        # Format p-value in scientific notation
        pval_formatted <- formatC(pval, format = "e", digits = 2)
        # Add stars based on significance level
        stars <- ""
        if (pval < 0.001) {
          stars <- '<span style="font-family: monospace;">***</span>'
        } else if (pval < 0.01) {
          stars <- '<span style="font-family: monospace;">** </span>'  # Use two non-breaking spaces
        } else if (pval < 0.05) {
          stars <- '<span style="font-family: monospace;"> * </span>'  # Use three non-breaking spaces
        }
        # Format pathway label with name, p-value, and subcategory, with asterisks in monospace
        pathway_label <- paste0(pathway_name, " - ", pid, ' - ', pathway_subcategory, " - p-val: ", pval_formatted, " ", stars)
        ridge_df <- rbind(ridge_df, data.frame(pathway = rep(pathway_label, length(log2fc_values)), 
                                               log2fc = log2fc_values, 
                                               gene_ratio = rep(gene_ratio, length(log2fc_values)), 
                                               category = rep(pathway_category, length(log2fc_values)),
                                               subcategory = rep(pathway_subcategory, length(log2fc_values))))
      }
    }
    return(ridge_df)
  }
  
  # Updated the ridge_data function call with p-value formatting and star flagging
  ridge_df_up <- ridge_data(upregulated_pathway_ids, upregulated_pvals, "Upregulated", cats)
  ridge_df_down <- ridge_data(downregulated_pathway_ids, downregulated_pvals, "Downregulated", cats)
  ridge_df <- rbind(ridge_df_up, ridge_df_down)
  
  # Create a combined category for sorting
  ridge_df$combined_category <- paste(ridge_df$category, ridge_df$subcategory, sep = " - ")
  
  # Order by category and subcategory
  ridge_df <- ridge_df %>%
    arrange(category, subcategory)
  
  # Create a factor to maintain the order in the plot
  ridge_df$pathway <- factor(ridge_df$pathway, levels = unique(ridge_df$pathway))
  
  # Function to save the ridge plot and ensure the Images directory exists
  save_ridge_plot <- function(plot, file_name, path = 'Images') {
    # Check if the directory exists, if not, create it
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    full_path <- file.path(path, file_name)
    ggsave(filename = full_path, plot = plot, device = "png", width = 20, height = 8, bg = "white")
  }
  
  # Plot ridge plots with binning
  if (nrow(ridge_df) > 0) {
    category_colors <- c("Metabolism" = "red1", "Genetic Information Processing" = "deepskyblue2", 
                         "Environmental Information Processing" = "green3", "Cellular Processes" = "darkorchid1", 
                         "Organismal Systems" = "orange", "Human Diseases" = "coral4", "Unknown" = "gray")
    
    ridge_plot <- ggplot(ridge_df, aes(x = log2fc, y = pathway, fill = gene_ratio, color = category)) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), name = "Gene Ratio") +
      scale_color_manual(values = category_colors) +
      theme_ridges() +
      theme(legend.position = "right",
            axis.text.y = ggtext::element_markdown()) +  # Use ggtext to enable markdown rendering
      labs(title = paste("Ridge Plots for", cell_type, "Pathways"),
           x = "abs.log2FC",
           y = "Pathways",
           fill = "Gene Ratio",
           color = "Category")
    
    print(ridge_plot)
    
    save_ridge_plot(ridge_plot, paste0(cell_type, "_KEGG_density_plot.png"))
    
  } else {
    cat("No data available for ridge plots.\n")
  }
  
  
  return(list(upregulated = upregulated_pathways, downregulated = downregulated_pathways, gobpres = gobpres, keggres = keggres))
}



##################### process multiple files function ###################

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

################################### Test Code ##################################

# Test single dataset
results <- analyze_pathways(df = df_gob, cell_type = 'GOB', p_val = 0.05, lfc = 0, export_pathway_files = F, cats = cats)
View(results$keggres$greater)

# Test all datasets in folder
crypt_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/crypt/differential_expression'
villi_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/villi/differential_expression'
pseudo_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/crypt_pb/differential_expression'
crp_out_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/r_studio/KEGG_Results/crypt'
vil_out_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/r_studio/KEGG_Results/villi'
sdo_out_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/r_studio/KEGG_Results/pseudo'

process_files(input_directory = crypt_path, output_directory = crp_out_path, p_val = 0.05, lfc = 0, export_pathway_files = F, cats = cats)

# Report details about specific pathway
pathway_report(df_gob, "mmu04010")

# Look at KEGG pathway graph of specific pathway
specific_pathway_analysis = function(df, pid, p_val = 0.05){
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
}

specific_pathway_analysis(df_enp, pid = 'mmu03013', p_val = 0.05)
pathway_report(df_isc, kegg_pathway = 'mmu04530')


count_kegg_pathways("Tuba1b")

################################### Misc. ######################################

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161