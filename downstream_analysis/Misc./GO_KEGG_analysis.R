# Load packages
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

# Load data
setwd("/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/python/crypt/differential_expression")
df_pan = read.csv('diff_exp_CD_vs_HFD_Paneth_cell.csv')
df_ent = read.csv('diff_exp_CD_vs_HFD_Enterocyte.csv')
df_eec = read.csv('diff_exp_CD_vs_HFD_EEC.csv')
df_isc = read.csv('diff_exp_CD_vs_HFD_ISC.csv')

################################## GO FUNC #####################################

perform_GO_enrichment <- function(df, ONTOLOGY = NULL, pval_cutoff = 0.05, orgdb = org.Mm.eg.db, showCategory = 25) {

  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(orgdb, keys = df$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Create a named vector of logFC values
  geneList <- df$logFC
  names(geneList) <- entrez_ids[df$X]
  geneList <- na.omit(geneList)
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform GSEA
  gsea_result <- gseGO(geneList = geneList, OrgDb = orgdb, ont = ONTOLOGY, pvalueCutoff = pval_cutoff, eps = 0)
  
  # Append GO term IDs to the descriptions
  gsea_result@result$Description <- paste0(gsea_result@result$Description, " (", gsea_result@result$ID, ")")
  
  # Assign dotplot to an object
  p <- enrichplot::dotplot(gsea_result, showCategory = showCategory, split = ".sign") + 
    facet_grid(. ~ .sign) +
    theme(axis.text.y = element_text(size = 8), 
          plot.title = element_text(size = 14), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))
  
  print(p)
  
  # Extract the gene sets for the emapplot
  geneSets2plot <- p$data$Description
  
  # Calculate pairwise term similarities for the emapplot
  gsea_result <- pairwise_termsim(gsea_result)
  
  # Create and print the emapplot with a limited number of categories
  p2 = emapplot(gsea_result, showCategory, group_category = T, group_legend = T, cex_label_group = 1, node_label = "category", nCluster = 6) #https://github.com/YuLab-SMU/clusterProfiler/issues/347
  print(p2)
}

perform_GO_enrichment(df_eec, ONTOLOGY = 'BP')

###################################Test KEGG FUNC###############################

analyze_pathways <- function(df, cell_type = 'Cell_Type', p_val = 0.05, export_pathway_files = TRUE, cats) {
  original_wd <- getwd()
  
  # Subset only the significantly DEGs
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
  
  # Extract and process downregulated KEGG pathways
  downregulated_pathways <- data.frame(id = rownames(keggres$less), keggres$less) %>%
    tibble::as_tibble() %>%
    filter(p.val < p_val) %>%
    .[['id']] %>%
    as.character()
  downregulated_pathway_ids <- substr(downregulated_pathways, start = 1, stop = 8)
  
  # Find intersection of upregulated and downregulated pathways
  intersected_pathways <- intersect(upregulated_pathway_ids, downregulated_pathway_ids)
  if (length(intersected_pathways) > 0) {
    cat("Intersection of upregulated and downregulated pathways:\n")
    cat(paste(intersected_pathways, collapse = ", "), "\n")
  } else {
    cat("No intersection between upregulated and downregulated pathways.\n")
  }
  
  # Create timestamp for result folder
  timestamp <- format(Sys.time(), "%H%M%S_%d%m%Y")
  result_folder <- paste0("KEGG_Results_", cell_type, '_', timestamp)
  dir.create(result_folder)  # Create main result folder
  
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
        pathview(gene.data = foldchanges, pathway.id = pid, species = 'mmu')
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
  
  # Prepare data for ridge plots
  ridge_data <- function(pathway_ids, pathway_direction, cats) {
    ridge_df <- data.frame()
    for (pid in pathway_ids) {
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
        # Format pathway label with name, code, and subcategory
        pathway_label <- paste0(pathway_name, " - ", pathway_subcategory)
        ridge_df <- rbind(ridge_df, data.frame(pathway = rep(pathway_label, length(log2fc_values)), 
                                               log2fc = log2fc_values, 
                                               gene_ratio = rep(gene_ratio, length(log2fc_values)), 
                                               category = rep(pathway_category, length(log2fc_values)),
                                               subcategory = rep(pathway_subcategory, length(log2fc_values))))
      }
    }
    return(ridge_df)
  }
  
  ridge_df_up <- ridge_data(upregulated_pathway_ids, "Upregulated", cats)
  ridge_df_down <- ridge_data(downregulated_pathway_ids, "Downregulated", cats)
  ridge_df <- rbind(ridge_df_up, ridge_df_down)
  
  # Create a combined category for sorting
  ridge_df$combined_category <- paste(ridge_df$category, ridge_df$subcategory, sep = " - ")
  
  # Order by category and subcategory
  ridge_df <- ridge_df %>%
    arrange(category, subcategory)
  
  # Create a factor to maintain the order in the plot
  ridge_df$pathway <- factor(ridge_df$pathway, levels = unique(ridge_df$pathway))
  
  # Plot ridge plots with binning
  if (nrow(ridge_df) > 0) {
    category_colors <- c("Metabolism" = "red", "Genetic Information Processing" = "blue", 
                         "Environmental Information Processing" = "green", "Cellular Processes" = "purple", 
                         "Organismal Systems" = "orange", "Human Diseases" = "brown", "Unknown" = "gray")
    
    ridge_plot <- ggplot(ridge_df, aes(x = log2fc, y = pathway, fill = gene_ratio, color = category)) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), name = "Gene Ratio") +
      scale_color_manual(values = category_colors) +
      theme_ridges() +
      theme(legend.position = "right") +
      labs(title = paste("Ridge Plots for", cell_type, "Pathways"),
           x = "abs.log2FC",
           y = "Pathways",
           fill = "Gene Ratio",
           color = "Category")
    
    print(ridge_plot)
  } else {
    cat("No data available for ridge plots.\n")
  }
  
  return(list(upregulated = upregulated_pathways, downregulated = downregulated_pathways, intersected = intersected_pathways))
}

results <- analyze_pathways(df = df_isc, cell_type = 'EEC', p_val = 0.05, export_pathway_files = F, cats = cats)

##################################### Misc. ####################################

# Produce list of genes in specific go/kegg term
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


###### TEST CODE ######
# Assuming 'foldchanges' is a named vector with Entrez IDs as names and logFC as values
# Assuming 'splice_list' is a vector of Entrez IDs

splice_list = kegg.sets.mm$`mmu03040 Spliceosome`

# Extract Entrez IDs from 'foldchanges'
entrez_ids_foldchanges <- names(foldchanges)

# Find intersection
common_entrez_ids <- intersect(entrez_ids_foldchanges, splice_list)

# Count the number of common Entrez IDs
num_common_ids <- length(common_entrez_ids)

# Output the result
print(paste("Number of common Entrez IDs:", num_common_ids))


# Assuming 'foldchanges' is a named vector with Entrez IDs as names
# Assuming 'splice_list' is a vector of Entrez IDs
# Assuming 'df' is a data frame with a column 'Symbol' for gene symbols and 'log2FC' for log2 fold changes

# Extract Entrez IDs from 'foldchanges'
entrez_ids_foldchanges <- names(foldchanges)

# Find intersection of Entrez IDs in 'foldchanges' and 'splice_list'
common_entrez_ids <- intersect(entrez_ids_foldchanges, splice_list)

# Map common Entrez IDs to gene symbols
common_symbols <- mapIds(org.Mm.eg.db, keys = common_entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Remove NA values (if any)
common_symbols <- common_symbols[!is.na(common_symbols)]

# Filter 'df' to keep only the rows corresponding to the common gene symbols
df_common <- df %>% filter(X %in% common_symbols)

# Calculate the average absolute log2 fold change
average_abs_log2FC <- mean(df_common$abs.log2FC)

# Output the result
print(paste("Average absolute log2FC for common Entrez IDs:", average_abs_log2FC))










################################# KEGG FUNC ####################################

perform_KEGG_enrichment <- function(df, pval_cutoff = 0.05, organism = 'mmu', showCategory = 25) {
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, keys = df$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Create a named vector of logFC values
  geneList <- df$logFC
  names(geneList) <- entrez_ids[df$X]
  geneList <- na.omit(geneList)
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform KEGG GSEA
  kegg_gsea_result <- gseKEGG(geneList = geneList, organism = organism, pvalueCutoff = pval_cutoff, eps = 0)
  
  # Append KEGG pathway IDs to the descriptions
  kegg_gsea_result@result$Description <- paste0(kegg_gsea_result@result$Description, " (", kegg_gsea_result@result$ID, ")")
  
  # Remove unnecessary text from the descriptions
  kegg_gsea_result@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_gsea_result@result$Description)
  
  # Assign dotplot to an object
  p <- enrichplot::dotplot(kegg_gsea_result, showCategory = showCategory, split = ".sign") + 
    facet_grid(. ~ .sign) +
    theme(axis.text.y = element_text(size = 8), 
          plot.title = element_text(size = 14), 
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))
  
  print(p)
  
  # # Extract the gene sets for the emapplot
  # geneSets2plot <- p$data$Description
  # 
  # # Calculate pairwise term similarities for the emapplot
  # kegg_gsea_result <- pairwise_termsim(kegg_gsea_result)
  # 
  # # Create and print the emapplot with a limited number of categories
  # p2 <- enrichplot::emapplot(kegg_gsea_result, showCategory = geneSets2plot[1:showCategory], cex_label_category = 1)
  # print(p2)
  # 
  # pathview(gene.data = geneList, pathway.id = "mmu04930", species = "mmu", gene.idtype = "entrez")
}
perform_KEGG_enrichment(df_isc)


################################## KEGG ISC ####################################

# KEGG over-representation analysis
kegg_ora_ISC <- enrichKEGG(gene = significant_entrez_ids_ISC, organism = 'mmu', pvalueCutoff = 0.05)

# Modify the labels to include the KEGG pathway codes
kegg_ora_ISC@result$Description <- paste0(kegg_ora_ISC@result$Description, " (", kegg_ora_ISC@result$ID, ")")
kegg_ora_ISC@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_ora_ISC@result$Description)

# View and visualize results with cleaned and adjusted text labels
dotplot(kegg_ora_ISC, 
        title = 'KEGG Pathway Enrichment Analysis', 
        showCategory = 50, 
        font.size = 5) + 
  theme(axis.text.y = element_text(size = 8), 
        plot.title = element_text(size = 14), 
        legend.text = element_text(size = 10)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 80))

############################### GO Enterocytes #################################

# Load data
df_ent = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_Enterocyte.csv')

# Convert gene symbols to Entrez IDs
entrez_ids_ent <- mapIds(org.Mm.eg.db, keys = df_ent$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids_ent <- na.omit(entrez_ids_ent)

# Get significant genes
significant_genes_ent <- df_ent$X[df_ent$P.Value < 0.05]
significant_entrez_ids_ent <- entrez_ids_ent[significant_genes_ent]
significant_entrez_ids_ent <- na.omit(significant_entrez_ids_ent)

# GO over-representation analysis for Biological Process
go_bp_ent <- enrichGO(gene = significant_entrez_ids_ent, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.05)

# View and visualize results
head(go_bp_ent)
dotplot(go_bp_ent, title = 'Gene Ontology Analysis Enterocytes - Cellular Component')

############################## KEGG Enterocytes ################################

# KEGG over-representation analysis
kegg_ora_ent <- enrichKEGG(gene = significant_entrez_ids_ent, organism = 'mmu', pvalueCutoff = 0.05)

# View and visualize results
head(kegg_ora_ent)
dotplot(kegg_ora_ent, title = 'KEGG Pathway Enrichment Analysis - Enterocytes')