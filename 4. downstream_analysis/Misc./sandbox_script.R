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
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  result_folder <- paste0("KEGG_Results_", cell_type, timestamp)
  dir.create(result_folder)  # Create main result folder
  
  # Function to export pathway files individually
  export_pathways <- function(pathway_ids, fullnames, pathway_direction) {
    pathway_folder <- file.path(result_folder, paste0("KEGG_", cell_type, pathway_direction)) # creates main folder
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
    export_pathways(downregulated_pathway_ids, downregulated_pathways, "_Decreased")
    
    # Export upregulated pathways
    setwd(original_wd)
    print(upregulated_pathway_ids)
    export_pathways(upregulated_pathway_ids, upregulated_pathways, "_Increased")
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
        if (length(pathway_category) == 0) {
          pathway_category <- "Unknown"
        }
        pathway_name <- get_kegg_pathway_name(pid)
        ridge_df <- rbind(ridge_df, data.frame(pathway = rep(paste0(pathway_name, " - ", pid, "_", pathway_direction), length(log2fc_values)), 
                                               log2fc = log2fc_values, 
                                               gene_ratio = rep(gene_ratio, length(log2fc_values)), 
                                               category = rep(pathway_category, length(log2fc_values))))
      }
    }
    return(ridge_df)
  }
  
  ridge_df_up <- ridge_data(upregulated_pathway_ids, "Upregulated", cats)
  ridge_df_down <- ridge_data(downregulated_pathway_ids, "Downregulated", cats)
  ridge_df <- rbind(ridge_df_up, ridge_df_down)
  
  # Plot ridge plots with binning
  if (nrow(ridge_df) > 0) {
    category_colors <- c("Metabolism" = "red", "Genetic Information Processing" = "blue", 
                         "Environmental Information Processing" = "green", "Cellular Processes" = "purple", 
                         "Organismal Systems" = "orange", "Human Diseases" = "brown", "Unknown" = "gray")
    
    ridge_plot <- ggplot(ridge_df, aes(x = log2fc, y = pathway, fill = gene_ratio, color = category)) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, linewidth = 1.5) +
      scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +  # Use plasma color palette with limits from 0 to 1
      scale_color_manual(values = category_colors) +
      theme_ridges() +
      theme(legend.position = "right") +
      labs(title = paste("Ridge Plots for", cell_type, "Pathways"),
           x = "abs.log2FC",
           y = "Pathways",
           fill = "Gene Ratio",
           color = "Category") +
      scale_y_discrete(labels = function(x) sapply(strsplit(x, "_"), function(y) paste(y[2], y[1], sep = "_")))
    
    print(ridge_plot)
  } else {
    cat("No data available for ridge plots.\n")
  }
  
  return(list(upregulated = upregulated_pathways, downregulated = downregulated_pathways, intersected = intersected_pathways))
}

results <- analyze_pathways(df = df_eec, cell_type = 'EEC', p_val = 0.05, export_pathway_files = FALSE, cats = cats)

########################## Pathway Names ##################################

get_kegg_pathway_name <- function(pid) {
  pathway_info <- keggGet(pid)
  pathway_name <- pathway_info[[1]]$NAME
  pathway_name <- gsub(" - Mus musculus \\(house mouse\\)", "", pathway_name)
  return(pathway_name)
}


########################## Extract Categories ##################################

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

cats = parse_kegg_categories('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/r_studio/kegg_categories')

###############################################################################

# Find number of genes associated with go term
library(org.Mm.eg.db)
retrieved <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keytype = "GOALL", 
                                   keys = "GO:0030864", 
                                   columns = c("ENSEMBL", "SYMBOL"))

length(retrieved[!duplicated(retrieved$SYMBOL), ][,1])


################################################################################

go_ids <- go_bp_ISC@result$ID

# Convert to comma-separated format
go_ids_csv <- paste(go_ids, collapse = ", ")

# Write to a text file
write(go_ids_csv, file = "go_ids.txt")

# Print to console for verification
cat(go_ids_csv)


################################################################################

perform_GO_enrichment <- function(df, ONTOLOGY = NULL, pval_cutoff = 0.05, orgdb = org.Mm.eg.db, showCategory = 25) {
  # Load necessary libraries
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(orgdb, keys = df$index, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Create a named vector of logFC values
  geneList <- df$logFC
  names(geneList) <- entrez_ids[df$index]
  geneList <- na.omit(geneList)
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Perform GSEA
  gsea_result <- gseGO(geneList = geneList, OrgDb = orgdb, ont = ONTOLOGY, pvalueCutoff = pval_cutoff, eps = 0)
  
  # Append GO term IDs to the descriptions
  gsea_result@result$Description <- paste0(gsea_result@result$Description, " (", gsea_result@result$ID, ")")
  
  # Wrap long text labels
  gsea_result@result$Description <- str_wrap(gsea_result@result$Description, width = 50)
  
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
  p2 <- enrichplot::emapplot(gsea_result, showCategory = geneSets2plot[1:showCategory], cex_label_category = 1)
  print(p2)
}

# Example usage
perform_GO_enrichment(df_isc, ONTOLOGY = 'BP')
perform_GO_enrichment(df_isc, ONTOLOGY = 'MF')
perform_GO_enrichment(df_isc, ONTOLOGY = 'CC')


#################################################################################

df_isc = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_ISC.csv')
df = df_isc

# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(dplyr)
library(stringr)

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = df$index, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# Create a named vector of logFC values
geneList <- df$logFC
names(geneList) <- entrez_ids[df$index]
geneList <- na.omit(geneList)

# Sort geneList in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Perform GSEA
gsea_result <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = 'CC', pvalueCutoff = 0.05, eps = 0)

# Append GO term IDs to the descriptions
gsea_result@result$Description <- paste0(gsea_result@result$Description, " (", gsea_result@result$ID, ")")

# Assign dotplot to an object
p <- enrichplot::dotplot(gsea_result, showCategory = 25, split = ".sign") + 
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

# Create and print the emapplot
set.seed(123)
p2 <- enrichplot::emapplot(gsea_result, showCategory = 40)
p2 = emapplot(gsea_result, showCategory = 20, group_category = T, group_legend = T, cex_label_group = 1, node_label = "category", nCluster = 10) #https://github.com/YuLab-SMU/clusterProfiler/issues/347
print(p2)

################################################################################

# Naming emap clusters

# Extract clusters and GO terms
clusters <- data.frame(
  Term = p2$data$name,
  Cluster = names(p2$data$color2),
  stringsAsFactors = FALSE
)

# Define custom cluster names
custom_cluster_names <- c(
  "Metabolism Pathways",
  "Cell Cycle Regulation",
  "Immune Response",
  "Enter Annotation",
  "Enter Annotation",
  "Enter Annotation",
  "Enter Annotation",
  "Enter Annotation",
  "Enter Annotation",
  "Signal Transduction"
)

# Map custom cluster names to clusters
cluster_mapping <- data.frame(
  Cluster = as.character(1:length(custom_cluster_names)),
  CustomName = custom_cluster_names
)

# Map custom cluster names to clusters
clusters_mapped <- clusters %>%
  left_join(cluster_mapping, by = "Cluster")

# Create a custom legend data frame
custom_legend_df <- clusters %>%
  group_by(CustomName) %>%
  summarize(GO_Terms = paste(Term, collapse = "\n")) %>%
  ungroup()
