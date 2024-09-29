setwd('/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/4. downstream_analysis/')

###################### Define the 'analyse pathways' function ####################

analyze_pathways <- function(df, cell_type = 'Cell_Type', p_val = 0.05, lfc = 0.1, export_pathway_files = TRUE, cats) {
  original_wd <- getwd()
  
  # Subset only the significantly DEGs
  df <- df[df$adj.P.Val <= p_val, ]
  df <- df[abs(df$abs.log2FC) >= lfc, ]
  
  # Extract gene symbols
  gene_symbols <- df$X
  
  # Map gene symbols to Entrez IDs
  gene_entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  gene_entrez_ids <- as.character(gene_entrez_ids)
  
  # Prepare fold changes data
  foldchanges <- df$abs.log2FC
  names(foldchanges) <- gene_entrez_ids
  foldchanges <- na.omit(foldchanges)
  
  # Remove duplicate gene names
  foldchanges <- foldchanges[!duplicated(names(foldchanges))]
  
  foldchanges <- sort(foldchanges, decreasing = TRUE)
  
  # Perform KEGG pathway analysis using gseKEGG
  keggres <- gseKEGG(geneList = foldchanges,
                     organism = 'mmu',
                     keyType = "kegg",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = p_val,
                     verbose = FALSE)
  
  # Process upregulated KEGG pathways
  upregulated_pathways <- keggres@result %>%
    filter(NES > 0, p.adjust < p_val) %>%
    mutate(id = ID)
  upregulated_pathway_ids <- upregulated_pathways$ID
  upregulated_pvals <- upregulated_pathways$p.adjust
  
  # Process downregulated KEGG pathways
  downregulated_pathways <- keggres@result %>%
    filter(NES < 0, p.adjust < p_val) %>%
    mutate(id = ID)
  downregulated_pathway_ids <- downregulated_pathways$ID
  downregulated_pvals <- downregulated_pathways$p.adjust
  
  # Print informative messages
  cat("Total number of pathways found:", nrow(keggres@result), "\n")
  cat("Number of upregulated pathways:", nrow(upregulated_pathways), "\n")
  cat("Number of downregulated pathways:", nrow(downregulated_pathways), "\n")
  
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
    pathway_folder <- file.path(result_folder, paste0("KEGG_", cell_type, '_', pathway_direction))
    dir.create(pathway_folder)
    setwd(pathway_folder)
    for (pid in pathway_ids) {
      full_pathway_name <- fullnames$Description[fullnames$ID == pid]
      print(full_pathway_name)
      tryCatch({
        dir.create(full_pathway_name)
        setwd(full_pathway_name)
        pathview(gene.data = foldchanges, 
                 pathway.id = pid, 
                 species = 'mmu', 
                 expand.node = FALSE,
                 kegg.native = TRUE,
                 low = list(gene = "red"), 
                 mid = list(gene = "gray"), 
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
  
  ridge_data <- function(pathway_ids, pathway_pvals, pathway_direction, cats) {
    ridge_df <- data.frame()
    for (i in seq_along(pathway_ids)) {
      pid <- pathway_ids[i]
      pval <- pathway_pvals[i]
      genes_in_pathway <- keggres@geneSets[[pid]]
      log2fc_values <- foldchanges[names(foldchanges) %in% genes_in_pathway]
      gene_ratio <- length(log2fc_values) / length(genes_in_pathway)
      if (length(log2fc_values) > 0) {
        pathway_category <- unique(cats[cats$Pathway == pid, "Major"])
        pathway_subcategory <- unique(cats[cats$Pathway == pid, "Sub"])
        if (length(pathway_category) == 0) {
          pathway_category <- "Unknown"
        }
        if (length(pathway_subcategory) == 0) {
          pathway_subcategory <- "Unknown"
        }
        pathway_name <- pathway_name <- get_kegg_pathway_name(pid)
        pval_formatted <- formatC(pval, format = "e", digits = 2)
        stars <- ""
        if (pval < 0.001) {
          stars <- '<span style="font-family: monospace;">***</span>'
        } else if (pval < 0.01) {
          stars <- '<span style="font-family: monospace;">** </span>'
        } else if (pval < 0.05) {
          stars <- '<span style="font-family: monospace;"> * </span>'
        }
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
  
  ridge_df_up <- ridge_data(upregulated_pathway_ids, upregulated_pvals, "Upregulated", cats)
  ridge_df_down <- ridge_data(downregulated_pathway_ids, downregulated_pvals, "Downregulated", cats)
  ridge_df <- rbind(ridge_df_up, ridge_df_down)
  
  ridge_df$combined_category <- paste(ridge_df$category, ridge_df$subcategory, sep = " - ")
  ridge_df <- ridge_df %>% arrange(category, subcategory)
  ridge_df$pathway <- factor(ridge_df$pathway, levels = unique(ridge_df$pathway))
  
  save_ridge_plot <- function(plot, file_name, path = 'Images') {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
    full_path <- file.path(path, file_name)
    ggsave(filename = full_path, plot = plot, device = "png", width = 20, height = 8, bg = "white")
  }
  
  if (nrow(ridge_df) > 0) {
    category_colors <- c("Metabolism" = "red1", "Genetic Information Processing" = "deepskyblue2", 
                         "Environmental Information Processing" = "green3", "Cellular Processes" = "darkorchid1", 
                         "Organismal Systems" = "orange", "Human Diseases" = "coral4", "Unknown" = "gray")
    
    ridge_plot <- ggplot(ridge_df, aes(x = log2fc, y = pathway, fill = gene_ratio, color = category)) +
      geom_density_ridges(scale = 3, rel_min_height = 0.005, linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), name = "Gene Ratio") +
      scale_color_manual(values = category_colors) +
      theme_ridges() +
      theme(legend.position = "right",
            axis.text.y = ggtext::element_markdown()) +
      labs(title = paste("KEGG Pathway Enrichment - Density plots - ", cell_type),
           x = "log2FC",
           y = "Pathways",
           fill = "Gene Ratio",
           color = "Category")
    
    print(ridge_plot)
    
    save_ridge_plot(ridge_plot, paste0(cell_type, "_KEGG_density_plot.png"))
    
  } else {
    cat("No data available for ridge plots.\n")
  }
  
  return(list(upregulated = upregulated_pathways, downregulated = downregulated_pathways, keggres = keggres))
}



################################### Test Code ##################################

# Test single dataset
results <- analyze_pathways(df = spy_isc, cell_type = 'ISC', p_val = 0.05, lfc = 0, export_pathway_files = F, cats = cats)
# View(results$keggres$greater)

# Test all datasets in folder
crypt_path = '/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/'
# villi_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/villi/differential_expression'
# pseudo_path = '/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/crypt_pb/differential_expression'
crp_out_path = '/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/5. results_repository/5. KEGG_Results/crypt/'
# vil_out_path = '/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/4. downstream_analysis/KEGG_Results/villi'
# sdo_out_path = '/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/4. downstream_analysis/KEGG_Results/pseudo'

# process_files(input_directory = crypt_path, output_directory = crp_out_path, p_val = 0.05, lfc = 0, export_pathway_files = T, cats = cats)

# Report details about specific pathway
# pathway_report(df_gob, "mmu04010")



specific_pathway_analysis(spy_ent, pid = 'mmu04540', p_val = 0.05)
# pathway_report(df_isc, kegg_pathway = 'mmu04530')


# count_kegg_pathways("Tuba1b")

################################### Misc. ######################################

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161 