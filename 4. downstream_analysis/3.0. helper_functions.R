# Load necessary data
data(go.sets.mm)
data(go.subs.mm)

get_gene_list <- function(identifier, type) {
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

perform_go_enrichment <- function(df, ont) {
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

extract_top_go_terms <- function(gores, n = 5) {
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

calculate_gene_ratio <- function(go_id, df) {
  # Create a list of genes in the dataframe where adj.P.Val is less than 0.05
  significant_genes <- df$X[df$adj.P.Val < 0.05]
  
  first_10_chars <- substr(go_id, 1, 10)
  print(first_10_chars)
  
  # Get the list of genes in the GO term
  go_genes <- get_gene_list(first_10_chars, "GO")
  go_gene_symbols <- unique(go_genes$SYMBOL)
  
  # Calculate the gene ratio
  intersection <- intersect(significant_genes, go_gene_symbols)
  gene_ratio <- length(intersection) / length(go_gene_symbols)
  
  print(gene_ratio)
  
  return(gene_ratio)
}

create_combined_df <- function(df, ont_list = c("BP", "CC", "MF"), total_terms = 30) {
  result_df <- data.frame()
  terms_per_ont <- ceiling(total_terms / length(ont_list) / 2)  # Divide by 2 for up and down regulation
  
  for (ont in ont_list) {
    gores <- perform_go_enrichment(df, ont)
    top_terms <- extract_top_go_terms(gores, n = terms_per_ont * 2)  # Extract more terms than needed
    
    valid_terms <- data.frame()
    for (i in 1:nrow(top_terms)) {
      gene_ratio <- tryCatch({
        calculate_gene_ratio(top_terms$id[i], df)
      }, error = function(e) {
        return(NA)
      })
      
      if (!is.na(gene_ratio)) {
        top_terms$gene_ratio[i] <- gene_ratio
        valid_terms <- rbind(valid_terms, top_terms[i, ])
      }
      
      if (nrow(valid_terms) == terms_per_ont * 2) break  # Stop if we have enough valid terms
    }
    
    if (nrow(valid_terms) < terms_per_ont * 2) {
      warning(paste("Not enough valid GO terms found for ontology", ont, ". Only", nrow(valid_terms), "terms included."))
    }
    
    valid_terms$ont <- ont
    
    # Get GO term names
    valid_terms$Term <- sapply(valid_terms$id, function(id) {
      tryCatch(
        Term(GOTERM[[id]]),
        error = function(e) id
      )
    })
    
    result_df <- rbind(result_df, valid_terms)
  }
  
  return(result_df)
}

plot_go_enrichment <- function(df) {
  ggplot(df, aes(x = gene_ratio, y = reorder(Term, gene_ratio))) +
    geom_point(aes(color = direction, size = -log10(p.val))) +
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    scale_size_continuous(range = c(2, 8)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = "Gene Ratio", y = "GO Term", color = "Direction", size = "-log10(p-value)")
}

# Run the analysis
result_df <- create_combined_df(spy_isc)
go_plot <- plot_go_enrichment(result_df)
print(go_plot)