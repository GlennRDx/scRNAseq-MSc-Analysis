# Create DGE_results folder
dir.create("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/", showWarnings = FALSE, recursive = TRUE)

# Function to run DESeq2 analysis for a single cell type
run_deseq2 <- function(cell_type) {
  # Set paths
  base_dir <- "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/DESeq2_data_files/"
  count_matrix_path <- file.path(base_dir, cell_type, "count_matrix.csv")
  metadata_path <- file.path(base_dir, cell_type, "metadata.csv")
  
  # Read data
  counts <- read.csv(count_matrix_path, row.names=1)
  counts <- t(counts)
  counts <- round(counts)
  metadata <- read.csv(metadata_path, row.names=1)
  
  # Ensure Sample is a factor
  metadata$Sample <- rownames(metadata)
  metadata$Sample <- factor(metadata$Sample)
  metadata$Diet <- factor(metadata$Diet)
  
  # Ensure the row names of metadata match the column names of counts
  all(rownames(metadata) %in% colnames(counts))
  all(rownames(metadata) == colnames(counts))
  
  cat(cell_type, "\n")
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ Diet)
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, name="Diet_HFHSD_vs_CD")
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  res_df$gene_symbol <- rownames(res_df)
  
  # Rename columns
  colnames(res_df)[colnames(res_df) == "gene_symbol"] <- "X"
  colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "LogFC"
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "P.Value"
  colnames(res_df)[colnames(res_df) == "padj"] <- "adj.P.Val"
  
  # Reorder columns
  res_df <- res_df[, c("X", "baseMean", "LogFC", "lfcSE", "stat", "P.Value", "adj.P.Val")]
  
  # Order by adjusted p-value
  res_df <- res_df[order(res_df$adj.P.Val), ]
  
  # Save results with new file naming convention
  output_filename <- paste0("diff_exp_CD_vs_HFD_", gsub(" ", "_", cell_type), ".csv")
  write.csv(res_df, file=file.path("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/", output_filename), row.names=FALSE)
  
  return(res_df)
}

# Get list of cell types
cell_types <- list.dirs("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/DESeq2_data_files/", full.names=FALSE, recursive=FALSE)

# Run DESeq2 for each cell type
results_list <- lapply(cell_types, run_deseq2)