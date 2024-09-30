# Loop through the list and create volcano plots
# Loop through the list and create volcano plots
df_list = spy_list
for (name in names(df_list)) {
  
  df = df_list[[name]]
  
  # Extract the last three characters and convert to uppercase
  plot_title = toupper(substr(name, nchar(name) - 2, nchar(name)))
  
  # Define thresholds for significant DEGs
  significance_threshold <- 0.05
  logFC_threshold <- 0.5
  
  # Filter significantly differentially expressed genes
  sig_genes <- df[df$adj.P.Val < significance_threshold & abs(df$logFC) > logFC_threshold, ]
  
  # Print the number of significantly differentially expressed genes
  print(paste("Number of significantly differentially expressed genes for", plot_title, ":", nrow(sig_genes)))
  
  # Sort by logFC to get top up-regulated and down-regulated genes
  top_up_genes <- head(sig_genes[order(-sig_genes$logFC), ], 10)  # Top 5 up-regulated
  top_down_genes <- head(sig_genes[order(sig_genes$logFC), ], 10)  # Top 5 down-regulated
  
  # Print top 5 up-regulated and down-regulated genes
  print("Top 5 up-regulated genes:")
  print(top_up_genes$X)
  
  print("Top 5 down-regulated genes:")
  print(top_down_genes$X)
  
  # Create the volcano plot
  plot = EnhancedVolcano(df,
                         lab = df$X,
                         x = 'logFC',
                         y = 'adj.P.Val',
                         title = paste("Volcano plot for", plot_title),
                         subtitle = "CD vs HFD",
                         FCcutoff = logFC_threshold)
  
  # Display the plot (if in an interactive R session)
  print(plot)
}

