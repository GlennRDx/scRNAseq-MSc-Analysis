# Loop through the list and create volcano plots
for (name in names(dsq_list)) {
  
  plot = EnhancedVolcano(df,
                         lab = df$X,
                         x = 'logFC',
                         y = 'adj.P.Val',
                         title = paste("Volcano plot for", name),
                         subtitle = "CD vs HFD EEC")
  
  # Display the plot (if in an interactive R session)
  print(plot)
}