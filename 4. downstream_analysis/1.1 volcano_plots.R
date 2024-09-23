# Loop through the list and create volcano plots
df_list = spy_list
for (name in names(df_list)) {
  
  df = df_list[[name]]
  
  # Extract the last three characters and convert to uppercase
  plot_title = toupper(substr(name, nchar(name) - 2, nchar(name)))
  
  plot = EnhancedVolcano(df,
                         lab = df$X,
                         x = 'logFC',
                         y = 'adj.P.Val',
                         title = paste("Volcano plot for", plot_title),
                         subtitle = "CD vs HFD",
                         FCcutoff = 0.5)
  
  # Display the plot (if in an interactive R session)
  print(plot)
}
