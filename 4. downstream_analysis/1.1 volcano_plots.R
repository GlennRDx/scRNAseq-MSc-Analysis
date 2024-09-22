# Loop through the list and create volcano plots
df_list = spy_eec_list

for (name in names(df_list)) {
  
  df = df_list[[name]]
  
  # df = df[df$adj.P.Val <= 0.9, ]
  
  plot = EnhancedVolcano(df,
                         lab = df$X,
                         x = 'logFC',
                         y = 'adj.P.Val',
                         title = paste("Volcano plot for", name),
                         subtitle = "CD vs HFD EEC",
                         FCcutoff = 1)
  
  # Display the plot (if in an interactive R session)
  print(plot)
}