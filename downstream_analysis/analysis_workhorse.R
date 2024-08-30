# Process all DGE tables - Ridgeplots + KEGG Graphs
process_files(input_directory = crypt_path, output_directory = crp_out_path, p_val = 0.05, lfc = 0, export_pathway_files = T, cats = cats)

# Pathway
pathway_list = c('mmu04110', 'mmu03050')
directory = "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/downstream_analysis/KEGG_Results/crypt/Individual_Pathway_Analysis"
dataset = df_isc
cell_type = 'ISC'


# Iterate through pathways for cell type
for (pathway in pathway_list) {
  
  # Individual specific pathway analysis - Gene Heatmap
  pathway_heatmap(df_list, pid = pathway, scale_to_one = TRUE, remove_na_rows = TRUE, order_by_sum = TRUE, filename = "my_heatmap.png", directory = directory)
  
  # Individual Pathway Analysis - KEGG Graph
  specific_pathway_analysis(dataset, output_directory = directory, pid = pathway, p_val = 0.05)
  
}