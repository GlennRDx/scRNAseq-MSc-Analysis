# Process all DGE tables - Ridgeplots + KEGG Graphs
process_files(input_directory = crypt_path, output_directory = crp_out_path, p_val = 0.05, lfc = 0, export_pathway_files = T, cats = cats)

# ISC details
pathway_list = c('mmu04110', 'mmu03030', 'mmu03050')
dataset = df_isc
cell_type = 'ISC'

# # ENP details
# pathway_list = c('mmu04110', 'mmu03050')
# dataset = df_enp
# cell_type = 'ENP'

# ENT details
pathway_list = c('mmu03050', 'mmu03320', 'mmu04974', 'mmu04973', 'mmu01040', 'mmu00071', 'mmu04146', 'mmu04810')
dataset = df_ent
cell_type = 'ENT'

directory = "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/downstream_analysis/KEGG_Results/crypt/Individual_Pathway_Analysis"

# GOB details
pathway_list = c('mmu04612', 'mmu03320')
dataset = df_gob
cell_type = 'GOB'


# Iterate through pathways for cell type
for (pathway in pathway_list) {
  
  # Construct the folder path
  folder_path = paste0(directory, "/", get_kegg_pathway_name(pathway), '_', cell_type)
  
  # Create the directory
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  # Individual specific pathway analysis - Gene Heatmap
  pathway_heatmap(df_list, 
                  pid = pathway, 
                  scale_to_one = TRUE, 
                  remove_na_rows = TRUE, 
                  order_by_sum = TRUE,
                  output_dir = folder_path,
                  file_name = "heatmap.png")
  
  # Individual Pathway Analysis - KEGG Graph
  specific_pathway_analysis(dataset, output_directory = folder_path, pid = pathway, p_val = 0.05)
  
}
