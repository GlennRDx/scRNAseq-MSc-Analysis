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

directory = "/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/5. results_repository/6. pathway_heatmaps"

# GOB details
pathway_list = c('mmu04612', 'mmu03320')
dataset = spy_gob
cell_type = 'GOB'

# Macronutrient Metabolism Analysis
pathway_list = c('mmu04910', 'mmu04973', 'mmu01040', 'mmu00910')
dataset = spy_ent
cell_type = 'Macronutrient_Metabolism'

# Protein Folding Analysis
pathway_list = c('mmu03050', 'mmu04141')
dataset = df_enp
cell_type = 'Protein_Folding'

# Barrier Function Analysis
pathway_list = c('mmu04520', 'mmu04530', 'mmu04540')
dataset = spy_ent
cell_type = 'Barrier_Function'

# Inflammation Analysis
pathway_list = c('mmu04620', 'mmu04060', 'mmu04064')
dataset = spy_enp
cell_type = 'Inflammation'

# Iterate through pathways for cell type
for (pathway in pathway_list) {
  
  # Construct the folder path
  folder_path = paste0(directory, "/", get_kegg_pathway_name(pathway), '_', cell_type)

  # Create the directory
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  # Individual specific pathway analysis - Gene Heatmap
  pathway_heatmap(spy_list, 
                  pid = pathway, 
                  scale_to_one = TRUE, 
                  remove_na_rows = TRUE, 
                  order_by_sum = TRUE,
                  output_dir = folder_path,
                  file_name = "heatmap1.png")
  
  # Individual Pathway Analysis - KEGG Graph
  specific_pathway_analysis(dataset, output_directory = folder_path, pid = pathway, p_val = 0.05)
  
}





desmosome_genes <- c("Atp2a2", "Cav1", "Cdh1", "Ckap4", "Ctnnd1", "Dsc2", "Dsc3", "Dsg1", "Dsg2", "Dsg3", "Dsp", "Egfr", "Jup", "Kif2a", "Ktn1", "Pkp1", "Pkp2", "Pkp3", "Pkp4", "Prkca", "Sfn")
pathway_heatmap(spy_list, pid = pid, custom_gene_list = desmosome_genes, scale_to_one = T, remove_na_rows = T, order_by_sum = T)
