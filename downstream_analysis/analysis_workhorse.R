# Process all DGE tables - Ridgeplots + KEGG Graphs
process_files(input_directory = crypt_path, output_directory = crp_out_path, p_val = 0.05, lfc = 0, export_pathway_files = T, cats = cats)

# Pathway
pathway = 'mmu04630'

# Individual specific pathway analysis - Gene Heatmap
pathway_heatmap(df_list, pid = pathway, scale_to_one = T, remove_na_rows = F, order_by_sum = T)

# Individual Pathway Analysis - KEGG Graph
specific_pathway_analysis(df_isc, pid = pathway, p_val = 0.05)