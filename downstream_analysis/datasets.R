# Limma DGE datasets
df_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EEC.csv")
df_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EE_progenitor.csv")
df_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte.csv")
df_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte_progen.csv")
df_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_cell.csv")
df_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_progenitor.csv")
df_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_ISC.csv")
df_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_cell.csv")
df_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_progenitor.csv")
df_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_cell.csv")
df_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_progenitor.csv")

#Scanpy DGE datasets
keep_first_half <- function(df) {
  rows_to_keep <- floor(nrow(df) / 2)
  return(df[1:rows_to_keep, ])
}

df_eec = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_EEC.csv"))
df_eep = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_EEC Progenitor.csv"))
df_ent = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Enterocyte.csv"))
df_enp = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Enterocyte Progenitor.csv"))
df_gob = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Goblet.csv"))
df_gop = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Goblet Progenitor.csv"))
df_isc = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_ISC.csv"))
df_pan = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Paneth.csv"))
df_pap = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Paneth Progenitor.csv"))
df_tuf = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Tuft.csv"))
df_tup = keep_first_half(read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/differential_expression_cluster_Tuft Progenitor.csv"))

df_list = list(df_isc = df_isc,
               df_ent = df_ent,
               df_gob = df_gob,
               df_tuf = df_tuf,
               df_eec = df_eec)
