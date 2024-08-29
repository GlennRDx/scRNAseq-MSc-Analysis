# # Limma DGE datasets
# df_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EEC.csv")
# df_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EE_progenitor.csv")
# df_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte.csv")
# df_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte_progen.csv")
# df_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_cell.csv")
# df_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_progenitor.csv")
# df_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_ISC.csv")
# df_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_cell.csv")
# df_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_progenitor.csv")
# df_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_cell.csv")
# df_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_progenitor.csv")

# Scanpy DGE datasets
df_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_EEC.csv")
df_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_EEC Progenitor.csv")
df_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Enterocyte.csv")
df_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Enterocyte Progenitor.csv")
df_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Goblet.csv")
df_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Goblet Progenitor.csv")
df_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_ISC.csv")
df_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Paneth.csv")
df_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Paneth Progenitor.csv")
df_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Tuft.csv")
df_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Tuft Progenitor.csv")

df_list = list(df_isc = df_isc,
               df_ent = df_ent,
               df_gob = df_gob,
               df_tuf = df_tuf,
               df_eec = df_eec,
               df_eep = df_eep,
               df_enp = df_enp,
               df_gop = df_gop,
               df_pan = df_pan,
               df_pap = df_pap,
               df_tup = df_tup)

