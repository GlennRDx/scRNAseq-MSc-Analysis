# datasets
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
df_sdo = read.csv("/home/glennrdx/Documents/Research_Project/RNA-seq_Analysis/python/crypt_pb/differential_expression/diff_exp_CD_vs_HFD_Intestinal_epithe.csv")

df_list = list(df_isc = df_isc,
               df_ent = df_ent,
               df_gob = df_gob,
               df_tuf = df_tuf,
               df_eec = df_eec)

