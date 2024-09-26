# # Limma DGE datasets
# lma_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EEC.csv")
# lma_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_EE_progenitor.csv")
# lma_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte.csv")
# lma_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Enterocyte_progen.csv")
# lma_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_cell.csv")
# lma_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Goblet_progenitor.csv")
# lma_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_ISC.csv")
# lma_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_cell.csv")
# lma_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Paneth_progenitor.csv")
# lma_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_cell.csv")
# lma_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_limma/diff_exp_CD_vs_HFD_Tuft_progenitor.csv")
# 
# lma_list = list(lma_isc = lma_isc,
#                lma_enp = lma_enp,
#                lma_ent = lma_ent,
#                lma_gob = lma_gob,
#                lma_tuf = lma_tuf,
#                lma_eec = lma_eec)

# DESeq2 DGE datasets
dsq_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_EEC.csv")
dsq_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_EEC_Progenitor.csv")
dsq_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Enterocyte.csv")
dsq_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Enterocyte_Progenitor.csv")
dsq_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Goblet.csv")
dsq_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Goblet_Progenitor.csv")
dsq_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_ISC.csv")
dsq_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Paneth.csv")
dsq_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Paneth_Progenitor.csv")
dsq_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Tuft.csv")
dsq_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Tuft_Progenitor.csv")

dsq_list = list(dsq_isc = dsq_isc,
               dsq_enp = dsq_enp,
               dsq_ent = dsq_ent,
               dsq_gob = dsq_gob,
               dsq_tuf = dsq_tuf,
               dsq_eec = dsq_eec)

# Scanpy DGE datasets
# spy_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_EEC.csv")
# spy_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_EEC Progenitor.csv")
# spy_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Enterocyte.csv")
# spy_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Enterocyte Progenitor.csv")
# spy_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Goblet.csv")
# spy_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Goblet Progenitor.csv")
# spy_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_ISC.csv")
# spy_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Paneth.csv")
# spy_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Paneth Progenitor.csv")
# spy_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Tuft.csv")
# spy_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/diff_exp_CD_vs_HFD_Tuft Progenitor.csv")
# 
# spy_list = list(spy_isc = spy_isc,
#                spy_enp = spy_enp,
#                spy_ent = spy_ent,
#                spy_gob = spy_gob,
#                spy_tuf = spy_tuf,
#                spy_eec = spy_eec)

# spy_eil = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/EEC_subpopulations/diff_exp_CD_vs_HFD_EEC (I_L).csv")
# spy_eex = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/EEC_subpopulations/diff_exp_CD_vs_HFD_EEC (X).csv")
# spy_eed = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/EEC_subpopulations/diff_exp_CD_vs_HFD_EEC (Delta).csv")
# spy_ech = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/EEC_subpopulations/diff_exp_CD_vs_HFD_EEC (EC).csv")
# spy_eek = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_scanpy/EEC_subpopulations/diff_exp_CD_vs_HFD_EEC (K).csv")
# 
# spy_eec_list = list(spy_eil = spy_eil,
#                     spy_eex = spy_eex,
#                     spy_eed = spy_eed,
#                     spy_ech = spy_ech,
#                     spy_eek = spy_eek)
