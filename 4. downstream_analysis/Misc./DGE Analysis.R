# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

df_ent = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_Enterocyte.csv')
df_isc = read.csv('/home/glennrd/Documents/Research_Project/RNA-seq_Analysis/Data/differential_expression/diff_exp_CD_VS_HFD_ISC.csv')

names(df_ent)[names(df_ent) == 'abs.log2FC'] <- 'log2FoldChange'
names(df_ent)[names(df_ent) == 'P.Value'] <- 'pvalue'

df_ent_top15 = head(df_ent[order(df_ent$pvalue),],10)$index

EnhancedVolcano(df_ent,
                lab = df_ent$index,
                selectLab = df_ent_top15,
                x = 'log2FoldChange',
                y = 'pvalue',
                labSize = 4,
                drawConnectors = T,
                FCcutoff = 0.1,
                title = 'CD vs HFD Enterocyte Differential Expression')

names(df_isc)[names(df_isc) == 'abs.log2FC'] <- 'log2FoldChange'
names(df_isc)[names(df_isc) == 'P.Value'] <- 'pvalue'

df_isc_top15 = head(df_isc[order(df_isc$pvalue),],10)$index

EnhancedVolcano(df_isc,
                lab = df_isc$index,
                selectLab = df_isc_top15,
                x = 'log2FoldChange',
                y = 'pvalue',
                labSize = 4,
                drawConnectors = T,
                FCcutoff = 0.1,
                title = 'CD vs HFD ISC Differential Expression')