
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

x = c('TAGLN', 'CSTA', 'HK2', 'H2AFZ', 'IGHG4', 'BAG1', 'SLC39A11', 'SREBF1', 'FARP1', 'TNFSF10',
      'COX6C', 'TTLL12', 'IFT122', 'DNAJB1', 'AL391001.1', 'RPL17', 'AGPAT2', 'COL4A1', 'IGFBP3',
      'PRR15', 'GPX3', 'CRNDE', 'WARS', 'NDRG1', 'MUC1', 'DCN', 'PVALB', 'CCL19', 'COL1A1', 'PERP',
      'LUM', 'CCND1', 'TMEM258', 'TOP2A', 'COL2A1', 'UBE2S', 'TSPAN1', 'IGHM', 'DNAJC1', 'S100G',
      'KPNA2', 'PFKFB3', 'APOD', 'ABHD2', 'TIMP3', 'CGA', 'FAM162A', 'SIVA1', 'TCIM', 'CNOT6',
      'SCGB1D2', 'CEBPD', 'CXCL9', 'DEGS1', 'REPS2', 'B4GALT1', 'TIMP1', 'PTGDS', 'COL4A2', 'COMMD3',
      'BMI1', 'SRD5A3', 'CERS6', 'FOS', 'BNIP3', 'IFIT1', 'THY1', 'ALB', 'S100P', 'APOL1', 'C6orf141',
      'HES1', 'ING1', 'DRAIC', 'CACYBP', 'IGKC', 'CTPS2', 'S100A9', 'MRPL45', 'RPS14', 'MMP11', 'FGFR1',
      'GFRA1', 'EXOC2', 'ISLR', 'TFF1', 'SERPINA3', 'CRISP3', 'GDF15', 'TAP1', 'COL18A1', 'SERHL2',
      'HMGB2', 'CPNE7', 'KCTD3', 'ITPR1', 'HLA-DQA1', 'THBS1', 'THBS2', 'SECTM1', 'MMP9', 'LGALS8',
      'VWF', 'IL32', 'NR4A1', 'HTRA1', 'CST1', 'IGHG1', 'SYAP1', 'RNASE1', 'TUBA1C', 'ACADSB', 'MZB1',
      'RPL41', 'HSPH1', 'COL5A1', 'CXCL14', 'CBWD5', 'COMP', 'FN1', 'RAMP1', 'SNRPF', 'TFF3', 'RAB30',
      'CKS2', 'AC087379.2', 'AQP3', 'SOCS2', 'NPEPPS', 'G3BP1', 'COL6A3', 'DEGS2', 'ADIPOR2', 'SHISA2',
      'EFHD1', 'KRT37', 'PLXNB1', 'POSTN', 'HEBP2', 'TRMT112', 'COL3A1', 'HEBP1', 'EEF1A2', 'CKB', 'YPEL3',
      'CPB1', 'VTCN1', 'NDUFC1', 'CISH', 'RERG', 'SPP1', 'IGHA1', 'IFITM1', 'IGLC2', 'MT1X', 'SFRP2',
      'TSPAN13', 'GNG5', 'AGR3', 'SLITRK6', 'MCCC2', 'AMIGO2', 'TRIB1', 'MT1E', 'SCGB2A2', 'KIAA1324',
      'CCN2', 'SNCG', 'IGLC1', 'MGST1', 'TACSTD2', 'AGR2', 'MUC5B', 'S100A6', 'ATP9A', 'COL1A2', 'ZFAS1',
      'UGDH', 'KRT17', 'LONP2', 'IGHG3', 'MESP1', 'FCGR3A', 'TGM2', 'LINC02224', 'NEBL', 'COL5A2', 'KLHDC7B',
      'MCCD1', 'IGLC3', 'IGHG2', 'TRAC', 'TGFBI', 'TMEM165', 'COL15A1', 'SERPINH1', 'ARMT1', 'SULF1')



trans = bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

kk <- enrichKEGG(gene = trans$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)


dotplot(kk, showCategory = 20)
ggsave("KEGG_enriched.pdf", width = 6.5, height = 6, dpi=300)
# barplot(kk, showCategory = 10)
#
# cnetplot(kk, showCategory = 5)
#
# head(kk)
#
# z = setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# head(z@result)
#
# dotplot(z, showCategory = 10)
#
# barplot(z, showCategory = 10)
#
# cnetplot(z, showCategory = 5)

result <- enrichGO(gene = trans$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)

dotplot(result,showCategory = 20)
ggsave("GO enriched.pdf", width = 6.5, height = 10, dpi=300)


