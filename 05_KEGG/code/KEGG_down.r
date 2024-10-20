library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)
library(org.Mm.eg.db)
library(biomaRt)

#导入需要进行GO的genes
sig.gene<-read.csv(file="../input/P16_DEG_down.csv")

#如果对应的基因没有ENTREZID，要先转换为ENTREZID，这里我们使用ensembl_gene_id转ENTREZID
head(sig.gene)
gene<-sig.gene$Row.names
gene <- gsub("\\..*", "", gene)
head(gene)
# 使用Biomart获取基因的其他ID,包括ENTREZID，ensembl数据库的版本可以调整
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
options(timeout = 4000000)
option_info<- getBM(
  attributes = c("ensembl_gene_id","mgi_symbol","entrezgene_id"),
  filters = 'ensembl_gene_id',
  values = gene,
  mart = mart)

#输出转换好的entrezgene_id
entrez_ids <- option_info$entrezgene_id


#KEGG分析
KEGG_down <- enrichKEGG(
  gene = entrez_ids,  
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

tiff("HSFko downregulate Enrichment KEGG.tif",res = 600, compression = "lzw", width=6, height=10, units="in")
barplot(KEGG_down ,title="Enrichment KEGG _down")
dev.off()

head(kk,2)
EKEGG_down <- setReadable(KEGG_down , 'org.Mm.eg.db', 'ENTREZID')
write.csv(EKEGG_down@result, file = "KEGG_down.csv")



