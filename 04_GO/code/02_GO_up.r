library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)
library(org.Mm.eg.db)
library(biomaRt)

#导入需要进行GO的genes
sig.gene<-read.csv(file="../input/P16_DEG_up.csv")

#如果对应的基因没有ENTREZID，要先转换为ENTREZID，这里我们使用ensembl_gene_id转ENTREZID
#先去除ensembl_gene_id中的版本号
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

GO_up<-enrichGO(gene       = entrez_ids,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)


#画图，这里提供了两个例子，根据需要调整
#矢量图、所有部分类别
svg(file="../output/GO_up.svg", width=6, height=10)
barplot(GO_up, showCategory=20, title=expression(paste(bold("upregulated genes in "), bolditalic("Hsf5"), bold(" KO"))))+facet_grid(ONTOLOGY~., scale="free")
dev.off()

#tif、MF
tiff("../output/GO_up_MF.tif",res = 600, compression = "lzw", width=6, height=10, units="in")
barplot(GO_up, font.size=9.5, showCategory=20,ONTOLOGY="MF",title=expression(paste(bold("upregulated genes in "), bolditalic("Hsf5"), bold(" KO"))))+facet_grid(ONTOLOGY~., scale="free")
dev.off()

#输出结果供画图和写文章使用s
#这行代码可以把GO结果中的ENTREZIDmap成Gene name
EGO_up <- setReadable(GO_up, 'org.Mm.eg.db', 'ENTREZID')

write.csv(EGO_up@result, file = "../output/GO_up.csv")
