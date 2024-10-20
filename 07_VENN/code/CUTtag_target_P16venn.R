

#确定要使用的集合

#获得P16HSF5_DEG_down_unique_mgi_symbol
P16HSF5_DEG_down <- read.csv("../input/P16_DEG_down_FC0.6.csv")
P16HSF5_DEG_down$Row.names <- gsub("\\..*", "", P16HSF5_DEG_down$Row.names)
head(P16HSF5_DEG_down)

#获取第二个集合
HSF5_CUTTAG_GENE <- read.csv("../input/HSF5_CUTTAG_GENE.CSV")
HSF5_CUTTAG_GENE_unique_ensembl_ids <- unique(na.omit(HSF5_CUTTAG_GENE$ENSEMBL))

head(HSF5_CUTTAG_GENE)


#画VENN图
library(VennDiagram)
venn.diagram(list(P16_Hsf5_KO_DWON =P16HSF5_DEG_down$Row.names,HSF5_CUTTAG = HSF5_CUTTAG_GENE_unique_ensembl_ids),
             fill = c("red2", "royalblue"), 
             cat.default.pos = "text", inverted = TRUE,
             cex = 1.5, 
             filename = "../output/HSF5_CUTTAG_P16_DOWN.tiff")
help(venn.diagram)


#输出VENN图中各部分的信息
#取交集
cuttag_intersect_P16_ensembl_gene_id <- intersect(P16HSF5_DEG_down$Row.names,HSF5_CUTTAG_GENE_unique_ensembl_ids)
cuttag_DEG_down_subset_intersection <- HSF5_CUTTAG_GENE[HSF5_CUTTAG_GENE$ENSEMBL %in% cuttag_intersect_P16_ensembl_gene_id,]
P16HSF5_DEG_down_subset_intersection <- P16HSF5_DEG_down[P16HSF5_DEG_down$Row.names%in% cuttag_intersect_P16_ensembl_gene_id,]

write.csv(cuttag_DEG_down_subset_intersection,file = "../output/cuttag_DEG_down_subset_intersection.csv",row.names = F)
write.csv(P16HSF5_DEG_down_subset_intersection,file = "../output/P16HSF5_DEG_down_subset_intersection(cuttag).csv",row.names = F)


#取补集
P16HSF5_not_cuttag_ensembl_gene_id <- setdiff(P16HSF5_DEG_down$Row.names,HSF5_CUTTAG_GENE_unique_ensembl_ids)
P16HSF5_DEG_down_subset_just_P16HSF5_not_cuttag <- P16HSF5_DEG_down[P16HSF5_DEG_down$Row.names %in% P16HSF5_not_cuttag_ensembl_gene_id,]

write.csv(P16HSF5_DEG_down_subset_just_P16HSF5_not_cuttag,file = "../output/P16HSF5_DEG_down_subset_just_P16HSF5_not_cuttag.csv",row.names = F)

