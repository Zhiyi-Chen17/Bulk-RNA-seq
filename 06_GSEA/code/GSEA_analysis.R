library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(tidyverse)

#定义我们进行GSEA分析需要使用的基因集，这里我们使用了几个自定义的基因集
chosed_gene_set <- read.csv(file = "../input/gsea_hsf5_gene_se.csv")
long_data <- gather(chosed_gene_set, key = 'variable', value = 'value')
chosed_genes<-na.omit(long_data)
colnames(chosed_genes) <- c("term","gene")

#导入差异分析的结果，包含基因名和FC
sig.gene<-read.csv(file="../input/P16_DEG.csv")

# we want the log2 fold change 
original_gene_list <- sig.gene$log2FoldChange

# name the vector
names(original_gene_list) <- sig.gene$symbol

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#进行GSEA分析，注意调整参数
P16_gsea <- GSEA(gene_list, TERM2GENE = chosed_genes,pvalueCutoff = 100, maxGSSize = 50000,nPermSimple=100000,
        
                          pAdjustMethod = "BH")

#这是一个输出图片的示例，获得P16_gsea对象后，可以使用GSEA_PICTURES.Rmd打印所有需要的图片
tiff("HSF5_P16_GSEA_all_noTABLE.tif",res = 600, compression = "lzw", width=14, height=8, units="in")
gseaplot2(P16_gsea,
          title = "HSF5_P16",  
          geneSetID = 1:7, #选择要绘制的基因集
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, #展示的内容
          pvalue_table = F) # 显示p值
dev.off()

#输出GSEA分析结果
p16_gsea_df = P16_gsea@result
write.csv(p16_gsea_df,file="../output/P16_gsea.csv",row.names=F)


