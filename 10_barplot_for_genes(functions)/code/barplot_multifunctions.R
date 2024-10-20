library(ggplot2)
#读取数据
library(readxl)

#从差异分析结果中提取一个多个子集，并作log2FC排序，画barplot
# 绘图
plot_category <- function(df, category_name) {      
  p <- ggplot(df, aes(x = symbol, y = log2FoldChange, fill = -log10(padj))) +      
    geom_col(position = "dodge") +  # 使用dodge避免重叠，如果symbol有重复的话  
    scale_fill_gradient(low = "royalblue", high = "red2") +      
    labs(x = "symbol", y = "log2FoldChange", fill = "-log10(p.adj)", title = category_name) + 
    # scale_x_discrete(position = "top")  +
    theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75,size=12)) + 
    geom_hline(yintercept = 0.6, color = "black", linetype = "dashed", size = 0.5)+
    geom_hline(yintercept = -0.6, color = "black", linetype = "dashed", size = 0.5)
  
  return(p)      
}


#导入差异分析结果和目的基因子集
P16_DEG <- read.csv("../input/P16_DEG.csv")
functions_and_genes <- read_excel("../input/hsf5_barplot_geneset.xlsx", 
                                  sheet = "Sheet1", col_names = T)



#############
spermatogonia_id <- na.omit(unique(functions_and_genes$spermatogonia))
meiosisI_id <- na.omit(unique(functions_and_genes$`meiosis prophase I`))
Sperm_dev_id <- na.omit(unique(functions_and_genes$`spermatid development`))
Sertoli_cell_id <- na.omit(unique(functions_and_genes$`Sertoli_cell_GOtermSertoli cell development`))


#提取基因信息
subset_spermatogonia <- P16_DEG[(P16_DEG$symbol) %in% spermatogonia_id,]
subset_meiosisI <- P16_DEG[(P16_DEG$symbol) %in% meiosisI_id,]
subset_Sperm_dev <- P16_DEG[(P16_DEG$symbol) %in% Sperm_dev_id,]
subset_Sertoli_cell <- P16_DEG[(P16_DEG$symbol) %in% Sertoli_cell_id,]



#整合数据到一张图未使用
# 给每个数据框添加类别列
subset_spermatogonia$category <- "spermatogonia"
subset_meiosisI$category <- "meiosisI"
subset_Sperm_dev$category <- "Sperm_dev"
subset_Sertoli_cell$category <- "Sertoli_cell"



# 合并数据框
df <- rbind(subset_spermatogonia, subset_meiosisI, subset_Sperm_dev, subset_Sertoli_cell)


#按类及log2FC的排序来排
subset_spermatogonia <- subset_spermatogonia[order(subset_spermatogonia$log2FoldChange), ]  
subset_meiosisI <- subset_meiosisI[order(subset_meiosisI$log2FoldChange), ]  
subset_Sperm_dev <- subset_Sperm_dev[order(subset_Sperm_dev$log2FoldChange), ]  
subset_Sertoli_cell <- subset_Sertoli_cell[order(subset_Sertoli_cell$log2FoldChange), ]  

#将symbol数据类型转化为因子
df$symbol = factor(df$symbol, levels = c(subset_spermatogonia$symbol,subset_meiosisI$symbol,
                                         subset_Sperm_dev$symbol,subset_Sertoli_cell$symbol
))

                                         
#输出基因排序，方便后续作图
gene_order <- c("spermatogonia",subset_spermatogonia$symbol,
       "meiosisI",subset_meiosisI$symbol,"spermatid development",subset_Sperm_dev$symbol,"Sertoli_cell",subset_Sertoli_cell$symbol)
write.table(gene_order,file = "../output/gene_order.txt")


#画图
p1 <- plot_category(df, "sperm")
#svg(file="Sperm_dev_90c.svg", width=12, height=7.5)
tiff("../output/Sperm_dev_90C.tif",res = 600, compression = "lzw", width=12, height=7.5, units="in")
p1
dev.off()



