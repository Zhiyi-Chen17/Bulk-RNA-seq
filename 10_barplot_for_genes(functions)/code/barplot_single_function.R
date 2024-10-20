library(ggplot2)
#读取数据
library(readxl)
#从差异分析结果中提取一个基因子集，并作log2FC排序
plot_category <- function(df, category_name) {    
  p <- ggplot(df, aes(x = log2FoldChange, y = reorder(symbol, -log2FoldChange), fill = -log10(padj))) +    
    geom_col() +    
    scale_fill_gradient(low = "royalblue", high = "red2") +    
    labs(x = "log2FoldChange", y = "symbol", fill = "-log10(p.adj)", title = category_name)+    
    theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(position = "top")+geom_vline(xintercept = 0)+geom_vline(xintercept = -0.6, color = "black", linetype = "dashed", size = 0.4)
  return(p)    
}

#导入差异分析结果和目的基因子集
P16_DEG <- read.csv("../input/P16_DEG.csv")
functions_and_genes <- read_excel("../input/hsf5_barplot_geneset.xlsx", 
                                  sheet = "Sheet1", col_names = T)



#example1 作piRNA_process图
piRNA_PROCESS_ID <- na.omit(unique(functions_and_genes$`piRNA process (GO+ARTICLE)`))
subset_piRNA_PROCESS <- P16_DEG[(P16_DEG$symbol) %in% piRNA_PROCESS_ID,]
p <- plot_category(subset_piRNA_PROCESS, "piRNA_process")
tiff("../output/piRNA_process.tif",res = 600, compression = "lzw", width=8, height=5, units="in")
p
dev.off()



#example2 作MSCI图
MSCI_id <- na.omit(unique(functions_and_genes$`ARTICLE:Meiotic Silencing in Mammals`))
subset_MSCI <- P16_DEG[(P16_DEG$symbol) %in% MSCI_id,]
p <- plot_category(subset_MSCI, "MSCI")
tiff("../output/subset_MSCI.tif",res = 600, compression = "lzw", width=8, height=5, units="in")
p
dev.off()


