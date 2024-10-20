# 加载必要的包  
library(ggplot2)  
library(dplyr)  
library(openxlsx)
#读取GO结果
GO_up <- read.csv(file = "../output/GO_UP.csv")  

#提取每个类别("BP", "MF", "CC"),单独对p.adjust排序并输出图像
for (i in c("BP", "MF", "CC")) {    
  # 子集数据，排序，并取前10个 
  subset_data <- GO_up %>%    
    filter(ONTOLOGY == i) %>%    
    mutate(neg_log_padj = -log10(p.adjust)) %>%    
    arrange(desc(neg_log_padj)) %>%    
    head(10)    
  
  # 创建并保存条形图  
  tiff(paste0('../output/',"UP_GO_", i, ".tif"), res = 600, compression = "lzw", width = 10, height = 6, units = "in")  
  
  # 使用print()输出
  print(ggplot(subset_data, aes(x = neg_log_padj, y = reorder(Description, neg_log_padj))) +  
          geom_bar(stat = "identity", fill = "red2") +  
          labs(title = paste("UPregulated GO -", i),  
               x = "-log10(p.adjust)",  
               y = "Description")   
        + theme_minimal() 
        + theme_bw()
        + theme(panel.grid = element_blank())
  )  
  dev.off()  
}

