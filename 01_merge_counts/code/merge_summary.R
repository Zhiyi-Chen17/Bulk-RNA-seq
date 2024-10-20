library(dplyr)  

#输入样本名称以构成完整的文件名
samples <- c('SC_Ythdc2_iKO_1','SC_Ythdc2_iKO_2','SC_Ythdc2_iKO_3','SC_Ythdc2_HT_1','SC_Ythdc2_HT_2','SC_Ythdc2_HT_3') # 假设有更多的样本  

# 初始化一个空列表来保存每个样本的数据
result_list <- list() 

# 使用for循环遍历samples向量中的每个字符串  
for (i in samples) {  
  
  # 使用paste0函数构建文件名  
  file_name <- paste0('../input/',i, ".summary.txt")  
  
  # 读取文件并选择第1列ensembl_id和第7列count数
  data <- read.table(file_name, sep="", header=T) %>%   
    dplyr::select(1, 7)  
  #更改colname,让所有数据可以合并
  colnames(data) <- c("gene_id", paste0(i))
  # 将结果添加到列表中  
  result_list[[length(result_list) + 1]] <- data  
}  

#按gene_id列将所有count合并成一个数据框
combined_data <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), result_list)  

# 查看结果  
head(combined_data)

#输出count martix
write.csv(combined_data,file="../output/gene_count_matrix.csv",row.names = F)

