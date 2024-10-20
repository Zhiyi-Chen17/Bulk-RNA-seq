up_data<-read.csv("../input/P16_DEG_up.csv")
down_data<-read.csv("../input/P16_DEG_down.csv")


up_table <- as.data.frame(table(up_data$chrom)) 
down_table <- as.data.frame(table(down_data$chrom)) 

colnames(up_table) <- c("Chromosome","Upregulated_genes")
colnames(down_table) <- c("Chromosome","Downregulated_genes")
data <- merge(up_table,down_table,all=TRUE)


# 已知小鼠所有染色体名称
mouse_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

# 检查并添加缺失的染色体
for (i in mouse_chromosomes) {
  if (!(i %in% data$Chromosome)) {
    new_row <- data.frame(Chromosome = i, Upregulated_genes = 0, Downregulated_genes = 0)
    data <- rbind(data, new_row)
  }
}

data[is.na(data)] <- 0
print(data)

write.csv(data,file="../output/DEGs_count_by_Chr.csv",row.names = F)


#画图
library(ggplot2)
# 在Chromosome列中去除chr字符
data$Chromosome <- gsub("chr", "", data$Chromosome)

# 将 Chromosome 列转换为 factor，并指定顺序
data$Chromosome <- factor(data$Chromosome, levels = paste0(c(1:19, "X", "Y")))

# 将数据从宽格式转换为长格式（tidy data）
data_long <- tidyr::gather(data, key = "Type", value = "Number_of_DEGs", -Chromosome)
# 绘制分组条形图

tiff("../output/DEG-conut-up-vs-down-by-chr.tif",res = 600, compression = "lzw", width=6, height=6, units="in")
p <- ggplot(data_long, aes(x = Chromosome, y = Number_of_DEGs, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Chromosome Frequencies",
       x = "Chromosome",
       y = "Number of DEGs",
       fill = "Legend")+scale_fill_manual(values = c("Downregulated_genes" = "royalblue","Upregulated_genes" = "red2" ))+ theme_bw() + theme(panel.grid=element_blank())+theme(legend.position = c(0.8, 0.85))
# 调整图例位置
# 打印图形
print(p)
dev.off()

