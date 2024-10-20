
GO_up <- read.csv(file = "../output/GO_UP.csv")  


# 不分GO类别，画全部下调的（前40个）  
goUP_sorted <- GO_up %>%  
  mutate(neg_log_padj = -log10(p.adjust)) %>%  
  arrange(desc(neg_log_padj)) %>%  
  head(40)  

tiff("../output/UP_GO_all.tif", res = 600, compression = "lzw", width = 10, height = 10, units = "in")  
ggplot(goUP_sorted, aes(x = neg_log_padj, y = reorder(Description, neg_log_padj))) +  
  geom_bar(stat = "identity", fill = "red2") +  #下调基因对应Term将fill改为 royalblue
  labs(title = "UPregulated GO Terms",  
       x = "-log10(p.adjust)",  
       y = "Description") +  
  theme_minimal() + theme(panel.grid = element_blank())  
dev.off()

