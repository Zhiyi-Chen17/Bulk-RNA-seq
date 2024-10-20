library(tidyr)

library(readxl)
library(ComplexHeatmap)
mmc2_2_ <- read_excel("mmc2 (2).xlsx")
chosed_gene_set <- read_excel("chosed_gene_set.xlsx", 
                              sheet = "Sheet2")



#############################
undif <- mmc2_2_[(mmc2_2_$gene)%in%chosed_gene_set$`spermatogonia`,]
undif  <- undif[,c(8,7,3)]
undif_wide<-spread(undif,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值undif_c1 <- undif[(undif$cluster)%in% "6",]
undif_wide <- as.data.frame(undif_wide)
rownames(undif_wide) <- undif_wide$gene
undif_wide <- undif_wide[,-1]


#############################
meio <- mmc2_2_[(mmc2_2_$gene)%in%chosed_gene_set$`meiosis prophase I`,]

meio  <- meio[,c(8,7,3)]
meio_wide<-spread(meio,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值meio_c1 <- meio[(meio$cluster)%in% "6",]
meio_wide <- as.data.frame(meio_wide)
rownames(meio_wide) <- meio_wide$gene
meio_wide <- meio_wide[,-1]



#############################
sperm <- mmc2_2_[(mmc2_2_$gene)%in%chosed_gene_set$`sperm_assembly`,]
#sperm$legend <- c(rep("sperm",601))

sperm  <- sperm[,c(8,7,3)]
sperm_wide<-spread(sperm,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值sperm_c1 <- sperm[(sperm$cluster)%in% "6",]
sperm_wide <- as.data.frame(sperm_wide)
rownames(sperm_wide) <- sperm_wide$gene
sperm_wide <- sperm_wide[,-1]


all_col_sperm <- sperm_wide
all_col_meio <- meio_wide
all_col_undif <- undif_wide


#######################
all_col_undif$'1' <-  c(rep(0,170))
all_col_undif$'2' <-  c(rep(0,170))
all_col_undif$'3' <-  c(rep(0,170))
all_col_undif$'4' <-  c(rep(0,170))
all_col_undif$'5' <-  c(rep(0,170))
all_col_undif$'6' <-  c(rep(0,170))
all_col_undif$'7' <-  c(rep(0,170))
all_col_undif$'9' <-  c(rep(0,170))
all_col_undif$'11' <-  c(rep(0,170))
all_col_undif$'12' <-  c(rep(0,170))
all_col_undif$'13' <-  c(rep(0,170))
all_col_undif$'14' <-  c(rep(0,170))




#####

all_col_meio$'1' <-  c(rep(0,66))
all_col_meio$'2' <-  c(rep(0,66))
all_col_meio$'3' <-  c(rep(0,66))
all_col_meio$'4' <-  c(rep(0,66))
all_col_meio$'5' <-  c(rep(0,66))
all_col_meio$'7' <-  c(rep(0,66))
all_col_meio$'11' <-  c(rep(0,66))
all_col_meio$'13' <-  c(rep(0,66))
all_col_meio$'14' <-  c(rep(0,66))
###
all_col_sperm$'15' <-  c(rep(0,152))
all_col_sperm$'16' <-  c(rep(0,152))
all_col_sperm$'17' <-  c(rep(0,152))



##################



d<-rbind(all_col_undif,all_col_meio)

all <-rbind(d,all_col_sperm)

all_mat<- as.matrix(all)
all_mat_scale = t(scale(t(all_mat)))
#############



column_order = c("16","10","8","6","9","12","2","7","13","11","5","3","1","4","14","15","17")
#
svg(file="confirm2.svg", width=10, height=65)

Heatmap(all_mat_scale, row_split = c(rep('A# spermatogonia', 170), rep('B# meiosis prophase I', 66),rep('C# spermatid development',152)),
        name = 'row_scaled_avg_logFC', 
        column_order = column_order,
        cluster_rows = F, 
        cluster_columns = F )
dev.off()
#
rownames(all_mat_scale) <- NULL

tiff("confirm1.tif", res = 600, compression = "lzw", width = 10, height = 10, units = "in")  

Heatmap(all_mat_scale, row_split = c(rep('A# spermatogonia', 170), rep('B# meiosis prophase I', 66),rep('C# spermatid development',152)),
        name = 'row_scaled_avg_logFC', 
        column_order = column_order,
        cluster_rows = F, 
        cluster_columns = F )
dev.off()


