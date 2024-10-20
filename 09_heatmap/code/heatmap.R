library(tidyr)
library(ComplexHeatmap)

#导入基因表达矩阵,以TPM为例
TPM <- as.data.frame(read_csv("../input/TPM.csv"))

#rownames(TPM) <- TPM$`Gene Name`  #这段代码可以加入行名

#构造矩阵
TPM <- TPM[,-1]
TPM_mat<- as.matrix(TPM)

##如果需要比较统一基因在不同样本中的表达，需要对行进行标准化
TPM_mat_scale = t(scale(t(TPM_mat)))

#作图,根据要求对heatmap进行定制
#svg(file="heatmap.svg", width=10, height=65)#矢量图
tiff("../output/heatmap.tif", res = 600, compression = "lzw", width = 10, height = 10, units = "in")  

Heatmap(TPM_mat_scale,
        name = 'row_scaled_data', 
        cluster_rows = F, 
        cluster_columns = T )
dev.off()


