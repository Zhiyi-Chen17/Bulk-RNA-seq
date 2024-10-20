if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')
#? https://github.com/kevinblighe/EnhancedVolcano
library(EnhancedVolcano)
res <- read.csv("../input/KO_VS_WT_DEG.csv", header = T)
dup_rows <- duplicated(res$ensembl_gene_id)
# 输出重复�
res[dup_rows, ]
head(res)

#定义配色方案
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('grey30', nrow(res))

# set the base name/label as 'Mid'
names(keyvals) <- rep('None', nrow(res))

# modify keyvals for transcripts with fold change > 1
keyvals[which(res$log2FoldChange > 0.6 & res$padj<=0.05)] <- 'red2'
names(keyvals)[which(res$log2FoldChange > 0.6 & res$padj<=0.05)] <- 'Log2FC>=0.6 & padj<=0.05'

# modify keyvals for transcripts with fold change < -1
keyvals[which(res$log2FoldChange < -0.6& res$padj<=0.05)] <- 'royalblue'
names(keyvals)[which(res$log2FoldChange < -0.6& res$padj<=0.05)] <- 'Log2FC<=-0.6 & padj<=0.05'

unique(names(keyvals))
unique(keyvals)

#绘图
help(EnhancedVolcano)

caption <- "1066 downregulated genes               84 upregulated genes"
#svg(file="../output/Volcano1.svg", width=8, height=10) #这段代码可以输出矢量图
tiff(file="../output/Volcano1.tif",res = 600, compression = "lzw", width=8, height=10, units="in")
EnhancedVolcano(res,
                lab = res$mgi_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'Fold change'), 
                ylab = bquote('-'~Log[10]~ '(padj)'),
                title = ~italic("YTHDC2")~'_iKO vs WT',
                subtitle = "Differential expression",
                #legendLabels = c('NS',bquote('abs('~Log[2]~ 'FC)>=1'),'padj<=0.05',
                              # bquote('abs('~Log[2]~ 'FC)>=1 & padj<=0.05')),
                legendLabSize = 12,
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 0.6,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                pointSize = 2.0,
                caption = '',
                labSize=3,
                xlim = c(-7.5, 7.5),
                ylim = c(-1, 20),
                #selectLab = c('Ddx4','Stk31','Tdrd1','Piwil1','Tdrd12','Piwil2','Tdrd7','Tdrd6','Mael','Pld6','Henmt1','Tdrd9','Tdrd5','Hsf5','Bbc3','Btg2','Vegfa','Hrk','Adnp','Jun','Phlda3','Bcl3'),
                selectLab = c(""),
            
                drawConnectors = TRUE) + 
                annotate(
                  geom = "text", x = -5, y = 20, 
                  label = caption, hjust = 0, vjust = 1, size = 5.2
                )

dev.off()

