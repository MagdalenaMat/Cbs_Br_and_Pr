###August experiment fraction analysis
data = read.csv( "~/Desktop/s7s_Breast/data/subset/August_2016_FFPE_LCM.raw_reads.csv")
data$gene.name = NULL

library(AnnotationDbi)
library(org.Hs.eg.db)

data$symbol = mapIds(org.Hs.eg.db,
                     keys = as.character(data$Ensembl.ID),
                     column = c("SYMBOL"),
                     keytype = "ENSEMBL",
                     multiVals = "first")

data1 = data[complete.cases(data),]
data_agr = aggregate(. ~ symbol, data = data1, max)
rownames(data_agr) = data_agr$symbol
data_agr$symbol = NULL
data_agr$Ensembl.ID = NULL

to_plot = data.frame(na = matrix(NA, dim(data_agr)[1],1))
for(i in c("normal","EN","DCIS","IDC","LN")){
  to_plot = cbind(to_plot, data_agr[,grepl(i,colnames(data_agr))])
}
to_plot$na = NULL
colnames(to_plot)
rownames(to_plot) = rownames(data_agr)

desine = data.frame(type = c(rep("normal",24), rep("EN",8), rep("DCIS",16), rep("IDC",26),rep("LN",16)))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = to_plot,
                              colData = desine,
                              design= ~ type)
vsd = vst(dds,blind=FALSE)
PCA = plotPCA(vsd,"type")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p

dds <- estimateSizeFactors(dds)
norm_data = counts(dds, normalized = TRUE)
colnames(norm_data)

to_fractios = data.frame(names = rownames(norm_data), norm_data)
write.table(to_fractios, "~/Desktop/s7s_Breast/to_fractions/to_fractions_august.txt", sep = "\t", row.names = F, col.names = T, quote = F)

###fraction analysis

fractions = read.csv("~/Desktop/s7s_Breast/to_fractions/august_fractions.csv")
rownames(fractions) = fractions$Mixture
fractions$Mixture = NULL
type = as.factor(c(rep("normal",24), rep("EN",8), rep("DCIS",16), rep("IDC",26),rep("LN",16)))

to_plot = as.matrix(fractions[1:22])
LN = to_plot[75:90,]
library(pheatmap)
library(RColorBrewer)

col1 <- brewer.pal(8, "Set2")


mydf <- data.frame(row.names = rownames(to_plot), category = type)

# add row annotations
myBreaks = unique(c(seq(0, 0.19, length=49), 0.2, seq(0.21,max(to_plot), length=50)))
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
         annotation_row = mydf,gaps_row = c(24,32,48,74),
         breaks = myBreaks, show_rownames = F)

