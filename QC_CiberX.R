## Qc dor CibersortX deconvolution

data = read.csv( "~/Desktop/s7s_Breast/data/subset/2017-05-24_single-cell_Smart-3SEQ.raw_reads.csv")
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
data_agr = data_agr[,!grepl("control|ablation",colnames(data_agr))]
colnames(data_agr)
desine = data.frame(type = c(rep("bulk_DCIS",6), rep("single_DCIS",10), rep("bulk_MO",6), rep("single_MO",10),rep("mix",2)))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = data_agr,
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

to_fractios = norm_data[,!grepl("single",colnames(norm_data))] 
colnames(to_fractios)
to_fractios = data.frame(names = rownames(to_fractios), to_fractios)

###fraction analysis
write.table(to_fractios, "~/Desktop/s7s_Breast/to_fractions/to_fractions.txt", sep = "\t", row.names = F, col.names = T, quote = F)
fractions = read.csv("~/Desktop/s7s_Breast/to_fractions/bulk_DCIS_Macrophage_mix.csv")
rownames(fractions) = fractions$Mixture
fractions$Mixture = NULL
type = as.factor(c(rep("bulk_DCIS",6),rep("bulk_MO",6),rep("mix",2)))

to_plot = as.matrix(fractions[1:22])

wLN = rbind(to_plot, LN)
type_w_LN = as.factor(c(rep("bulk_DCIS",6),rep("bulk_MO",6),rep("mix",2),rep("LN",16)))

library(pheatmap)
library(RColorBrewer)
col1 <- brewer.pal(8, "Set2")

mydf <- data.frame(row.names = rownames(to_plot), category = type)
mydfLN = data.frame(row.names = rownames(wLN), category = type_w_LN)

# add row annotations
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, annotation_row = mydf,gaps_row = c(6, 12) )
pheatmap(wLN, cluster_cols = F, cluster_rows = F, annotation_row = mydfLN,gaps_row = c(6, 12, 14) )

######


mix_sml = norm_data[,33:34]
bulk_MO = norm_data[,17:22]
single_MO = norm_data[,23:32]
bulk_DCIS = norm_data[,1:6] 

setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/Filtered/")
DCIS = read.table("CIBERSORTxGEP_20180314_not_batch_corrected_DCIS.txt_GEPs_Filtered.txt", sep = "\t", header = TRUE)
DCIS_MO = DCIS[, c(1,6)]

DCIS_MO1 = DCIS_MO[complete.cases(DCIS_MO),]
DCIS_MO1$GeneSymbol = sub("[.]","-",DCIS_MO1$GeneSymbol)
all(DCIS_MO1$GeneSymbol %in% rownames(mix_sml))

mix_sml = mix_sml[rownames(mix_sml) %in% DCIS_MO1$GeneSymbol,]
all(DCIS_MO1$GeneSymbol == rownames(mix_sml))

bulk_MO = bulk_MO[rownames(bulk_MO) %in% DCIS_MO1$GeneSymbol,]
single_MO = single_MO[rownames(single_MO) %in% DCIS_MO1$GeneSymbol,]
bulk_DCIS = bulk_DCIS[rownames(bulk_DCIS) %in% DCIS_MO1$GeneSymbol,]

cor(DCIS_MO1$Monocytes,bulk_MO,use="na.or.complete",method="spearman")
cor(DCIS_MO1$Monocytes,bulk_MO,use="na.or.complete",method="pearson")

par(mfrow = c(3,3))
for(name in colnames(bulk_DCIS)){
  plot(bulk_DCIS[,name], DCIS_MO1$Monocytes, ylim = c(0,15000), xlim = c(0,15000))
}

par(mfrow = c(3,3))
for(name in colnames(bulk_MO)){
  plot(bulk_MO[,name], DCIS_MO1$Monocytes, ylim = c(0,15000), xlim = c(0,15000))
}

par(mfrow = c(3,3))
for(name in colnames(bulk_DCIS)){
  plot(log2(bulk_DCIS[,name]+1), log2(DCIS_MO1$Monocytes+1), ylim = c(0,14), xlim = c(0,14))
}

par(mfrow = c(3,3))
for(name in colnames(bulk_MO)){
  plot(log2(bulk_MO[,name]+1), log2(DCIS_MO1$Monocytes+1), ylim = c(0,14), xlim = c(0,14))
}
