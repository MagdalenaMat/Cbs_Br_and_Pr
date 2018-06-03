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

to_plot = data_agr

pdata = read.csv("~/Desktop/s7s_Breast/data/pdata/design_August_2016_FFPE_LCM.csv")
pdata = data.frame(pdata, names = paste(pdata$patient,pdata$stage,pdata$biol.rep, sep = "_"))
table(pdata$Batch, pdata$stage)
table(pdata$stage)

pdata = pdata[pdata$names %in% colnames(to_plot),]
colnames(to_plot) = gsub("[.]","-", colnames(to_plot))
pdata = pdata[match(colnames(to_plot),pdata$names),]
colnames(to_plot) == pdata$names
pdata$stage = gsub("-","_", pdata$stage)
pdata$stage = factor(pdata$stage)
pdata$Batch = factor(pdata$Batch)

library(DESeq2)
library(plotly)
dds <- DESeqDataSetFromMatrix(countData = to_plot,
                              colData = pdata,
                              design= ~ stage)
vsd = vst(dds,blind=FALSE)
PCA = plotPCA(vsd,"Batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p

dds <- estimateSizeFactors(dds)
norm_data = counts(dds, normalized = TRUE)

keep = rowSums(norm_data > 1) >= 4
norm_counts = norm_data[keep,]
dim(norm_data)
dim(norm_counts)

pdata1 = pdata[,c(1,3,5)]
table(pdata$Batch,pdata$stage)
table(pdata$stage)
library(sva)
batch = pdata1$Batch
modcombat = model.matrix(~as.factor(stage), data=pdata1) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(norm_counts+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

combat_edata = round(combat_edata)
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData = pdata1,
                              design= ~ stage)
vsd = vst(dds,blind=FALSE)
PCA = plotPCA(vsd,"stage")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p

to_hires = data.frame(names = rownames(combat_edata), combat_edata)
pdata$names = gsub("-","_", pdata$names)
colnames(to_hires) = gsub("[.]","_",colnames(to_hires))
write.table(to_hires, "~/Desktop/s7s_Breast/data/August_batch_corrected_data/August_BC_full.txt", sep = "\t", row.names = F, col.names = T, quote = F)

all(colnames(to_hires)[-1] == pdata$names)
pdata$stage = factor(pdata$stage, levels = c("normal","EN","DCIS","IDC","AVL","met_no_ECE","met_ECE","ECE","LN"))

to_write = list()
for(i in levels(pdata$stage)){
  to_write[[i]] = combat_edata[,pdata1$stage == i]
}

setwd("~/Desktop/s7s_Breast/data/August_batch_corrected_data/")
for(i in names(to_write)){
  to_write[[i]] = data.frame(genenames = rownames(combat_edata), to_write[[i]])
  write.table(to_write[[i]], paste0("August_BC_", i, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}

write.table(pdata1, "~/Desktop/s7s_Breast/data/August_batch_corrected_data/pdata_August.txt", sep = "\t", row.names = F, col.names = T, quote = F)


to_fractios = data.frame(names = rownames(norm_data), norm_data)
write.table(to_fractios, "~/Desktop/s7s_Breast/to_fractions/to_fractions_august.txt", sep = "\t", row.names = F, col.names = T, quote = F)

###fraction analysis

#fractions = read.csv("~/Desktop/s7s_Breast/to_fractions/august_fractions.csv")
fractions = read.table("~/Desktop/s7s_Breast/data/August_batch_corrected_data/absolute_mode_august_full/CIBERSORTx_August_BC_full.txt_Adjusted.txt",
                       header = T, row.names = 1, sep = "\t")

fractions = fractions[,1:22]
pdata = read.table("~/Desktop/s7s_Breast/data/August_batch_corrected_data/pdata_August.txt",
                   header = T, sep = "\t")


all(as.character(pdata$names) == rownames(fractions))

######
sorted_fr = data.frame(matrix(NA, 1, dim(fractions)[2]))
colnames(sorted_fr) = colnames(fractions)
for(i in levels(pdata$stage)){
  sorted_fr = rbind(sorted_fr, fractions[pdata$stage == i,])
}
sorted_fr = sorted_fr[-1,]
pdata_s = pdata[match(rownames(sorted_fr),pdata$names),]
all(rownames(sorted_fr) == pdata_s$names)
######
type = factor(pdata_s$stage)

to_plot = as.matrix(sorted_fr)
LN = to_plot[75:90,]
library(pheatmap)
library(RColorBrewer)

col1 <- brewer.pal(8, "Set2")


mydf <- data.frame(row.names = rownames(to_plot), category = type)

# add row annotations
#myBreaks = unique(c(seq(0, 0.19, length=49), 0.2, seq(0.21,max(to_plot), length=50)))
to_boxplot = to_plot
to_plot[to_plot >= 0.3] = 0.3
#to_plot = to_plot[1:74,]
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
         annotation_row = mydf,gaps_row = c(24,32,48,74,90,104,108,112),
         show_rownames = F)
max(to_plot)

#####boxplots
library(ggpubr)
#to_boxplot = to_boxplot[-c(grep("LN",rownames(to_plot))),]
#type = type[-c(grep("LN",rownames(to_plot)))] 
to_ggplots = data.frame(to_boxplot, type)
#colnames(to_ggplots)[9] = "Treg"
my_comparisons <- list( c("normal","EN"),c("EN", "DCIS"),c("DCIS","IDC"), 
                        c("normal", "DCIS"), c("EN","IDC"), c("normal", "IDC"),
                        c("IDC","AVL"),c("DCIS","AVL"),c("AVL","met_no_ECE"),
                        c("met_no_ECE", "met_ECE"), c("met_ECE","ECE"))
path = "~/Desktop/s7s_Breast/data/August_batch_corrected_data/absolute_mode_august_full/plots/boxplots/"
for( i in colnames(to_ggplots)[1:22]){
  png(paste0(path,i,"400.png"), width = 600, height = 600)
  plot(ggboxplot(to_ggplots, x = "type", y = i,  outlier.colour = NA,
                 color = "type", palette = "jco", add = "jitter") + 
         theme(text = element_text(size=20)) + rotate_x_text() +
         stat_compare_means(comparisons = my_comparisons))
  dev.off()
}

