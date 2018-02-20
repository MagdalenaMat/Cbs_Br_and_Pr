#######Magdalena Matusiak
#######17 feb 2018
#######prostate cancer batches data pooled and analysed by PCA


path = "~/Desktop/s3s_Prostate/data/"
setwd(path)
files <- list.files(path=path)

#check if gene names are same and in the same order in all files  
names = list()
for(i in 1:length(files)){
  print(i)
  namb = as.character(read.csv(files[i])[,2])
  names[[i]] = namb
}
for(i in 1:(length(names)-1) ){
  print(all(names[[i]][1] == names[[c(i+1)]][1]))
}

for (i in 1: length(files)){
  r = read.csv(files[i])
  print(colnames(r))
}

gene.names = as.character(read.csv(files[1])[,1])

names <- data.frame(matrix(gene.names,58037,1))
dim(names)

results = list()
desine = data.frame(name = NA, batch = NA, subtype = NA)
for(j in c("normal|Normal","BPH", "HGPIN","LGPIN","T3$|T3_","T34","T4$|T4_","T45","T5")){
  for (i in 1: length(files)){
    r = read.csv(files[i])
    cols = as.matrix(r[,grepl(j,colnames(r)), drop = FALSE])
    if(ncol(cols) != 0){
      results[[j]] = cbind(results[[j]],cols)
      des = data.frame(name = colnames(cols), batch = files[i], subtype = j)
      desine = rbind(desine, des)
    }
    print(paste(j,files[i],ncol(cols), colnames(cols)))
  }
  results[[j]] = data.frame(GeneNames = gene.names, results[[j]])
  print(ncol(results[[j]]))
}




desine = desine[-1,]
dim(desine)
sum(lengths(results))
desine$subtype = gsub("normal|Normal","normal",desine$subtype)
desine$subtype = gsub("T3$|T3_","T3",desine$subtype)
desine$subtype = gsub("T4$|T4_","T4",desine$subtype)


names(results) = c("normal","BPH","HGPIN","LGPIN", "T3", "T34", "T4", "T45", "T5")
desine$subtype = factor(desine$subtype)
levels(desine$subtype) = c("BPH","HGPIN","LGPIN", "normal", "T3", "T34", "T4", "T45", "T5")

results_df = data.frame(names = results$normal$GeneNames)
for(i in names(results)){
  results_df = cbind(results_df, results[[i]][,c(2:length(colnames(results[[i]])))])
}
dim(results_df)

#deduplicate gene names names, leave the row with the max value per gene  
results_df <- aggregate(. ~ names, data = results_df, max)
rownames(results_df) = results_df$names
results_df$names = NULL

desine$name[!(desine$name == colnames(results_df))] = "V60.T34.1"
all(desine$name == colnames(results_df))

desine$totalCount = colSums(results_df)

desine$size_f = NA

for(i in c(300000,100000,20000)){
  desine[desine$totalCount < i,"size_f"] = paste0("smaller_than_",i)
}  

desine[is.na(desine$size_f), "size_f"] = "bigger_than_300K"

desine$size_f = factor(desine$size_f)
levels(desine$size_f) = c("bigger_than_300K",
                           "between_20K_and_100K",
                           "between_100K_and_300K")
table(desine$size_f)
desine[desine$size_f == "between_20K_and_100K",c("name","batch","totalCount")]

desine$batch = sub("\\_.*","", desine$batch)


library(ggplot2)
setwd("~/Desktop/s3s_Prostate/data/")
png("total_reads.png", width = 1200, height = 600)
ggplot(desine, aes(x = name, y = totalCount, 
                    fill = size_f)) + 
  facet_grid(. ~ batch, labeller = label_value, scales = "free", space = "free") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank()) +
  scale_fill_brewer(palette="Dark2")
dev.off()

sorted.res = results_df[rowSums(results_df) > 0 , ]

library(DESeq2)
library(plotly)
dds <- DESeqDataSetFromMatrix(countData = sorted.res,
                              colData = desine,
                              design= ~ subtype)

vsd <- vst(dds, blind=FALSE) #blind=FALSE means that the group means are taken into consideration when estimating the dispersion-mean relationship

t = list(size = 18)
m <- list(
  l = 50,
  r = 50,
  b = 50,
  t = 50,
  pad = 4)


PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
p
export(p, "PCA_stage.png")
#color by librery size
PCA = plotPCA(vsd, "totalCount")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
            colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
p
export(p, "PCA_total_count.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
p
export(p, "PCA_batch.png")

###############normalize and export data for analysis
#library size normalize the data
dds <- estimateSizeFactors(dds)
DESEq2_norm_counts = counts(dds, normalized = TRUE)

# remove potential batch effect
keep = rowSums(DESEq2_norm_counts > 1) >= 4
norm_counts = DESEq2_norm_counts[keep,]
dim(DESEq2_norm_counts)
dim(norm_counts)

library(sva)
batch = desine$batch
modcombat = model.matrix(~as.factor(subtype), data=desine) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(norm_counts+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

combat_edata = round(combat_edata)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData = desine,
                              design= ~ batch + subtype)
levels(dds$subtype)

vsd <- vst(dds, blind=FALSE)

library(vsn)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd), ranks = FALSE)$gg +ggtitle("batch corrected NormTransformed") + ylim(0,3)
meanSdPlot(assay(vsd), ranks = FALSE)$gg +ggtitle("batch corrected VST_Transformed") + ylim(0,3)


setwd("~/Desktop/s3s_Prostate/data/")
PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
p
export(p, "PCA_stage_batchcorr_onDESeq_norm_counts.png")
#color by librery size
PCA = plotPCA(vsd, "totalCount")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
            colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
p
export(p, "PCA_total_count_batchcarr_onDESeq_norm_counts.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
p
export(p, "PCA_batch_batchcorr_onDESeq_norm_counts.png")

#annotate gene names to HUGO symbols
mixture = data.frame(GeneNames = rownames(combat_edata), combat_edata)

library(AnnotationDbi)
library(org.Hs.eg.db)

genes = as.character(mixture$GeneNames)
genes.list <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
list1 = data.frame(lapply(genes.list, type.convert), stringsAsFactors=FALSE)
list2 = t(list1)
list2 = data.frame(list2)
list2$ens_id = rownames(list2)
ex1 = merge(list2, mixture, by.x="ens_id", by.y="GeneNames")
ex1$ens_id = NULL
ex2 = ex1[complete.cases(ex1),]
colnames(ex2)[1] = "genes"
table(duplicated(ex2[,1]))
ex3 = aggregate(. ~ genes, data = ex2, max)
rownames(ex3) = ex3$genes
ex3$genes = NULL

##################
all(colnames(ex3) == desine$name)
normal = ex3[,desine$subtype == "normal"]
BPH_LG_HGPIN = ex3[,desine$subtype %in% c("BPH","LGPIN","HGPIN")]
T3T34 = ex3[,desine$subtype %in% c("T3","T34")]
T4T45T5 = ex3[,desine$subtype %in% c("T4","T45","T5")]


setwd("~/Desktop/s3s_Prostate/data_for_deco/")
dataPr = list(normal, BPH_LG_HGPIN, T3T34, T4T45T5)
names(dataPr) = c("normal", "BPH_LG_HGPIN", "T3T34", "T4T45T5")
for(i in names(dataPr)){
  dataPr[[i]] = data.frame(genenames = rownames(ex3), dataPr[[i]])
  write.table(dataPr[[i]], paste0("pr_", i, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}


test =data.frame(genes = rownames(normal),normal)
colnames(normal)
