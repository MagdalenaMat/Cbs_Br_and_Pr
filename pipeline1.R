#select the stages from files

path = "~/Desktop/s7s_Breast/data/subset/"
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


gene.names = as.character(read.csv(files[1])[,1])

names <- data.frame(matrix(gene.names,58037,1))
dim(names)

results = list()
desine = data.frame(name = NA, batch = NA, subtype = NA)
for(j in c("DCIS|_CIS|LCIS", "EN|ALH|ADH|FEA|UDH|LobNeo|_75_L", "IDC|invasive|idc|ILC","normal|NL")){
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

desine$subtype = gsub("DCIS\\|_CIS\\|LCIS","DCIS",desine$subtype)
desine$subtype = gsub("EN\\|ALH\\|ADH\\|FEA\\|UDH\\|LobNeo\\|_75_L","EN",desine$subtype)
desine$subtype = gsub("IDC\\|invasive\\|idc\\|ILC","IDC",desine$subtype)
desine$subtype = gsub("normal\\|NL","normal",desine$subtype)

names(results) = c("DCIS","EN","IDC","normal")
lengths(results)
#delete scDCIS
dim(results[["DCIS"]])
results[["DCIS"]] = results[["DCIS"]][,!grepl("DCIS_single|DCIS_ablation", colnames(results[["DCIS"]]))]

desine1 = desine[!grepl("DCIS_single|DCIS_ablation",desine$name),]
dim(desine1)

results_df = data.frame(names = results$DCIS$GeneNames)
for(i in names(results)){
  results_df = cbind(results_df, results[[i]][,c(2:length(colnames(results[[i]])))])
}
dim(results_df)

#deduplicate gene names names, leave the row with the max value per gene  
results_df <- aggregate(. ~ names, data = results_df, max)
rownames(results_df) = results_df$names
results_df$names = NULL
all(desine1$name == colnames(results_df))
desine1$batch = sub("\\_.*","", desine1$batch)
#desine1$batch = sub("\\..*","", desine1$batch)

desine1$batch[desine1$batch == "August"] = "160801" 

###remove normal prostate
results_df = results_df[,!(desine1$subtype == "normal" & desine1$batch == "170713")]
desine1 = desine1[!(desine1$subtype == "normal" & desine1$batch == "170713"),]
all(desine1$name == colnames(results_df))
###

#calculate total read counts after duplicate gene names were removed
desine1 = cbind(desine1, 
                sizes = colSums(results_df))

desine1$size_f = NA

for(i in c(300000,100000,20000)){
  desine1[desine1$sizes < i,"size_f"] = paste0("smaller_than_",i)
}  

desine1[is.na(desine1$size_f), "size_f"] = "bigger_than_300K"

desine1$size_f = factor(desine1$size_f)
levels(desine1$size_f) = c("bigger_than_300K",
                           "between_20K_and_100K",
                           "smaller_than_20K",
                           "between_100K_and_300K")
table(desine1$size_f)
desine1[desine1$size_f == "smaller_than_20K",c("name","batch", "sizes")]

library(ggplot2)
setwd("~/Desktop/s7s_Breast/data/pictures/20180314_9batches")
png("total_reads_all_batches.png", width = 1200, height = 600)
ggplot(desine1, aes(x = name, y = sizes, 
                    fill = size_f)) + 
  facet_grid(. ~ batch, labeller = label_value, scales = "free", space = "free") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank()) +
  scale_fill_brewer(palette="Dark2")
dev.off()

setwd("~/Desktop/s7s_Breast/data/")
write.table(desine1, "desine_20180314_9batches.txt", sep = "\t")
write.table(results_df, "breast.data.combined.agregated.gene.names_20180314_9_batches.txt", sep = "\t", row.names = TRUE)

setwd("~/Desktop/s7s_Breast/data/")
desine1 = read.table("desine_20180314_9batches.txt", sep = "\t")
results_df = read.table( "breast.data.combined.agregated.gene.names_20180314_9_batches.txt", sep = "\t", row.names = 1)

#filter out samples that have colSums < 10000
dim(results_df)
sorted.res = results_df[rowSums(results_df) > 0 , desine1$sizes > 20000]
dim(sorted.res)
desine2 = desine1[desine1$sizes > 20000,]
all(desine2$name == colnames(sorted.res))
desine2$batch = factor(desine2$batch)
desine2$batch = sub("MeganBatch","Megan",desine2$batch)
desine2$subtype = factor(desine2$subtype)
desine2$total_count = desine2$sizes > 300000

desine2$Megan = NA
desine2$Megan[grep("Megan", desine2$batch)] = "Megan"
desine2$Megan[is.na(desine2$Megan)] = "rest"

desine2$batch18 = NA
desine2$batch18[grep("^18", desine2$batch)] = "batches_2018"
desine2$batch18[is.na(desine2$batch18)] = "rest"

sorted.res = as.matrix(sorted.res)
rownames(desine2) = desine2$name

norm_sorted.res = sorted.res[,desine2$subtype == "normal"] 
norm_desine2 = desine2[desine2$subtype == "normal",]

DCIS_sorted.res = sorted.res[,desine2$subtype == "DCIS"] 
DCIS_desine2 = desine2[desine2$subtype == "DCIS",]

#make dds object
library(DESeq2)
library(plotly)
dds <- DESeqDataSetFromMatrix(countData = sorted.res,
                              colData = desine2,
                              design= ~ subtype)

vsd <- vst(dds, blind=FALSE) #blind=FALSE means that the group means are taken into consideration when estimating the dispersion-mean relationship

t = list(size = 18)
m <- list(
  l = 50,
  r = 50,
  b = 50,
  t = 50,
  pad = 4)

setwd("~/Desktop/s7s_Breast/data/pictures/20180314_9batches/DCIS_only/")
PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
p
export(p, "PCA_9_batches.png")
#color by librery size
PCA = plotPCA(vsd, "sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
        colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
p
export(p, "PCA_total_9_batches.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
p
export(p, "PCA_batch_9_batches.png")
#color by total_count
PCA = plotPCA(vsd, "total_count")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "total count bigget than 300K reads",legend = list(orientation = 'h'))
p
export(p, "PCA_count_300K_9_batches.png")
#color by Megan

PCA = plotPCA(vsd, "Megan")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "batches from Megan",legend = list(orientation = 'h'))
p
export(p, "PCA_Megan_all_batches.png")
#color by sizes
PCA = plotPCA(vsd, "size_f")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "size factor colored",legend = list(orientation = 'h'))
p
export(p, "PCA_size_factor_normal_9_batches.png")

#color by betches from 2018
PCA = plotPCA(vsd, "batch18")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "size factor colored",legend = list(orientation = 'h'))
p
export(p, "PCA_batches_2018_9_batches.png")

table(desine2$subtype, desine2$batch18)
table(desine2$batch, desine2$batch18)
table(desine2$batch, desine2$subtype)
table(desine2$subtype)
#library size normalize the data
dds <- estimateSizeFactors(dds)
DESEq2_norm_counts = counts(dds, normalized = TRUE)


# remove potential batch effect

#TPM = 1E6 * sweep(sorted.res, 2, colSums(sorted.res), "/")
keep = rowSums(DESEq2_norm_counts > 1) >= 4
norm_counts = DESEq2_norm_counts[keep,]
dim(DESEq2_norm_counts)
dim(norm_counts)

library(sva)
batch = desine2$batch
modcombat = model.matrix(~as.factor(subtype), data=desine2) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(norm_counts+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

combat_edata = round(combat_edata)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData = desine2,
                              design= ~ batch + subtype)
levels(dds$subtype)

vsd <- vst(dds, blind=FALSE)

library(vsn)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd), ranks = FALSE)$gg +ggtitle("batch corrected NormTransformed") + ylim(0,3)
meanSdPlot(assay(vsd), ranks = FALSE)$gg +ggtitle("batch corrected VST_Transformed") + ylim(0,3)

setwd("~/Desktop/s7s_Breast/data/pictures/8batches/batch_corrected/")
PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
p
export(p, "PCA_stage_batch_corr_all_batches.png")
#color by librery size
PCA = plotPCA(vsd, "sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
            colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
p
export(p, "PCA_total_count_batch_corr_all_batches.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
p
export(p, "PCA_batch_batch_corr_all_batches.png")

#color by betches from 2018
PCA = plotPCA(vsd, "batch18")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "size factor colored",legend = list(orientation = 'h'))
p
export(p, "PCA_batches_2018_batch_corr_all_batches.png")

#annotate gene names to HUGO symbols
#mixture = data.frame(GeneNames = rownames(combat_edata), combat_edata)
mixture = data.frame(GeneNames = rownames(DESEq2_norm_counts), DESEq2_norm_counts)

library(AnnotationDbi)
library(org.Hs.eg.db)

mixture$genes = mapIds(org.Hs.eg.db,
                       keys = as.character(mixture$GeneNames),
                       column = c("SYMBOL"),
                       keytype = "ENSEMBL",
                       multiVals = "first")

mixture$GeneNames = NULL
mixture = mixture[complete.cases(mixture),]


ex3 = aggregate(. ~ genes, data = mixture, max)
rownames(ex3) = ex3$genes
ex3$genes = NULL

##########prepare data for fraction analysis
pdata = c()
to_fractions = data.frame(na = matrix(NA, dim(ex3)[1],1))
for(i in c("normal","EN","DCIS","IDC","LN")){
  to_fractions = cbind(to_fractions, ex3[,grepl(i,desine2$subtype)])
  pdata = c(pdata, rep(i, sum(grepl(i,desine2$subtype))))
}
to_fractions$na = NULL

pheno_data = data.frame(names = colnames(to_fractions), subtypes = pdata)
to_fractions = data.frame(names = rownames(to_fractions), to_fractions)
write.table(to_fractions, "~/Desktop/s7s_Breast/to_fractions/all_data_to_fractions.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(pheno_data, "~/Desktop/s7s_Breast/to_fractions/pdata_all_samples.txt", sep = "\t", row.names = F, col.names = T, quote = F)

fractions = read.csv("~/Desktop/s7s_Breast/to_fractions/all_data_fractions.csv")
rownames(fractions) = fractions$Mixture
fractions$Mixture = NULL
type = as.factor(pheno_data$subtypes)

to_plot = as.matrix(fractions[1:22])
library(pheatmap)
library(RColorBrewer)

col1 <- brewer.pal(8, "Set2")

mydf <- data.frame(row.names = rownames(to_plot), category = type)

table(type)
# add row annotations
myBreaks = unique(c(seq(0, 0.29, length=49), 0.3, seq(0.31,max(to_plot), length=50)))
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
         annotation_row = mydf,gaps_row = c(44,71,242),
         breaks = myBreaks, show_rownames = F)


############

all(colnames(ex3) == desine2$name)
normal = ex3[,desine2$subtype == "normal"]
EN = ex3[,desine2$subtype == "EN"]
DCIS = ex3[,desine2$subtype == "DCIS"]
IDC = ex3[,desine2$subtype == "IDC"]

table(desine2$subtype)

dataPr = list(normal, EN, DCIS, IDC, ex3)
names(dataPr) = c("normal", "EN", "DCIS", "IDC", "full_data")

setwd("~/Desktop/s7s_Breast/data_for_decon/20180314_8batches_not_corrected_for_batch/")
for(i in names(dataPr)){
  dataPr[[i]] = data.frame(genenames = rownames(ex3), dataPr[[i]])
  write.table(dataPr[[i]], paste0("20180314_not_batch_corrected_", i, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}


