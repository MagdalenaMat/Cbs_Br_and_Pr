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

#delete scDCIS
dim(results[["DCIS"]])
results[["DCIS"]] = results[["DCIS"]][,!grepl("DCIS_single|DCIS_ablation", colnames(results[["DCIS"]]))]

desine1 = desine[!grepl("DCIS_single|DCIS_ablation",desine$name),]
dim(desine1)

results_df = data.frame(names = results$DCIS$GeneNames)
for(i in names(results)){
  results_df = cbind(results_df, results[[i]][,c(2:length(colnames(results[[i]])))])
}

#deduplicate gene names names, leave the row with the max value per gene  
results_df <- aggregate(. ~ names, data = results_df, max)
rownames(results_df) = results_df$names
results_df$names = NULL
all(desine1$name == colnames(results_df))
desine1$batch = sub("\\_.*","", desine1$batch)
desine1$batch = sub("\\..*","", desine1$batch)

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
desine1[desine1$size_f == "smaller_than_20K",c("name","batch")]

library(ggplot2)
png("total_reads.png", width = 1200, height = 600)
ggplot(desine1, aes(x = name, y = sizes, 
                    fill = size_f)) + 
  facet_grid(. ~ batch, labeller = label_value, scales = "free", space = "free") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank()) +
  scale_fill_brewer(palette="Dark2")
dev.off()

setwd("~/Desktop/s7s_Breast/")
write.table(desine1, "desine.BC.txt", sep = "\t")
write.table(results_df, "breast.data.combined.agregated.gene.names.txt", sep = "\t", row.names = TRUE)

desine1 = read.table("desine.BC.txt", sep = "\t")
results_df = read.table("breast.data.combined.agregated.gene.names.txt", sep = "\t", row.names = 1)

#filter out samples that have colSums < 10000
dim(results_df)
sorted.res = results_df[rowSums(results_df) > 0 , desine1$sizes > 10000]
dim(sorted.res)
desine2 = desine1[desine1$sizes > 10000,]
all(desine2$name == colnames(sorted.res))
desine2$batch = factor(desine2$batch)
desine2$batch = sub("MeganBatch","Megan",desine2$batch)
desine2$subtype = factor(desine2$subtype)
desine2$total_count = desine2$sizes > 300000

desine2$Megan = NA
desine2$Megan[grep("Megan", desine2$batch)] = "Megan"
desine2$Megan[is.na(desine2$Megan)] = "rest"

sorted.res = as.matrix(sorted.res)
rownames(desine2) = desine2$name
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

setwd("~/Desktop/s7s_Breast/data")
PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
export(p, "PCA_stage.png")
#color by librery size
PCA = plotPCA(vsd, "sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
        colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
export(p, "PCA_total_count.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
export(p, "PCA_batch.png")
#color by total_count
PCA = plotPCA(vsd, "total_count")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 10),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "total count bigget than 300K reads",legend = list(orientation = 'h'))
export(p, "PCA_count_300K.png")
#color by Megan

PCA = plotPCA(vsd, "Megan")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 10),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "batches from Megan",legend = list(orientation = 'h'))
export(p, "PCA_Megan.png")
#color by sizes
PCA = plotPCA(vsd, "size_f")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set2", 
        marker = list(size = 10),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "size factor colored",legend = list(orientation = 'h'))
export(p, "PCA_size_factor.png")


#library size normalize the data
dds <- estimateSizeFactors(dds)
DESEq2_norm_counts = counts(dds, normalized = TRUE)

# remove potential batch effect

TPM = 1E6 * sweep(sorted.res, 2, colSums(sorted.res), "/")
keep = rowSums(TPM > 1) >= 4
TPM1 = TPM[keep,]

library(sva)
batch = desine2$batch
modcombat = model.matrix(~as.factor(subtype), data=desine2) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(TPM1+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

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


setwd("~/Desktop/s7s_Breast/data")
PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by stage", 
         legend = list(orientation = 'h'), 
         margin = m)
export(p, "PCA_stage_batchcorr.png")
#color by librery size
PCA = plotPCA(vsd, "sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
            colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by total count in every sample",legend = list(orientation = 'h'))
export(p, "PCA_total_count_batchcarr.png")
#color by batch
PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, colors = "Set2", 
            marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>%
  layout(title = "colored by batch",legend = list(orientation = 'h'))
export(p, "PCA_batch_batchcorr.png")


PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20))%>%
  layout(title = "ComBat batch corrected")

PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20))%>%
  layout(title = "ComBat batch corrected")


#annotate gene names to HUGO symbols

