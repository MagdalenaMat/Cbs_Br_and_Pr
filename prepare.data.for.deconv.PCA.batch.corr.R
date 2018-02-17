#prepare samples accross diferent experiments for PCA


path = "~/Desktop/s7s_Breast/data/subset/"
#path = "~/Desktop/breast_cancer_progression/data" #old data plus megan data where megan data are very different in comparison to other batches on PCA

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

results_df <- aggregate(. ~ names, data = results_df, max)
rownames(results_df) = results_df$names
results_df$names = NULL
all(desine1$name == colnames(results_df))
desine1$batch = sub("\\_.*","", desine1$batch)
desine1$batch = sub("\\..*","", desine1$batch)
desine1 = cbind(desine1, 
                sizes = colSums(results_df))

desine1$size_f = NA

for(i in c(300000,100000,10000)){
    desine1[desine1$sizes < i,"size_f"] = paste0("smaller_than_",i)
}  

desine1[is.na(desine1$size_f), "size_f"] = "bigger_than_0.3M"
desine1$size_f =factor(desine1$size_f)

levels(desine1$size_f) = c("bigger_than_0.3M","smaller_than_0.3M")
table(desine1$size_f)

library(ggplot2)
setwd("~/Desktop/s7s_Breast/data")
png("total_reads_5batches.png", width = 1200, height = 800)
ggplot(desine1, aes(x = name, y = sizes, 
                    fill = size_f)) + 
  facet_grid(. ~ batch, labeller = label_value, scales = "free", space = "free") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Dark2")
dev.off()

################################################################################################################
setwd("~/Desktop/s7s_Breast/")
write.table(desine1, "desine.BC_5batches.txt", sep = "\t")
write.table(results_df, "breast.data.combined.agregated.gene.names_5batches.txt", sep = "\t", row.names = TRUE)
################################################################################################################

dim(results_df)
sorted.res = results_df[rowSums(results_df) > 0 , desine1$sizes > 10000]
dim(sorted.res)
#boxplot(sorted.res)

#annotate
library(AnnotationDbi)
library(org.Hs.eg.db)
anot.counts = list()
for(i in names(results)){
  mixture = results[[i]]
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
  ex3 = ex2[!duplicated(ex2[,1]),]
  anot.counts[[i]] = ex3
}

#genes.list <- mapIds(org.Hs.eg.db, keys = "ENSG00000232810", column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

##############do PCA on results

for(i in 1:(length(anot.counts)-1) ){
  print(all(anot.counts[[i]][1] == anot.counts[[c(i+1)]][1]))
}

#nam.fac = as.character()
breast.data = data.frame(gene.names = anot.counts[[1]][,1])
for(name in names(anot.counts)){
  bdat = anot.counts[[name]]
  bdat = bdat[,2:length(colnames(bdat))]
  breast.data = cbind(breast.data, bdat)
  #nam.fac = c(nam.fac, rep(name, length(colnames(bdat))))
}


rownames(breast.data) = breast.data[,1]
breast.data = breast.data[,-1]

TPM = 1E6 * sweep(breast.data, 2, colSums(breast.data), "/")
keep = rowSums(TPM > 1) >= 4
TPM1 = TPM[keep,]

#breast.data = breast.data[rowSums(breast.data) > 300,]
breast.data = as.matrix(breast.data)

colnames(breast.data)[is.na(pmatch(colnames(breast.data), desine1$name))]

desine2 = desine1[match(colnames(breast.data), desine1$name),]
desine2$name == colnames(breast.data)
desine2$batch = factor(desine2$batch)
desine2$batch = sub("\\_.*","", desine2$batch) 
desine2$subtype = factor(desine2$subtype, levels = c("normal","EN","DCIS","IDC"))

setwd("~/Desktop/s7s_Breast/")
write.table(desine2, "desine.BC.txt", sep = "\t")
write.table(breast.data, "breast.data.combined.txt", sep = "\t", row.names = TRUE)
# nam.fac = factor(nam.fac, levels = c("normal","EN","DCIS","IDC","ECE","met_ECE", "LN", "AVL"))
# design1 = data.frame(samples = colnames(breast.data),subtype = nam.fac)
# design2 = merge(design1, design.df, by.x = "samples", by.y = "name" )
# design1$samples[is.na(pmatch(design1$samples,design.df$name))]
# design$Id[is.na(pmatch(design$Id,colnames(data.ordered)))]

library(DESeq2)
library(plotly)

total.read = data.frame(size = colSums(breast.data))
total.read$names = rownames(total.read)
total.read$OK = total.read$size > 20000
library(ggplot2)
ggplot(total.read, aes(x = names, y = size)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Dark2")

#check librery size
# scale to 1 to 100 range
# libSize <- colSums(breast.data > 0) # number of genes detected per cell as proxy
# libSize <- libSize - min(libSize) + 1
# libSize <- libSize / max(libSize)
# libSize <- round(libSize*100)

all(desine2$name == total.read$names)
desine2 = cbind(desine2, Lib_size = total.read$size, OK = total.read$OK)
desine2$megan = grepl("Megan",desine2$batch)

dds <- DESeqDataSetFromMatrix(countData = breast.data,
                              colData = desine2,
                              design= ~ subtype)


vsd <- vst(dds, blind=FALSE) #blind=FALSE means that the group means are taken into consideration when estimating the dispersion-mean relationship

library(vsn)
library(ggplot2)
meanSdPlot(assay(vsd), ranks=FALSE)$gg + ggtitle("VST normalized") + ylim(0,3)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd), ranks = FALSE)$gg + ggtitle("normTransformed") + ylim(0,3)
PCA = plotPCA(vsd, "subtype")

dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20)) %>%
  layout(title = "DESeq2::vst transformed original data")

#color by librery size
PCA = plotPCA(vsd, "OK")

dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = colorRampPalette(c('blue', 'yellow', 'red'))(100), marker = list(size = 20)) %>%
  layout(title = "DESeq2::vst transformed original data, samples > 20000 reads")


table(desine2$batch, desine2$subtype)

PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
#dataPCA = cbind(dataPCA, subtype = dds$subtype)
#symbol = ~subtype, symbols = list("15","16","17","1")
p <- plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, 
             colors = "Set2", marker = list(size = 20)) %>%
  layout(title = "DESeq2::vst transformed original data")
p


normal = TPM[,desine2$subtype == "normal"]
EN = TPM[,desine2$subtype == "EN"]
DCIS = TPM[,desine2$subtype == "DCIS"]
IDC = TPM[,desine2$subtype == "IDC"]

setwd("~/Desktop/s7s_Breast/data_for_decon/6_batches_and_Megan/not_batch_corr/")

data = list(normal,EN,DCIS,IDC)
names(data) = c("normal","EN","DCIS","IDC")
for( i in names(data)){
  new = data.frame(genes = rownames(data[[i]]), data[[i]])
  write.table(new, paste0("6b_with_Megan_not_corr",i,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}


# remove potential batch effect
library(sva)
batch = desine2$batch
modcombat = model.matrix(~as.factor(subtype), data=desine2) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(TPM1+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

combat_edata = round(combat_edata)
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData = desine2,
                              design= ~ batch + subtype)
levels(dds$subtype)

vsd <- vst(dds, blind=FALSE)

library(vsn)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd), ranks = FALSE)$gg +ggtitle("batch corrected NormTransformed") + ylim(0,3)
meanSdPlot(assay(vsd), ranks = FALSE)$gg +ggtitle("batch corrected VST_Transformed") + ylim(0,3)


PCA = plotPCA(vsd, "subtype")
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20))%>%
  layout(title = "ComBat batch corrected")

PCA = plotPCA(vsd, "batch")
dataPCA = PCA$data
plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, text = ~ name, colors = "Set2", marker = list(size = 20))%>%
  layout(title = "ComBat batch corrected")

#####
normal = combat_edata[,desine2$subtype == "normal"]
EN = combat_edata[,desine2$subtype == "EN"]
DCIS = combat_edata[,desine2$subtype == "DCIS"]
IDC = combat_edata[,desine2$subtype == "IDC"]

setwd("~/Desktop/s7s_Breast/data_for_decon/6_batches_and_Megan/batch_corr/")

data = list(normal,EN,DCIS,IDC)
names(data) = c("normal","EN","DCIS","IDC")
for( i in names(data)){
  new = data.frame(genes = rownames(data[[i]]), data[[i]])
  write.table(new, paste0("bulk_corr_n_",i,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

