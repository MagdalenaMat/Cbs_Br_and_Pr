#analyse only experiment form 2016 and 2017

path = "~/Desktop/s7s_Breast/data/subset/"
setwd(path)
files <- list.files(path=path)
files = files[c(1:3,9)]

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

names <- data.frame(matrix(gene.names,length(gene.names),1))
dim(names)

results = list()
desine = data.frame(name = NA, batch = NA, subtype = NA)
for(j in c("normal|NL","EN|ALH|ADH|FEA|UDH|LobNeo|_75_L","DCIS|_CIS|LCIS", "IDC|invasive|idc|ILC")){
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

names(results) = c("normal","EN","DCIS","IDC")
lengths(results)

results_df = data.frame(names = results$DCIS$GeneNames)
for(i in names(results)){
  results_df = cbind(results_df, results[[i]][,c(2:length(colnames(results[[i]])))])
}
dim(results_df)
colnames(results_df)

all(colnames(results_df)[-1] == desine$name)
sizes = colSums(results_df[,-1])

library(AnnotationDbi)
library(org.Hs.eg.db)

results_df$symbol = mapIds(org.Hs.eg.db,
                             keys = as.character(results_df$names),
                             column = c("SYMBOL"),
                             keytype = "ENSEMBL",
                             multiVals = "first")
results_df$names = NULL
results_df = results_df[complete.cases(results_df),]
results_df = data.frame(gene.names = results_df$symbol, results_df)
results_df$symbol = NULL

#deduplicate gene names names, leave the row with the max value per gene  
results_df <- aggregate(. ~ gene.names, data = results_df, max)
rownames(results_df) = results_df$gene.names
results_df$gene.names = NULL
all(desine$name == colnames(results_df))
desine$batch = sub("\\_.*","", desine$batch)


desine$batch[desine$batch == "August"] = "160801" 

####august_pdata = read.table("~/Desktop/s7s_Breast/data/August_batch_corrected_data/pdata_August.txt", sep = "\t", header = T)
august_pdata = read.csv("~/Desktop/s7s_Breast/data/pdata/design_August_2016_FFPE_LCM.csv")
august_pdata = data.frame(august_pdata, 
                          names = paste(august_pdata$patient,august_pdata$stage,august_pdata$biol.rep, sep = "_"))


desine[desine$name %in% august_pdata$names,]$batch = paste("160801",
                                                           august_pdata[match(desine[desine$name %in% august_pdata$names,]$name, 
                                                                              august_pdata$names),]$Batch,sep = "_")

desine1 = desine
desine1$batch_ver = desine$batch
desine1[grep("160801",desine1$batch_ver),]$batch_ver = paste("160801",
                                                             august_pdata[match(desine[grep("160801",desine1$batch),]$name, 
                                                                                august_pdata$names),]$Batch,sep = "_")

desine$batch == desine1$batch_ver

desine$subtype = factor(desine$subtype, levels = c("normal","EN","DCIS","IDC"))
desine$batch = factor(desine$batch)

levels(desine$batch) = c(paste0(rep("batch",7),1:7))

names(sizes) == desine$name
desine$sizes = sizes

desine$year = grepl("16",desine$batch)
desine$year = sub(TRUE,"2016",desine$year)
desine$year = sub(FALSE,"2017",desine$year)

table(desine$batch,desine$subtype)
library(DESeq2)
library(plotly)
dds <- DESeqDataSetFromMatrix(countData = results_df,
                              colData = desine,
                              design= ~ subtype)
vsd = vst(dds,blind=FALSE)
PCA = plotPCA(vsd,"batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>% 
  layout(yaxis = list(tickfont = list(size = 20), 
                      title = "PC2", 
                      titlefont = list(size = 35)),
         xaxis = list(tickfont = list(size = 20), 
                      title = "<br><br><br><br><br><br> PC1", 
                      titlefont = list(size = 35)),
         legend = list(x = 0.1 , y = 100, orientation = 'h',font = list(size = 20)))
p
setwd("~/Desktop/toll_plots/PCA/")
export(p, "not_corrected_batch.pdf")


PCA = plotPCA(vsd,"subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600) %>% 
  layout(yaxis = list(tickfont = list(size = 20), 
                      title = "PC2", 
                      titlefont = list(size = 30)),
         xaxis = list(tickfont = list(size = 20), 
                      title = "PC1", 
                      titlefont = list(size = 30)),
         legend = list(x = 0.1 , y = 100 ,orientation = 'h', font = list(size = 35)))
p
setwd("~/Desktop/toll_plots/PCA/")
export(p, "not_corrected_subtype.pdf")

PCA = plotPCA(vsd,"year")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "not_corrected_year.png")

PCA = plotPCA(vsd,"sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = colorRampPalette(c('blue', 'yellow', 'red'))(100),
            marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "not_corrected_size.png")

table(desine$batch, desine$subtype)

dds <- estimateSizeFactors(dds)
norm_data = counts(dds, normalized = TRUE)

keep = rowSums(norm_data > 1) >= 4
norm_counts = norm_data[keep,]
dim(norm_data)
dim(norm_counts)

to_fractions = data.frame(names = rownames(norm_counts),norm_counts)
#########write norm counts without batch correction
write.table(to_fractions, "~/Desktop/s7s_Breast/data/4stages/normalized_counts_2016_2017.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)


norm_counts = read.table("~/Desktop/s7s_Breast/data/4stages/normalized_counts_2016_2017.txt", sep = "\t", quote = "")

all(colnames(norm_counts) == desine$names)

to_write = list()
for(i in levels(desine$subtype)){
  to_write[[i]] = norm_counts[,desine$subtype == i]
}

setwd("~/Desktop/s7s_Breast/data/4stages/")
for(i in names(to_write)){
  to_write[[i]] = data.frame(genenames = rownames(norm_counts), to_write[[i]])
  write.table(to_write[[i]], paste0("4_stages_2016_2017_NOT_BC_", i, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}
#########



library(sva)
batch = desine$batch
modcombat = model.matrix(~as.factor(subtype), data=desine) #variable of interest shuld be included
combat_edata = 2^ComBat(dat=log2(norm_counts+1), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE) 

combat_edata = round(combat_edata)
dds <- DESeqDataSetFromMatrix(countData = combat_edata,
                              colData = desine,
                              design= ~ subtype)
vsd = vst(dds,blind=FALSE)
PCA = plotPCA(vsd,"batch")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "batch_corrected_batch.png")

PCA = plotPCA(vsd,"subtype")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "batch_corrected_subtype.png")


PCA = plotPCA(vsd,"year")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = "Set2", marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "batch_corrected_year.png")

PCA = plotPCA(vsd,"sizes")
dataPCA = PCA$data
p = plot_ly(data = dataPCA, x = ~PC1, y = ~PC2, color = ~group, 
            text = ~ name, font = t, colors = colorRampPalette(c('blue', 'yellow', 'red'))(100),
            marker = list(size = 20),
            autosize = F, width = 900, height = 600)
p
setwd("~/Desktop/s7s_Breast/data/4stages/plots/PCA/")
export(p, "batch_corrected_size.png")



combat_edata = data.frame(gene.names = rownames(combat_edata), combat_edata)


write.table(combat_edata, "~/Desktop/s7s_Breast/data/4stages/combat_batch_corrected_2016_2017.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)

all(colnames(combat_edata) == desine$names)

to_write = list()
for(i in levels(desine$subtype)){
  to_write[[i]] = combat_edata[,desine$subtype == i]
}

setwd("~/Desktop/s7s_Breast/data/4stages/")
for(i in names(to_write)){
  to_write[[i]] = data.frame(genenames = rownames(combat_edata), to_write[[i]])
  write.table(to_write[[i]], paste0("2016_2017_BC_", i, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}

write.table(desine, "~/Desktop/s7s_Breast/data/4stages/2016_2017_pdata.txt", sep = "\t", row.names = F, col.names = T, quote = F)


