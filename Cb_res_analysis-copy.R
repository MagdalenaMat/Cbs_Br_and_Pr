mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


#############
path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/6_batches/Filtered/" #6b_
path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/batch_corr/Filtered/" #bulk_corr_n_
path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/not_batch_corr/Filtered/" #6b_with_Megan_not_corr

setwd(path)
files <- list.files(path=path)
files
#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_bulk_corr_n_",".txt_GEPs_Filtered.txt"),c("",""),i)
  namesGEP = c(namesGEP, nam)
}

namesGEP

#put all files in the list
filteredGEP = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  filteredGEP[[i]] = r
}
names(filteredGEP) = namesGEP

lapply(filteredGEP, dim) ##dim differ between diferent stages (runs of CibersortX)

gene_count = data.frame()
for(i in namesGEP){
  gene_count = rbind(gene_count, apply(filteredGEP[[i]], 2, function(x){sum(!is.na(x))}))
}

rownames(gene_count) = namesGEP
colnames(gene_count) =colnames(filteredGEP[[1]])

gene_count
#extract all gene names from the files and put them in the list to find intersect 
#of genes acros all the histological subtypes
results = list()
for (i in namesGEP){
  results[[i]] = as.character(filteredGEP[[i]][,1])
}

universe = Reduce(intersect, results)

#extract genes in monocyte signature common for all samples 
Mono_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Mono = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,2])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = c("DCIS","EN","IDC","normal")
#colnames(Mono_ex) = namesGEP

Mono_ex[Mono_ex<=1] <- NA
predicted_genes_Monocytes = apply(Mono_ex, 2, function(x){sum(!is.na(x))})

#write.table(predicted_genes_Monocytes, "predicted_genes_Monocytes_in_stages.txt", col.names = FALSE, quote = FALSE)

#Mono_ex_DCIS = Mono_ex[!is.na(Mono_ex$DCIS_Mono),] 
#cor(Mono_ex,use="pairwise.complete.obs",method="pearson") # can't use spearman here
al2_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 2,]
al3_mono = Mono_ex[rowSums(!is.na(Mono_ex)) >= 3,]
# library(heatmaply)
# heatmaply_na(al3_mono, showticklabels = c(T,F))


##############
#StErrors
#############

path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/6_batches/StEr/"
path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/batch_corr/StEr/"
path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/not_batch_corr/StEr/"
setwd(path)
files <- list.files(path=path)
files
#extract names
namesStEr = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_bulk_corr_n_",".txt_GEPs_StdErrs.txt"),c("",""),i)
  namesStEr = c(namesStEr, nam)
}

namesStEr

#put all files in the list
StEr = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  StEr[[i]] = r
}
names(StEr) = namesStEr


#extract genes in monocyte signature common for all samples 
Mono_ster = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesStEr){
  Error = StEr[[i]][StEr[[i]][,1] %in% universe,c(1,6)] #select gene_names and Monocytes
  Mono_ster = cbind(Mono_ster,Error)
}

all(as.character(Mono_ster[,2])==as.character(Mono_ster[,4]))
rownames(Mono_ster) = as.character(Mono_ster[,2])
Mono_ster = Mono_ster[, seq(3, length(Mono_ster), 2)]
colnames(Mono_ster) = c("DCIS","EN","IDC","normal")

all(rownames(Mono_ex) == rownames(Mono_ster))

# Mono_ster = data.frame(gene = rownames(Mono_ster), Mono_ster)
#Mono_ster = Mono_ster[,c(1,5,3,2,4)]
# Mono_ster$gene = sub("HLA.","HLA-", Mono_ster$gene) 

###########DE Chloe

Mono_ex = Mono_ex[,c(4,2,1,3)]
Mono_ster = Mono_ster[,c(4,2,1,3)]

Geps = Mono_ex
StdErr = Mono_ster

Zqvals = matrix(NA,nrow(Geps),1)
for(i in 1:(length(Geps)-1)){
  vBetaZ <- sapply(1:nrow(Geps), function(j) (Geps[j,i]-Geps[j,i+1])/sqrt(StdErr[j,i]^2+StdErr[j,i+1]^2))
  ZPs <- 2*pnorm(-abs(vBetaZ))
  Zqvals <- cbind(Zqvals,p.adjust(ZPs, method="BH"))
}

Zqvals = Zqvals[,-1]
colnames(Zqvals) = c("normal_EN", "EN_DCIS", "DCIS_IDC")
rownames(Zqvals) = rownames(Geps)

normal_EN = rownames(Zqvals)[which(Zqvals[,"normal_EN"]<= 0.1)]
EN_DCIS = rownames(Zqvals)[which(Zqvals[,"EN_DCIS"]<= 0.05)]
DCIS_IDC = rownames(Zqvals)[which(Zqvals[,"DCIS_IDC"]<= 0.05)]

ImGenes = list()
ImGenes[["IL"]] = c("IL10","IL21R","IL24","IL25")
ImGenes[["IFN"]] = c("IFNAR2","LYZ","NEK2")
genes = c("IL","IFN")

for(i in genes){
  Immuno = Mono_ex[rownames(Mono_ex) %in% ImGenes[[i]],]
  Immuno = data.frame(gene.names = rownames(Immuno), Immuno)
  #Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.expression")
  er = Mono_ster[rownames(Mono_ster) %in% Immuno$gene.names,]
  er = data.frame(gene = rownames(er), er)
  #er$gene = as.factor(as.character(er$gene))
  er_long = melt(er, id.vars = "gene")
  long_Imm = data.frame(long_Imm, sd = er_long$value)
  long_Imm$sd[is.na(long_Imm$gene.expression)] = NA
  print(ggplot(long_Imm,aes(gene.names,gene.expression, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+
          theme(text = element_text(size=20))+ 
          geom_errorbar(aes(ymin=gene.expression-sd, ymax=gene.expression+sd), width=.2, position=position_dodge(.9)) +
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red")) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1))+
          ggtitle("6_batches_with_Megan_NOT_batch_corrected"))
}

####Tcells
Tcell_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Tcell = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c(1,4)] #select gene_names and Tcell
  Tcell_ex = cbind(Tcell_ex,Tcell)
}

for(i in seq(4,length(Tcell_ex)-1,2)){
  print(all(as.character(Tcell_ex[,i]) == as.character(Tcell_ex[,c(i+2)])))
}

rownames(Tcell_ex) = as.character(Tcell_ex[,2])
Tcell_ex = Tcell_ex[, seq(3, length(Tcell_ex), 2)]
colnames(Tcell_ex) = c("AVL","DCIS","EN","IDC","LN","MetECE","normal")

Tcell_ex[Tcell_ex<=1] <- NA
al2_Tcell = Tcell_ex[rowSums(!is.na(Tcell_ex)) >= 2,]

to_reshape = data.frame(gene.names = rownames(al2_Tcell), al2_Tcell)
##########
#plot immune related gene values in different stages 
###########
library(plotly)
library(reshape2)
library(dplyr)
library(ggplot2)
to_reshape = data.frame(gene.names = rownames(al2_mono), al2_mono)
colnames(to_reshape)
to_reshape = to_reshape[,c(1,5,3,2,4)]

#compute log2 Fold Changes
change = mutate(to_reshape, EN_DCIS = log2(to_reshape$DCIS/to_reshape$EN), 
                DCIS_IDC = log2(to_reshape$IDC/to_reshape$DCIS),
                EN_IDC = log2(to_reshape$IDC/to_reshape$EN))
change$gene.names = sub("HLA.","HLA-", change$gene.names) 
change$gene.names = sub("NKX2.","NKX2-", change$gene.names) 
change$gene.names = sub("WT1.","WT1-", change$gene.names)

ImGenes = list()
ImGenes[["TLR_NFkB"]] = c("TLR3","TLR4","TLR1","TLR6","IL10RA","IL12A", "CYLD", "NFKBIA","NFKBIE", "MYD88")
ImGenes[["NLRs"]] = c("CASP1","PYCARD","GSDMD", "TLR3", "DHX9")
ImGenes[["HLA"]] = c("HLA-DMA","HLA-DPB1","HLA-DRB1","HLA-DQB1")
ImGenes[["TNF"]] = c("TNFSF12","TNFRSF1A","TRADD")
ImGenes[["NFkB"]] = c("NFKBIA","MYD88")
ImGenes[["B2M"]] = c("B2M","HLA-A","HLA-B","HLA-C")
ImGenes[["CD"]] = c("CD14","CD163")
ImGenes[["IFN"]] = c("IFNGR1","IFNGR2","NFKBIA")
ImGenes[["MRT"]] = c("UQCRB","UQCRC1","UQCRC2","UQCRQ")
ImGenes[["selected"]] = c("CYLD", "CASP1","TRADD")
ImGenes[["selected1"]] = c("HLA-B")
ImGenes[["selected2"]] = c("B2M")
ImGenes[["selected3"]] = c("NFKBIA", "PYCARD")
genes = c("TLR_NFkB","NLRs","CD","HLA","IFN","TNF","NFkB", "B2M","MRT","selected","selected1", "selected2", "selected3")  
#genes = c("CD", "HLA","IFN", "selected")

for(i in genes){
  Immuno = change[change$gene.names %in% ImGenes[[i]],]
  Immuno = Immuno[,c(1:5)]
  Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.expression")
  er = Mono_ster[Mono_ster$gene %in% Immuno$gene.names,]
  er$gene = as.factor(as.character(er$gene))
  er_long = melt(er, id.vars = "gene")
  long_Imm = data.frame(long_Imm, sd = er_long$value)
  long_Imm$sd[is.na(long_Imm$gene.expression)] = NA
  print(ggplot(long_Imm,aes(gene.names,gene.expression, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+
          theme(text = element_text(size=20))+ 
          geom_errorbar(aes(ymin=gene.expression-sd, ymax=gene.expression+sd), width=.2, position=position_dodge(.9)) +
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red")) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1))+
          ggtitle("6_batches_with_Megan_NOT_batch_corrected"))
}

  plot(rnorm(100))

png('selectedv2.png',width = 1200, height = 1200, res = 120)
print(ggplot(long_Imm,aes(gene.names,gene.expression, group = stage, fill =stage)) + 
        geom_col(position = "dodge")+ theme_classic()+
        theme(text = element_text(size=25))+
        scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red")))
dev.off()


########for bootstraped

gnames = unlist(ImGenes)
to_reshape = to_reshape[to_reshape[,1] %in% gnames,]
to_plot = list()
for(i in c("EN","DCIS","IDC")){
  to_plot[[i]] = as.matrix(to_reshape[,grepl(i,colnames(to_reshape))])
}

#stack dataframes together to have it in good format for plotting 
res.boot = data.frame()
for(i in 1:length(to_plot)){
  nstag = data.frame(t(to_plot[[i]]), stage = names(to_plot)[i])
  res.boot = rbind(res.boot,nstag)
}

for(gene in gnames){
  p = ggboxplot(res.boot, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

for(gene in gnames){
  p = ggboxplot(res.boot, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list(c("normal", "EN"), c("EN", "DCIS"), c("DCIS","IDC"),c("normal", "DCIS"),c("EN","IDC"),c("normal", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}


#####for data from 100x bootstrapp
library(reshape2)
library(dplyr)
library(ggplot2)

boot.immuno = read.table("/home/magda/Desktop/bootstrap100x_plots/long_res_boot_100_ImGenes.txt", sep = "\t", header = TRUE, row.names = 1)
boot.immuno = boot.immuno[boot.immuno$stage %in% c("EN","DCIS","IDC"),]
boot.immuno = boot.immuno[,colnames(boot.immuno) %in% c("CD14","CD163","NFKBIA","PYCARD","IL10RA","stage")]
long_Bimmuno = melt(boot.immuno, id.vars = "stage")
names(long_Bimmuno) = c("stage","gene.name","expression")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
df1 <- data_summary(long_Bimmuno, varname="expression", 
                    groupnames=c("stage", "gene.name"))
df1$gene.name=as.factor(df1$gene.name)
df1$stage = factor(df1$stage, levels = c("EN","DCIS","IDC"))
p <- ggplot(df1, aes(x=gene.name, y=expression, fill=stage)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=expression-sd, ymax=expression+sd), width=.2,
                position=position_dodge(.9))

p + scale_fill_brewer(palette="Paired") + 
  scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+
  theme_classic() +theme(axis.text=element_text(size=24, colour = "black"), axis.ticks = element_line(size = 2, colour = "black"),
                         axis.line = element_line(size = 2, colour = "black"), 
                         plot.title = element_text(hjust = 0.5), 
                         axis.title=element_text(size=22,face="bold"), 
                         legend.text = element_text(size=20), legend.title = element_text(size = 20)) 


#####################boxplots
library(ggpubr)
boot.immuno$stage = factor(boot.immuno$stage, levels = c("EN","DCIS","IDC"))
gnames = c("CD14","CD163","NFKBIA","PYCARD","IL10RA")
for(gene in gnames){
  p = ggboxplot(boot.immuno, x = "stage", y = gene,outlier.colour = NA,
                color = "stage",palette =c("dodgerblue3","orange","chartreuse4","red"))+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list( c("EN", "DCIS"), c("DCIS","IDC"),c("EN", "IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

