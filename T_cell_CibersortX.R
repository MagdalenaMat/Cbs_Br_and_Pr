#### Magdalena Matusiak
#### 20 JAN 2018
#### analysis of immune landscape changes in breast cancer progression
#### LCM bulk RNASeq deconvoluted with CibersortX
#### 6 batches and Megan's data batch corrected
#Tcells


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


path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/batch_corr/Filtered/"
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

gene_count = data.frame()
for(i in namesGEP){
  gene_count = rbind(gene_count, apply(filteredGEP[[i]], 2, function(x){sum(!is.na(x))}))
}

rownames(gene_count) = namesGEP
colnames(gene_count) =colnames(filteredGEP[[1]])

gene_count

lapply(filteredGEP, dim) ##dim differ between diferent stages (runs of CibersortX)

#extract all gene names from the files and put them in the list to find intersect 
#of genes acros all the histological subtypes
results = list()
for (i in namesGEP){
  results[[i]] = as.character(filteredGEP[[i]][,1])
}

universe = Reduce(intersect, results)

#extract genes in T_cell signature common for all samples 
Tcell_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Tcell = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c("GeneSymbol","T.cells")] #select gene_names and T.cells
  Tcell_ex = cbind(Tcell_ex,Tcell)
}

all(as.character(Tcell_ex[,2])==as.character(Tcell_ex[,4]))
rownames(Tcell_ex) = as.character(Tcell_ex[,2])
Tcell_ex = Tcell_ex[, seq(3, length(Tcell_ex), 2)]
colnames(Tcell_ex) = c("DCIS","EN","IDC","normal")


Tcell_ex[Tcell_ex<=1] <- NA
predicted_genes_Tcells = apply(Tcell_ex, 2, function(x){sum(!is.na(x))})

##############
#StErrors
#############

path = "~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/batch_corr/StEr/"

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


#extract genes in Tcell signature common for all samples 
Tcell_ster = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesStEr){
  Error = StEr[[i]][StEr[[i]][,1] %in% universe,c("GeneSymbol","T.cells")] #select gene_names and Tcells
  Tcell_ster = cbind(Tcell_ster,Error)
}

all(as.character(Tcell_ster[,2])==as.character(Tcell_ster[,4]))
rownames(Tcell_ster) = as.character(Tcell_ster[,2])
Tcell_ster = Tcell_ster[, seq(3, length(Tcell_ster), 2)]
colnames(Tcell_ster) = c("DCIS","EN","IDC","normal")

all(rownames(Tcell_ex) == rownames(Tcell_ster))

# Mono_ster$gene = sub("HLA.","HLA-", Mono_ster$gene) 

###########DE Chloe

Tcell_ex = Tcell_ex[,c(4,2,1,3)]
Tcell_ster = Tcell_ster[,c(4,2,1,3)]

Geps = Tcell_ex
StdErr = Tcell_ster

Zqvals = matrix(NA,nrow(Geps),1)
cnames = c()
for(i in names(Geps)){
  for(j in names(Geps)){
    if(i <= j){
      next
    }
    cnames= c(cnames, paste0(i,"_",j))
    vBetaZ <- sapply(1:nrow(Geps), function(x) (Geps[x,i]-Geps[x,j])/sqrt(StdErr[x,i]^2+StdErr[x,j]^2))
    ZPs <- 2*pnorm(-abs(vBetaZ))
    Zqvals <- cbind(Zqvals,p.adjust(ZPs, method="BH"))
    
  }
} 

Zqvals = Zqvals[,-1]
colnames(Zqvals) = cnames
rownames(Zqvals) = rownames(Geps)
dim(Zqvals)
head(Zqvals)
Zqvals = data.frame(Zqvals)

names(Zqvals)[5:6] = c("EN_IDC", "DCIS_IDC")

sPv = list()
for(i in colnames(Zqvals)){
  sPv[[i]] = data.frame(genes = rownames(Zqvals)[which(Zqvals[,i]<= 0.1)], Zqv = Zqvals[which(Zqvals[,i]<= 0.1),i, drop = FALSE])
}

LFC = list()
for(i in names(Zqvals)){
  fr = strsplit(i, "[_]")[[1]][1]
  sd = strsplit(i, "[_]")[[1]][2]
  LFC[[i]] = data.frame(names = rownames(Geps), FC = log2(Geps[,sd]/Geps[,fr]))
}

#extract only sig FC
LFC_s = list()
for(i in names(LFC)){
  LFC_s[[i]] = LFC[[i]][LFC[[i]][,"names"] %in% sPv[[i]][,"genes"],]
}

for(i in names(LFC_s)){
  print(all(as.character(LFC_s[[i]][,"names"]) == as.character(sPv[[i]][,"genes"])))
}

#merge LFC and p.adj for downstream GO/GSEA analysis
to_enr = list()
for(i in names(LFC_s)){
  to_enr[[i]] = cbind(LFC_s[[i]], Zqv = sPv[[i]][,i])
}

##############
#plot interesting genes
##############

ImGenes = list()
ImGenes[["IL"]] = c("FOXP3","CXCL2","IL18R1","IL25","IL24","CD207","CD70","CXCL6","CD25")
ImGenes[["IFN"]] = c("IFNGR2","LYZ","NFKBIA")
genes = c("IL","IFN")

library(reshape2)
for(i in genes){
  Immuno = Tcell_ex[rownames(Tcell_ex) %in% ImGenes[[i]],]
  Immuno = data.frame(gene.names = rownames(Immuno), Immuno)
  #Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.expression")
  er = Tcell_ster[rownames(Tcell_ster) %in% Immuno$gene.names,]
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
          ggtitle("batch_corrected"))
}


####BC GO enrichment

to_GO = list()
for(nm in names(to_enr)){
  to_GO[[paste0("DE_",nm,"_up")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] > 0), "names"]
  to_GO[[paste0("DE_",nm,"_down")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] < 0), "names"]
}

lengths(to_GO)

library("Rgraphviz")
library("DOSE")
library("clusterProfiler")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

#translate gene symbols to entrez ID

for(element in names(to_GO)){
  to_GO[[element]] = sub("[.]","-", to_GO[[element]]) 
}

to_GO_EI = list()
for(element in names(to_GO)){
  to_GO_EI[[element]] = bitr(to_GO[[element]], 
                             fromType = "SYMBOL",
                             toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)
}

###run GO analysis, I am running ont= BP (biological processes) /MF/CC
GO_res = list()
for(element in names(to_GO_EI)){
  GO_res[[element]] = enrichGO(gene         = to_GO_EI[[element]]$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable = TRUE)
}

#path = setwd("~/Desktop/breast_cancer_progression/plots/BC_30_10_17/enrichment_analysis")
for(element in names(GO_res)){
  #png(paste0(path,"/",element,".png"), width = 900, height = 1100)
  plot(dotplot(GO_res[[element]],showCategory=50)+labs(title= paste0("Tcell_",element)))
  #dev.off()
}

res = list()
for(element in names(GO_res)){
  if(nrow(GO_res[[element]]@result) == 0){
    next
  }
  res[[element]] = GO_res[[element]]@result[,c("ID","p.adjust", "Description","geneID","GeneRatio")]
}

#subset ID and p.adj for revigo analysis
to_REV1 = GO_res$DE_DCIS_IDC_up@result[,c("ID","p.adjust")]
write.table(to_REV1, "~/Desktop/breast_cancer_progression/plots/BC_30_10_17/GO_res_DCIS_IDC_up", sep = "\t", row.names = F, quote = FALSE)

######disect genes for ploting

to_plot = list()
for(i in c("T |T-","MHC|antigen|present","neutroph","lip|LDL","horm","vir")){
  print(i)
  for(j in names(res)){
    print(j)
    print(res[[j]][grepl(i,res[[j]][,"Description"]),"Description"])
    gr = res[[j]][grepl(i,res[[j]][,"Description"]),"geneID"]
    gr = unique(unlist(strsplit(gr,"/")))
    to_plot[[i]] = unique(c(to_plot[[i]],gr))
  }
  print(to_plot[[i]])
}

to_plot

library(heatmaply)
library(gplots)
library(RColorBrewer)
library(genefilter)



rownames(Geps) = sub("[.]","-", rownames(Geps)) 


genes = to_plot
names(genes) = c("T_cell_responce","Antigen_presentation","neutrophil_responce","Lipid_metabolism","Horm_responce","Viral_responce")
for(element in names(genes)){
  to_heatmap = Geps[rownames(Geps) %in% genes[[element]],]
  # to_heatmap1 = log2(to_heatmap+1)
  # to_heatmap1 = to_heatmap1 - rowMeans(to_heatmap1, na.rm = T)
  # minv = min(to_heatmap1, na.rm = T)
  # maxv = max(to_heatmap1, na.rm = T)
  # print(heatmaply(to_heatmap1, dendrogram = "row", main = element, fontsize_row = 15, fontsize_col = 15, margins = c(80,120,NA,NA),
  #                 scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", na.value = "black" ,midpoint = minv + ((maxv-minv)/2), 
  #                                                                         limits = c(minv, maxv))))
  # 
  # print(heatmaply(to_heatmap1, dendrogram = "row", main = element, fontsize_row = 15, fontsize_col = 15, margins = c(80,120,NA,NA),
  #                 scale_fill_gradient_fun = ggplot2::scale_fill_gradient(low = "cornflowerblue", high = "red",na.value = "black",midpoint = minv + ((maxv-minv)/2))))
  
  to_heatmap = as.matrix(to_heatmap)
  heatmap.2(to_heatmap, scale="row", na.color = "gray",
            trace="none", Rowv = TRUE, Colv = NA, margins=c(6,8), 
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
            main= element)
  
  
  
  to_barplot = data.frame(gene.names = rownames(to_heatmap), to_heatmap)
  long_barplot = melt(to_barplot, id.vars = "gene.names")
  colnames(long_barplot) = c("gene.names","stage","gene.count")
  print(ggplot(long_barplot,aes(gene.names,gene.count, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+ labs(y = "normalized expression")+
          theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+ggtitle(element))
}



####KEGG enrichment
KEGG_res = list()
for(element in names(to_GO_EI)){
  KEGG_res[[element]] = enrichKEGG(gene         = to_GO_EI[[element]]$ENTREZID,
                                   organism     = 'hsa',
                                   pvalueCutoff = 0.05)
}

#path = setwd("~/Desktop/breast_cancer_progression/plots/BC_30_10_17/enrichment_analysis")
for(element in names(KEGG_res)){
  #png(paste0(path,"/","KEGG_MO_",element,".png"), width = 900, height = 1100)
  plot(dotplot(KEGG_res[[element]],showCategory=50)+labs(title= paste0("KEGG_Macrophage_",element)))
  #dev.off()
}
