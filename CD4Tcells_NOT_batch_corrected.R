#### Magdalena Matusiak
#### May 15 2108 
#### analysis of immune landscape changes in breast cancer progression
#### LCM bulk RNASeq  deconvoluted with CibersortX
#### 2016 and 2017 batches data DESeq2 size factor normalized and NOT batch corrected cause we concluded analysing distrigution of notmal cances that we should not correct 

#Tcells separately CD8 and CD4


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


path = "~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/2016_2017_group_level_CBx_NOT_BC/Filtered/"
setwd(path)
files <- list.files(path=path)
files
#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_2016_2017_NOT_BC_",".txt_GEPs_Filtered.txt"),c("",""),i)
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
length(universe)
#extract genes in monocyte signature common for all samples 
colnames(filteredGEP$DCIS)
Tcell_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Tcell = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c("GeneSymbol","CD4Tcells")] #select gene_names and Tcells
  Tcell_ex = cbind(Tcell_ex,Tcell)
}

all(as.character(Tcell_ex[,2])==as.character(Tcell_ex[,4]))
rownames(Tcell_ex) = as.character(Tcell_ex[,2])
Tcell_ex = Tcell_ex[, seq(3, length(Tcell_ex), 2)]
colnames(Tcell_ex) = c("AVL", "DCIS","EN","IDC","normal")


Tcell_ex[Tcell_ex<=1] <- NA
predicted_genes_Tcells = apply(Tcell_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_Tcells
##############
#StErrors
#############

path = "~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/2016_2017_group_level_CBx_NOT_BC/StErr/"

setwd(path)
files <- list.files(path=path)
files
#extract names
namesStEr = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_2016_2017_NOT_BC_",".txt_GEPs_StdErrs.txt"),c("",""),i)
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
Tcell_ster = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesStEr){
  Error = StEr[[i]][StEr[[i]][,1] %in% universe,c("GeneSymbol","CD4Tcells")] #select gene_names and Monocytes
  Tcell_ster = cbind(Tcell_ster,Error)
}

all(as.character(Tcell_ster[,2])==as.character(Tcell_ster[,4]))
rownames(Tcell_ster) = as.character(Tcell_ster[,2])
Tcell_ster = Tcell_ster[, seq(3, length(Tcell_ster), 2)]
colnames(Tcell_ster) = c("AVL","DCIS","EN","IDC","normal")

all(rownames(Tcell_ex) == rownames(Tcell_ster))

###########DE Chloe

Tcell_ex = Tcell_ex[,c(5,3,2,4,1)]
Tcell_ster = Tcell_ster[,c(5,3,2,4,1)]

Geps = Tcell_ex
dim(Geps)
Geps = Geps[rowSums(!is.na(Geps)) >= 2,]
StdErr = Tcell_ster[rownames(Tcell_ster) %in% rownames(Geps),]
all(rownames(Geps) == rownames(StdErr))

Zqvals = matrix(NA,nrow(Geps),1)
cnames = c()
for(i in seq(names(Geps))){
  for(j in seq(names(Geps))){
    if(i >= j){
      next
    }
    cnames= c(cnames, paste0(names(Geps)[i],"_",names(Geps)[j]))
    vBetaZ <- sapply(1:nrow(Geps), function(x) (Geps[x,i]-Geps[x,j])/sqrt(StdErr[x,i]^2+StdErr[x,j]^2)) #order does not matter
    ZPs <- 2*pnorm(-abs(vBetaZ))
    Zqvals <- cbind(Zqvals,p.adjust(ZPs, method="BH")) #### shouldn't I take NA into account? cause there is not so many comparisons as genes?!
    
  }
} 

Zqvals = Zqvals[,-1]
colnames(Zqvals) = cnames
rownames(Zqvals) = rownames(Geps)
dim(Zqvals)
head(Zqvals)
Zqvals = data.frame(Zqvals)

#extract only significant pValues
sPv = list()
for(i in colnames(Zqvals)){
  sPv[[i]] = data.frame(genes = rownames(Zqvals)[which(Zqvals[,i]<= 0.25)], Zqv = Zqvals[which(Zqvals[,i]<= 0.25),i, drop = FALSE])
}

#calculate all log2FoldChange
LFC = list()
for(i in names(Zqvals)){
  fr = strsplit(i, "[_]")[[1]][1]
  sd = strsplit(i, "[_]")[[1]][2]
  LFC[[i]] = data.frame(names = rownames(Geps), FC = log2(Geps[,sd]/Geps[,fr]))
}


#merge LFC and p.adj for downstream GSEA ALL GENES
for(i in names(LFC)){
  print(all(as.character(LFC[[i]][,"names"]) == rownames(Zqvals)))
}

to_GSEA = list()
for(i in names(LFC)){
  to_GSEA[[i]] = cbind(LFC[[i]], Zqv = Zqvals[,i])
}

for(i in names(to_GSEA)){
  to_GSEA[[i]] = to_GSEA[[i]][complete.cases(to_GSEA[[i]]),]
}

#log10 plots to visualize gene expression 
for(i in names(to_GSEA)){
  fr = strsplit(i, "[_]")[[1]][1]
  sd = strsplit(i, "[_]")[[1]][2]
  plot(log10(Geps[to_GSEA[[i]][,1],sd]),log10(Geps[to_GSEA[[i]][,1],fr]), xlim=c(0, 5), ylim=c(0, 5),
       xlab = "log10(normalized expression)", ylab = "log10(normalized expression)")
  title(i)
  abline(0,1)
}

#----------------------------------------------------------------------
#vulcano plot
for(i in names(to_GSEA)){
  to_GSEA[[i]] = data.frame(to_GSEA[[i]], ml10pv = -log10(to_GSEA[[i]][,3]), ord = 1:length(to_GSEA[[i]][,3]))
}

for(i in names(to_GSEA)){
  plot(to_GSEA[[i]][, 2], to_GSEA[[i]][, 4],
       xlab = paste0("log2(FC_",i,")"), ylab = "-log10(Zqv)")
}


#plot interactive vulcano plots-----------------------------------------------------------------------------------------------
library(reshape2)
library(plotly)
library(magrittr)
setwd("~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/Tcells_plots/")

# "GZMA","GZMB","ICOS","CD69","CTLA4","ITGAE",
# "FOXP3","IL10","TIGIT","IDO1","PDCD1","LAG3",
# "IL19","IL18RAP","CD8A","CD8B","CD151","CD9",
# "IL10","CTSB","CLIC1","CLIC3","ZAP70","TIMP4"
# 
# "CD40LG","CD69","CD8A","FASLG","ICOS","IFNG","GZMK","TBX21","TNFRSF9","TNF","CD151","IL18RAP","Il19","ZAP70"
# "CTLA4","ITGAE","IL10","LAG3","IL12A"

for(i in names(to_GSEA)){
  plotd = data.frame(LogFC = to_GSEA[[i]][, 2], LZqv = to_GSEA[[i]][, 4], gnames = to_GSEA[[i]][[1]])
  m = plotd[plotd$gnames %in% c("SOX4","RAG1","BATF3","CXCL14","CD40LG","CD69","CD8A","FASLG","ICOS","IFNG","GZMK", 
                                "CD3E","CXCR1","CASP8","IL7R","IL11RA", "TNFSF8",
                                "TBX21","TNFRSF9","TNF","CD151","IL18RAP","IL19","ZAP70",
                                "CTLA4","ITGAE","IL10","LAG3","IL12A", "CD2","TNFRSF18","TNFRSF4",
                                "GZMA","GZMB","ICOS","CD69","CTLA4","ITGAE",
                                "FOXP3","IL10","TIGIT","IDO1","PDCD1","LAG3",
                                "IL19","IL18RAP","CD8A","CD8B","CD151","CD9",
                                "IL10","CTSB","CLIC1","CLIC3","ZAP70","TIMP4",
                                "CD40LG","CD69","CD8A","FASLG","ICOS","IFNG","GZMK","TBX21","TNFRSF9","TNF","CD151","IL18RAP","Il19","ZAP70",
                                "CTLA4","ITGAE","IL10","LAG3","IL12A"),]
  a = list(x = m$LogFC, y = m$LZqv, text = m$gnames)
  p = plot_ly(data = plotd, x = ~LogFC, y = ~LZqv, text = ~ gnames, marker = list(size = 20) ) %>%
    layout(title = i, xaxis = list(title = paste0("log2(FC_",i,")")), yaxis = list(title = "-log10(Zqv)"),
           annotations = a)
  print(p)
  #export(p, paste0(i,"new_selected_10FDR_vulcano.png"))
}

library(ggplot2)
library(ggrepel)

for(i in names(to_GSEA)){
  to_GSEA[[i]]$Selected <- ifelse(to_GSEA[[i]]$Zqv < 0.05 , "FDR < 0.05", "Not Sig")
  path = "~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/Tcells_plots/CD4/vulcano_repeled/"
  png(paste0(path,i,".png"), width = 1200, height = 1200)
  p= ggplot(data.frame(to_GSEA[[i]]), aes(x = FC, y = ml10pv)) +
    geom_point(aes(color = Selected)) +
    scale_color_manual(values = c("red", "grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = subset(data.frame(to_GSEA[[i]]), to_GSEA[[i]]$Zqv < 0.05),
      aes(label = names),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    ggtitle(i) + labs(x = "LFC", y = "-log10 qvalue")
  plot(p)
  dev.off()
}



subd = Geps[rownames(Geps) %in% c("IL7R","CD25","CCR8","IDO1","TIM3","CD39","PDCD1LG2","CD274","CXCL13","CD204","CD205",
                                  "CXCL14","GZMA","GZMB","ICOS","CD69","CTLA4","ITGAE","IL7R",
                                  "FOXP3","IL10","TIGIT","IDO1","PDCD1","LAG3",
                                  "IL19","IL18RAP","CD8A","CD8B","CD151","CD9",
                                  "IL10","CTSB","CLIC1","CLIC3","ZAP70","TIMP4"),]
subd = data.frame(names = rownames(subd), subd)
l_subd = melt(subd, id.vars = "names")
names(l_subd) = c("gname","stage","expr")
p = plot_ly(data = l_subd, x = ~stage, y = ~expr, color = ~gname, text = ~gname) %>% 
  add_lines() %>%  layout(showlegend = FALSE)  
p


for(i in names(to_GSEA)){
  write.table(to_GSEA[[i]][,c(1,4)], paste0(i,".rnk"), sep = "\t", quote = FALSE, row.names = FALSE)
}

for(i in names(to_GSEA)){
  print(plot(-log10(to_GSEA[[i]][,3])))
}

#extract only sig FC
LFC_s = list()
for(i in names(LFC)){
  LFC_s[[i]] = LFC[[i]][LFC[[i]][,"names"] %in% sPv[[i]][,"genes"],]
}

for(i in names(LFC_s)){
  print(all(as.character(LFC_s[[i]][,"names"]) == as.character(sPv[[i]][,"genes"])))
}

#merge LFC and p.adj for downstream GO/GSEA analysis SIGNIFICANT
to_enr = list()
for(i in names(LFC_s)){
  to_enr[[i]] = cbind(LFC_s[[i]], Zqv = sPv[[i]][,i])
}



##############
#plot interesting genes
##############


ImGenes = list()
ImGenes[["Treg_recruitment"]] = c("CCL5","CXCL14") #"CCR5"
ImGenes[["memory_survival"]] = c("IL7R")
ImGenes[["term_diff_inh_rud"]] = c("CTLA4") #,"TIGIT","FOXP3","IL2RA","CD39"
ImGenes[["term_diff_stm_rud"]] = c("CD2") #,"TNFRSF9","TNFRSF4", "TNFRSF9","TNFRS18"
ImGenes[["neg_sig"]] = c("CTLA4","LAG3","ITGAE")
ImGenes[["act_sig1"]] = c("CD2","CD8A","CD3E")
ImGenes[["act_sig2"]] = c("FASLG") #"CD151","IFNG","TBX21","BATF3"
ImGenes[["fate"]] = c("SOX4","RAG1","BATF3","CXCL14")
ImGenes[["act"]] = c("GZMA","GZMB","GZMK","ICOS","CD69","IL19","IL18RAP","CD8A","CD8B",
                     "CD151","CD9","CD32","ZAP70", "CCL4","CCL5","CXCR4","CCR5")
ImGenes[["migration"]] = c("CTSB","CLIC1","CLIC3","TIMP4")
ImGenes[["inh"]] = c("CTLA4","ITGAE","FOXP3","IL10","TIGIT","IDO1","PDCD1","LAG3","HAVCR2") 
ImGenes[["Treg"]] = c("CTLA4","ITGAE","IL10","LAG3","HAVCR2")
ImGenes[["diff"]] = c("TNFRSF9","TNFRSF18","CTLA4","CD2")
ImGenes[["Tcell_activation"]] = c("CD40LG","CD8A","FASLG","ICOS","IFNG","GZMK","TBX21","TNFRSF9","TNF","CD151","IL18RAP","Il19","ZAP70")
genes = c("Treg_recruitment","memory_survival","term_diff_inh_rud","term_diff_stm_rud","neg_sig","act_sig1","act_sig2",
          "fate","act","migration","inh","Treg","Tcell_activation","diff")


library(reshape2)
for(i in genes){
  Immuno = Tcell_ex[rownames(Tcell_ex) %in% ImGenes[[i]],]
  Immuno = data.frame(gene.names = rownames(Immuno), Immuno)
  #Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  Immuno = Immuno[,-6]
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.expression")
  er = Tcell_ster[rownames(Tcell_ster) %in% Immuno$gene.names,]
  er = data.frame(gene = rownames(er), er)
  er = er[,-6]
  #er$gene = as.factor(as.character(er$gene))
  er_long = melt(er, id.vars = "gene")
  long_Imm = data.frame(long_Imm, sd = er_long$value)
  long_Imm$sd[is.na(long_Imm$gene.expression)] = NA
  print(ggplot(long_Imm,aes(gene.names,gene.expression, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+
          theme(text = element_text(size=20))+ 
          geom_errorbar(aes(ymin=gene.expression-sd, ymax=gene.expression+sd), width=.2, position=position_dodge(.9)) +
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red","palevioletred3")) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1))+
          ggtitle(i))
}


####BC GO enrichment

to_GO = list()
for(nm in names(to_enr)){
  to_GO[[paste0("DE_",nm,"_up")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] > 0), "names"]
  to_GO[[paste0("DE_",nm,"_down")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] < 0), "names"]
}

lengths(to_GO)[-c(7,8,13,14,17,18,19,20)]

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


path = "~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/Tcells_plots/CD8/GO_enrichment/"
for(element in names(GO_res)){
  png(paste0(path,element,".png"), width = 900, height = 1100)
  plot(dotplot(GO_res[[element]],showCategory=50)+labs(title= paste0("CD8Tcell_",element)))
  dev.off()
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
write.table(to_REV1, "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/DCIS_IDC_up.txt", sep = "\t", row.names = F, quote = FALSE)

######disect genes for ploting

to_plot = list()
for(i in c("oxygen|hypoxia|oxydative","MHC|antigen|present","neutroph","extracellular|collagen",
           "folding|incorrect|unfolded","p53|apoptosis","cytokine|interferon","actin|Wnt|polarity",
           "stress|apoptosis","catabolic", "metabolic|lipoprotein|lipid|cholesterol|sterol")){
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
names(genes) = c("responce_to_hypoxia","Antigen_presentation","neutrophil_responce","Extracellular_matrix_organization",
                 "Unfolded_protein_responce","DNA_damage_apoptosis","IFNI_cytokine_regulation","Cell_polarity_Wnt",
                 "stress_responce","catabolic_processes","metabolic_lipoprotein")

####-----------------------------------------------------------------------
for(element in names(genes)){
  to_heatmap = Geps[rownames(Geps) %in% genes[[element]],]
  to_heatmap = as.matrix(to_heatmap)
  path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/heatmaps/"
  png(paste0(path,element,"_scalled_no_dendro.png"), width = 1000, height = 1200)
  heatmap.2(to_heatmap, scale = "row", na.color = "gray",
            trace="none", Rowv = F, Colv = F, dendrogram = "none", margins=c(6,8),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
            main= element)
  dev.off()
  
  path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/barplots/"
  to_barplot = data.frame(gene.names = rownames(to_heatmap), to_heatmap)
  long_barplot = melt(to_barplot, id.vars = "gene.names")
  colnames(long_barplot) = c("gene.names","stage","gene.count")
  to_barplot$gene.names = sub("-",".",to_barplot$gene.names)
  er = StdErr[rownames(StdErr) %in% to_barplot$gene.names,]
  er = data.frame(gene = rownames(er), er)
  er_long = melt(er, id.vars = "gene")
  long_barplot = data.frame(long_barplot, sd = er_long$value)
  long_barplot$sd[is.na(long_barplot$gene.count)] = NA
  png(paste0(path,element,"_woB2M_10FDR.png"), width = 1200, height = 900)
  print(ggplot(long_barplot,aes(gene.names,gene.count, group = stage, fill =stage)) +
          geom_col(position = "dodge")+ theme_classic()+ labs(y = "normalized expression")+
          geom_errorbar(aes(ymin=gene.count-sd, ymax=gene.count+sd), width=.2, position=position_dodge(.9)) +
          theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+ggtitle(element))
  dev.off()
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
