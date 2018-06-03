#macrophages test if same terms are enriched when I changed megaclasses

#### Magdalena Matusiak
#### NOT batch corrected
#### analysis of immune landscape changes in breast cancer progression
#### LCM bulk RNASeq deconvoluted with CibersortX
#### 8 batches data DESeq2 size factor normalized



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


path = "~/Desktop/s7s_Breast/data_after_decon/CD8CD4_separate_run_not_b_c/Filtered/"
setwd(path)
files <- list.files(path=path)
files
#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_20180314_not_batch_corrected_",".txt_GEPs_Filtered.txt"),c("",""),i)
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
genecount = gene_count[c(4,2,1,3),]

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
Mono_ex = data.frame(matrix(NA, nrow = length(universe), ncol = 1))
for (i in namesGEP){
  Mono = filteredGEP[[i]][filteredGEP[[i]][,1] %in% universe,c("GeneSymbol","Monocytes")] #select gene_names and Monocytes
  Mono_ex = cbind(Mono_ex,Mono)
}

all(as.character(Mono_ex[,2])==as.character(Mono_ex[,4]))
rownames(Mono_ex) = as.character(Mono_ex[,2])
Mono_ex = Mono_ex[, seq(3, length(Mono_ex), 2)]
colnames(Mono_ex) = c("DCIS","EN","IDC","normal")


Mono_ex[Mono_ex<=1] <- NA
predicted_genes_Monocytes = apply(Mono_ex, 2, function(x){sum(!is.na(x))})
predicted_genes_Monocytes
##############
#StErrors
#############

path = "~/Desktop/s7s_Breast/data_after_decon/CD8CD4_separate_run_not_b_c/StErr/"

setwd(path)
files <- list.files(path=path)
files
#extract names
namesStEr = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_20180314_not_batch_corrected_",".txt_GEPs_StdErrs.txt"),c("",""),i)
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
  Error = StEr[[i]][StEr[[i]][,1] %in% universe,c("GeneSymbol","Monocytes")] #select gene_names and Monocytes
  Mono_ster = cbind(Mono_ster,Error)
}

all(as.character(Mono_ster[,2])==as.character(Mono_ster[,4]))
rownames(Mono_ster) = as.character(Mono_ster[,2])
Mono_ster = Mono_ster[, seq(3, length(Mono_ster), 2)]
colnames(Mono_ster) = c("DCIS","EN","IDC","normal")

all(rownames(Mono_ex) == rownames(Mono_ster))

# Mono_ster$gene = sub("HLA.","HLA-", Mono_ster$gene) 

###########DE Chloe

Mono_ex = Mono_ex[,c(4,2,1,3)]
Mono_ster = Mono_ster[,c(4,2,1,3)]

Geps = Mono_ex
dim(Geps)
Geps = Geps[rowSums(!is.na(Geps)) >= 2,]
StdErr = Mono_ster[rownames(Mono_ster) %in% rownames(Geps),]
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


#ranked list of genes by FC
FC_rnk = list()
for(i in names(to_GSEA)){
  FC_rnk[[i]] = to_GSEA[[i]][order(to_GSEA[[i]][,"FC"],decreasing = T),]
}


vec_GSEA = list()
for(i in names(FC_rnk)){
  vec_GSEA[[i]] = FC_rnk[[i]][,"FC"]
  names(vec_GSEA[[i]]) = FC_rnk[[i]][,"names"]
}

####GSEA for GO terms check BP MF
resGSEA = list()
for(element in names(vec_GSEA)){
  resGSEA[[element]] = gseGO(geneList     = vec_GSEA[[element]],
                             OrgDb        = org.Hs.eg.db,
                             ont          = "BP",
                             keyType = "SYMBOL",
                             nPerm        = 1000,
                             minGSSize    = 20,
                             maxGSSize    = 500,
                             pvalueCutoff = 0.1,
                             verbose      = FALSE)
}

library(ggplot2)
for(element in names(resGSEA)){
  if(nrow(resGSEA[[element]]@result) > 0){
    plot(dotplot(resGSEA[[element]],showCategory=50) + labs(title= paste0("GO_MF_GSEA_",element)))
  }
}


########draw GSEA plot
gseaplot(resGSEA$DCIS_IDC, geneSetID = "GO:0001568",title = "blood vessel development")
gseaplot(resGSEA$DCIS_IDC, geneSetID = "GO:0002683",title = "negative regulation of immune system process")


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
setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/volcano/")
#"ENO1","P4HB", "CLTB","TAP2","HLA.C",
#"B2M","GAPDH","HK3", "BAX","C4BPA",
#"TREM2", "IDO1","CTSD","UBC","UBB","ANXA2"

#"ATM","PML","DDX5","DROSHA","BAX","FOXO3A","FOXO3","CASP8", "TP53I3"

for(i in names(to_GSEA)){
  plotd = data.frame(LogFC = to_GSEA[[i]][, 2], LZqv = to_GSEA[[i]][, 4], gnames = to_GSEA[[i]][[1]])
  m = plotd[plotd$gnames %in% c("ENO1","P4HB", "CLTB","TAP2","HLA.C","MMP9","MMP13","MMP14","MMP1","MMP11", "CASP8",
                                "B2M","GAPDH", "BAX","C4BPA", "STAT6", "TLR5",
                                "IDO1","CTSD","UBC","UBB","ANXA2",
                                "TREM2", "APOE", "CD68","CHIT1", #activation
                                "MARCO","FN1","NRP2","SPP1","CD276", #M2
                                "CD64","FN1","MSR1","MMP14","CTSD", "CTSA","CTSC","CTSD","CTSK","CTSO","MARCO","VEGFB", #co-vary TAM in rudenski
                                "IL4R", "ITGAL", "CD74", "HLA.DRA", "ISG15", "CD14",#monocytes
                                "FGF1", "FGF12", "FGF13", "FGF14", "FGF17", "FGF18", "FGF2", "FGF20", "FGF21", "FGF22","FGF23","FGF7","FGF9","EGF","EGFR", "VEGFA","VEGFB","VEGFC",   
                                "TGFA", "TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBI", "TGFBR1", 
                                "TGFBR2", "TGFBR3","TGFBRAP1"),]
  a = list(x = m$LogFC, y = m$LZqv, text = m$gnames, font = list(size = 18))
  p = plot_ly(data = plotd, x = ~LogFC, y = ~LZqv, text = ~ gnames, marker = list(size = 20) ) %>%
    layout(title = i, xaxis = list(title = paste0("log2(FC_",i,")")), yaxis = list(title = "-log10(Zqv)"),
           annotations = a)
  print(p)
  #export(p, paste0(i,"_10FDR_vulcano.png"))
}

subd = Geps[rownames(Geps) %in% c("MMP9","MMP13","NFKBIA","CD163", "CD68","STAT3","DROSHA","CXCL10","IL12","IL13RA1","ENO1","P4HB","C4BPA"),]
subd = data.frame(names = rownames(subd), subd)
l_subd = melt(subd, id.vars = "names")
names(l_subd) = c("gname","stage","expr")
p = plot_ly(data = l_subd, x = ~stage, y = ~expr, color = ~gname, text = ~gname) %>% 
  add_lines() %>%  layout(showlegend = FALSE)  
p

#setwd("~/Desktop/s7s_Breast/data_after_decon/9_megaclasses/with_Megan_data/batch_corr/GSEA/input/")

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

############plot volcano of only significant LFC
for(i in names(to_enr)){
  to_enr[[i]] = data.frame(to_enr[[i]], ml10pv = -log10(to_enr[[i]][,"Zqv"]))
}

setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/volcano/")

for(i in names(to_enr)){
  plotd = data.frame(LogFC = to_enr[[i]][, "FC"], LZqv = to_enr[[i]][, "ml10pv"], gnames = to_enr[[i]][,"names"])
  m = plotd[plotd$gnames %in% c("ENO1","P4HB", "CLTB","TAP2","HLA.C","MMP9","MMP13","MMP14","MMP1","MMP11", "CASP8",
                                "B2M","GAPDH", "BAX","C4BPA","TREM2", "APOE", "STAT6", "TLR5",
                                "TREM2", "IDO1","CTSD","UBC","UBB","ANXA2",
                                "TGFA", "TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBI", "TGFBR1", 
                                "TGFBR2", "TGFBR3","TGFBRAP1"),]
  a = list(x = m$LogFC, y = m$LZqv, text = m$gnames, font = list(size = 18))
  p = plot_ly(data = plotd, x = ~LogFC, y = ~LZqv, text = ~ gnames, marker = list(size = 20) ) %>%
    layout(title = i, xaxis = list(title = paste0("log2(FC_",i,")")), yaxis = list(title = "-log10(Zqv)"),
           annotations = a)
  print(p)
  #export(p, paste0(i,"significant_10FDR_vulcano.png"))
}


##############
#plot interesting genes
##############
# "MMP9","MMP13","MMP11","MMP10","NFKBIA","CD163", "CD68","STAT3",
# "DROSHA","CXCL10","IL13RA1","ENO1","P4HB", "CLTB","TAP2",
# "ERAP1", "ERAP2","B2M","GAPDH","HK3",
# "TREM2","DAP12"

ImGenes = list()
ImGenes[["TAMactivation"]] = c("TREM2", "CHIT1", "TLR5") #activation
ImGenes[["M2_pol1"]] = c("TGFB1","IL4R","STAT6") #M2 polarizing
ImGenes[["M2_pol2"]] = c("NFKBIA","TGFBI")
ImGenes[["M2_rud"]] = c("FN1","NRP2") #M2
ImGenes[["angio1"]] = c("ANXA2","FN1")
ImGenes[["angio2"]] = c("MMP14","MMP2","VEGFB","NRP2")
ImGenes[["TAM_rud"]] = c("CD64","FN1","MSR1","MMP14","CTSD", "CTSA","CTSC","CTSD","CTSK","CTSO","MARCO","VEGFB") #co-vary TAM in rudenski
ImGenes[["Mono_rud"]] = c("HLA.DRA", "CD14") #monocytes
ImGenes[["E_G_F"]] = c("FGF","VEGFB","VEGFC")
ImGenes[["ER_stress"]] = c("P4HB","UBB","UBC","CTSD","PSMB4","PSMD11","PSMD14","PSMD2","PSME1")
ImGenes[["antigen"]] = c("B2M","HLA.A","HLA.B","HLA.C")
ImGenes[["TAP2"]] = c("TAP2","CD207","NOS3")
ImGenes[["M2"]] = c("TGFB1", "TGFBI","NFKBIA")
ImGenes[["phagocytosis"]] = c("TREM2","CD47")
ImGenes[["phagocytosis2"]] =c("CLTB","CTSD", "CTSA","CTSC","CTSD","CTSK","CTSO")
ImGenes[["glykolysis"]] = c("GAPDH","ENO1")
ImGenes[["STAT"]] = c("STAT3","STAT6")
ImGenes[["IFN1"]] = c("B2M","HLA.A","HLA.B","HLA.C","HLA.F")
ImGenes[["IFN2"]] = c("IFI35","IFIT3","IRF7","STAT1")
ImGenes[["apoptosis1"]] = c("ATM","PML","BAX")
ImGenes[["apaptosis2"]] = c("CASP8","TP53I3")
ImGenes[["MMPs"]] = c("MMP13","TIMP1","TIMP2","MMP14","MMP11") #"MMP9"
ImGenes[["oxidatie_stress"]] = c("FOXO3","BAX")
genes = c("TAMactivation","M2_rud","M2_pol1","M2_pol2",
          "angio1","angio2","Mono_rud","E_G_F","ER_stress",
          "antigen","TAP2","M2","phagocytosis","phagocytosis2", 
          "glykolysis","STAT","IFN1","IFN2","apoptosis1",
          "apaptosis2","MMPs","oxidatie_stress")

library(reshape2)
for(i in genes){
  Immuno = Mono_ex[rownames(Mono_ex) %in% ImGenes[[i]],]
  Immuno = data.frame(gene.names = rownames(Immuno), Immuno)
  #Immuno$gene.names = as.factor(as.character(Immuno$gene.names))
  print(Immuno)
  Immuno = Immuno[,c(1,3,4,5)]
  #Immuno[,c(2,3)] =  log2(Immuno[,c(2,3)]+1)
  long_Imm = melt(Immuno, id.vars = "gene.names")
  colnames(long_Imm) = c("gene.names","stage","gene.expression")
  er = Mono_ster[rownames(Mono_ster) %in% Immuno$gene.names,]
  er = data.frame(gene = rownames(er), er)
  er = er[,c(1,3,4,5)]
  #er[,c(2,3)] =  log2(er[,c(2,3)]+1)
  #er$gene = as.factor(as.character(er$gene))
  er_long = melt(er, id.vars = "gene")
  long_Imm = data.frame(long_Imm, sd = er_long$value)
  long_Imm$sd[is.na(long_Imm$gene.expression)] = NA
  
  print(ggplot(long_Imm,aes(gene.names,gene.expression, group = stage, fill =stage)) + 
          geom_col(position = "dodge")+ theme_classic()+
          theme(text = element_text(size=20))+ 
          geom_errorbar(aes(ymin=gene.expression-sd, ymax=gene.expression+sd), width=.2, position=position_dodge(.9)) +
          scale_fill_manual(values=c("orange","chartreuse4","red")) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1))+
          ggtitle(i))
}

#"dodgerblue3","orange",

test = to_enr$DCIS_IDC
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


path = setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/dotplots/GO/")
for(element in names(GO_res)){
  #png(paste0(path,"/",element,".png"), width = 900, height = 1100)
  plot(dotplot(GO_res[[element]],showCategory=50)+labs(title= paste0("Macrophage_",element)))
  #dev.off()
}

dev.off()

res = list()
for(element in names(GO_res)){
  if(nrow(GO_res[[element]]@result) == 0){
    next
  }
  res[[element]] = GO_res[[element]]@result[,c("ID","p.adjust", "Description","geneID","GeneRatio")]
}

#subset ID and p.adj for revigo analysis
to_REV1 = GO_res$DE_DCIS_IDC_up@result[,c("ID","p.adjust")]
write.table(to_REV1, "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/DCIS_IDC_up.txt", sep = "\t", row.names = F, quote = FALSE)

######disect genes for ploting

to_plot = list()
for(i in c("extracellular|movement|location|motility|migration","MHC|antigen|present","neutroph","collagen",
           "folding|unfolded","virus|interferon","axon|neuron|projection",
           "catabolic", "autophag","angiogen","ossification|skeletal","transform","ephrin")){
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
names(genes) = c("migration","Antigen_presentation","neutrophil_responce","Extracellular_matrix_organization",
                 "Unfolded_protein_responce","IFNI_cytokine_regulation","neuron",
                 "catabolic_processes","autophagy","angiogenesis","bone_rormation","TNFb_responce","ephrin")

####-----------------------------------------------------------------------
for(element in names(genes)){
  to_heatmap = Geps[rownames(Geps) %in% genes[[element]],]
  to_heatmap = as.matrix(to_heatmap)
  path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/heatmaps/"
  png(paste0(path,element,"_scalled_no_dendro.png"), width = 1000, height = 1200)
  heatmap.2(to_heatmap, scale = "row", na.color = "gray",
            trace="none", Rowv = F, Colv = F, dendrogram = "none", margins=c(6,8),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
            main= element)
  dev.off()
  
  path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/barplots/"
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

#path = setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/dotplots/KEGG/")
for(element in names(KEGG_res)){
  #png(paste0(path,"/","KEGG_MO_",element,".png"), width = 900, height = 1100)
  plot(dotplot(KEGG_res[[element]],showCategory=50)+labs(title= paste0("KEGG_Macrophage_",element)))
  #dev.off()
}

# pt_code = c()
# for(element in names(KEGG_res)){
#   pt_code = c(pt_code,KEGG_res[[element]]$ID)
# }

KEGG_res$DE_DCIS_IDC_up@result[,c("ID","p.adjust", "Description","geneID","GeneRatio")]
pathway = KEGG_res$DE_DCIS_IDC_up@result[,"ID"]
pathway = c(pathway,"hsa04612","hsa03050","hsa04140","hsa04612")

#geneList = merge(to_GO_EI$DE_DCIS_IDC_up,to_enr$DCIS_IDC,by.x = "SYMBOL", by.y = "names")

dcis_idc = to_enr$DCIS_IDC
dcis_idc$names = sub("[.]","-", dcis_idc$names)

library(AnnotationDbi)
library(org.Hs.eg.db)

dcis_idc$entrez = mapIds(org.Hs.eg.db,
                         keys=dcis_idc$names,
                         column=c("ENTREZID"),
                         keytype = "SYMBOL",
                         multiVals = "first")


pathList = dcis_idc$FC
names(pathList) = dcis_idc$entrez
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")

library(pathview)
setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/KEGG_pathways/run2/")
pathview(gene.data  = pathList,
         pathway.id = pathway,
         species    = "hsa",
         limit      = list(gene=max(abs(pathList)), cpd=1))

###KEGG Module

KEGG_mod = list()
for(element in names(to_GO_EI)){
  KEGG_mod[[element]] = enrichMKEGG(gene         = to_GO_EI[[element]]$ENTREZID,
                                   organism     = 'hsa')
}
to_GO_EI$DE_normal_EN_up$ENTREZID

#######MSigDB c7 immunological processes enrichment
library(GSEABase)
library("clusterProfiler")
c7 <- read.gmt("~/Desktop/s7s_Breast/gene_sets/c7.all.v6.1.entrez.gmt")

c7enr_res = list()
for(element in names(to_GO_EI)){
  c7enr_res[[element]] = enricher(gene         = to_GO_EI[[element]]$ENTREZID,
                                  TERM2GENE=c7)
}

for(element in names(c7enr_res)){
  #png(paste0(path,"/","KEGG_MO_",element,".png"), width = 900, height = 1100)
  plot(dotplot(c7enr_res[[element]],showCategory=50))
  #dev.off()
}


##GSEA

for(element in names(FC_rnk)){
  FC_rnk[[element]][,"entrez"] = mapIds(org.Hs.eg.db,
                         keys= as.character(FC_rnk[[element]][,"names"]),
                         column=c("ENTREZID"),
                         keytype = "SYMBOL",
                         multiVals = "first")
}

vec_GSEA_ez = list()
for(i in names(FC_rnk)){
  vec_GSEA_ez[[i]] = FC_rnk[[i]][,"FC"]
  names(vec_GSEA_ez[[i]]) = FC_rnk[[i]][,"entrez"]
}


c7_GSEA = list()
for(element in  names(vec_GSEA_ez)){
  c7_GSEA[[element]] = GSEA(vec_GSEA_ez[[element]],
        TERM2GENE=c7, 
        verbose=FALSE,
        nPerm        = 1000,
        minGSSize    = 10,
        maxGSSize    = 500,
        pvalueCutoff = 0.1)
}

library(ggplot2)
for(element in names(c7_GSEA)){
  if(nrow(c7_GSEA[[element]]@result) > 0){
    plot(dotplot(c7_GSEA[[element]],showCategory=50) + labs(title= paste0("MSigDB_c7_GSEA_",element)))
  }
}

