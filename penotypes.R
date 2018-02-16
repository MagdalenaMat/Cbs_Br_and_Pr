phenotypes = read.csv("~/Desktop/breast_cancer_progression/table5.csv")
phenotypes[phenotypes == ""] = NA 
phenotypes = as.list(phenotypes, na.rm = T)
phenotypes = lapply(phenotypes, na.omit)
phenotypes = lapply(phenotypes, toupper)

signatures = Geps[rownames(Geps) %in% unlist(phenotypes), ]
signatures = data.frame(gene.names = rownames(signatures), signatures)

pheno = as.character()
for(name in names(phenotypes)){
  ph = rep(name, length(phenotypes[[name]]))
  pheno = c(pheno, ph)
}
all_phen = data.frame(genes = unlist(phenotypes), pheno = as.factor(pheno))
sig = merge(signatures, all_phen, by.x = "gene.names", by.y = "genes")

sig = sig[!duplicated(sig$gene.names),]
sig = sig[rowSums(!is.na(sig[,2:5])) >= 2,]

for(phen in levels(sig$pheno)){
  to_heatmap = sig[sig$pheno == phen, c(2:5)]
  to_heatmap = as.matrix(to_heatmap)
  if(nrow(to_heatmap) >=1){
    heatmap.2(to_heatmap, scale="row", na.color = "gray",
              trace="none", Rowv = FALSE, Colv = FALSE, margins=c(6,8), 
              col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
              main= element)
  }
}

pheno2 = as.character()
for(name in names(phenotypes)[6:7]){
  ph = rep(name, length(phenotypes[[name]]))
  pheno2 = c(pheno2, ph)
}

pheno1 = as.character()
for(name in names(phenotypes)[1:5]){
  ph = rep(name, length(phenotypes[[name]]))
  pheno1 = c(pheno1, ph)
}


all_phen2 = data.frame(genes = unlist(phenotypes[6:7]), pheno = as.factor(pheno2))
all_phen1 = data.frame(genes = unlist(phenotypes[1:5]), pheno = as.factor(pheno1))


library(heatmaply)
for(element in names(phenotypes)){
  to_heatmap = sig[sig$gene.names %in% phenotypes[[element]], c("gene.names","normal","EN","DCIS","IDC")]
  if(nrow(to_heatmap) >=1){
    rownames(to_heatmap) = to_heatmap$gene.names
    to_heatmap$gene.names = NULL
    to_heatmap1 = log2(to_heatmap+1)
    to_heatmap1 = to_heatmap1 - rowMeans(to_heatmap1, na.rm = T)
    minv = min(to_heatmap1, na.rm = T)
    maxv = max(to_heatmap1, na.rm = T)
    print(heatmaply(to_heatmap1, margins = c(50,80,NA,0), dendrogram = "none", main = element, #scale = "row",
                    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = minv + ((maxv-minv)/2),
                                                                            limits = c(minv, maxv))))
    
  }
}


to_heatmap = change[change$gene.names %in% all_phen2$genes, c("gene.names","normal","EN","DCIS","IDC")]
to_heatmap = merge(all_phen2, to_heatmap, by.x = "genes", by.y = "gene.names")
to_heatmap = to_heatmap[order(to_heatmap$pheno),]
to_heatmap =  to_heatmap[!duplicated(to_heatmap[,1]),]
rownames(to_heatmap) = to_heatmap$genes
to_heatmap$genes = NULL
draw = log2(to_heatmap[,2:5]+1)
#draw = apply(draw, 1, function(x){scale(x, center = F, scale = T)})
#draw1 = t(scale(t(draw), center = FALSE, scale = apply(draw, 1, min, na.rm = TRUE)))
draw = draw - rowMeans(draw, na.rm = T)
minv = min(draw, na.rm = T)
maxv = max(draw, na.rm = T)

print(heatmaply(draw, margins = c(50,80,NA,0), dendrogram = "none", main = element, 
                row_side_colors = to_heatmap$pheno, 
                scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                                                        limits = c(minv, maxv)), na.value = "black"))


########Macrophage markers

markers = c("SIGLEC1","APOE","TREM2","CEBPB","CEBPD","CD81","SMARCA2")
plot_MO_m = change[change$gene.names %in% markers, c("gene.names","normal","EN","DCIS","IDC")]
plot_MO_m$gene.names = as.factor(as.character(plot_MO_m$gene.names))
long_barplot = melt(plot_MO_m, id.vars = "gene.names")
colnames(long_barplot) = c("gene.names","stage","gene.count")
print(ggplot(long_barplot,aes(gene.names,gene.count, group = stage, fill =stage)) + 
        geom_col(position = "dodge")+ theme_classic()+
        theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+ggtitle("MO_markets"))



test= change[!is.na(change$DCIS_IDC),c("gene.names","DCIS","IDC","DCIS_IDC")]
test2 = data.frame(gene.names = test$gene.names,scale(test[,2:3], center = TRUE, scale = TRUE), DCIS_IDC = test$DCIS_IDC)
test3 = mutate(test2, FC = abs(test2$DCIS_IDC) > 0.5)
plot_ly(data = test3, x= ~log10(test3$DCIS+1), y= ~log10(test3$IDC+1), type = 'scatter', text= ~test3$gene.names, color = ~test3$FC )
