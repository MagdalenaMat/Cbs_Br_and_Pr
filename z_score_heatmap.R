library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )

mono_DCIS = read.table("~/Desktop/results/CIBERSORTxHiRes_Job4_Monocytes_Window40.txt", header = T)
rownames(mono_DCIS) = mono_DCIS[,1]
mono_DCIS = mono_DCIS[,-1]
M1_M2 = phenotypes[6:7]
mon_filter = unlist(M1_M2)

mono_DCIS = mono_DCIS[rownames(mono_DCIS) %in% mon_filter,]
mono_DCIS = mono_DCIS[complete.cases(mono_DCIS),]

ph = data.frame(gene = NA, pheno = NA)
for(i in names(M1_M2)){
  for(j in M1_M2[[i]]){
    ph = rbind(ph, c(gene = j, pheno =  i)) 
  }
}
ph = ph[-1,]

pdata = data.frame(gene = NA, pheno = NA)
for(i in rownames(mono_DCIS)){
  pdata = rbind(pdata, c(gene = i, pheno = ph[ph$gene == i, 2]))
}
pdata = pdata[-1,]


mono_DCIS = round(mono_DCIS, 2)
mono.DCIS = as.matrix(mono_DCIS)
filt = apply(mono_DCIS - rowMeans(mono_DCIS),1,sum) != 0
to_hm = mono.DCIS[filt,]
rownames(mono.DCIS) == pdata$gene
pdata = pdata[filt,]
topVarGenes <- head( order( rowVars( to_hm ), decreasing=TRUE ), 50 )


hc.g = hclust(as.dist(1-cor(t(to_hm))))
hc.p = hclust(as.dist(1-cor(to_hm)))
plot(hc.p)

heatmap.2( to_hm, scale="row", 
           trace="none", Rowv = as.dendrogram(hc.g), Colv = NA, 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           RowSideColors = c(M1 ="gray", M2 ="darkgreen")[
            pdata$pheno])



library(heatmaply)
heatmaply(mono.DCIS[filt,], margins = c(50,80,NA,0), scale = "row")

ncol(to_hm)

