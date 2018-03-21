#rudenski_immune_phenotypes
#read.csv("~/Desktop/BC_ssRNAseq_rudenski/pheno_magda.csv")
phenotypes = read.csv("~/Desktop/signatures/Polyak/Polyak_phenotypes.csv", sep = "\t")
phenotypes[phenotypes == ""] = NA 
phenotypes = as.list(phenotypes, na.rm = T)
phenotypes = lapply(phenotypes, na.omit)
phenotypes = lapply(phenotypes, toupper)

for(element in names(phenotypes)){
  phenotypes[[element]] = sub("-",".",phenotypes[[element]])
}

#for all bartches analysis run with for(element in names(phenotypes)[-14]
#for analysis from 20180219 run with for(element in names(phenotypes)[-c(3,18,20)]
#for T cells
library(gplots)
for(element in names(phenotypes)[6:7]){
  to_heatmap = Geps[rownames(Geps) %in% phenotypes[[element]],]
  to_heatmap = as.matrix(to_heatmap)
  keep = complete.cases(to_heatmap[,3:4])
  to_heatmap = to_heatmap[keep,]
  #path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/phenotypes/"
  #png(paste0(path,element,"_log2_cor_ward.D_scaled.png"), width = 900, height = 1100)
  cor_to_heatmap = as.dist(1-cor(t(to_heatmap),use="na.or.complete",method = "spearman"))
  hcl_row = hclust(cor_to_heatmap, "ward.D")
  heatmap.2(log2(to_heatmap+1), trace='none', na.color = "gray",margins=c(6,8),main= element,
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
            Rowv=as.dendrogram(hcl_row),
            Colv=NA,
            scale = "row",
            dendrogram = "row")
  #dev.off()
  to_barplot = data.frame(gene.names = rownames(to_heatmap), to_heatmap)
  long_barplot = melt(to_barplot, id.vars = "gene.names")
  colnames(long_barplot) = c("gene.names","stage","gene.count")
  to_barplot$gene.names = sub("-",".",to_barplot$gene.names)
  er = StdErr[rownames(StdErr) %in% to_barplot$gene.names,]
  er = data.frame(gene = rownames(er), er)
  er_long = melt(er, id.vars = "gene")
  long_barplot = data.frame(long_barplot, sd = er_long$value)
  long_barplot$sd[is.na(long_barplot$gene.count)] = NA
  #png(paste0(path,element,"_woB2M_10FDR.png"), width = 1200, height = 900)
  print(ggplot(long_barplot,aes(gene.names,gene.count, group = stage, fill =stage)) +
          geom_col(position = "dodge")+ theme_classic()+ labs(y = "normalized expression")+
          geom_errorbar(aes(ymin=gene.count-sd, ymax=gene.count+sd), width=.2, position=position_dodge(.9)) +
          theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+ggtitle(element))
  #dev.off()
}
dev.off()

#for T cells
for(element in names(phenotypes)[c(1,2,3,4,5,14,16,17,18,19,20)]){
  to_heatmap = Geps[rownames(Geps) %in% phenotypes[[element]],]
  to_heatmap = as.matrix(to_heatmap)
  keep = complete.cases(to_heatmap[,4])
  to_heatmap = to_heatmap[keep,]
  #path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/phenotypes/"
  #png(paste0(path,element,"_log2_cor_ward.D_scaled.png"), width = 900, height = 1100)
  heatmap.2(log2(to_heatmap+1), trace='none', na.color = "gray",margins=c(6,8),main= element,
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),Colv = F, scale = "row",
            dendrogram = "none")
  #dev.off()
  to_barplot = data.frame(gene.names = rownames(to_heatmap), to_heatmap)
  long_barplot = melt(to_barplot, id.vars = "gene.names")
  colnames(long_barplot) = c("gene.names","stage","gene.count")
  to_barplot$gene.names = sub("-",".",to_barplot$gene.names)
  er = StdErr[rownames(StdErr) %in% to_barplot$gene.names,]
  er = data.frame(gene = rownames(er), er)
  er_long = melt(er, id.vars = "gene")
  long_barplot = data.frame(long_barplot, sd = er_long$value)
  long_barplot$sd[is.na(long_barplot$gene.count)] = NA
  #png(paste0(path,element,"_woB2M_10FDR.png"), width = 1200, height = 900)
  print(ggplot(long_barplot,aes(gene.names,gene.count, group = stage, fill =stage)) +
          geom_col(position = "dodge")+ theme_classic()+ labs(y = "normalized expression")+
          geom_errorbar(aes(ymin=gene.count-sd, ymax=gene.count+sd), width=.2, position=position_dodge(.9)) +
          theme(text = element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values=c("dodgerblue3","orange","chartreuse4","red"))+ggtitle(element))
  #dev.off()
}
