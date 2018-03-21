######heatmap
cols1<- palette(brewer.pal(8, "Dark2"))
cols = cols1[desine2$subtype]
cols = cols1[cellType]
colnames(conv_dat)

to_heatmap = conv_dat[rownames(conv_dat) %in% c("TNF","IL1B","IL1RN","TGFA","TGFB1","TGFB2"),]

to_heatmap = as.matrix(to_heatmap)
row_cor_to_heatmap = as.dist(1-cor(t(to_heatmap),use="na.or.complete",method = "spearman"))
hcl_row = hclust(row_cor_to_heatmap, "ward.D")
heatmap.2(log2(to_heatmap+1), trace='none', na.color = "gray",margins=c(8,8),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          Rowv=NA,
          Colv = NA, 
          ColSideColors = cols,
          scale = "row",
          dendrogram = "row")

heatmap.2(log2(to_heatmap+1), trace='none', na.color = "gray",margins=c(8,8),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          distfun=function(x) dist(x, method='manhattan'), hclustfun=function(x)hclust(x,method="ward.D"),
          Colv = T, 
          ColSideColors = cols,
          scale = "row",
          dendrogram = "row")


######
rownames(ex3)[grep("TGF", rownames(ex3))]

desine2$subtype = factor(desine2$subtype, levels = c("normal","EN","DCIS","IDC"))
desine2$name
