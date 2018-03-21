to_heatmap = Geps[rownames(Geps) %in% phenotypes$M2.Macrophage.Polarization,]
to_heatmap = as.matrix(to_heatmap)
keep = complete.cases(to_heatmap[,3:4])
to_heatmap = to_heatmap[keep,]

# depending on the analysis, the data can be centered and scaled by row or column. 
# default parameters correspond to the ones in the heatplot function. 
distCor <- function(x) as.dist(1-cor(t(x)))
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  if (scale=="row") z <- t(scale(t(x)))
  if (scale=="col") z <- scale(x)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  hcl_row <- hclust(distCor(t(z)), method=method)
  return(list(data=z, Rowv=as.dendrogram(hcl_row)))
}

z <- zClust(to_heatmap)

library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.2(z$data, trace='none', col=rev(cols), Rowv=z$Rowv, Colv=NA)
#tu jestem

cor_to_heatmap = as.dist(1-cor(t(to_heatmap),use="na.or.complete",method = "spearman"))
hcl_row = hclust(cor_to_heatmap, "ward.D")
heatmap.2(scale(t(to_heatmap)), trace='none', na.color = "gray",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          Rowv=T, 
          Colv=NA,
          scale = "row",
          dendrogram = "row")
dev.off()


library(massageR)
library(gplots)

z <- heat.clust(to_heatmap,
                scaledim="row",
                zlim=c(-3,3),
                zlim_select = c("dend","outdata"),
                reorder=c("row"),
                distfun  = function(x) as.dist(1-cor(t(x))), 
                hclustfun= function(x) hclust(x, method="average"),
                scalefun = scale)


heatmap.2(z$data,
          Rowv=z$Rowv, 
          Colv=z$Colv,
          trace="none",
          scale="none",
          symbreaks = TRUE,
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          margins=c(4,7)
)

#range normalize function
rg_nl = function(data,x,y){
  m = min(data, na.rm = T)
  range = max(data, na.rm = T) - m
  data = (data - m)/range
  range2 = y - x
  normalized = (data * range2) + x 
  print(normalized)
}
rg_nl(to_heatmap, -1, 1)

########
cdist <- as.dist(1 - cor(t(to_heatmap), use="na.or.complete", method = "spearman"))
hc <- hclust(cdist, "average")

# draw a heatmap

heatmap.2(to_heatmap, na.color = "gray",
          Rowv=as.dendrogram(hc),
          Colv=NA,
          scale = "row",
          dendrogram="row",
          trace="none",
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          main = "scaled, hclust(1-cor(t(data))")

heatmap.2(log2(to_heatmap+1), na.color = "gray",
          trace="none", Rowv = T, Colv=NA, dendrogram = "row",margins=c(6,8),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          main= "log2(data+1),dist(), not scaled")

zlim=c(-3,3)
z <- pmin(pmax(log2(to_heatmap+1), zlim[1]), zlim[2])
heatmap.2(log2(to_heatmap+1), na.color = "gray",
          Rowv=as.dendrogram(hc),
          Colv=NA,
          dendrogram="row",
          trace="none",
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          main = "log2(data+1), hclust(1-cor(t(data))")


###############
to_heatmap = Geps[rownames(Geps) %in% genes$responce_to_hypoxia,]
to_heatmap = as.matrix(to_heatmap)
cdist <- as.dist(1 - cor(t(to_heatmap), use="na.or.complete", method = "spearman"))
hc <- hclust(cdist, "average")
heatmap.2(log2(to_heatmap+1), na.color = "gray",
          Rowv = T, Colv=NA,
          dendrogram="row",
          margins = c(6,6),
          trace="none",
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          main = "scaled/log2/cor/average/")
test = scale(t(to_heatmap), center = T, scale = F)

