plotdat = Geps[1:20,]
colnames(plotdat) = c("MO_normal","MO_EN","MO_DCIS","MO_IDC")
plotnorm = normal[1:20,]
library(heatmaply)
heatmaply(plotnorm,dendrogram = "none", margins = c(150,100,NA,0))
dev.off()
install.packages("heatmaply", )
