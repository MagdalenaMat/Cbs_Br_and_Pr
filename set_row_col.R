mtcars = mtcars

library(gplots)
library(RColorBrewer)
cols1<- palette(brewer.pal(8, "Dark2"))
cols = cols1[mtcars$cyl]


to_plot = as.matrix(mtcars)
heatmap.2(to_plot, 
          trace='none',
          RowSideColors = cols)

library(heatmaply)
heatmaply(to_plot, RowSideColors = cols)

library(pheatmap)
mydf <- data.frame(row.names = rownames(to_plot), category = mtcars$cyl)
pheatmap(to_plot, 
         annotation_row = mydf)
