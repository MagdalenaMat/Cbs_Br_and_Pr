#immune cell fractions analysis form Cibersort

data = read.table("~/Desktop/s7s_Breast/data_after_decon/20180314_8batches/CIBERSORT.Output_Job13.txt",
                  quote = "", row.names = 1, header = T, sep = "\t")

colnames(data)[1:23]

data = data[,1:23]
data = data.frame(data, stage = NA)

for(i in c("DCIS|_CIS|LCIS", "EN|ALH|ADH|FEA|UDH|LobNeo|_75_L", "IDC|invasive|idc|ILC","normal|NL")){
  data$stage[grep(i,rownames(data))] = i
}

data$stage = gsub("DCIS\\|_CIS\\|LCIS","DCIS",data$stage)
data$stage = gsub("EN\\|ALH\\|ADH\\|FEA\\|UDH\\|LobNeo\\|_75_L","EN",data$stage)
data$stage = gsub("IDC\\|invasive\\|idc\\|ILC","IDC",data$stage)
data$stage = gsub("normal\\|NL","normal",data$stage)
data$P.value

library(ggplot2)

data$stage = factor(data$stage, levels = c("normal", "EN", "DCIS", "IDC"))

ggplot(data,aes(stage, P.value, fill = stage)) + geom_boxplot() + 
  geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
  geom_hline(yintercept=0.05, color = "red")

for(cell in cells){
  p = ggboxplot(res2, x = "stage", y = cell,outlier.colour = NA,
                color = "stage")+
    geom_jitter(position=position_jitter(width=.2, height=0),size = 1) +
    theme(panel.grid.major = element_blank(), legend.position = "none")
  my_comparisons <- list(c("EN","DCIS"),c("DCIS", "IDC"),c("EN","IDC"))
  print(p + stat_compare_means(comparisons = my_comparisons, size = 1))
}

sum(data$P.value <= 0.05)
subset = data[data$P.value <= 0.05,]

library(gplots)
library(RColorBrewer)
cols1<- palette(brewer.pal(8, "Dark2"))
cols = cols1[subset$stage]


to_plot = as.matrix(subset[1:22])
heatmap.2(to_plot, dendrogram = "none", 
          trace='none', 
          margins=c(12,10),
          col = colorRampPalette(c("darkblue", "brown1"))(100),
          Rowv=NA,
          Colv = NA, 
          RowSideColors = cols)

library(heatmaply)
heatmaply(to_plot, margins = c(200, 170), RowSideColors = cols,
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.1, limits = c(0, 1)))
