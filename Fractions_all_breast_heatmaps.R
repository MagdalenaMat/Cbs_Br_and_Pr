#analyse which samples for CibersortX run have detectable immune cell infiltrates
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


path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/Franctions/"
setwd(path)
files <- list.files(path=path)
files
#extract names
namesGEP = as.character()
for (i in files){
  nam = mgsub(c("CIBERSORTxGEP_20180314_not_batch_corrected_",".txt_Fractions-Adjusted.txt"),c("",""),i)
  namesGEP = c(namesGEP, nam)
}

namesGEP

#put all files in the list
Fractions = list()
for (i in 1: length(files)){
  r = read.table(files[i], sep = "\t", header = TRUE)
  Fractions[[i]] = r
}
names(Fractions) = namesGEP

for(element in names(Fractions)){
  rownames(Fractions[[element]]) = Fractions[[element]][,"Mixture"]
  Fractions[[element]] = Fractions[[element]][,2:23] 
}


#library(gplots)

# for(element in names(Fractions)){
#   to_plot = as.matrix(t(Fractions[[element]]))
#   heatmap.2(to_plot, margins = c(8,11),#scale = "row",
#             trace="none",distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x)hclust(x,method="ward.D"),
#             col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#             main= element)
# }


# mydf <- data.frame(row.names = rownames(to_plot), category = type)
# 
# # add row annotations
# myBreaks = unique(c(seq(0, 0.19, length=49), 0.2, seq(0.21,max(to_plot), length=50)))
# pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
#          annotation_row = mydf,gaps_row = c(24,32,48,74),
#          breaks = myBreaks, show_rownames = F)

library(pheatmap)
library(RColorBrewer)

all_fractions = matrix(NA, 1, 22)
colnames(all_fractions) = colnames(Fractions$DCIS)
for(element in c("normal","EN","DCIS","IDC")){
  all_fractions = rbind(all_fractions, Fractions[[element]])
}
all_fractions = all_fractions[-1,]
class(all_fractions)

to_plot = as.matrix(all_fractions)

stages = c()
for(i in c("normal","EN","DCIS","IDC")){
  stages = c(stages, rep(i,nrow(Fractions[[i]])))
  print(nrow(Fractions[[i]]))
}


stages = factor(stages, levels = c("normal","EN","DCIS","IDC"))

cols1<- palette(brewer.pal(8, "Dark2"))
cols = cols1[as.factor(stages)]

all_fractions = data.frame(all_fractions, m1m2 = all_fractions$Macrophages.M1+all_fractions$Macrophages.M2)
# to_plot = as.matrix(all_fractions)
# heatmap.2(to_plot,
#           margins = c(11,8), 
#           RowSideColors = cols, 
#           dendrogram = "col", 
#           #scale = "col", 
#           #Rowv = NA,
#           trace="none",
#           distfun=function(x) as.dist(1-cor(t(x))), 
#           hclustfun=function(x)hclust(x,method="ward.D"),
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#           main= "all", labRow = F)

mydf <- data.frame(row.names = rownames(to_plot), category = stages)

# add row annotations
#myBreaks = unique(c(seq(0, 0.19, length=49), 0.2, seq(0.21,max(to_plot), length=50)))
to_plot[to_plot >= 0.5] = 0.5
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
         annotation_row = mydf,gaps_row = c(44,71,242),
         #cellwidth = 4, cellheight = 4,
         show_rownames = F)

#########boxplots
class(all_fractions)
to_boxplot = as.matrix(all_fractions)
class(to_boxplot)
library(beeswarm)
boxplot(to_boxplot[,6] ~ stages,
        outline = FALSE, main = "T cells memory resting")    
beeswarm(to_boxplot[,6] ~ stages,
         pch = 20, 
         col = "lightseagreen",
         add = TRUE)

library(ggpubr)
to_ggplots = data.frame(to_boxplot, stages)
#colnames(to_ggplots)[9] = "Treg"
my_comparisons <- list( c("normal","EN"),c("EN", "DCIS"),c("DCIS","IDC"), 
                        c("normal", "DCIS"), c("EN","IDC"), c("normal", "IDC") )
path = "~/Desktop/s7s_Breast/to_fractions/boxplots/all_samples/"
for( i in colnames(to_ggplots)[1:23]){
  png(paste0(path,i,"400.png"), width = 400, height = 400)
  plot(ggboxplot(to_ggplots, x = "stages", y = i,  outlier.colour = NA,
                 color = "stages", palette = "jco", add = "jitter") + 
         theme(text = element_text(size=20)) +
         stat_compare_means(comparisons = my_comparisons))
  dev.off()
}


