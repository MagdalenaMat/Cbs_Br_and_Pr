#CIDR imputes data and runs dimentionality reduction
path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/Filtered/"
setwd(path)
files <- list.files(path=path)
files
stages = c("DCIS", "EN","IDC","normal")
stag = rep(stages,each = 9)
stag = factor(stag, levels = c("normal","EN","DCIS","IDC"))
class(stag)
sub_to_SI = list()
for(i in 1: length(files)){
  r = read.csv(files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
  stage = r[r[,1] %in% universe,] 
  colnames(stage) = paste0(colnames(stage),rep(stages[i],9))
  sub_to_SI[[i]] = stage
}

for(i in 1:(length(sub_to_SI)-1) ){
  print(all(sub_to_SI[[i]][1] == sub_to_SI[[c(i+1)]][1]))
}

conv_dat = data.frame(names = sub_to_SI[[1]][1])
for(i in seq(length(sub_to_SI))){
  conv_dat = data.frame(conv_dat, sub_to_SI[[i]][,2:10])
}
rownames(conv_dat) = conv_dat$GeneSymbolDCIS
conv_dat$GeneSymbolDCIS = NULL

#conv_dat2 = conv_dat[rowSums(conv_dat, na.rm = T) > 10,]
conv_data = conv_dat[rowSums(!is.na(conv_dat)) > 8, ]

library(cidr)
cellType = colnames(read.csv(files[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""))[2:10]
cellType = rep(cellType,4)
cellType <- factor(cellType)
types <- levels(cellType)

## Assign each cell type a color
scols <- c("red","blue","green","brown","pink","purple","darkgreen","grey","yellow")
cols <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
}


conv_data[is.na(conv_data)] = 0
sum(is.na(conv_data))
sData <- scDataConstructor(as.matrix(conv_data))
sData <- determineDropoutCandidates(sData)
sData <- wThreshold(sData)
sData <- scDissim(sData)
sData <- scPCA(sData)
sData <- nPC(sData)
nCluster(sData)
sData <- scCluster(sData)
plot(sData@PC[,c(1,2)], col=cols,
     pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")

to_plot = sData@PC[,c(1,2)]
colnames(to_plot) = c("PC1","PC2")
to_plot = data.frame(to_plot, group = cellType, name = colnames(conv_data), stage = stag, size = colSums(conv_data))

library(ggplot2)
ggplot(to_plot, aes(PC1, PC2, color = group, shape = stage)) +
  geom_point(size = 7) + scale_color_brewer(palette="Paired") +
  scale_shape_manual(values=c(15:18)) +
  ggtitle("CIDR dim reduced") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20))


ggplot(to_plot, aes(PC1, PC2, color = size, shape = stage)) +
  geom_point(size = 7) + 
  scale_color_gradient(low="blue", high="red") +
  scale_shape_manual(values=c(15:18)) +
  ggtitle("CIDR dim reduced") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20))


##########Monocytes only
#conv_data = conv_dat[rowSums(!is.na(conv_dat2)) > 8, ]
Mono = conv_data[,grep("Monocytes",colnames(conv_data))]
colnames(Mono) = stages


sum(is.na(Mono))
## Assign each cell type a color
cols <- c("red","blue","green","brown")

sum(is.na(Mono))
sData <- scDataConstructor(as.matrix(Mono))
sData <- determineDropoutCandidates(sData)
sData <- wThreshold(sData)
sData <- scDissim(sData)
sData <- scPCA(sData)
sData <- nPC(sData)
nCluster(sData)
sData <- scCluster(sData)
plot(sData@PC[,c(1,2)], col=cols,
     pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")

to_plot = sData@PC[,c(1,2)]
colnames(to_plot) = c("PC1","PC2")
to_plot = data.frame(to_plot, group = colnames(Mono), name = colnames(Mono), size = colSums(Mono))
library(plotly)
plot_ly(data = to_plot, x = ~PC1, y = ~PC2, color = ~group, 
        text = ~ name, colors = "Set3", 
        marker = list(size = 20),
        autosize = F, width = 900, height = 600) %>%
  layout(title = "Monocytes CIDR dimentionality reduction",legend = list(orientation = 'h'))

ggplot(to_plot, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 7) + scale_color_brewer(palette="Set2") +
  scale_shape_manual(values=c(15:18)) +
  ggtitle("CIDR dim reduced") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20))

########tsne

library(parallel)
library(ggplot2)
library(Rtsne)
library(DESeq2)

options(mc.cores = detectCores())
conditions <- data.frame(group = cellType, name = colnames(conv_data), stage = stag, size = colSums(conv_data, na.rm = T))
ldata <- log2(conv_data+1) #directly put to tsne
dist.mat <- dist(t(ldata))
tsne.iter <- mclapply(5:10, function(p) Rtsne(dist.mat,
                                              perplexity = p,
                                              check_duplicates = F
))
tsne.frame.iter <- lapply(tsne.iter, function(tsne.result) {
  coords <- tsne.result$Y
  colnames(coords) <- c("x", "y")
  data.frame(coords, conditions)	
}
)
setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/")
pdf("tsne_after_deconv.pdf", width = 6, height = 5)
for (p in 1:length(tsne.iter)) {
  print(
    ggplot(tsne.frame.iter[[p]], aes(x, y, color = group, shape = stage)) + 
      geom_point(size = 5) +
      ggtitle(paste0("perplexity = ", p, ", cost =", sum(tsne.iter[[p]]$costs))) +
      scale_color_gradient(low="blue", high="red") +
      scale_shape_manual(values=c(15:18)) +
      ggtitle("tsne dim reduced") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
dev.off()

for (p in 1:length(tsne.iter)) {
  print(
    ggplot(tsne.frame.iter[[p]], aes(x, y, color = group, shape = stage)) + 
      geom_point(size = 5) +
      ggtitle(paste0("tsne ","perplexity = ", p, ", cost =", sum(tsne.iter[[p]]$costs))) +
      scale_color_brewer(palette="Paired") +
      scale_shape_manual(values=c(15:18)) +
      theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))
  )
}

for (p in 1:length(tsne.iter)) {
  print(
    ggplot(tsne.frame.iter[[p]], aes(x, y, color = size, shape = stage)) + 
      geom_point(size = 5) +
      ggtitle(paste0("tsne ","perplexity = ", p, ", cost =", sum(tsne.iter[[p]]$costs))) +
      scale_color_gradient(low="blue", high="red") +
      scale_shape_manual(values=c(15:18)) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20))
  )
}

