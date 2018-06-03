##plot fractions
fractions = read.table("~/Desktop/s7s_Breast/data/2016_2017_exp_batch_corrected/2016_2017/CIBERSORTx_to_hires_bc_2016_2017_data.txt_Adjusted.txt",
                       header = T, row.names = 1, sep = "\t")

fractions = fractions[,1:22]
pdata = read.table("~/Desktop/s7s_Breast/data/2016_2017_exp_batch_corrected/2016_2017_pdata.txt",
                   header = T, sep = "\t")

#------
fractions = read.table("~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/fractions_not_BC/CIBERSORTx_Job2_Results.csv",
                       header = T, row.names = 1, sep = ",")

fractions = fractions[,1:22]
pdata = read.table("~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/fractions_not_BC/2016_2017_pdata.txt",
                   header = T, sep = "\t")
pdata = pdata[1:176,]

all(pdata$name == rownames(fractions))

pdata$subtype =factor(pdata$subtype, levels = c("normal","EN","DCIS","IDC"))
type = factor(pdata$subtype)

to_plot = as.matrix(fractions)

library(pheatmap)
library(RColorBrewer)

col1 <- brewer.pal(8, "Set2")

mydf <- data.frame(row.names = rownames(to_plot), category = type)

# add row annotations
#myBreaks = unique(c(seq(0, 0.19, length=49), 0.2, seq(0.21,max(to_plot), length=50)))
to_boxplot = to_plot
to_plot[to_plot >= 0.5] = 0.5

png(paste0("~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/fractions_not_BC/","BCx_franctions_hm_0_5.png"), 
    width = 1000, height = 1200)
pheatmap(to_plot, cluster_cols = F, cluster_rows = F, 
         annotation_row = mydf,gaps_row = c(37,64,133,176),
         show_rownames = F)
dev.off()
max(to_plot)

#####boxplots
library(ggpubr)
to_ggplots = data.frame(to_boxplot, type)
to_ggplots = cbind(Macrophages = apply(to_ggplots[,c(14:16)],1,sum), to_ggplots)
colnames(to_ggplots)[10] = "Treg"
colnames(to_ggplots)[7] = "Trm"

my_comparisons <- list( c("normal","EN"),c("EN", "DCIS"),c("DCIS","IDC"), 
                        c("normal", "DCIS"), c("EN","IDC"), c("normal", "IDC"))
path = "~/Desktop/toll_plots/fractions/"
#s7s_Breast/data/2016_2017_NOT_batch_corrected/fractions_not_BC/
for( i in colnames(to_ggplots)[1:23]){
  tiff(paste0(path,i,"600.tiff"))
  plot(ggboxplot(to_ggplots, x = "type", y = i,  outlier.colour = NA,
                 color = "type", palette = c("dodgerblue3","orange","chartreuse4","red"), add = "jitter", size = 1) + 
         theme(text = element_text(size=30)) + rotate_x_text() +
         stat_compare_means(comparisons = my_comparisons) + labs(x = ""))
  dev.off()
}

