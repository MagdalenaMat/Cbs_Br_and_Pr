ex = read.table("~/Desktop/s7s_Breast/to_fractions/all_data_to_fractions.txt", sep = "\t",  header = T, row.names = 1)
pdata = read.table("~/Desktop/s7s_Breast/to_fractions/pdata_all_samples.txt", sep = "\t",  header = T, row.names = 1)


Tcell = c("PDCD1","CD274","IL10","CCL5","CCR5","CXCL14","IL7R","CTLA4","TIGIT","FOXP3","IL2RA","CD39","CD2","TNFRSF9","TNFRSF4", "TNFRSF9","TNFRS18",
          "CTLA4","LAG3","ITGAE","CD3G","PFKP","FOXP3","GATA3","B3GNT2","CD2","CD8A","CD3E","CD151","FASLG","IFNG","TBX21","BATF3","PRF1",
          "SOX4","RAG1","BATF3","CXCL14","GZMA","GZMB","GZMK","ICOS","CD69","IL19","IL18RAP","CD8A","CD8B",
          "CD151","CD9","CD32","ZAP70", "CCL4","CCL5","CXCR4","CCR5")
MAC = c("TREM2", "CHIT1", "TLR5","TGFB1","IL4R","STAT6","NFKBIA","TGFBI","ANXA2","FN1","MMP14","MMP2","VEGFB","NRP2",
        "MMP13","TIMP1","TIMP2","MMP14","MMP11","B2M","HLA.A","HLA.B","HLA.C","HLA.F","IFI35","IFIT3","IRF7","STAT1",
        "CLTB","CTSD", "CTSA","CTSC","CTSD","CTSK","CTSO","TAP2","NOS3","GAPDH","ENO1","ATM","PML","BAX",
        "P4HB","UBB","UBC","CTSD","PSMB4","PSMD11","PSMD14","PSMD2","PSME1")
to_boxplot = t(ex[rownames(ex) %in% Tcell,])
all(colnames(ex) == rownames(pdata))
pdata$subtypes = factor(pdata$subtypes, levels = c("normal","EN","DCIS","IDC"))
stages = pdata$subtypes
to_ggplots = data.frame(to_boxplot, stages)

my_comparisons <- list( c("normal","EN"),c("EN", "DCIS"),c("DCIS","IDC"), 
                        c("normal", "DCIS"), c("EN","IDC"), c("normal", "IDC") )
path = "~/Desktop/s7s_Breast/bulk_plots/T_cells/"
for( i in colnames(to_ggplots)[1:42]){
  png(paste0(path,i,".png"), width = 1200, height = 1200)
  plot(ggboxplot(to_ggplots, x = "stages", y = i,  outlier.colour = NA,
                 color = "stages", palette = "jco", add = "jitter") + 
         theme(text = element_text(size=20)) +
         stat_compare_means(comparisons = my_comparisons))
  dev.off()
}
