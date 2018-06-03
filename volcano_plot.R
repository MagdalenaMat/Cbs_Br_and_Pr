TAMact = c("TREM2", "CHIT1", "TGFB1","IL4R","TGFBI","STAT6","NFKBIA") #"TLR5",
M2pol = c("TGFB1","IL4R")
angio = c("ANXA2","FN1","MMP14","MMP2","VEGFB","NRP2")
phago_ant = c("CTSA","CTSO","CTSC","CTSK","CTLB","TAP2","HLA.C","B2M") #"CTSD",
stress = c("P4HB","BAX","UBB","UBC","ENO1","GAPDH")

processes = list(TAMact,M2pol,angio,phago_ant,stress)
names(processes) = c("TAMact","M2pol","angio","phago_ant","stress")

setwd("~/Desktop/s7s_Breast/data/4stages/NOT_BC_CBx_group/plots/volcano_selected/")
for(i in names(to_GSEA)[6]){
  plotd = data.frame(LogFC = to_GSEA[[i]][, 2], LZqv = to_GSEA[[i]][, 4], gnames = to_GSEA[[i]][[1]])
  for(selected in names(processes)) {
    m = plotd[plotd$gnames %in% processes[[selected]],]
    a = list(x = m$LogFC, y = m$LZqv, text = m$gnames,font = list(size = 18)) #text = m$gnames,
    p = plot_ly(data = plotd, x = ~LogFC, y = ~LZqv, text = ~ gnames, marker = list(size = 20) ) %>%
      layout(title = i, xaxis = list(title = paste0("log2(FC_",i,")")), yaxis = list(title = "-log10(Zqv)"),
             annotations = a)
    print(p)
  export(p, paste0(selected,"_",i,"_25FDR_vulcano.png"))
  }
}
dev.off()


#######
TAMact = c("PYCARD","CASP1","IFI35","STAT1","IFIT5") #"EIF2AK2"
endocytosis_antigen_presentation = c("CLTB", "CTSK","AP1S1","AP1S2","ARF1","TAP2","B2M","HLA.B","HLA.C","RHOA") 
M2pol = c("TGFB1","TGFBI","HCK","MMP14","VEGFB","TIMP2","ZC3H12A") #,"VEGFC","TIMP1",
phagocytosis = c("MARCO","SPP1","TREM2","CXCR4")
TNF = c("TNFRSF1A","TRAF5","RIPK1","RBCK1","TNFAIP3") 
ER_stress = c("P4HB","UBB","UBC","CTSD","PSMB4","PSMD11","PSMD2")

processes = list(TAMact,endocytosis_antigen_presentation, M2pol, phagocytosis, TNF, ER_stress)
names(processes) = c("TAM_activation","Endocytosis_antigen_presentation","M2pol","Foam_cell_phenotype","TNF_signaling", "ER_stress")


for(i in names(to_GSEA)[8]){
  for(selected in names(processes)){
    print(selected)
    plotd = data.frame(to_GSEA[[i]])
    plotd$Selected <- ifelse(plotd$names %in% processes[[selected]],"T","F")
    path = "~/Desktop/toll_plots/vulcano/DCIS_IDC/"
    png(paste0(path,i,selected,".png"), width = 800, height = 800)
    plot(ggplot(plotd, aes(x = FC, y = ml10pv)) +
      geom_point(aes(color = Selected)) +
      scale_color_manual(values = c("gray", "red")) +
      theme_bw(base_size = 20) + theme(legend.position = "bottom") +
      geom_text_repel(
        data = subset(plotd, plotd$names %in% processes[[selected]]),
        aes(label = names),
        size = 10,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")) +
      ggtitle(paste0(i,"_", selected)) + labs(x = "LFC", y = "-log10 qvalue") +
      theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  }
}
