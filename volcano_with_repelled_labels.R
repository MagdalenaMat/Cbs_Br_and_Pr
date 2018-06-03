library(ggplot2)
library(ggrepel)

for(i in names(to_GSEA)){
  to_GSEA[[i]]$Significant <- ifelse(to_GSEA[[i]]$Zqv < 0.05, "FDR < 0.05", "Not Sig")
  path = "~/Desktop/s7s_Breast/data"
  png(paste0(path,element,"_scalled_no_dendro.png"), width = 1000, height = 1200)
  plot(ggplot(data.frame(to_GSEA[[i]]), aes(x = FC, y = ml10pv)) +
         geom_point(aes(color = Significant)) +
         scale_color_manual(values = c("red", "grey")) +
         theme_bw(base_size = 12) + theme(legend.position = "bottom") +
         geom_text_repel(
           data = subset(data.frame(to_GSEA[[i]]), Zqv < 0.05),
           aes(label = names),
           size = 5,
           box.padding = unit(0.35, "lines"),
           point.padding = unit(0.3, "lines")) +
         ggtitle(i)
  )
  dev.off()
}