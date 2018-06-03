##for all data analysis
library(ggplot2)
library(forcats)

for(i in c("oxygen|hypoxia|oxydative","MHC|antigen|present","neutroph","extracellular|collagen",
           "folding|incorrect|unfolded","p53|apoptosis","cytokine|interferon","actin|Wnt|polarity",
           "stress|apoptosis","catabolic", "metabolic|lipoprotein|lipid|cholesterol|sterol")){
  print(i)
  test = GO_res$DE_DCIS_IDC_up@result[grepl(i,GO_res$DE_DCIS_IDC_up@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_size_continuous(range = c(4, 6)) +
    theme_bw() + theme(text = element_text(size=14),
                       axis.text.y = element_text(size=20)) +
    scale_colour_gradient(limits=c(0, 0.10), low="red") +
    ylab(NULL) + ggtitle(i))
}

for(i in c("neuron")){
  print(i)
  test = GO_res$DE_EN_DCIS_down@result[grepl(i,GO_res$DE_EN_DCIS_down@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_size_continuous(range = c(4, 6)) +
          theme_bw(base_size = 14) + theme(text = element_text(size=15)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL) + ggtitle(i))
}
dev.off()

####for reanalysis of 20190219 on March 5th

library(ggplot2)
library(forcats)

for(i in c("oxygen|hypoxia|oxydative","MHC|antigen|present","neutroph","extracellular|collagen",
           "folding|incorrect|unfolded","p53|apoptosis","cytokine|interferon","actin|Wnt|polarity",
           "stress|apoptosis","catabolic", "metabolic|lipoprotein|lipid|cholesterol|sterol")){
  print(i)
  test = GO_res$DE_DCIS_IDC_up@result[grepl(i,GO_res$DE_DCIS_IDC_up@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_with_EN_batch_corr/plots/dotplots/"
  png(paste0(path,i,".png"), width = 800, height = 800)
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_size_continuous(range = c(4, 6)) +
          theme_bw(base_size = 14) + theme(text = element_text(size=15)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL) + ggtitle(i))
  dev.off()
}

###for 8 batches no batch correction


library(ggplot2)
library(forcats)

for(i in c("extracellular|movement|location|motility|migration","MHC|antigen|present","collagen",
           "folding|unfolded","virus|interferon","axon|neuron|projection",
           "catabolic", "autophag","angiogen","ossification|skeletal","transform","ephrin")){
  print(i)
  test = GO_res$DE_DCIS_IDC_up@result[grepl(i,GO_res$DE_DCIS_IDC_up@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  #path = "~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/dotplots/GO/DCIS_IDC/"
  #png(paste0(path,i,".png"), width = 800, height = 800)
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_size_continuous(range = c(4, 6)) +
          theme_bw(base_size = 14) + theme(text = element_text(size=15)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL) + ggtitle(i))
  #dev.off()
}


###4batches  19 May 2018


library(ggplot2)
library(forcats)

for(i in c("MHC|antigen|present","tumor necrosis","oxygen|hypoxia","neutroph")){
  print(i)
  test = GO_res$DE_DCIS_IDC_up@result[grepl(i,GO_res$DE_DCIS_IDC_up@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  path = "~/Desktop/s7s_Breast/data/4stages/NOT_BC_CBx_group/plots/GO_enrichment/selected_dotplots/"
  png(paste0(path,i,".png"), width = 800, height = 800)
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_size_continuous(range = c(4, 6)) +
          theme_bw(base_size = 14) + theme(text = element_text(size=25)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL) + ggtitle(i))
  dev.off()
}


#2016_2017 batches magda megaclasses

library(ggplot2)
library(forcats)

for(i in c("MHC|antigen|present","neutroph","T cell","innate","defense","oxygen|hypoxia",
           "NF-kappa","necrosis factor|necrotic","virus")){
  print(i)
  test = GO_res$DE_DCIS_IDC_up@result[grepl(i,GO_res$DE_DCIS_IDC_up@result$Description),]
  #setSize = as.numeric(strsplit(test$GeneRatio,"/")[[1]][2])
  nodeSize = as.numeric(sapply(strsplit(test$BgRatio,"/"), function(x) x[1]))
  
  #test$GeneRatio = sapply(test$Count, function(x) x/setSize)
  test = data.frame(test, annotRatio = test$Count/nodeSize)
  path = "~/Desktop/s7s_Breast/data/2016_2017_NOT_batch_corrected/Mac_plots/selected_dotplots/"
  png(paste0(path,i,".png"), width = 1200, height = 600)
  print(ggplot(test, aes(x = annotRatio, y = fct_reorder(Description, annotRatio))) + 
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_size_continuous(range = c(4, 6)) +
          theme_bw(base_size = 14) + theme(text = element_text(size=25)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL))
  dev.off()
}
