
to_GO = list()
for(nm in names(to_enr)){
  to_GO[[paste0("DE_",nm,"_up")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] > 0), "names"]
  to_GO[[paste0("DE_",nm,"_down")]] = to_enr[[nm]][which(to_enr[[nm]][,"FC"] < 0), "names"]
}

lengths(to_GO)

library("Rgraphviz")
library("DOSE")
library("clusterProfiler")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

#translate gene symbols to entrez ID

for(element in names(to_GO)){
  to_GO[[element]] = sub("[.]","-", to_GO[[element]]) 
}

to_GO_EI = list()
for(element in names(to_GO)){
  to_GO_EI[[element]] = bitr(to_GO[[element]], 
                             fromType = "SYMBOL",
                             toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)
}

###run GO analysis, I am running ont= BP (biological processes) /MF/CC
GO_res = list()
for(element in names(to_GO_EI)){
  GO_res[[element]] = enrichGO(gene         = to_GO_EI[[element]]$ENTREZID,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable = TRUE)
}


#path = setwd("~/Desktop/s7s_Breast/data_after_decon/8_batches_not_batch_corr/plots/dotplots/GO/")
for(element in names(GO_res)){
  #png(paste0(path,"/",element,".png"), width = 900, height = 1100)
  plot(dotplot(GO_res[[element]],showCategory=50)+labs(title= paste0("Macrophage_",element)))
  #dev.off()
}

dev.off()


#ranked list of genes by FC
FC_rnk = list()
for(i in names(to_GSEA)){
  FC_rnk[[i]] = to_GSEA[[i]][order(to_GSEA[[i]][,"FC"],decreasing = T),]
}


vec_GSEA = list()
for(i in names(FC_rnk)){
  vec_GSEA[[i]] = FC_rnk[[i]][,"FC"]
  names(vec_GSEA[[i]]) = FC_rnk[[i]][,"names"]
}

####GSEA for GO terms check BP MF
resGSEA = list()
for(element in names(vec_GSEA)){
  resGSEA[[element]] = gseGO(geneList     = vec_GSEA[[element]],
                             OrgDb        = org.Hs.eg.db,
                             ont          = "BP",
                             keyType = "SYMBOL",
                             nPerm        = 1000,
                             minGSSize    = 20,
                             maxGSSize    = 500,
                             pvalueCutoff = 0.1,
                             verbose      = FALSE)
}

library(ggplot2)
for(element in names(resGSEA)){
  if(nrow(resGSEA[[element]]@result) > 0){
    plot(dotplot(resGSEA[[element]],showCategory=50) + labs(title= paste0("GO_BP_GSEA_",element)))
  }
}


########draw GSEA plot
gseaplot(resGSEA$DCIS_IDC, geneSetID = "GO:0001568",title = "blood vessel development")
gseaplot(resGSEA$DCIS_IDC, geneSetID = "GO:0002683",title = "negative regulation of immune system process")

#######MSigDB c7 immunological processes enrichment
library(GSEABase)
library("clusterProfiler")
c7 <- read.gmt("~/Desktop/s7s_Breast/gene_sets/c7.all.v6.1.entrez.gmt")

c7enr_res = list()
for(element in names(to_GO_EI)){
  c7enr_res[[element]] = enricher(gene         = to_GO_EI[[element]]$ENTREZID,
                                  TERM2GENE=c7)
}

for(element in  names(c7enr_res)){
  plot(dotplot(c7enr_res[[element]],showCategory=50)+labs(title= paste0("MSigDB_c7_enrich_",element)))
}


##GSEA

for(element in names(FC_rnk)){
  FC_rnk[[element]][,"entrez"] = mapIds(org.Hs.eg.db,
                                        keys= as.character(FC_rnk[[element]][,"names"]),
                                        column=c("ENTREZID"),
                                        keytype = "SYMBOL",
                                        multiVals = "first")
}

vec_GSEA_ez = list()
for(i in names(FC_rnk)){
  vec_GSEA_ez[[i]] = FC_rnk[[i]][,"FC"]
  names(vec_GSEA_ez[[i]]) = FC_rnk[[i]][,"entrez"]
}


c7_GSEA = list()
for(element in  names(vec_GSEA_ez)){
  c7_GSEA[[element]] = GSEA(vec_GSEA_ez[[element]],
                            TERM2GENE=c7, 
                            verbose=FALSE,
                            nPerm        = 1000,
                            minGSSize    = 10,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.1)
}

library(ggplot2)
for(element in names(c7_GSEA)){
  if(nrow(c7_GSEA[[element]]@result) > 0){
    plot(dotplot(c7_GSEA[[element]],showCategory=50) + labs(title= paste0("MSigDB_c7_GSEA_",element)))
  }
}

###########c2.cp.v6.1.entrez.gmt c2 canonical pthways
c2 <- read.gmt("~/Desktop/s7s_Breast/gene_sets/c2.cp.v6.1.entrez.gmt")

c2enr_res = list()
for(element in names(to_GO_EI)[11]){
  c7enr_res[[element]] = enricher(gene         = to_GO_EI[[element]]$ENTREZID,
                                  TERM2GENE=c2)
}

for(element in names(c2enr_res)){
  plot(dotplot(c2enr_res[[element]],showCategory=50))
}


##GSEA


c2_GSEA = list()
for(element in  names(vec_GSEA_ez)){
  c2_GSEA[[element]] = GSEA(vec_GSEA_ez[[element]],
                            TERM2GENE=c2, 
                            verbose=FALSE,
                            nPerm        = 1000,
                            minGSSize    = 10,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.1)
}

library(ggplot2)
for(element in names(c2_GSEA)){
  if(nrow(c2_GSEA[[element]]@result) > 0){
    plot(dotplot(c2_GSEA[[element]],showCategory=50) + labs(title= paste0("MSigDB_c2_GSEA_",element)))
  }
}
