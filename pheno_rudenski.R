#rudenski_immune_phenotypes
phenotypes = read.csv("~/Desktop/BC_ssRNAseq_rudenski/pheno_magda.csv")
phenotypes[phenotypes == ""] = NA 
phenotypes = as.list(phenotypes, na.rm = T)
phenotypes = lapply(phenotypes, na.omit)
phenotypes = lapply(phenotypes, toupper)

phenotypes = phenotypes[c(1:5,12:18)]
names(phenotypes) =c("Treg", "CD8_act", "anti_infl","anergy","pro_infl", "M1","M2","Teff_cytotoxic","IFN_I","IFN_II","hypoxia","Tcell_term_diff")
pheno = as.character()
for(name in names(phenotypes)){
  ph = rep(name, length(phenotypes[[name]]))
  pheno = c(pheno, ph)
}
all_phen = data.frame(genes = unlist(phenotypes), pheno = as.factor(pheno))


data = read.csv("~/Desktop/s7s_Breast/s7s/170624_NS500418_0623_AH2F27BGX3.reads.csv", header = T)
library(dplyr)
data = mutate(gen = sub("\\..*", "",data$gene.name.nodup), data) #removes all the numbers after the dot
data = data[,c(2,40)]
genes = data[,1, drop = F]

to_HR = all_phen$genes[!duplicated(all_phen$genes)]
to_HR[is.na(pmatch(to_HR, data$gen))]

to_HR = data.frame(genes = as.character(to_HR))
setwd("~/Desktop/BC_ssRNAseq_rudenski")
write.table(to_HR, "small_rudenski_pheno.txt", row.names = F, col.names = F, quote = F)


data_1 = read.table("~/Desktop/s7s_Breast/data_for_deconvolution/bulkDCIS.txt", header = T, sep = "\t")
to_HR$genes[is.na(pmatch(to_HR$genes, data_1$genes))]

to_HR$genes[is.na(pmatch(to_HR$genes, rownames(TPM)))]
