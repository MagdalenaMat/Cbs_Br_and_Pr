path = "~/Desktop/after_decon/full/full_not_corr/"
setwd(path)
files = Sys.glob("*Monocytes.txt")

results = list()
for(i in files){
  results[[i]] = as.character(read.table(i, sep = "\t", header = T)[,1])
}

names(results) = c("DCIS","EN","IDC","normal")

universe = Reduce(intersect, results)

to_heat = data.frame(genes = universe)
for(i in files){
  dat = read.table("DCIS_Monocytes.txt", sep = "\t", header = T)
  dat = dat[dat[,1] %in% universe,]
  to_heat = 
}
