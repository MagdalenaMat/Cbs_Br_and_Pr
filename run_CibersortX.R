setwd("~/Desktop/Chloe2/")

# GROUP MODE

sigfile = "LM22.txt"
truth = "Mono_corect.txt"
#mergedclasses = "GEPs_LM9.txt"
sourceGEP = "LM22_source_GEPs.txt"

mergedclasses = "9_mergedclasses_newformat.txt"

# m = read.table(mixturefile, header=T, sep = "\t", check.names = F)
# s = read.table(sigfile, header=T, sep = "\t", check.names = F)
#t = read.table(truth, header=T, sep = "\t", check.names = F)
#mc = read.table(mergedclasses, header=T, sep = "\t", check.names = F)
# sg = read.table(sourceGEP, header=T, sep = "\t", check.names = F)

source("./CIBERSORTxGEP.R")
for(i in Sys.glob("6b*")){
  mixturefile = i
  CIBERSORTxGEP(mixture = mixturefile, 
                sigmatrix = sigfile,
                classes = mergedclasses,
                groundtruth = truth,
                label = i,
                rmbatch=T,
                sourceGEPs=sourceGEP,
                redocibersort = TRUE,
                redo = T, # You probably don't need to run CIBERSORT everytime for your analysis, so you can set that to FALSE
                QN=F, 
                plots=F,
                outdir="not_batch_corr") # Since you are not plotting, you should not get a warning here
}


gep_macrophages = read.table("CIBERSORTx_output/CIBERSORTxGEP_1_GEPs_Filtered.txt", header = T, sep = "\t")

# CIBERSORTxGEP(mixture, sigmatrix, classes, cibresults = NA, label= i, groundtruth= groundtruth, threads = 1, nsampling = 30,
#              maxsamples= "NA" , degclasses = "", rmbatch = TRUE,  sourceGEPs = "LM22_source_GEPs.txt", 
#              redocibersort=TRUE,plots=FALSE)

mixture=mixturefile 
sigmatrix=sigfile
classes = mergedclasses
groundtruth="Monocytes.txt"
rmbatch=T
sourceGEPs=sourceGEP
QN=T
cibresults = "NA"
label = ""
maxsamples="NA"
degclasses = "NA"
redo = TRUE
threads = "NA"
plots = TRUE
nsampling = 30
outdir=outdirglobal
redocibersort = TRUE


predicted = paste(outdir,"/CIBERSORTxGEP_",label,"GEPs.txt",sep="")
dolog= TRUE
CVs=paste(outdir,"/CIBERSORTxGEP_",label,"GEPs_CVs.txt",sep="")
CVthresh = -1
label=origlabel
cibresults=paste(outdir,"/CIBERSORTxGEP_",label,"Weights.txt",sep="")
outdir=outdir
subsetgenes="NA"
maxnum="All"

  