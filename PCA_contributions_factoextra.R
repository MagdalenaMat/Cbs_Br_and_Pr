library("FactoMineR")
library("devtools")
library("factoextra")

edata <- breast.data
design.file = desine2

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = breast.data,
                              colData = desine2,
                              design= ~ subtype)
vsd <- vst(dds, blind=FALSE)

#subset 3000 most variable genes
rv = rowVars(assay(vsd)) 
select = order(rv, decreasing=TRUE)[1:3000]
vsd = assay(vsd)[select,]


centered = vsd - rowMeans(vsd)
dim(centered)

# We transpose in order to run PCA
t_centered = t(centered)
dim(t_centered)

# The categorical variables must be merged with the data table.
t_centered = cbind.data.frame(desine2$subtype, desine2$batch, t_centered)


# tell PCA the first two columns are categorical variables
res.pca <- PCA(t_centered, quali.sup = 1:2 , graph = FALSE)
fviz_pca_ind(res.pca, habillage = 1, axes = c(3,4))
plot(res.pca,habillage=1,axes=c(1,3),cex=0.7) 


# Make the scree plot using the package factoextra:
fviz_screeplot(res.pca, ncp=10)

# We can visualize the most important variables associated with a given PC, you can decide to show only the top contributing variables.
fviz_pca_contrib(res.pca, choice = "var", axes = 1, top = 10)

# We can  carry out the PCA and build the graph of individuals by colouring the individuals according to the variable case where habillage defines variable to colour (here the ???rst variable in the table). The size of the text can be changed using the setting cex ("cex = 0.7" instead of 1 by default):
plot(res.pca,habillage=1,cex=0.7)
plot(res.pca,habillage=2,cex=0.7)
plot(res.pca,habillage=2,axes=c(3,4),cex=0.7)
plot(res.pca,habillage=2,axes=c(4,5),cex=0.7)

# It is also possible to visualise other components of the PCA
plot(res.pca,habillage=1,axes=c(1,3),cex=0.7) 
plot(res.pca,habillage=2,axes=c(1,3),cex=0.7)


# Extract genes that most linked to the dimension (proba = 0.05 by default) by dimdesc function:
gene_list = dimdesc(res.pca)
gene_list_PC1 = gene_list$Dim.1$quanti
gene_list_PC2 = gene_list$Dim.2$quanti
gene_list_PC3 = gene_list$Dim.3$quanti
