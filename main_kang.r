library(CellMix)
library(GEOquery)
library(R.matlab) # install.packages("R.matlab")

######## https://support.bioconductor.org/p/58585/ ######

#  load kang data (normally requires an internet connection to GEO)
# kangdata <- getGEO('GSE25219', destdir="/home/lior/R-data/CellMix/")


########## load kang data matrix from mat file #############
kangdata <- readMat('/cortex/data/microarray/human/Kang2011/kang_samples_with_ontology.mat')
expression_matrix = kangdata[["data"]]
gene_names = kangdata[["gene.names"]]
gene_names <- unlist(gene_names,recursive=FALSE)
rownames(expression_matrix) <- gene_names


######### load cell-type markers #####################
brainmarkers <- MarkerList(file = "brain_markers.txt", header = TRUE)
summary(brainmarkers)

# check that the markers are matched

geneIds(brainmarkers)$Neuron %in% gene_names
geneIds(brainmarkers)$Astrocytes %in% gene_names
geneIds(brainmarkers)$Oligodendrocytes %in% gene_names


################# Use NMF with markers ###############


# deconvolve using KL divergence metric
kl <-ged(expression_matrix, brainmarkers,"ssKL", rng = 12345, nrun = 20)



# plot against known proportions
profplot(kangdata, kl)
# check consistency of most expressing cell types in signatures
g <-MarkerList(predict(x,"features"), names =names(m))
basismarkermap(g, kl)
# correlation with known signaturesbasiscor(x, kl)


