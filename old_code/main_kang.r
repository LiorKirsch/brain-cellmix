library(CellMix)
library(GEOquery)
library(R.matlab)
setwd('~/workspace/cell_mix/')

########## load kang data matrix from mat file #############
kangdata <- readMat('kang_samples_adults.mat')
outputFileNMF <- 'cellmix_output/cellmix_kang_nmf.mat'
outputFileDeconf <- 'cellmix_output/cellmix_kang_Deconf.mat'

# kangdata <- readMat('kang_cortex_samples_adults.mat')
# outputFileNMF <- 'cellmix_output/cellmix_kang_cortex_nmf.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_kang_cortex_Deconf.mat'

expression_matrix = kangdata[["data"]]
gene_names = kangdata[["gene.names"]]
gene_names <- unlist(gene_names,recursive=FALSE)
rownames(expression_matrix) <- gene_names


######### load cell-type markers #####################
brainmarkers <- MarkerList(file = "brain_markers_human.txt", header = TRUE)
summary(brainmarkers)

# check that the markers are matched

geneIds(brainmarkers)$Neuron %in% gene_names
geneIds(brainmarkers)$Astrocytes %in% gene_names
geneIds(brainmarkers)$Oligodendrocytes %in% gene_names

############ checks on subset ############
# samples_subset_size = 500
# gene_subset_size = 100
# all_markers_symbols = c( geneIds(brainmarkers)$Neurons, geneIds(brainmarkers)$Astrocytes, geneIds(brainmarkers)$Oligodendrocytes)
# is_member = gene_names %in% all_markers_symbols
# 
# is_member[sample(1:length(gene_names), gene_subset_size)] = TRUE
# num_samples = dim(expression_matrix)[2]
# sample_subset = sample(1:num_samples, samples_subset_size) 
# matrix_subset = expression_matrix[is_member, sample_subset]
# kl_subset <-ged(matrix_subset, brainmarkers,"ssKL", rng = 12345, nrun = 20)
# 
# celltype_profile_test = kl_subset@fit@W
# proportions_test = kl_subset@fit@H


################# Use NMF with markers ###############

# deconvolve using KL divergence metric
kl <-ged(expression_matrix, brainmarkers,"ssKL", rng = 12345, nrun = 20)
celltype_profile = kl@fit@W
proportions = kl@fit@H
cell_types = brainmarkers@names
writeMat(outputFileNMF, celltype_profile=celltype_profile, proportions=proportions, cell_types=cell_types, gene_names=gene_names)



########## use DECONF ############
#  deconvolve using KL divergence metric
dec <-ged(expression_matrix, brainmarkers,"deconf", rng = 12345, nrun = 20)

celltype_profile = dec@fit@W
proportions = dec@fit@H
cell_types = brainmarkers@names
writeMat(outputFileDeconf, celltype_profile=celltype_profile, proportions=proportions, cell_types=cell_types, gene_names=gene_names)

hist(proportions["Astrocytes",])
hist(proportions["Neuron",])
hist(proportions["Oligodendrocytes",])