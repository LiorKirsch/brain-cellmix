library(CellMix)
library(GEOquery)
library(R.matlab)
setwd('~/workspace/cell_mix/')

########## load kang data matrix from mat file #############
inputFile <- 'kang_samples_adults.mat'
outputFileNMF <- 'cellmix_output/cellmix_kang_nmf.mat'
outputFileDeconf <- 'cellmix_output/cellmix_kang_Deconf.mat'

# inputFile <- 'kang_cortex_samples_adults.mat'
# outputFileNMF <- 'cellmix_output/cellmix_kang_cortex_nmf.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_kang_cortex_Deconf.mat'

markersFile <- "brain_markers_human.txt"
do_cellmix(inputFile, outputFileNMF, markersFile ,"NMF")
do_cellmix(inputFile, outputFileDeconf, markersFile ,"DECONF")

