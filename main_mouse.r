library(CellMix)
library(GEOquery)
library(R.matlab) # install.packages("R.matlab")

########## load data matrix from mat file #############
inputFile <- 'deconv_for_zapala.mat'
outputFileNMF <- 'cellmix_output/cellmix_zapala_nmf.mat'
outputFileDeconf <- 'cellmix_output/cellmix_zapala_Deconf.mat'
# inputFile <- 'deconv_for_akahoshi.mat'
# outputFileNMF <- 'cellmix_output/cellmix_akahoshi_nmf.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_akahoshi_Deconf.mat'
# inputFile <- 'deconv_for_allen_cortex.mat'
# outputFileNMF <- 'cellmix_output/cellmix_allen_cortex_nmf.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_allen_cortex_Deconf.mat'

markersFile <- "brain_markers_mouse.txt"
# markersFile <- "brain_markers_mouse_allen.txt"
do_cellmix(inputFile, outputFileNMF, markersFile ,"NMF")
do_cellmix(inputFile, outputFileDeconf, markersFile ,"DECONF")




