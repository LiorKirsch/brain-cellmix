library(CellMix)
library(GEOquery)
library(R.matlab) # install.packages("R.matlab")
source("cellmix_func.R")


########## load data matrix from mat file #############
# inputFile <- 'deconv_for_bernard.mat'
# outputFileNMF <- 'cellmix_output/cellmix_bernard_nmf.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_bernard_Deconf.mat'

inputFile <- 'deconv_blueprint_monkey_micro.mat'
outputFileNMF <- 'cellmix_output/cellmix_blueprint_monkey_micro.mat'
outputFileDeconf <- 'cellmix_output/cellmix_blueprint_monkey_micro_deconf.mat'

# inputFile <- 'deconv_blueprint_monkey_macro.mat'
# outputFileNMF <- 'cellmix_output/cellmix_blueprint_monkey_macro.mat'
# outputFileDeconf <- 'cellmix_output/cellmix_blueprint_monkey_macro_deconf.mat'

markersFile <- "brain_markers_monkey.txt"
do_cellmix(inputFile, outputFileNMF, markersFile ,"NMF")
do_cellmix(inputFile, outputFileDeconf, markersFile ,"DECONF")

