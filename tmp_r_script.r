setwd("/home/lior/workspace/brain-cellmix")
source("cellmix_func.R")
do_cellmix("deconv input/deconv_for_zapala.mat","markers subsets/cellmix_zapala_nmf.mat12","markers subsets/brain_markers_mouse.txt12","NMF") 
