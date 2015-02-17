setwd("/home/lab/lior/Projects/cellmix")
source("cellmix_func.R")
do_cellmix("deconv input/deconv_for_zapala.mat","cellmix results/cellmix_zapala_nmf.mat","markers subsets/brain_markers_mouse.txt2","NMF") 
