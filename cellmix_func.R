library(CellMix)
library(GEOquery)
library(R.matlab)


do_cellmix <- function(input_MAT_filename,output_MAT_filename,marker_file_name,cellmix_method) {
  
  ########## load data matrix from mat file #############
  data <- readMat(input_MAT_filename)
  
  expression_matrix = data[["expression"]]
  gene_names = data[["gene.symbols"]]
  gene_names <- unlist(gene_names,recursive=FALSE)
  rownames(expression_matrix) <- gene_names
  
  
  ######### load cell-type markers #####################
  brainmarkers <- MarkerList(file = marker_file_name, header = TRUE)
  summary(brainmarkers)
  
  # check that the markers are matched
  
  geneIds(brainmarkers)$Neurons %in% gene_names
  geneIds(brainmarkers)$Astrocytes %in% gene_names
  geneIds(brainmarkers)$Oligodendrocytes %in% gene_names
  
  print('applying cellmix')
  if (cellmix_method == "NMF") {
    ################# Use NMF with markers ###############
    decnvoloved_data <-ged(expression_matrix, brainmarkers,"ssKL", rng = 123456, nrun = 20, log = FALSE) # ,verbose = TRUE)
  }
  else if (cellmix_method == "NMFforb") {
    decnvoloved_data <-ged(expression_matrix, brainmarkers,"ssFrobenius", rng = 123456, nrun = 20)
  }
  else if (cellmix_method == "meanProfile") {
    decnvoloved_data <-ged(expression_matrix, brainmarkers,"meanProfile", rng = 123456, nrun = 20)
  }
  else if (cellmix_method == "DSA") {
    decnvoloved_data <-ged(expression_matrix, brainmarkers,"DSA", rng = 123456, nrun = 20)
  }
  else if (cellmix_method == "DECONF") {
    ########## use DECONF ############
    decnvoloved_data <-ged(expression_matrix, brainmarkers,"deconf", rng = 123456, nrun = 20)
  }
  else {
    print('unknown cellmix option')
  }
  
  
  celltype_profile = decnvoloved_data@fit@W
  proportions = decnvoloved_data@fit@H
  cell_types = brainmarkers@names
  print('saving MAT file')
  writeMat(output_MAT_filename, celltype_profile=celltype_profile, proportions=proportions, cell_types=cell_types, gene_names=gene_names)
  
}
