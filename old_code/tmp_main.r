library(CellMix)


#  load data (normally requires an internet connection to GEO)
acr <-ExpressionMix("GSE20300", verbose = 2)

# estimate proportions using signatures from Abbas et al. (2009)
res <-gedBlood(acr, verbose = TRUE)


# aggregate into CBC
cbc <-asCBC(res)dim(cbc)
# plot against actual CBC
profplot(acr, cbc)
# plot cell proportion differences between groups
boxplotBy(res, acr$Status, main ="Cell proportions vs Transplant status")





################# Use NMF with markers ###############
                    
# generate random data with 5 markers per cell type
x <-rmix(3, 20, 200, markers = 5) 
m <-getMarkers(x)
# deconvolve using KL divergence metric
kl <-ged(x, m,"ssKL", rng = 12345, nrun = 20)



# plot against known proportions
profplot(x, kl)
# check consistency of most expressing cell types in signatures
g <-MarkerList(predict(x,"features"), names =names(m))
basismarkermap(g, kl)
# correlation with known signaturesbasiscor(x, kl)


########## use DECONF ############
#  deconvolve using KL divergence metric
dec <-ged(x, m,"deconf", rng = 12345, nrun = 20)
# plot against known proportions
profplot(x, dec)
# check consistency of most expressing cell types in signatures
basismarkermap(g, dec)
# correlation with known signatures
basiscor(x, dec)
