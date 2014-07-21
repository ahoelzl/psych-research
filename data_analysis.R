### read data set (SPSS data file) ### 

library(foreign)
IST.ges <- read.spss("IST Grundmodul.sav",use.value.labels=F,to.data.frame=T)

## only sum-variables
IST <- IST.ges[,c(4:12)]

facs <- IST




corM <- cor(facs, use="pairwise.complete.obs", method="pearson")

fa.ges <- fa(facs, nfactors=3, max.iter=100, fm="ml", rotate="promax", method="pearson")
comparing <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))

zuordnung.ges <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))
