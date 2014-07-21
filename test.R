maptsim <- function(nobs,nrep) {
  library(psych)
  
  
  numfac <- 1:nrep
  
  set.seed(1000)
  for (i in 1:nrep){
    
    sp <- sim.structure(fx=fx,Phi=Phi,n=nobs,uniq=fa.ges$uniquenesses,items=T,cat=5)  
    
    
    corM <- sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=nobs, items=T, cat=5)$r
    
    map.sp <- VSS(sp$r, n = 10, rotate = "promax", diagonal = FALSE, fm = "mle", n.obs=nobs,plot=TRUE,title="Anzahl der Faktoren")
    
    numfac[i] <- which.min(map.sp$map)
    # sqres[i] <- map.sp$vss.stats$sqresid
    results <- numfac
    
  }
  results
}

no<-100


RMSE <- rep(NA,length(no))
sd <- rep(NA,length(no))

for (i in 1:length(no)){
  
  out1 <- maptsim(no[i],1000) 
  RMSE[i] <- sqrt((sum((5-out1[1])^2))/1000)
  sd[i] <- sd(out1)
}
RMSEmap <- RMSE
sdmap <- sd