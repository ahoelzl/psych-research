library(psych)
source("faktorensimulation.R")
source("cmdsolve.R")
source("GesamtDatenAnalysen.R")
source("Vergleichsverfahren.R")
source("compareClusterings.R")

compareClusterings <- function(NLmean,NLsd,phimean,phisd, comparing=1,toSimulate, nrep=1, addError=F) {
  #  jpeg(paste("/home/andreas/Desktop/Bachelorarbeit/",NL,phi,comparing,".jpg"), width = 600, height = 400)
  
  
  
  results <- c()
  sumresults <- c(rep(0,length(toSimulate)))

    print(paste(NLmean,NLsd,phimean,phisd))
   corM <- setCorrelationMatrix(NLmean,NLsd,phimean,phisd, addError)
 #   cat("corMatrix set! ", NLmean,NLsd,phimean,phisd)
#    completecor <- completeCor(corM)
#    completecorcor <- completeCorCor(corM)
 #   averagecor <- averageCor(corM)
#    averagecorcor <- averageCorCor(corM)
#    kmeans <- cmdsolve(corM)
    

    results <- getResults(corM, toSimulate, comparing)

    results
    
    
    
   # completecorresult <- completecorresult + sum(vergleich(zuordnung.ges,completecor, comparing))
  #  completecorcorresult <- completecorcorresult +  sum(vergleich(zuordnung.ges,completecorcor,comparing))
  #  averagecorresult <- averagecorresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
  #  averagecorcorresult <- averagecorcorresult + sum(vergleich(zuordnung.ges,averagecorcor, comparing))
  #  kmeansresult <- kmeansresult + sum(vergleich(zuordnung.ges,kmeans, comparing))

  #results <<- c(completecorresult/nobs, completecorcorresult/nobs,averagecorresult/nobs,averagecorcorresult/nobs,  kmeansresult/nobs)
  
  
 # barplot(results, beside=T,main="Erkennen von Faktorenstrukturen",
#          xlab="Clusterverfahren bei", ylab=mainText, ylim=c(0,1.0), col=gray.colors(length(results)), sub=optionstext,cex.names=1.3,font.names=2,cex.axis=,cex.main=1.5,cex.lab=1.5,cex.sub=1.5)
  
 # legend(x="topleft",col=c("green","red"), legend= resultnames, fill=gray.colors(5), inset=c(0.01, 0), ce=0.4, bty="n")
}


getResults <- function(corM, toSimulate, zuordnung.ges, comparing) {
  completecorresult <- 0
  completecorresult <- 0
  completecorcorresult <- 0
  averagecorresult <- 0
  averagecorcorresult <- 0
  kmeansresult <- 0
  kmeansCorResult<-0
  kmeansNewResult<-0
  completecorcornometricresult <-0
  completecornometricresult <-0
  averagecornometricresult <- 0
  averagecorcornometricresult <- 0
  completeCorCorCorResult <- 0
  averageCorCorCorResult <- 0
  kmeanscorcorResult <- 0
  cmdsolveCorResult <- 0
  kMeansOnDistancesCorResult <- 0
  results <- c()
  resultnames <- toSimulate
  nrep <- 1
  for(sim in toSimulate) {
    print(sim)
    if(sim=="averagecor") {
      averagecor <- averageCor(corM)
      
      averagecorresult <- averagecorresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
      results <-append(results,averagecorresult/nrep)
    } else if(sim=="averagecorcor") {
      
      averagecorcor <- averageCorCor(corM)
      averagecorcorresult <- averagecorcorresult  + sum(vergleich(zuordnung.ges,averagecorcor, comparing))
      
      results <-append(results,averagecorcorresult/nrep)
    } else if(sim=="completecor") {
      
      completecor <- completeCor(corM)
      completecorresult <- completecorresult  + sum(vergleich(zuordnung.ges,completecor, comparing))
      results <-append(results,completecorresult/nrep)
    } else if(sim=="completecorcor") {
      
      completecorcor <- completeCorCor(corM)
      completecorcorresult <- completecorcorresult  + sum(vergleich(zuordnung.ges,completecorcor, comparing))
      results <-append(results,completecorcorresult/nrep)
    }  else if(sim=="kmeansmds") {   
      kmeans <- cmdsolve(corM)
      kmeansresult <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult/nrep)
    } else if(sim=="Dim1") {   
      kmeans <- cmdsolve(corM,dim=1)
      kmeansresult1 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult1/nrep)
    }  else if(sim=="Dim2") {   
      kmeans <- cmdsolve(corM,dim=2)
      kmeansresult2 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult2/nrep)
    }  else if(sim=="Dim3") {   
      kmeans <- cmdsolve(corM,dim=3)
      kmeansresult3 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult3/nrep)
    }  else if(sim=="Dim4") {   
      kmeans <- cmdsolve(corM,dim=4)
      kmeansresult4 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult4/nrep)
    } else if(sim=="Dim10") {   
      kmeans <- cmdsolve(corM,dim=10)
      kmeansresult10 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult10/nrep)
    } else if(sim=="Dim20") {   
      kmeans <- cmdsolve(corM,dim=20)
      kmeansresult20 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult20/nrep)
    } else if(sim=="Dim30") {   
      kmeans <- cmdsolve(corM,dim=30)
      kmeansresult30 <- kmeansresult  + sum(vergleich(zuordnung.ges,kmeans, comparing))
      results <-append(results,kmeansresult30/nrep)
    }  else if(sim=="kmeanskoord") {
      kmeansNew <- kMeansOnDistances(corM)
      kmeansNewResult <- kmeansNewResult  + sum(vergleich(zuordnung.ges,kmeansNew, comparing))
      results <-append(results,kmeansNewResult/nrep)
    }  else if(sim=="kmeanscor") {
      kmeansCor <- kmeansCor(corM)
      kmeansCorResult <- kmeansCorResult  + sum(vergleich(zuordnung.ges,kmeansCor, comparing))
      results <-append(results,kmeansCorResult/nrep)
    } else if(sim=="completecorcorcor") {
      completeCorCorCor <- completeCorCorCor(corM)
      completeCorCorCorResult <- completeCorCorCorResult  + sum(vergleich(zuordnung.ges,completeCorCorCor, comparing))
      results <-append(results,completeCorCorCorResult/nrep)
    } else if(sim=="averagecorcorcor") {
      averageCorCorCor <- averageCorCorCor(corM)
      averageCorCorCorResult <- averageCorCorCorResult  + sum(vergleich(zuordnung.ges,averageCorCorCor, comparing))
      results <-append(results,averageCorCorCorResult/nrep)
    } else if(sim=="kmeanscorcor") {
      print("kmeanscorcor")
      kmeanscorcor<- kmeansCorCor(corM)
      print(kmeanscorcor)
      kmeanscorcorResult <- kmeanscorcorResult  + sum(vergleich(zuordnung.ges,kmeanscorcor, comparing))
      print(print(kmeanscorcorResult))
      results <-append(results,kmeanscorcorResult/nrep)
    } else if(sim=="kmeansmdscor") {
      cmdsolveCor <- cmdsolveCor(corM)
      cmdsolveCorResult <-  cmdsolveCorResult   + sum(vergleich(zuordnung.ges, cmdsolveCor, comparing))
      results <-append(results,cmdsolveCorResult/nrep)
    }  else if(sim=="kMeansondistancescor ") {
      kMeansOnDistancesCor  <- kMeansOnDistancesCor(corM)
      kMeansOnDistancesCorResult <-  kMeansOnDistancesCorResult   + sum(vergleich(zuordnung.ges,kMeansOnDistancesCor, comparing))
      results <-append(results,kMeansOnDistancesCorResult/nrep)
    }  else if(sim=="averagecornom") {
      averagecor <- averageCorNoMetric(corM)
      averagecornometricresult <- averagecornometricresult  + sum(vergleich(zuordnung.ges,averagecor, comparing))
      results <-append(results,averagecornometricresult/nrep)
    } else if(sim=="averageccnom") {
      averagecorcor <- averageCorCorNoMetric(corM)
      averagecorcornometricresult <- averagecorcornometricresult + sum(vergleich(zuordnung.ges,averagecorcor, comparing))
      results <-append(results,averagecorcornometricresult/nrep)
    } else if(sim=="completecornom") {
      completecor <- completeCorNoMetric(corM)
      completecornometricresult <- completecornometricresult  + sum(vergleich(zuordnung.ges,completecor, comparing))
      results <-append(results,completecornometricresult/nrep)
    } else if(sim=="completeccnom") {
      completecorcor <- completeCorCorNoMetric(corM)
      completecorcornometricresult <- completecorcornometricresult  + sum(vergleich(zuordnung.ges,completecorcor, comparing))
      results <-append(results,completecorcornometricresult /nrep)
    }  else if(sim=="faclust") {
      completecorcor <- fclustering(corM)
      completecorcornometricresult <- completecorcornometricresult  + sum(vergleich(zuordnung.ges,completecorcor, comparing))
      results <-append(results,completecorcornometricresult /nrep)
    } 
  }
  results
  
}



# getClusterSimiliarity.simulation <- function(NL.mus, Kor.mus, toSimulate, comparing) {
#   #compareClusterings <- function(NLmean,NLsd,phimean,phisd, comparing,toSimulate, nrep=1, addError=F) 
#   ##Die Bedingungen in den R's werden hier erzeugt
#   
#   r.names <- c()
#   descriptions <- ""
#   for(i in 1:length(NL.mus)) {
#     r.names[i] <- paste0("constr ", i)
#     descriptions <- paste0(descriptions,  " und " , r.names[i] , " mit NL von ",
#                            NL.mus[i], " und Faktorkorrelation von ", Kor.mus[i], " \n  ")
#   }
#   
#   rs <- matrix(nrow= length(NL.mus), ncol=length(toSimulate)) 
#   rownames(rs) <- r.names
#   colnames(rs) <- toSimulate
#   for(i in 1:length(NL.mus))  {
#     r1 <- compareClusterings(NL.mus[i],0,Kor.mus[i],0,comparing=1,toSimulate, addError=F)
#     print("now")
#     print(r1)
#     rs[i,] <- r1
#   }
#   
#   paintTable(rs, "Clusterübereinstimmung bei EFA-Simulation", paste0("type ",type, "\n" , descriptions))
#   rs
# }



###method 1 : NL alle gleich, ergeben zusammen Kommunalität
####### nur eine NL; entspricht Kommunalität
####### zwei, entsprechen zusammen Kommunalität
getClusterSimiliarity.simulation.methods <- function(methods, zuordnung.ges,  toSimulate, fa.ges) {
  #compareClusterings <- function(NLmean,NLsd,phimean,phisd, comparing,toSimulate, nrep=1, addError=F) 
  ##Die Bedingungen in den R's werden hier erzeugt
  
  r.names <- c()
  descriptions <- ""
  #for(i in 1:length(NL.mus)) {
  #  r.names[i] <- paste0("constr ", i)
  #  descriptions <- paste0(descriptions,  " und " , r.names[i] , " mit NL von ",
  #                         NL.mus[i], " und Faktorkorrelation von ", Kor.mus[i], " \n  ")
  #}
  descriptions<- "bei allen NL gleich und entsprechen Kommunalität (NL.equal)"
  
  rs <- matrix(nrow= 3, ncol=length(toSimulate)) 
  
  for( i in 1:length(methods)) {
    
  method <- methods[i]  
    
  if(method==2) {
    descriptions<- "bei einer NL und entsprechen Kommunalität (NL.one)"
  } else if(method==3) {
    descriptions<- "bei zwei NL und entsprechen Kommunalität (NL.two)"
  }
  

  colnames(rs) <- toSimulate
  #for(i in 1:length(NL.mus))  {
    #r1 <- compareClusterings(NL.mus[i],0,Kor.mus[i],0,1,toSimulates, addError=addError)
    
 
 loads <- NL.equal(fa.ges$loadings)
  
 # loads <- NL.fixed(fa.ges$loadings, 0.2)
  
  if(method==2) {
    loads <- NL.one(fa.ges$loadings)
  } else if(method == 3) {
    loads <- NL.two(fa.ges$loadings)
  }

  
  Phi <- fa.ges$Phi
  
  corM <- sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=0)$model
  
 
  rs[i,] <- getResults(corM, toSimulate,zuordnung.ges, comparing=1)
  
  }
  
  rownames(rs) <- c("Sim1", "Sim2", "Sim3")
  
  
  
  
  
  
  
  paintTable(rs, "Clusteruebereinstimmung bei EFA", paste0("\n" , descriptions))
  rs
}


###method 1 : NL alle gleich, ergeben zusammen Kommunalität
####### nur eine NL; entspricht Kommunalität
####### zwei, entsprechen zusammen Kommunalität
getClusterSimiliarity.simulation.samples.methods <- function(methods, zuordnung.ges,  toSimulate, fa.ges, nobs, nrep) {
  #compareClusterings <- function(NLmean,NLsd,phimean,phisd, comparing,toSimulate, nrep=1, addError=F) 
  ##Die Bedingungen in den R's werden hier erzeugt
  
  rs.list <- list()
  
  for(c in 1:nrep) {
  
  r.names <- c()
  descriptions <- ""
  #for(i in 1:length(NL.mus)) {
  #  r.names[i] <- paste0("constr ", i)
  #  descriptions <- paste0(descriptions,  " und " , r.names[i] , " mit NL von ",
  #                         NL.mus[i], " und Faktorkorrelation von ", Kor.mus[i], " \n  ")
  #}
  descriptions<- "bei allen NL gleich und entsprechen Kommunalität (NL.equal)"
  
  rs <- matrix(nrow= 3, ncol=length(toSimulate)) 


  for( i in 1:length(methods)) {
    
    method <- methods[i]  
    
    if(method==2) {
      descriptions<- "bei einer NL und entsprechen Kommunalität (NL.one)"
    } else if(method==3) {
      descriptions<- "bei zwei NL und entsprechen Kommunalität (NL.two)"
    }
    
    
    colnames(rs) <- toSimulate
    #for(i in 1:length(NL.mus))  {
    #r1 <- compareClusterings(NL.mus[i],0,Kor.mus[i],0,1,toSimulates, addError=addError)
    
    
 
    loads <- NL.equal(fa.ges$loadings)
    Phi <- fa.ges$Phi
    # loads <- NL.fixed(fa.ges$loadings, 0.2)
    
    if(method==2) {
      loads <- NL.one(fa.ges$loadings)
    } else if(method == 3) {
      loads <- NL.two(fa.ges$loadings)
    }
    
    corM <- sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=nobs, raw=T)$r
 
    
    rs[i,] <- getResults(corM, toSimulate,zuordnung.ges, comparing=1)
  }
  
  rownames(rs) <- c("Sim1", "Sim2", "Sim3")
  rs.list[[c]] <- rs
  }
  
  rs <- rs.list[[1]]
  if(length(rs.list) > 1) {
  for(l in 2:length(rs.list)) {
  rs <- rs + rs.list[[l]]
  }
  }
  rs <- rs/length(rs.list)
  
  
  paintTable(rs, "Clusteruebereinstimmung bei EFA-Stichproben", paste0("\n" , descriptions, " nrep ", nrep, " nobs ", nobs))
  rs
}
