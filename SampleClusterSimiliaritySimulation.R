source("compareClusterings.R")

simulate <- function(cor.sp,nobs,nrep,k,method,wholecut,compareMethod,data,...) {
  
  print(method)
  v1 <- 1:nrep
  v2 <- 1:nrep
  w1 <- list()
  w2 <- list()
  l1  <-  list()
  l2  <-  list()
 
  
  
  cut.sp1 <- c()
  
  vergleich1 <- 0
  n11 <- 0
  n00 <- 0
    if(method=="completecor")  {
      cut.sp1 <- completeCor(cor.sp,k)
    } else if(method=="completecorcor") {
      cut.sp1 <- completeCorCor(cor.sp,k)
    } else if(method=="averagecor") {
      cut.sp1 <- averageCor(cor.sp,k)
    } else if(method=="averagecorcor") {
      cut.sp1 <- averageCorCor(cor.sp,k)
    } else if(method=="kmeansmds") {
      cut.sp1 <- cmdsolve(cor.sp,k)
    }else if(method=="Dim1") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=1)
    }else if(method=="Dim2") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=2)
    }else if(method=="Dim3") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=3)
    }  else if(method=="Dim4") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=4)
    }  else if(method=="Dim10") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=10)
    }  else if(method=="Dim20") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=20)
    }   else if(method=="Dim30") {
      cut.sp1 <- cmdsolve(cor.sp,k,dim=30)
    } else if(method=="kmeanskoord") {
      cut.sp1 <- kMeansOnDistances(cor.sp,k)
    } else if(method=="kmeanscor") {
      cut.sp1 <- kmeansCor(cor.sp,k)
    } else if(method=="completecornom")  {
      cut.sp1 <- completeCorNoMetric(cor.sp,k)
    } else if(method=="completecorcornom") {
      cut.sp1 <- completeCorCorNoMetric(cor.sp,k)
    } else if(method=="averagecornom") {
      cut.sp1 <- averageCorNoMetric(cor.sp,k)
    } else if(method=="averagecorcornom") {
      cut.sp1 <- averageCorCorNoMetric(cor.sp,k)
    }  else if(method=="totalclust") {
      cut.sp1  <- kmeans(t(data),centers=5,nstart=100)$cluster
    } else if(method=="averagecorcorcor") {
      cut.sp1  <- averageCorCorCor(cor.sp,k)
    } else if(method=="completecorcorcor") {
      cut.sp1  <- completeCorCorCor(cor.sp,k)
    }  else if(method=="kmeansmdscor") {
      cut.sp1  <- cmdsolveCor(cor.sp)
    } else if(method=="kmeanscorcor") {
      cut.sp1  <- kmeansCorCor(cor.sp,k)
    } else if(method=="kmeansneucor") {
      cut.sp1  <- kMeansOnDistancesCor(cor.sp,k)
    } else if(method=="faclust") {
      cut.sp1 <- fclustering(cor.sp,k)
    } else if(method=="clustofvar") {
      cut.sp1 <- varClust(data,k)
    }

  print(cut.sp1)
  print(wholecut)
    n00 <- n00 + vergleich(wholecut,cut.sp1, compareMethod)[1]
    n11 <- n11 + vergleich(wholecut,cut.sp1, compareMethod)[2]
print(n00)

  result <- c(n00/nrep, n11/nrep)
  print(result)
  result
}


simulateClusterSamplesComparison <- function(facs,toSimulate,compareWith, allnobs, nrep,numbercluster,compareMethod ) {
  simulationresults <- c()
  counter <- 1
for(nobs in allnobs) {
 
  #numbercluster4
sumvalues <- c(rep(0,2*length(toSimulate)))
resultnames <- c()
resultvalues <- c()
  for (i in 1:nrep){
    
    daten.sp <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    resultvalues <- c()

    cor.sp <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
   for(sim in toSimulate) {
     cat(sim)
    newresult <-  simulate(cor.sp,nobs,nrep, k=numbercluster,method=sim,compareWith[[sim]], compareMethod,daten.sp)
    cat("rrrsult " , newresult)
    resultvalues <- append(resultvalues,newresult)
    resultnames <- sim
  }
    sumvalues  <- sumvalues + resultvalues
  }
  resultvalues <- sumvalues
  print(resultvalues)
  #resultnames <- c("Average,  Korrelation", "Complete, Korrelation",
   #                "Averag, Korrelation der Korrelation", "Complete, Korrelation der Korrelation", "kmeans cmd")
  
  simulationinfo <- list()
  
  simulationinfo[["resultvalues"]] <- resultvalues
  
  simulationinfo[["resultnames"]] <- resultnames
  
  simulationinfo[["samplesize"]] <- nobs
  
  simulationinfo[["nreps"]] <- nrep
  
#  simulationinfo[["numbercluster"]] <- numbercluster
  
  simulationresults[[counter]] <-  simulationinfo

  counter <- counter + 1
}

  print(simulationresults)

##bringt diese Daten in ein richtiges Format für die grafische Darstellung
dataasvector11 <- c()
dataasvector00 <- c()
samplesizes <- c()


for(i in 1:length(simulationresults)) {
  print(simulationresults[[i]][["resultvalues"]])
  for(j in 1:length(simulationresults[[i]][["resultvalues"]])) { 
    dataasvector00 <- as.vector(na.omit(append(dataasvector00, simulationresults[[i]][["resultvalues"]][2*j-1])))
    samplesizes <- append(samplesizes,paste(simulationresults[[i]][["samplesize"]]) )
  }
}


for(i in 1:length(simulationresults)) {
  print(simulationresults[[i]][["resultvalues"]])
  for(j in 1:length(simulationresults[[i]][["resultvalues"]])) { 
    dataasvector11 <- as.vector(na.omit(append(dataasvector11, simulationresults[[i]][["resultvalues"]][2*j])))
    samplesizes <- append(samplesizes,paste(simulationresults[[i]][["samplesize"]]) )
  }
}


#dataasmatrix00 <- matrix(c(dataasvector00), ncol=length(toSimulate), byrow=T)
#dataasmatrix11 <- matrix(c(dataasvector11), ncol=length(toSimulate), byrow=T)
dataasmatrixall <- matrix(c(dataasvector00+dataasvector11), ncol=length(toSimulate), byrow=T)


#relevant für den Rand-Index in dem nur die wahr positiven oder wahr negativen
#und für den Rand-Index Spezifität und Sensitivität
#ansonsten immer dataasmatrixall verwenden 
dataasmatrix <-  dataasmatrixall
mainText <- "Rand-Index mit Clusteranzahl=5"
#dataasmatrix <- dataasmatrix11
#mainText <- "Rand-Index Sensitivität"
#dataasmatrix <- dataasmatrix00
#mainText <- "Rand-Index Spezifität"
if(compareMethod==2) {
  mainText <- "Meila-Heckerman-Measure"
} 


colnames(dataasmatrix) <- toSimulate
rownames(dataasmatrix) <- allnobs

print(t(dataasmatrix))

  
  dataasmatrix
}


#corcorM <- cor(corM, use="pairwise.complete.obs", method="pearson")
#if(!exists("averagecor")) {
  k<-5
averagecor <- averageCor(corM,k) 
averagecorcor <- averageCorCor(corM,k) 
completecor <- completeCor(corM,k)
completecorcor <- completeCorCor(corM,k)
kmeansmds <- cmdsolve(corM,k)


varclust <- varClust(facs, k)
kmeanscor <- kmeansCor(corM,k)
#averagecornometric <- averageCorNoMetric(corM,k) 
#averagecorcornometric <- averageCorCorNoMetric(corM,k) 
#completecornometric <- completeCorNoMetric(corM,k)
#completecorcornometric <- completeCorCorNoMetric(corM,k)
#completecorcornometric <- completeCorCorNoMetric(corM,k)
#totalclust <- kmeans(t(facs),centers=5,nstart=100)$cluster
fclust <- fclustering(corM,k)
#averagecorcorcor <- averageCorCorCor(corM)
#completecorcorcor <- completeCorCorCor(corM)
#kmeansmdscor <- cmdsolveCor(corM,k)
#kmeanscorcor <- kmeansCorCor(corM,k)
#kmeansneucor <- kMeansOnDistancesCor(corM,k)
#}


getClusterSimiliarity.samples <- function(nrep, numbercluster) {
  
  r <- simulateClusterSamplesComparison(facs,toSimulate,compareWith, allnobs, nrep=nrep,numbercluster=numbercluster,compareMethod=1 )
  paintTable(r, "Clusterübereinstimmung bei Simulation von Stichprobendaten", paste0("nrep ", nrep, " clusternumber ", numbercluster))
}




compareWith <- list("averagecor" =averagecor, "completecor" = completecor, "averagecorcor" = averagecorcor,
                    "completecorcor"= completecorcor,"kmeansmds" = kmeansmds, "faclust" = fclust, "kmeanscor" = kmeanscor, "clustofvar"=varclust)



