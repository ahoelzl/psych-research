library(NbClust)
library(clValid)




#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluadvanced <- function (data,nobs,nrep,type="kmeans") {
  
  v <- c()
  for (i in 1:nrep){
    samps <- sample((x=1:nrow(data)), size=nobs, replace=T)
    
    print(length(samps))
    daten.sp <- data[samps,]
    
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    print(cor[1,2])
    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    result <- getClusterNumbers(points=points, type=type)
    number.cluster <- as.numeric(as.character(optimalScores(result)[,3]))
    print( number.cluster)
    v[[i]] <- number.cluster 
  }
  v
}



#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluadvanced.whole <- function (data,type="kmeans") {
  
    
    daten.sp <- data
    
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")

    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    result <- getClusterNumbers(points=points, type=type)
    cat("result: ", result)

  result
}

getClusterNumbers <- function(points, type="kmeans") {
  result <- 0
  if(type=="kmeans" || type == "kmeansmds") {
    result <- clValid(obj=points, nClust=2:15,clMethods="kmeans", validation=c("internal","stability"))
  } else if(type=="complete" || type=="completecor" || type=="completecorcor") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
  } else if(type=="average" || type=="averagecor" || type=="averagecorcor") {
    result <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="average")
  }else {
   cat("not found! ", type)
   result <- 0
  }
  
  number.cluster <- as.numeric(as.character(optimalScores(result)[,3]))
  names(number.cluster) <- method.names.normal
  number.cluster
}


EFA.Cluster.number <- function(daten.sp,n) {
  
  data <- daten.sp
  map.sp <- VSS(daten.sp, rotate = "promax", fm = "mle", title="Anzahl der Faktoren", n.obs=n)
  map <-    which.min(map.sp$map)
  
  
  pa.sp <- fa.parallel(daten.sp, fm="ml",n.iter=100, n.obs=n)
  
  paralell.ncomp <- pa.sp$ncomp
  paralell.nfact <- pa.sp$nfact
  
  result = tryCatch({
    aic <- 1:14
    for (j in 1:14) {
      fa.sp <- fa(daten.sp, nfactors=j, max.iter=100, fm="ml", rotate="promax", method="pearson", n.obs=n)
      aic[j] <- (fa.sp$STATISTIC)-(2*(ncol(data)*(ncol(data)-1)/2-(ncol(data)*j+(j*(j-1)/2))))
    }
    aicmin <- which.min(aic)
},  error = function(e) {
  aicmin <- paralell.ncomp
})
 
  
  
  clusternumber.EFA.methods <<- c("MAP", "Paralell-mcomp", "Paralell-nfact", "AIC")
  clusternumbers <- c(map, paralell.ncomp, paralell.nfact, aicmin )
  
}

drawNumberClusterAdvanced<- function(facs,nrep, type="kmeans") {

  whole.cluster.number <- c()
if(type=="kmeans")   {
  whole.cluster.number <- whole.cluster.number.kmeans
} else if(type=="average") {
  whole.cluster.number <- whole.cluster.number.average
} else if(type=="complete") {
  whole.cluster.number <- whole.cluster.number.complete
}
  
method.names <- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")
result.names <- c("whole", "Var", "Bias")  
  nobs <- c(500)
  
  
  resultsmatrix <- matrix(nrow=length(result.names), ncol=length(method.names))
  rownames(resultsmatrix) <- result.names
  colnames(resultsmatrix) <- method.names

  for(o in 1:length(nobs)) {
  vector <- as.vector(numcluadvanced(data=facs, nobs=nobs[o], nrep=nrep, type=type))

    m <- matrix( ncol=length(vector), nrow=length(vector[[1]]))
    
    for(i in 1:length(vector)) {
      for(j in 1:length(vector[[1]])) {
        m[j,i] <- vector[[i]][j] - whole.cluster.number[j]
      }
    }
    
  
  
   # par(mfrow=c(3,3))
    ###aufteilen auf vektoren der einzelnen Methoden, die dann geplottet werden
    for(i in 1:dim(m)[1]) {
  #    drawBarplot(m[i,],ylab=paste("Faktorenanalyse"),nob=nob,type=type, cex.lab=1.5, method= measNames(result)[i])
      method= measNames(result)[i]
      method.var <- var(m[i,])
      method.bias <- mean(m[i,])
      method.whole <-  whole.cluster.number[i]
       
      resultsmatrix[1,i] <- method.whole
      resultsmatrix[2,i] <- method.var
      resultsmatrix[3,i] <- method.bias
      
    }
  
  }
  
  resultsmatrix
  
}



if(!exists("whole.cluster.number.kmeans")) {
whole.cluster.number.kmeans <- numcluadvanced.whole(facs, type="kmeans")
whole.cluster.number.average <- numcluadvanced.whole(facs, type="average")
whole.cluster.number.complete <- numcluadvanced.whole(facs, type="complete")
}

cluster.numbers.whole.complete
#r1 <- drawNumberClusterAdvanced(facs,nrep=10, type="kmeans")

#completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs,nrep=30, type="average")

#completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs,nrep=30, type="complete")

#t <- drawKmeans(1,type="average") 

#test <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
