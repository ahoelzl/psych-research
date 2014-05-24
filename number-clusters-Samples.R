library(NbClust)
library(clValid)

source("faktorensimulation.R")
source("cmdsolve.R")
source("Vergleichsverfahren.R")

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
    number.cluster <- getClusterNumbers(points=points,cor.sp = cor, type=type, n=nobs)
    
    
    print( number.cluster)
    v[[i]] <- number.cluster 
  }
  v
}



#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluadvanced.whole <- function (data,type="kmeans", numbermethod = "internal", n) {
  
    
    daten.sp <- data
    
    cor <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")
    print(cor[1,2])
    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)

    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points
    cat("2")
    
    points.save <<- points
    type.save <<- type

    number.cluster <- getClusterNumbers(points=points, cor.sp=cor, type=type, numbermethod, n=n)
 
    print( number.cluster)
    v <- number.cluster 
 
    v
}



getClusterNumbers <- function(points,cor.sp, type="kmeans", numbermethod="internal", n) {
  
  validationtype <- c("internal")
  if(numbermethod %in% c("APN", "ADM", "AD", "FOM")) {
    validationtype <- c("stability")
  }
  
  if(numbermethod == "both") {
    validationtype <- c("stability", "internal")
  }
  
  dims <- dim(cor.sp)[1]
  minv <- min(dims-1, 10)
  
  if(type=="kmeans" || type=="kmeansmds") {
    result <- clValid(obj=points, nClust=2:minv,clMethods="kmeans", validation=c(validationtype), verbose=T)
    result.names <-  measNames(result)
     #result <-   as.numeric(as.character(optimalScores(result)[,3]))
    measures <- measures(result)
    mins <-  apply(measures,1, function(x)  which(x == min(x)))
    maxs <- apply(measures,1, function(x)  which(x == max(x)))
    
    result <- c(unlist(mins[1])[1] + 1 , unlist(maxs[2])[1] + 1 ,unlist(maxs[3])[1]+1)
    if(is.na(result[1])) {
    result[1] <- 2
    }
    names(result) <- result.names
  } else if(type=="complete" || type=="completecor" || type=="completecorcor") {
    result <- clValid(obj=points, nClust=2:minv,clMethods="hierarchical", validation=c(validationtype), method="complete")
    result.names <-  measNames(result)
    
    #result <-   as.numeric(as.character(optimalScores(result)[,3]))
    measures <- measures(result)
    mins <-  apply(measures,1, function(x)  which(x == min(x)))
    maxs <- apply(measures,1, function(x)  which(x == max(x)))
    
    result <- c(unlist(mins[1])[1] + 1 , unlist(maxs[2])[1] + 1 ,unlist(maxs[3])[1]+1)
    if(is.na(result[1])) {
      result[1] <- 2
    }
    names(result) <- result.names
  } else if(type=="average" || type=="averagecor" || type=="averagecorcor") {
    result <- clValid(obj=points, nClust=2:minv,clMethods="hierarchical", validation=c(validationtype), method="average")
    result.names <-  measNames(result)
    #result <-   as.numeric(as.character(optimalScores(result)[,3]))
    measures <- measures(result)
    mins <-  apply(measures,1, function(x)  which(x == min(x)))
    maxs <- apply(measures,1, function(x)  which(x == max(x)))
    
    result <- c(unlist(mins[1])[1] + 1 , unlist(maxs[2])[1] + 1 ,unlist(maxs[3])[1]+1)
    if(is.na(result[1])) {
      result[1] <- 2
    }
    names(result) <- result.names
  } else if(type=="varclust") {
    
    dist.m <- dist(points, method="euclidean")
    dunns <- c()
    conns<- c()
    sills <- c()
    for(t in 2:(minv-1)) {
      cluster <- varClust(points, k=t)
      dist.m = dist(t(points), method="euclidean")
      names(cluster) <- rownames(as.matrix(dist.m))
      conns <- append(conns,connectivity(dist.m, cluster))
      dunns <- append(dunns,dunn(dist.m, cluster))
      sil <- silhouette(cluster, dist.m)
      sills <- append(sills,  summary(sil)$avg.width)
    }
    
    con <- which.min(conns)+1
    if(length(con)==0) {
      con <- 2
    }
    result <- c(con, which.max(dunns)+1,which.max(sills)+1)
    
    } else if(type=="varclust2"){
      dist.m <- dist(points, method="euclidean")
      dunns <- c()
      conns<- c()
      sills <- c()
      for(t in 2:(minv-1)) {
        cluster <- varClust2(points, k=t)
        dist.m = dist(t(points), method="euclidean")
        names(cluster) <- rownames(as.matrix(dist.m))
        conns <- append(conns,connectivity(dist.m, cluster))
        dunns <- append(dunns,dunn(dist.m, cluster))
        sil <- silhouette(cluster, dist.m)
        sills <- append(sills,  summary(sil)$avg.width)
      }
      
      con <- which.min(conns)+1
      if(length(con)==0) {
        con <- 2
      }
      result <- c(con, which.max(dunns)+1,which.max(sills)+1)
    }  else if(type=="faclust") {
    result <- EFA.Cluster.number(cor.sp = cor.sp, n=n)
    ##(type=="kmeanscor")
  } else  {
  
    result <- clValid(obj=cor.sp, nClust=2:minv,clMethods="kmeans", validation=c(validationtype))
    result.names <-  measNames(result)
    #result <-   as.numeric(as.character(optimalScores(result)[,3]))
    measures <- measures(result)
    mins <-  apply(measures,1, function(x)  which(x == min(x)))
    maxs <- apply(measures,1, function(x)  which(x == max(x)))
    
    result <- c(unlist(mins[1])[1] + 1 , unlist(maxs[2])[1] + 1 ,unlist(maxs[3])[1]+1)
    if(is.na(result[1])) {
      result[1] <- 2
    }
    names(result) <- result.names
  }

  result
}



EFA.Cluster.number <- function(cor.sp, n) {
  
  cat("n ", n)
  data <- cor.sp
  daten.sp <- cor.sp
  map.sp <- VSS(daten.sp, rotate = "promax", fm = "mle", title="Anzahl der Faktoren", plot=F, n.obs=n)
  map <-    which.min(map.sp$map)
  
  
  pa.sp <- fa.parallel(daten.sp, fm="ml",n.iter=100, n.obs=n)

  paralell.ncomp <- pa.sp$ncomp
  paralell.nfact <- pa.sp$nfact
  
 aicmin <-  tryCatch(expr={
  aic <- 1:14
  for (j in 1:14) {
    fa.sp <- fa(daten.sp, nfactors=j, max.iter=100, fm="ml", rotate="promax", method="pearson", n.obs=n)
    aic[j] <- (fa.sp$STATISTIC)-(2*(ncol(data)*(ncol(data)-1)/2-(ncol(data)*j+(j*(j-1)/2))))
  }
  aicmin <- which.min(aic)
  }, error=function(cond) {
    aicmin <-  paralell.ncomp
  })
  
  

  clusternumbers <- c(map, paralell.ncomp, paralell.nfact, aicmin )
  
}




drawNumberClusterAdvanced<- function(facs,nrep, type="kmeans", nobs = 200) {

  whole.cluster.number <- c()
if(type=="kmeans")   {
  whole.cluster.number <- whole.cluster.number.kmeans
} else if(type=="average") {
  whole.cluster.number <- whole.cluster.number.average
} else if(type=="complete") {
  whole.cluster.number <- whole.cluster.number.complete
} else if(type=="faclust") {
  whole.cluster.number <-  whole.cluster.number.faclust
}  else if(type=="kmeanscor") {
  whole.cluster.number <-  whole.cluster.number.kmeanscor
} else if(type=="kmeansmds") {
  whole.cluster.number <- whole.cluster.number.kmeans
} else if(type=="varclust"){
  whole.cluster.number <- whole.cluster.number.varclust
} else if(type=="varclust2"){
  whole.cluster.number <- whole.cluster.number.varclust2
}

  
method.names <- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")
result.names <- c("whole", "Var", "Bias")  

  
  
  resultsmatrix <- matrix(nrow=length(result.names), ncol=length(method.names))
  rownames(resultsmatrix) <- result.names
  colnames(resultsmatrix) <- method.names

  for(o in 1:length(nobs)) {
  vector <- as.vector(numcluadvanced(data=facs, nobs=nobs, nrep=nrep, type=type))

    m <- matrix( ncol=length(vector), nrow=length(vector[[1]]))
    
    for(i in 1:length(vector)) {
      for(j in 1:length(vector[[1]])) {
        whole.number <-  whole.cluster.number[j]
        if(class(whole.number)=="list") {
          whole.number <- unlist(whole.number)
        }
        
        v <-  vector[[i]][j]
        
        if(class(v)=="list") {
          v <- unlist(v)
        }
        m[j,i] <- v- whole.number
      }
    }
    
  
  cat("m:", m)
  mglobal <<- m
  
   # par(mfrow=c(3,3))
    ###aufteilen auf vektoren der einzelnen Methoden, die dann geplottet werden
    for(i in 1:dim(m)[1]) {
  #    drawBarplot(m[i,],ylab=paste("Faktorenanalyse"),nob=nob,type=type, cex.lab=1.5, method= measNames(result)[i])
      if(type=="faclust") {
        method = method.names.EFA[i]
      } else {
      method= method.names[i]
      }
      
      
      method.var <- var(m[i,])
      method.bias <- mean(m[i,])
      method.whole <-  whole.cluster.number[i]
       
      if(class(method.whole)=="list") {
        method.whole <- unlist(method.whole)
      }
      resultsmatrix[1,i] <- round(method.whole, digits=4)
      resultsmatrix[2,i] <- round(method.var, digits=4)
      resultsmatrix[3,i] <- round(method.bias, digits=4)
    }
  }
  
  cat("resultsmatrix   " , resultsmatrix)
  
  global.resultsmatrix <<- resultsmatrix
  
  resultsmatrix
  
}

#nur laufen lassen wenn noch nicht vorhanden
#if(!exists("whole.cluster.number.kmeans")) {
whole.cluster.number.kmeans <- numcluadvanced.whole(facs, type="kmeans")
whole.cluster.number.average <- numcluadvanced.whole(facs, type="average")
whole.cluster.number.complete <- numcluadvanced.whole(facs, type="complete")
whole.cluster.number.faclust <- numcluadvanced.whole(facs, type="faclust", n=dim(facs)[1])
whole.cluster.number.kmeanscor <- numcluadvanced.whole(facs, type="kmeanscor")
whole.cluster.number.varclust <- numcluadvanced.whole(facs, type="varclust")
whole.cluster.number.varclust2 <- numcluadvanced.whole(facs, type="varclust2")
#}


getClusterNumberBiasVariance.samples <- function(nrep, types, nobs) {
  
  rs <- matrix(nrow = 5, ncol= (length(types) - 1) * length(method.names) + length( method.names.EFA) + 20 )
  rs[1, ] <- ""
  
  globalsave <<- c()
for(i in 1:(length(types))) {
  
  cat("i", i)
  type <- types[i]
r1 <- drawNumberClusterAdvanced(facs,nrep=nrep, type=type, nobs=nobs)
 # globalsave[[i]] <<- r1
  #r1 <- globalsave[[i]]
  cat(paste0("type : ", type, " ", r1))
rs[1, (i-1) * length(method.names) + 1] <- types[i]
  if(type=="faclust") {
    rs[2, (i-1) * length(method.names) + 1:length(method.names.EFA)] <- method.names.EFA
    rs[3:5, (i-1) * length(method.names) + 1:length(method.names.EFA)] <- r1[,1:length(method.names.EFA)]
  } else {
    rs[2, (i-1) * length(method.names) + 1:length(method.names)] <- method.names
    rs[3:5, (i-1) * length(method.names) + 1:length(method.names)] <- r1[,1:length(method.names)]
  }

#rs 

}
  
  rownames(rs) <- c("clustertype", "clusternumber", "whole", "Var", "Bias")
  
  
  names <-  rs[2,] %in% c("APN","AD","ADM","FOM",NA) 
  numbers <-  rs[5,] %in% c(NA) 
  
  clustertype <- rep(types,each=3)
  clustertype <- c(clustertype,"faclust")
  
  clusternumber <- rs[2,][!names]
  whole <- rs[3,][!numbers]
  var <- rs[4,][!numbers]
  bias <- rs[5,][!numbers]
  
  result <- cbind(clustertype, clusternumber, whole,var,bias)
  
result
  
  paintTable(result, "Clusteranzahlsgenauigkeit bei Samples", paste0("type ",type, " nobs ", nobs))
}

#r1 <- getClusterNumberBiasVariance.samples(50,"kmeans")
#r2 <- getClusterNumberBiasVariance.samples(50,"average")
#r3 <- getClusterNumberBiasVariance.samples(50,"complete")
#completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs,nrep=30, type="average")

#completeCCorrelationCorrelation <- drawNumberComparisonKmeans(facs,nrep=30, type="complete")

#t <- drawKmeans(1,type="average") 

#test <- clValid(obj=points, nClust=2:15,clMethods="hierarchical", validation=c("internal","stability"), method="complete")
