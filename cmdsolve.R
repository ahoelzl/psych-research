library(rootSolve)
library(GLDEX)
model <- function(x) {
  
  returnVector <- c()
  
  size <- length(distmatrix[1,])
  
  diffMatrix <- matrix(0,nrow=size,ncol=(size-1))
  counter <- 1
  for(i in 2:size) {
    for(j in 1:(i-1)) {
      diffMatrix[i,j] <- counter
      counter <-counter+1
    }
  }
  
  for(i in 1:(size-1)) {
    for(j in (i+1):size) {
      sumResult <- 0
      for(k in 1:(size-1)) {
        #   cat("i" , i, " j: ",  j , " k: ", k);
        #  cat("first " , (diffMatrix[i,k]), " second: ",  diffMatrix[j,k] , "\n");
        if(diffMatrix[i,k] != 0 && diffMatrix[j,k] != 0) {
          #   cat("add1 " , diffMatrix[i,k], diffMatrix[j,k]," ^2 ", "\n");
          sumResult <- sumResult + (x[diffMatrix[i,k]] - x[diffMatrix[j,k]])^2
        }
        if(diffMatrix[i,k] == 0 && diffMatrix[j,k] != 0) {
          #  cat("add2" , diffMatrix[j,k] , "^2)" , "\n");
          sumResult <- sumResult + x[diffMatrix[j,k]]^2
        }
        if(diffMatrix[i,k] != 0 && diffMatrix[j,k] == 0) {
          #     cat("add3 " , diffMatrix[i,k] , "^2" , "\n");
          sumResult <- sumResult+ (x[diffMatrix[i,k]])^2
        }
        
      }
      sumResult <- sumResult - distmatrix[i,j]
      # cat("eine Zeile: " , sumResult, "\n")
      returnVector <- append(returnVector, sumResult)
    }
  }
  
  returnVector
}


oldModel <- function(x) {
  spass <- 1
  cat("x1 " , x[spass] , "\n");
  F1 <- x[spass]^2 - 1
  F2 <- x[2]^2 + x[3]^2  - 5
  F3 <- x[4]^2 + x[5]^2 + x[6]^2 - 14
  F4 <- (x[1] - x[4])^2 + x[5]^2 + x[6]^2 - 13
  F5 <- (x[1] - x[2])^2 + x[3]^2 - 4
  F6 <- (x[2] - x[4])^2 + (x[3] - x[5])^2 + x[6]^2 - 9
  
  
  F7<-sum(F1,F2,F3)
  cat("F7 " , F7 , "\n");
  c(F1, F2, F3, F4, F5, F7)
}

callModel2 <- function(distMatrix) {
  
  size <- length(distMatrix[1,])
  
  diffMatrix <- matrix(0,nrow=size,ncol=(size-1))
  counter <- 1
  for(i in 2:size) {
    for(j in 1:(i-1)) {
      diffMatrix[i,j] <- counter
      counter <-counter+1
    }
  }
  
  returnvector <- c()
  
  
  coordMatrix <<- matrix(0,size, (size-1))
  
  counter <- 0
  for(i in 2:size) {
    # cat("i" , i)
    globali<<-i
    distMatrix <<-distMatrix
    solution <<- multiroot(f = model2, start = rep(1, (i-1)), verbose=F)
    #cat("solution ", solution$root)
    coordMatrix[i,] <<- c(solution$root,rep(0,length(coordMatrix[i,]) - length(solution$root)))
    
  }
  
}

model2 <- function(x,..) {
  # cat("globali: ", globali)
  i <- globali
  returnvector <- c()
  distMatrix
  for(j in 1:(i-1)) {
    sumResult <- 0
    for(k in 1:(i-1)) {
      #   cat("i: " , i, " j: ",  j , " k: ", k, "\n");
      sumResult <- sumResult + (x[k] - coordMatrix[j,k])^2
      
      #  cat("add1 ",  k ," add2: ", coordMatrix[j,k], "\n")
      
      
      
    }
    #  cat("distance: ", distMatrix[j,i], "\n")
    sumResult <- sumResult - distMatrix[j,i]
    #   cat("append!", "\n")
    returnvector <- append(returnvector, sumResult)
  }
  returnvector
}

createDistMatrix <- function(coords) {
  dimension <- dim(coords)
  dist <- matrix(0,nrow= dimension, ncol= dimension)
  
  for(i in 1:dimension) {
    for(j in i:dimension) {
      dist[i,j] <-sqrt(dist2(coords[i,], coords[j,]))
      dist[j,i] <-sqrt(dist2(coords[i,], coords[j,]))
    }
  }
  dist
}

dist2 <- function(a,b) {
  sum((a-b)^2)
}

vectorlength <- function(size) {
  sum <- 0
  for(i in 1:size) {
    sum <- sum + i
  }
  sum
}

#testpoint:

#die Koordinaten der Punkte

kMeansOnDistances <- function(corM,k=5) {
  distmatrix <- getDist(corM,F)
  distMatrix <<- getDist(corM,F)
  
  size <- length(distmatrix[1,])
  
  #print("start kmeans!");
  startvector <- rep(0.1, vectorlength(size-1))
  
  callModel2(distmatrix)
  #  print(coordMatrix)
  clustering <- kmeans(coordMatrix, centers=k,nstart=300)
  kmeans <- clustering$cluster
  names(kmeans) <- rownames(corM)
  kmeans
}



cmdsolve <- function(corM,k=5, dim=0) {
  d <- getDist(corM,T)
  dim.max <-  dim(corM)[1] - 1
  if(dim == 0) {
    dim <- dim.max
  }
 # dcomp <- getDist(corM,F)
  fit <- cmdscale(d,eig=TRUE, k=dim) # k is the number of dim
  fit # view results
  clustering <- kmeans(fit$points, centers=k,nstart=300)
  kmeans <- clustering$cluster
  
  distmatrix <- matrix(0, nrow=k,ncol=dim(fit$points)[1])
  
  for(i in 1:k) {
    for(j in 1:dim(fit$points)[1]) {
      distance <- sqrt(sum((clustering$centers[i,] - fit$points[j,])^2))
      distmatrix[i,j] <- distance
    }
  }
  
  distpoints <- dist(fit$points, upper=T, diag=T)
  corPoints <- as.matrix(1 - 2*distpoints^2)

  centerCors <- as.matrix(1 - 2*distmatrix^2)
  centerCors <- t(centerCors)
  varimax(centerCors)
  
  rownames(centerCors) <-  rownames(fit$points)
  colnames(distmatrix) <- rownames(fit$points)
  distmatrix <- t(distmatrix)
  
  names(kmeans) <- rownames(corM)
  kmeans
}


cmdsolve.loading <- function(corM,k=5, dim=0) {
  d <- getDist(corM,T)
  
  sign <- sign(as.matrix(d)*corM)
  diag(sign) <- 1
  dim.max <-  dim(corM)[1] - 1
  if(dim == 0) {
    dim <- dim.max
  }
  # dcomp <- getDist(corM,F)
  fit <- cmdscale(d,eig=TRUE, k=dim) # k is the number of dim
  fit # view results
  clustering <- kmeans(fit$points, centers=k,nstart=300)
  kmeans <- clustering$cluster
  
  distmatrix <- matrix(0, nrow=k,ncol=dim(fit$points)[1])
  
  for(i in 1:k) {
    for(j in 1:dim(fit$points)[1]) {
      distance <- sqrt(sum((clustering$centers[i,] - fit$points[j,])^2))
      distmatrix[i,j] <- distance
    }
  }
  
  distpoints <- dist(fit$points, upper=T, diag=T)
  corPoints <- as.matrix(1 - 2*distpoints^2)
  
  centerCors <- as.matrix(1 - 2*distmatrix^2)
  centerCors <- t(centerCors)
  varimax(centerCors)
  
  rownames(centerCors) <-  rownames(fit$points)
  colnames(distmatrix) <- rownames(fit$points)
  distmatrix <- t(distmatrix)
  
  names(kmeans) <- rownames(corM)
  centerCors <- centerCors * sign
  centerCors
}


kmeansCor <- function(corM,k=5) {
  
  kmeans(corM, centers=k,nstart=300)$cluster
  
}


kmeansCorCor <- function(corM,k=5) {
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  kmeans(corCorM, centers=k,nstart=300)$cluster
  
}



cmdsolveCor <- function(corM,k=5, dim=39) {
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  d <- getDist( corCorM,T)
  dcomp <- getDist( corCorM,F)
  fit <- cmdscale(d,eig=TRUE, k=dim) # k is the number of dim
  fit # view results
  clustering <- kmeans(fit$points, centers=k,nstart=300)
  kmeans <- clustering$cluster
  names(kmeans) <- rownames(corM)
  kmeans
}

kMeansOnDistancesCor <- function(corM,k=5) {
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  distmatrix <- getDist(corCorM,F)
  distMatrix <<- getDist(corCorM,F)
  
  size <- length(distmatrix[1,])
  
  #print("start kmeans!");
  startvector <- rep(0.1, vectorlength(size-1))
  
  callModel2(distmatrix)
  #  print(coordMatrix)
  clustering <- kmeans(coordMatrix, centers=k,nstart=300)
  kmeans <- clustering$cluster
  names(kmeans) <- rownames(corM)
  kmeans
}


varClust <- function(cor.sp,k) {

  clust2 <- kmeansvar(cor.sp, init = k, nstart=300)
  clust2$cluster
}



varClust2 <- function(cor.sp,k) {
  tree <- hclustvar(cor.sp)
  
  P3<-cutreevar(tree,k)
  
  P3$cluster
}
