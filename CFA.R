library(lavaan)
library(clValid)




getClustering <- function(cor.sp, k, method, daten.sp=NULL) {
  cut.sp1 <- c()
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
  } else if(method=="faclust" || method=="efa") {
    cut.sp1 <- fclustering(cor.sp,k)
  } else if(method=="clustofvar") {
    cut.sp1 <- varClust(daten.sp,k)
  } else if(method=="clustofvar2") {
    cut.sp1 <- varClust2(daten.sp,k)
  }
  cut.sp1
}


getLoadingMatrix <- function(cor.sp, k, method, daten.sp=NULL) {
  cut.sp1 <- c()
 
    if(method=="kmeansmds") {
    cut.sp1 <- cmdsolve.loading(cor.sp,k)
    } else if(method=="faclust" || method=="efa") {
    cut.sp1 <- fclustering.loading(cor.sp,k)
    }
  
  cut.sp1
}



getCFASimiliarity <- function(facs, nobs,methods, daten.sp1, daten.sp2,  efa=F) {
  
  results <- matrix(0,nrow=length(method.names.normal)+4, ncol=7)
  #  daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
  cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
  
  efa.cluster.number <-  EFA.Cluster.number(cor.sp1, n=nobs)
  for(h  in 1:length(methods)) {
    method <- methods[h]
    
  method.names.normal <- c( "Connectivity", "Dunn" ,"Silhouette")
  
  measures <- c()
    
    
  #  daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp2 <- cor(as.matrix(daten.sp2), use="pairwise.complete.obs", method="pearson")
    
    
    if(method=="faclust" || method=="efa") {
      number.clusters <-  efa.cluster.number
      names(number.clusters) <- method.names.EFA
    } else {
      number.clusters <- numcluadvanced.whole(daten.sp1, type=method, n=nobs)
      number.clusters <- append(number.clusters,  efa.cluster.number)
      names(number.clusters) <- c(method.names.normal, method.names.EFA)
    }
    
measures <- c()

    for(u in 1:length(number.clusters)) {

      number.cluster <- number.clusters[u]
    #names(number.cluster) <- method.names
      cat("------------------------------", method, " ----------- ", names(number.clusters)[u], "\n")  

      
      
    k <- number.cluster
    
   
    
      if(class(k)=="list") {
        k <- unlist(k)[1]
      }
      
      if(k == 0) {
        k <- 1
      }
    cat("k :", k)
    clustering <- getClustering(cor.sp1, k, method, daten.sp=daten.sp1)
    cat("clustering :", clustering)
    latent.zuweisung <- ''
    
    for(i in 1:k) {
      clustername = as.character(paste0("class",i))
      cluster.name <- names(which(clustering==i))
      if(length(cluster.name) > 0) {
      latent.zuweisung <- paste(latent.zuweisung, clustername, " =~ " )
      
      for(name in cluster.name) {
        latent.zuweisung <- paste(latent.zuweisung, name, "+ ")
      }
      
      #removes the , at the end
      
      latent.zuweisung <- substring(latent.zuweisung, 1, nchar(latent.zuweisung) - 2)
      
      latent.zuweisung <- paste(latent.zuweisung,  " ",  sep='\n')
      }
    }
    
      print(latent.zuweisung)
    
    frame <- as.data.frame(daten.sp2)
    
    fit <- cfa(latent.zuweisung,data = frame, control=list(
      eval.max = 20000L, 
      iter.max = 10000L, 
      trace = 0L, 
      abs.tol = (.Machine$double.eps * 10) ,
      rel.tol = 1e-5, # 1e-10 seems 'too strict' 
      x.tol = 1.5e-8 ,
      step.min = 2.2e-14 ))
  
    cat("beforelogLk")
    
    test <-  try(logLik(fit))
      
    if(  class(test) == "try-error")  {

      
      results <- matrix(0,nrow=length(method.names.normal)+4, ncol=5)
      return(results)
     # meas <- 0
    } else {
      cat("before meas")
    meas <- fitMeasures(fit, c("BIC"))
    cat("meas: ", meas)
    
    cat(latent.zuweisung, "\n")
    cat(" Fit: ", meas,  "\n")
    cat("------------", "\n")
    }
    
      w <- 0
  if(method == "efa") {
    w <- 3
  }
results[w+u,h] <- meas
    
    }
  }

 results
  
  
}




getCFASimiliarity.loading <- function(facs, nobs,methods, daten.sp1, daten.sp2,  efa=F) {
  
  results <- matrix(0,nrow=length(method.names.normal)+4, ncol=7)
  #  daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
  cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
  
  efa.cluster.number <-  EFA.Cluster.number(cor.sp1, n=nobs)
  for(h  in 1:length(methods)) {
    method <- methods[h]
    
    method.names.normal <- c("Connectivity", "Dunn" ,"Silhouette")
    
    measures <- c()
    
    
    
    #  daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp2 <- cor(as.matrix(daten.sp2), use="pairwise.complete.obs", method="pearson")
    
    
    if(method=="faclust" || method=="efa") {
      number.clusters <-  efa.cluster.number
      names(number.clusters) <- method.names.EFA
    } else {
      number.clusters <- numcluadvanced.whole(daten.sp1, type=method, n=nobs)
      number.clusters <- append(number.clusters,  efa.cluster.number)
      names(number.clusters) <- c(method.names.normal, method.names.EFA)
    }
    
    measures <- c()
    
    for(u in 1:length(number.clusters)) {
      
      number.cluster <- number.clusters[u]
      
#names(number.cluster) <- method.names
      cat("------------------------------", method, " ----------- ", names(number.clusters)[u], "\n")  
      
      k <- number.cluster
      
      if(k == 0) {
        k <- 1
      }
      
      cat("k :", k)
      clustering <- getClustering(cor.sp1, k, method, daten.sp=daten.sp1)
      
      loads <- getLoadingMatrix(cor.sp1, k, method, daten.sp=daten.sp1)
      
      setzero <- function(x) {
        bins <- x == max(x)
        x[!bins] <- 0
        x
      }
      
      loads <- t(apply(loads, 1, setzero))
      
      cat("clustering :", clustering)
      latent.zuweisung <- ''
      
      for(i in 1:k) {
        clustername = as.character(paste0("class",i))
        cluster.name <- names(which(clustering==i))
        loadings <- loads
        if(length(cluster.name) > 0) {
          latent.zuweisung <- paste(latent.zuweisung, clustername, " =~ " )
          
          for(r in 1:length(cluster.name)) {
            name <- cluster.name[r]
            load <- sum(loadings[rownames(loadings)==name,])
            latent.zuweisung <- paste(latent.zuweisung,load," * ", name, "+ ")
          }
          
          #removes the , at the end
          
          latent.zuweisung <- substring(latent.zuweisung, 1, nchar(latent.zuweisung) - 2)
          
          latent.zuweisung <- paste(latent.zuweisung,  " ",  sep='\n')
        }
      }
      
      
      frame <- as.data.frame(cor.sp2)
      
      fit <- cfa(latent.zuweisung,data = frame, control=list(rel.tol =  1e-5))
      
      cat("beforelogLk")
      
      test <-  try(logLik(fit))
      
      if(  class(test) == "try-error")  {
        
        
        results <- matrix(0,nrow=length(method.names.normal)+4, ncol=5)
        return(results)
      } else {
        cat("before meas")
        meas <- fitMeasures(fit, c("BIC"))
        cat("meas: ", meas)
        
        cat(latent.zuweisung, "\n")
        cat(" Fit: ", meas,  "\n")
        cat("------------", "\n")
      }
      
      w <- 0
      if(method == "efa") {
        w <- 3
      }
      results[w+u,h] <- meas
      
    }
  }
  
  results
  
  
}





runCFR <- function(nrep, nobs, methods = c("efa", "averagecor","completecor", "kmeansmds", "kmeanscor", "clustofvar", "clustofvar2")) {

results <- c()

results.matrix <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
colnames(results.matrix) <- methods

results.matrix.var <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
colnames(results.matrix.var) <- methods
####simulate the data

samples <- list()
for(i in 1:nrep) {
  samples[[i]] <- sample(x=1:nrow(facs), size=nobs, replace=T)
}

resultlist <- list()

s <- matrix(0,nrow=length(method.names.normal)+4, ncol=7)
for(z in 1:nrep) {
    daten.sp1 <- facs[samples[[z]],]
  #  daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    
    ##das zweite keine Stichprobe sondern alle Daten
    
   daten.sp2 <- facs
    
    save.value <-  getCFASimiliarity(facs, nobs=nobs, methods=methods, efa=efa, daten.sp1 = daten.sp1,
                                daten.sp2 = daten.sp2)

    
    save.value
    rownames(save.value) <- c(method.names.normal, method.names.EFA)
    colnames(save.value) <- methods
    
    if(save.value[4,4]==0) {
      z <- z - 1
      next;
    }
    
    s <- s + save.value/nrep
    
    print(s)
   # resultlist[[z]] <- save.value
    
}

s <- s

avgsum <- sum(s)/(49-3) 

s <- s - avgsum
s[c(1,2,3),1] <- 0


paintTable(s, "Varianz des BIC bei konfirmatorischer CFA", paste0("nrep ", nrep))
#results.m <- t(as.matrix(results, 1)
paintTable(s, "BIC bei konfirmatorischer CFA", paste0("nrep ", nrep, " nobs ", nobs))
}



runCFR.loading <- function(nrep, nobs, methods = c("efa","kmeansmds")) {
  
  
  
  results <- c()
  
  results.matrix <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
 # colnames(results.matrix) <- methods
  
  results.matrix.var <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
 # colnames(results.matrix.var) <- methods
  ####simulate the data
  
  samples <- list()
  for(i in 1:nrep) {
    samples[[i]] <- sample(x=1:nrow(facs), size=nobs, replace=T)
  }
  
  resultlist <- list()
  
  s <- matrix(0,nrow=length(method.names.normal)+4, ncol=7)
  for(z in 1:nrep) {
    daten.sp1 <- facs[samples[[z]],]
    #  daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    
    ##das zweite keine Stichprobe sondern alle Daten
    
    daten.sp2 <- facs
    
    save.value <-  getCFASimiliarity.loading(facs, nobs=nobs, methods=methods, efa=efa, daten.sp1 = daten.sp1,
                                     daten.sp2 = daten.sp2)
    
    
    save.value
    rownames(save.value) <- c(method.names.normal, method.names.EFA)
  #  colnames(save.value) <- methods
    
    if(save.value[2,2]==0) {
      z <- z - 1
      next;
    }
    
    s <- s + save.value
    # resultlist[[z]] <- save.value
    
  }
  
  s <- s/nrep
  
  
  
  paintTable(s, "Varianz des BIC bei konfirmatorischer CFA - loadings", paste0("nrep ", nrep))
  #results.m <- t(as.matrix(results, 1)
  paintTable(s, "BIC bei konfirmatorischer CFA - loadings", paste0("nrep ", nrep, " nobs ", nobs))
}




runCFR_random <- function(nrep, nobs, methods = c("efa" , "averagecor","completecor", "kmeansmds")) {
  
  
  
  results <- c()
  
  results.matrix <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
  colnames(results.matrix) <- methods
  
  results.matrix.var <- matrix(0, ncol=length(methods), nrow=length(method.names.normal)+4)
  colnames(results.matrix.var) <- methods
  
  
  for(i  in 1:length(methods)) {
    results <- matrix(0,ncol=(length(method.names.normal)+4), nrow=nrep)
    for(z in 1:nrep) {
      daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
        daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
      
      ##das zweite keine Stichprobe sondern alle Daten
      
      #daten.sp2 <- facs
      
      save.value <-  getCFASimiliarity(facs, nobs=nobs, method=methods[i], efa=efa, daten.sp1 = daten.sp1,
                                       daten.sp2 = daten.sp2)
      
      if(save.value[5,5] == 0) {
        next
      }
      
      if(length(save.value) == 4) {
        results[z,] <- c(rep(0,length(method.names.normal)), save.value)
      } else {
        results[z,] <-save.value
      }
    }
    
    result <- apply(results, MARGIN=2, FUN=mean)
    result.var <- apply(results, MARGIN=2, FUN=var)
    
    results.matrix[,i] <- result
    
    results.matrix.var[,i] <- result.var
  }
  
  rownames(results.matrix) <- c(method.names.normal, method.names.EFA)
  rownames(results.matrix.var) <- c(method.names.normal, method.names.EFA)
  
  paintTable(results.matrix.var, "Varianz des BIC bei konfirmatorischer CFA_random", paste0("nrep ", nrep))
  #results.m <- t(as.matrix(results, 1)
  paintTable(results.matrix, "BIC bei konfirmatorischer CFA_random", paste0("nrep ", nrep, " nobs ", nobs))
}



output.cor.matrices <- function(nrep = 100, size=c(100,200,500,1000)) {
  
  par(mfrow=c(1,length(size)))
  for(s in size) {
  cor.total.vector <- c()
  
  for(i in 1:nrep) {
  daten.sp1 <- facs[sample(x=1:nrow(facs), size=s, replace=T),]
  
  daten.sp1 <- facs[sample(x=1:nrow(facs), size=s, replace=T),]
  
  cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
  
  cor.total <- cor(as.matrix(daten.sp2), use="pairwise.complete.obs", method="pearson")
  
  cor.total <- cor(as.matrix(facs), use="pairwise.complete.obs", method="pearson")

  changed <- cor.sp1 - cor.total
  
  
  cor.total.vector <- append(cor.total.vector, as.vector(changed))
 
  }
  
  hist(changed, main=paste0("size:", s))
  
}
}




