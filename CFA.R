library(lavaan)
library(clValid)




getClustering <- function(cor.sp, k, method) {
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
  } else if(method=="faclust") {
    cut.sp1 <- fclustering(cor.sp,k)
  }
  cut.sp1
}



getCFASimiliarity <- function(facs, nobs,method, daten.sp1, daten.sp2,  efa=F) {
  
  
  method.names.normal <- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")
  
  measures <- c()

  #  daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
    
  #  daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    cor.sp2 <- cor(as.matrix(daten.sp2), use="pairwise.complete.obs", method="pearson")

    

    number.clusters <- numcluadvanced.whole(daten.sp1, type=method)

measures <- c()
  
    for(number.cluster in number.clusters) {

    #names(number.cluster) <- method.names
    cat("------------------------------", "\n")    

    k <- number.cluster
    
    if(k == 0) {
      k <- 1
    }
    
    cat("k :", k)
    clustering <- getClustering(cor.sp1, k, method)
    cat("clustering :", clustering)
    latent.zuweisung <- ''
    
    for(i in 1:k) {
      clustername = as.character(paste0("class",i))
      cluster.name <- names(which(clustering==i))
      latent.zuweisung <- paste(latent.zuweisung, clustername, " =~ " )
      
      for(name in cluster.name) {
        latent.zuweisung <- paste(latent.zuweisung, name, "+ ")
      }
      
      #removes the , at the end
      
      latent.zuweisung <- substring(latent.zuweisung, 1, nchar(latent.zuweisung) - 2)
      
      latent.zuweisung <- paste(latent.zuweisung,  " ",  sep='\n')
      
    }
    
    cat("latent.zuweiung: ", latent.zuweisung)
    
    frame <- as.data.frame(cor.sp2)
    
    fit <- cfa(latent.zuweisung,data = frame)
  
    cat("beforelogLk")
    
    test <-  try(logLik(fit))
      
    if(  class(test) == "try-error")  {

      print("error thrown")
      meas <- -1
    } else {
      cat("before meas")
    meas <- fitMeasures(fit, c("BIC"))
    cat("meas: ", meas)
    
    cat(latent.zuweisung, "\n")
    cat(" Fit: ", meas,  "\n")
    cat("------------", "\n")
    }
    
  
measures <- append(measures, meas)
    
    }
  
  names(measures) <- names(number.clusters)
  
  measures
  
  
}

runCFR <- function(nrep, nobs, efa=F) {

methods <- c("averagecor","completecor","averagecorcor", "completecorcor", "kmeansmds")
clusternumber.names <- c("APN", "AD" ,"ADM", "FOM","Connectivity", "Dunn"  ,       "Silhouette")  

if(efa) {
  methods <- c("faclust")
  clusternumber.names <- c("MAP", "Paralell-mcomp", "Paralell-nfact", "AIC")
}
results <- c()

results.matrix <- matrix(0, ncol=length(methods), nrow=length(clusternumber.names))
colnames(results.matrix) <- methods
rownames(results.matrix) <- clusternumber.names


for(i  in 1:length(methods)) {
result <- 0
for(z in 1:nrep) {
    daten.sp1 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    ##daten.sp2 <- facs[sample(x=1:nrow(facs), size=nobs, replace=T),]
    
    ##das zweite keine Stichprobe sondern alle Daten
    
    daten.sp2 <- facs
    
    result <- result + getCFASimiliarity(facs, nobs=nobs, method=methods[i], efa=efa, daten.sp1 = daten.sp1,
                                               daten.sp2 = daten.sp2)
}
result <- result/nrep
    results.matrix[,i] <- result
  }





#results.m <- t(as.matrix(results, 1)
paintTable(results.matrix, "BIC bei konfirmatorischer CFA", paste0("nrep ", nrep))
}



output.cor.matrices <- function(nrep = 100, size=c(100,200,500,1000)) {
  
  par(mfrow=c(1,length(size)))
  for(s in size) {
  cor.total.vector <- c()
  
  for(i in 1:nrep) {
  daten.sp1 <- facs[sample(x=1:nrow(facs), size=s, replace=T),]
  
  cor.sp1 <- cor(as.matrix(daten.sp1), use="pairwise.complete.obs", method="pearson")
  
  cor.total <- cor(as.matrix(facs), use="pairwise.complete.obs", method="pearson")

  changed <- cor.sp1 - cor.total
  
  
  cor.total.vector <- append(cor.total.vector, as.vector(changed))
 
  }
  
  hist(changed, main=paste0("size:", s))
  
}
  
  
}




