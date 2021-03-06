library(gridExtra)

source("faktorensimulation.R")
source("cmdsolve.R")
source("Vergleichsverfahren.R")


#bestimmt für Stichproben der Größe nobs die Anzahl der Cluster
#simuliert das ganze nrep map
numcluadvanced.simulation <- function (cor,nrep,type="kmeans", n, data=NULL) {
  
  v <- c()
  for (i in 1:nrep){

    ####überführe in das Koordinatensystem
    dim <- dim(cor)[1] - 1
    dist <- getDist(cor, F)
    fit <- cmdscale(d=dist,eig=TRUE, k=dim) # k is the number of dim
    points <- fit$points

    number.cluster  <<- getClusterNumbers(points=points, cor.sp=cor, type=type, n=n, data=data)

    
    print( number.cluster)
    v[[i]] <- number.cluster 
  }
  v
}



drawNumberClusterAdvanced.simulation <- function(cor,nrep, type="kmeans",n, whole.number, data) {
  
  whole.cluster.number <- c()

    whole.cluster.number <- c( whole.number, whole.number, whole.number,
                               whole.number, whole.number, whole.number, whole.number)
  
  if(type %in% c("efa", "faclust")) {
    method.names <- method.names.EFA
  } else {
    method.names <- method.names.normal
  }
  
  result.names <- c("whole", "Bias")  

  
  
  resultsmatrix <- matrix(nrow=length(result.names), ncol=length(method.names))
  rownames(resultsmatrix) <- result.names
  colnames(resultsmatrix) <- method.names
  
  for(o in 1:length(nobs)) {
    vector <- as.vector(numcluadvanced.simulation(cor, nrep=nrep, type=type, n=n, data))
    
    m <- matrix( ncol=length(vector), nrow=length(vector[[1]]))
    
    for(i in 1:length(vector)) {
      for(j in 1:length(vector[[1]])) {
        
        v <- vector[[i]][j] 
        
        if(class(v)=="list") {
          v <- unlist(v)[1]
        }
        whole.number <-  whole.cluster.number[j]
        
        if(class(whole.number) == "list") {
          whole.number <- unlist(whole.number)
        }
        
        m[j,i] <- v - whole.number
      }
    }
    
    
    
    # par(mfrow=c(3,3))
    ###aufteilen auf vektoren der einzelnen Methoden, die dann geplottet werden
    for(i in 1:dim(m)[1]) {
      #    drawBarplot(m[i,],ylab=paste("Faktorenanalyse"),nob=nob,type=type, cex.lab=1.5, method= measNames(result)[i])
      #method= measNames(result)[i]
      method.var <- var(m[i,])
      method.bias <- mean(m[i,])
      method.whole <-  whole.cluster.number[i]
      
      resultsmatrix[1,i] <- method.whole
      resultsmatrix[2,i] <- method.bias
  
    }
    
  }
  
  resultsmatrix
  
}

method.names <<- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")

getClusterNumberBias.simulation <- function(NL.mus,  Kor.mus,  type="kmeans") {
  
  r.names <- c()
  descriptions <- ""
  for(i in 1:length(NL.mus)) {
    r.names[i] <- paste0("constr ", i)
    descriptions <- paste0(descriptions,  " und " , r.names[i] , " mit NL von ",
                           NL.mus[i], " und Faktorkorrelation von ", Kor.mus[i], " \n  ")
  }
  
 
  
  rs <- matrix(nrow=length(NL.mus), ncol=length(method.names))
  rownames(rs) <- r.names
  colnames(rs) <- method.names
  for(i in 1:length(NL.mus)) {
    corM1 <- setCorrelationMatrix(NL.mus[i],0,Kor.mus[i],0, F)
    r1 <- drawNumberClusterAdvanced.simulation(corM1,1,type="kmeans")
    rs[i,] <- r1[2,]
  }
  paintTable(rs, "Clusteranzahlsabweichung bei EFA-Similation", 
             paste0( "type=", type, " \n ", descriptions)) 
}





getClusterNumberBias.simulation.methods.original <- function(method, fa.ges) {
  
  r.names <- c()
  descriptions <- ""
  
  
  descriptions<- "bei allen NL gleich und entsprechen Kommunalität (NL.equal)"
  
  if(method==2) {
    descriptions<- "bei einer NL und entsprechen Kommunalität (NL.one)"
  } else if(method==3) {
    descriptions<- "bei zwei NL und entsprechen Kommunalität (NL.two)"
  }
  
  
  
  loads <- NL.equal(fa.ges$loadings)
  
  # loads <- NL.fixed(fa.ges$loadings, 0.2)
  
  if(method==2) {
    loads <- NL.one(fa.ges$loadings)
  } else if(method == 3) {
    loads <- NL.two(fa.ges$loadings)
  }
  
  zuordnung.ges <- apply(loads,1,function(x) which.max(abs(x)))
  
  Phi <- fa.ges$Phi
  set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  corM <- sim.structure(fx=loads,Phi=Phi,n=0)$model
  
  types <- c("kmeans", "average", "complete")
  
  rs <- matrix(nrow=length(types), ncol=length(method.names))
  rownames(rs) <- types
  colnames(rs) <- method.names
  
  
  for(i in 1:length(types)) {
    r <<- drawNumberClusterAdvanced.simulation(corM,1,types[i])
    rs[i,] <- r[2,]
  }
  
  paintTable(rs, "Clusteranzahlsabweichung bei EFA-Similation mit 5 Faktoren",
             paste0(" \n ", descriptions))
}

getClusterNumberBias.simulation.methods <- function(types, methods, fa.ges) {
 
  r.names <- c()
  descriptions <- ""
 
  rs <- matrix(nrow=(length(types)-1)*length(method.names) + length(method.names.EFA)+100, ncol=5)
  colnames(rs) <- c("clustermethod", "clusternumber" , "Sim1", "Sim2", "Sim3")

  colnames.rs <- c()
  
  for(m  in 1:length(methods)) {
    
  method <- methods[m]
  descriptions<- paste0("Sim1 mit allen NL gleich und entsprechen Kommunalität (NL.equal)", "\n", 
                        "Sim 2 bei einer NL und entsprechen Kommunalität (NL.one)", "\n",
                    "Sim3 bei zwei NL und entsprechen Kommunalität (NL.two)")
  
  loads <- NL.equal(fa.ges$loadings)
  
  # loads <- NL.fixed(fa.ges$loadings, 0.2)
  
  if(method==2) {
    loads <- NL.one(fa.ges$loadings)
  } else if(method == 3) {
    loads <- NL.two(fa.ges$loadings)
  }
  
  
  zuordnung.ges <- apply(loads,1,function(x) which.max(abs(x)))
  
  Phi <- fa.ges$Phi
    
   # for(z in 1:nrep) {
  set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
  sim <- sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=nobs, raw=F, items=T, cat=5)
  corM <- sim$r
  data <- sim$observed
  
  
 

  
  for(i in 1:(length(types))) {
  r <- drawNumberClusterAdvanced.simulation(cor,1,types[i], n=100, whole.number= dim(fa.ges$loadings)[2], data=data)

  cat("r: ", r)
  method.names.size <- 3
  
  if(types[i]=="faclust") {
    method.names.size <- length(method.names.EFA)
  }
  
  for(mn in 1:method.names.size) {
    rs[(i-1)*length(method.names) + mn,2+m] <- r[2,mn]
    method.name <- method.names[mn]
    if(mn == 1) {
    rs[(i-1)*length(method.names) + mn,1] <- types[i]
    } else {
      rs[(i-1)*length(method.names) + mn,1] <- ""
    }
    if(types[i] == "faclust") {
      rs[(i-1)*length(method.names) + mn,2] <- method.names.EFA[mn]
    } else {
    rs[(i-1)*length(method.names) + mn,2] <- method.names[mn]
  }
  }
  }

 # }
  }
  
  rs <- t(rs)
  names <-  rs[2,] %in% c("APN","AD","ADM","FOM",NA) 
  numbers <-  rs[5,] %in% c(NA) 
  
  clustertype <- rep(types,each=3)
  clustertype <- c(clustertype,"faclust")
  
  clusternumber <- rs[2,][!names]
  clusternumber <- rep(c("Connectivity","Dunn","Silhouette"),6)
  clusternumber <- c(clusternumber, "MAP", "Parallel-ncomp", "Parallel-nfact", "AIC")
  Sim1 <- rs[3,][!numbers]
  Sim2 <- rs[4,][!numbers]
  Sim3 <- rs[5,][!numbers]
  
  result <- cbind(clustertype, clusternumber, Sim1,Sim2,Sim3)
  
  result
  
  
  paintTable(result, "Clusteranzahlsabweichung bei EFA-Similation mit 5 Faktoren", 
             paste0(" \n ", descriptions)) 
}

#NL.mus <- c(0,0.1,0.2,0.1)
#Kor.mus <- c(0,0,0,0.4)

#getClusterNumberBias.simulation(NL.mus,Kor.mus,type="kmeans")

getClusterNumberBias.simulation.methods.samples <- function(types, methods, fa.ges, nobs=200, nrep=1000) {
  
  number.matrix1 <- matrix(nrow=nrep, ncol=46)
  number.matrix2 <- matrix(nrow=nrep, ncol=46)
  number.matrix3 <- matrix(nrow=nrep, ncol=46)
  number.matrix4 <- matrix(nrow=nrep, ncol=46)
  number.matrix5 <- matrix(nrow=nrep, ncol=46)
  
  items <- T

  
  for(c in 1:nrep) {
  
    print("-----------------")
    print(c)
    print("-----------------")
    
  r.names <- c()
  descriptions <- ""
  
  rs <- matrix(nrow=46, ncol=7)
  colnames(rs) <- c("clustermethod", "clusternumber" , "Sim1", "Sim2", "Sim3", "Sim4", "Sim5")
  
  colnames.rs <- c()
  
  for(m  in 1:length(methods)) {
    
    method <- methods[m]
    descriptions<- paste0("Sim1 mit allen NL gleich und entsprechen Kommunalität (NL.equal)", "\n", 
                          "Sim 2 bei einer NL und entsprechen Kommunalität (NL.one)", "\n",
                          "Sim3 bei zwei NL und entsprechen Kommunalität (NL.two)")
    
    loads <- NL.equal(fa.ges$loadings)
    
    # loads <- NL.fixed(fa.ges$loadings, 0.2)
    
    if(method==2) {
      loads <- NL.one(fa.ges$loadings)
    } else if(method == 3) {
      loads <- NL.two(fa.ges$loadings)
    } else if(method==4) {
      loads <- NL.original(fa.ges$loadings)
    } else if(method==5) {
      loads <- NL.original(fa.ges$loadings)
    }
    
    zuordnung.ges <- apply(loads,1,function(x) which.max(abs(x)))
    
    Phi <- fa.ges$Phi
    
    # for(z in 1:nrep) {
    set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
    if(method==5) {
      sim <-  sim.structure.stella(fx=loads,Phi=Phi, uniq=fa.ges$residual, n=nobs,  items=T, cat=5)
    } else {
    sim <-  sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=nobs,  items=T, cat=5)
    }
    corM <- sim$r
    data <-  sim$observed
    
    cors <- corM
    if(method==2) {
    print("##############################")
    print("##############################")
    print(cors[1,2])
    cat("method ", m)
    print("###############################")
    print("##############################")
    }
    
    
    for(i in 1:(length(types))) {
    
      r <- drawNumberClusterAdvanced.simulation(cor=cors,1,types[i], n=nobs, whole.number=max(zuordnung.ges), data=data)
      if(types[i] == "kmeans") {
      print(types[i])
      print( r)
      }
      method.names.size <- 3
      method.names <- c("Connectivity","Dunn","Silhouette")
      
      if(types[i]=="faclust") {
        method.names.size <- length(method.names.EFA)
        method.names <- c("MAP","Paralell-ncomp","Paralell-nfact", "AIC")
      }
      
      for(mn in 1:method.names.size) {
        rs[(i-1)*length(method.names) + mn,2+m] <- r[2,mn]
        method.name <- method.names[mn]
        if(mn == 1) {
          rs[(i-1)*length(method.names) + mn,1] <- types[i]
        } else {
          rs[(i-1)*length(method.names) + mn,1] <- ""
        }
        if(types[i] == "faclust") {
          rs[(i-1)*length(method.names) + mn,2] <- method.names.EFA[mn]
        } else {
          rs[(i-1)*length(method.names) + mn,2] <- method.names[mn]
        }
      }
    }
    
    # }
  }
 
    
  number.matrix1[c,] <- as.numeric(rs[,3])
  number.matrix2[c,] <- as.numeric(rs[,4])
  number.matrix3[c,] <- as.numeric(rs[,5])
  number.matrix4[c,] <- as.numeric(rs[,6])
  number.matrix5[c,] <- as.numeric(rs[,7])
  
  rownames(rs) <- colnames.rs
  
  
  }
  
  rs.clustermethod <- na.omit(rs[,1])
  rs.clusternumber<- na.omit(rs[,2])
  
  
  mean1 <-   apply(number.matrix1, FUN=mean, MARGIN=2)
  mean2 <-   apply(number.matrix2, FUN=mean, MARGIN=2)
  mean3 <-   apply(number.matrix3, FUN=mean, MARGIN=2)
  mean4 <-   apply(number.matrix4, FUN=mean, MARGIN=2)
  mean5 <-   apply(number.matrix5, FUN=mean, MARGIN=2)
  
  mean1 <- round(na.omit(mean1),2)
  mean2 <- round(na.omit(mean2),2)
  mean3 <- round(na.omit(mean3),2)
  mean4 <- round(na.omit(mean4),2)
  mean5 <- round(na.omit(mean5),2)
  
  whole <- rep(dim(fa.ges$loadings)[2],  length(mean1))
  
  
  var1 <-   round(sqrt(na.omit(apply(number.matrix1, FUN=var, MARGIN=2))),2)
  var2 <-   round(sqrt(na.omit(apply(number.matrix2, FUN=var, MARGIN=2))),2)
  var3 <-   round(sqrt(na.omit(apply(number.matrix3, FUN=var, MARGIN=2))),2)
  var4 <-   round(sqrt(na.omit(apply(number.matrix4, FUN=var, MARGIN=2))),2)
  var5 <-   round(sqrt(na.omit(apply(number.matrix5, FUN=var, MARGIN=2))),2)
  

  rs.clustermethod <- c(rep(types, each=3),"faclust")
  rs.clusternumber <- c(rep(c("Connectivity", "Dunn","Silhouette"),6), "MAP","Paralell-mcomp","Paralell-nfact","AIC")
  
  mean.matrix <- (cbind(rs.clustermethod, rs.clusternumber,whole, mean1,mean2,mean3,mean4, mean5))
  
  var.matrix <- (cbind(rs.clustermethod, rs.clusternumber, var1,var2,var3,var4, var5))
  

  
  
  paintTable(mean.matrix, "durchsch. Clusteranzahlsabweichung bei EFA-Similation mit 5 Faktoren", 
             paste0(" \n ", descriptions)) 
  
  paintTable(var.matrix, "Std.-abweichung der Clusteranzahlsabweichung bei EFA-Similation mit 5 Faktoren", 
             paste0(" \n ", descriptions, " nrep ", nrep, " nobs ", nobs)) 
}

