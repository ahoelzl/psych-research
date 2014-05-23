
#berechnet die Clusterung für den gesamten Datensatz für alle Clusterarten mit denen die 
#Cluster der Stichproben dann verglichen werden.

#numbercluster1 <- 13
#numbercluster2 <- 4
#numbercluster3 <- 11
#numbercluster4 <- 8
#numberclusterfactor <- 11



#die Korrelation und Korrelation der Korrelation
#cor <- cor(facs, use="pairwise.complete.obs", method="pearson")
#corcor <- cor(cor, use="pairwise.complete.obs", method="pearson")

#
#berechnet alle Clusterungen für den Gesamtdatensatz!
#

getDist <- function(corM, asdist) {
  if(asdist) {
  as.dist(sqrt(0.5-0.5*abs(corM)), upper=T, diag=T)
  } else {
    sqrt(0.5-0.5*abs(corM))
  }
}

#getDist <- function(corM, asdist) {
#  if(asdist) {
#    as.dist(1-corM)
#  } else {
#    1-corM
#  }
#}


getDistNoMetric <- function(corM, asdist) {
  if(asdist) {
    as.dist(1-corM)
  } else {
    1-corM
  }
}



completeCor <- function(corM,k=5,show=F) {
#Methode complete und Korrelation
d.ges1 <- getDist(corM,T)
hca.ges1 <- hclust(d.ges1, method="complete", members=NULL)
completecor5 <- cutree(hca.ges1,k)
if(show) {
plot(hca.ges1, hang=-1, main="Dendrogramm mit Complete und Korrelation ",sub=optionstext) # graphisch am sinnvollsten w?ren hier 4 Cluster, die bei einer minimalen Korrelation von 0 entstehen (also h=1)
rect.hclust(hca.ges1, k, border="red")
}
completecor5
}


completeCorCor <- function(corM,k=5,show=F) {
#Methode complete und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
d.ges3 <- getDist(corCorM,T)
hca.ges3 <- hclust(d.ges3, method="complete", members=NULL)
completecorcor5 <- cutree(hca.ges3,k)
if(show) {
plot(hca.ges3, hang=-1, main="Dendrogramm mit Complete und Korrelation der Korrelation",sub=optionstext)
rect.hclust(hca.ges3, k, border="red")
}
completecorcor5
}

averageCor <- function(corM,k=5,show=F) {
#Methode average und Korrelation
  globalcor <<-corM
d.ges5 <- getDist(corM, T)
hca.ges5 <- hclust(d.ges5, method="average", members=NULL)
averagecor5 <- cutree(hca.ges5,k)
if(show) {
plot(hca.ges5, hang=-1, main="Dendrogramm mit Average und Korrelation ", sub=optionstext)
rect.hclust(hca.ges5, k, border="red")
}

averagecor5
}

averageCorCor <- function(corM,k=5,show=F) {
#Methode average und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
d.ges7 <- getDist(corCorM, T)
hca.ges7 <- hclust(d.ges7, method="average", members=NULL)
averagecorcor5 <- cutree(hca.ges7,k)
if(show) {
plot(hca.ges7, hang=-1, main="Dendrogramm mit Average und Korrelation der Korrelation",sub=optionstext)
rect.hclust(hca.ges7, k, border="red")
}
averagecorcor5
}

completeCorNoMetric <- function(corM,k=5,show=F) {
  #Methode complete und Korrelation
  d.ges1 <- getDistNoMetric(corM,T)
  hca.ges1 <- hclust(d.ges1, method="complete", members=NULL)
  completecor5 <- cutree(hca.ges1,k)
  if(show) {
    plot(hca.ges1, hang=-1, main="Dendrogramm mit Complete und Korrelation ",sub=optionstext) # graphisch am sinnvollsten w?ren hier 4 Cluster, die bei einer minimalen Korrelation von 0 entstehen (also h=1)
    rect.hclust(hca.ges1, k, border="red")
  }
  completecor5
}

completeCorCorNoMetric <- function(corM,k=5,show=F) {
  #Methode complete und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  d.ges3 <- getDistNoMetric(corCorM,T)
  hca.ges3 <- hclust(d.ges3, method="complete", members=NULL)
  completecorcor5 <- cutree(hca.ges3,k)
  if(show) {
    plot(hca.ges3, hang=-1, main="Dendrogramm mit Complete und Korrelation der Korrelation",sub=optionstext)
    rect.hclust(hca.ges3, k, border="red")
  }
  completecorcor5
}

averageCorNoMetric <- function(corM,k=5,show=F) {
  #Methode average und Korrelation
  d.ges5 <- getDistNoMetric(corM, T)
  hca.ges5 <- hclust(d.ges5, method="average", members=NULL)
  averagecor5 <- cutree(hca.ges5,k)
  if(show) {
    plot(hca.ges5, hang=-1, main="Dendrogramm mit Average und Korrelation ", sub=optionstext)
    rect.hclust(hca.ges5, k, border="red")
  }
  
  averagecor5
}
averageCorCorNoMetric <- function(corM,k=5,show=F) {
  #Methode average und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  d.ges7 <- getDistNoMetric(corCorM, T)
  hca.ges7 <- hclust(d.ges7, method="average", members=NULL)
  averagecorcor5 <- cutree(hca.ges7,k)
  if(show) {
    plot(hca.ges7, hang=-1, main="Dendrogramm mit Average und Korrelation der Korrelation",sub=optionstext)
    rect.hclust(hca.ges7, k, border="red")
  }
  averagecorcor5
}




completeCorCorCor <- function(corM,k=5,show=F) {
  #Methode complete und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  corCorM  <-  cor(corCorM, use="pairwise.complete.obs", method="pearson")
  d.ges3 <- getDist(corCorM,T)
  hca.ges3 <- hclust(d.ges3, method="complete", members=NULL)
  completecorcor5 <- cutree(hca.ges3,k)
  if(show) {
    plot(hca.ges3, hang=-1, main="Dendrogramm mit Complete und Korrelation der Korrelation",sub=optionstext)
    rect.hclust(hca.ges3, k, border="red")
  }
  completecorcor5
}

averageCorCorCor <- function(corM,k=5,show=F) {
  #Methode average und Korrelation der Korrelation
  corCorM  <-  cor(corM, use="pairwise.complete.obs", method="pearson")
  corCorM  <-  cor(corCorM, use="pairwise.complete.obs", method="pearson")
  d.ges7 <- getDist(corCorM, T)
  hca.ges7 <- hclust(d.ges7, method="average", members=NULL)
  averagecorcor5 <- cutree(hca.ges7,k)
  if(show) {
    plot(hca.ges7, hang=-1, main="Dendrogramm mit Average und Korrelation der Korrelation",sub=optionstext)
    rect.hclust(hca.ges7, k, border="red")
  }
  averagecorcor5
}




fclustering <- function(corM,k=5,show=F) {
  fa.ges6 <- fa(corM, nfactors=k, max.iter=100, fm="ml", rotate="promax", method="pearson")
  zuordnung.ges6 <- apply(fa.ges6$loadings,1,function(x) which.max(abs(x)))
  zuordnung.ges6
}


fclustering.loading <- function(corM,k=5,show=F) {
  fa.ges6 <- fa(corM, nfactors=k, max.iter=100, fm="ml", rotate="promax", method="pearson")
  zuordnung.ges6 <- apply(fa.ges6$loadings,1,function(x) which.max(abs(x)))
  fa.ges6$loadings
}


##Faktoranalyse
#fa.ges <- fa(facs, nfactors=numberclusterfactor, max.iter=100, fm="ml", rotate="promax", method="pearson")
#fa.ges6 <- fa(facs, nfactors=6, max.iter=100, fm="ml", rotate="promax", method="pearson")
#zuordnung.ges <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))
#zuordnung.ges6 <- apply(fa.ges6$loadings,1,function(x) which.max(abs(x)))
#plot(hca.ges7, hang=-1, main="Dendrogramm mit Average und Korrelation der Korrelation")
#abline(h=0.8, col="red")