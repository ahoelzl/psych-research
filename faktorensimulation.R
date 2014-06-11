library(psych)

## es sollen aus der Faktorenanalyse (Nebenladungen, Faktorkorrelationen, Einzigartigkeiten) 
## des Gesamtdatensatzes Korrelationsmatrizen simuliert werden, mit denen dann wieder eine EFA
## gerechnet wird

## dazu erst eine kurze Analyse der Nebenladungen des Hauptdatensatzes:
#facs <- na.omit(facs)
#fa.ges <- fa(facs, nfactors=5, max.iter=100, fm="ml", rotate="promax", method="pearson")
#zuordnung.ges <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))
## Funktion für Histrogramm der Nebenladungen 

NL.hist <- function(x) {
  result=matrix(rep(1,(nrow(x)*(ncol(x)-1))),nrow=nrow(x),ncol=(ncol(x)-1))
  for (i in 1:nrow(x)) {
    result[i,] <- x[i,-which.max(x[i,])]
  }
  hist(result,breaks=20,main="Verteilung der Nebenladungen") # hier kann man verschiedene Dinge berechnen
}

# Zur Analyse der Faktorinterkorrelationen:

Phi.sum <- function(x) {
  result=matrix(rep(1,(nrow(x)*(ncol(x)-1))),nrow=nrow(x),ncol=(ncol(x)-1))
  for (i in 1:nrow(x)) {
    result[i,]  <- x[i,-i]
  }
  max(as.vector(result))
}




# Ergebnis: MW der NL: 0.0068; min der NL: -0.4403; max der NL: 0.4749, sd: 0.164

## ausgehend von diesen Analysen werden folgende Entscheidungen getroffen: 
## NL und Phi werden auf 0.0, 0.01, 0.05 und 0.1 festgesetzt (alle gleich hoch) und in einem zweiten Schritt
## werden sie nicht alle gleich gesetzt, sondern aus einer NV gezogen 


# Ladungsmatrix ohne Nebenladungen 
keineNL <- function(x){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=0
  }
  x
}

# Korrelationsmatrix ohne Korrelationen 
Phi.0 <- function(x){
  for (i in 1:nrow(x)) {
    x[i,-i]=0
  }
  x
}



# Ladungsmatrix mit Nebenladung 0.1
NL.fixed <- function(x,value){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=value
  }
  x
}

# Korrelationsmatrix mit Korrelationen 0.1 
Phi.fixed <- function(x, value){
  for (i in 1:nrow(x)) {
    x[i,-i]=value
  }
  x
}


# Ladungsmatrix mit Nebenladungen aus NV gezogen
NL.NV <- function(x,mean,sd){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=rnorm(ncol(x)-1,mean=mean,sd=sd)
  }
  x
}

# Korrelationsmatrix mit Korrelationen aus NV 
Phi.NV <- function(x,mean, sd){
  for (i in 1:nrow(x)) {
    x[i,-i]=rnorm(ncol(x)-1,mean=mean,sd=sd)
  }
  x
}

# Nebenladungen aus Gleichverteilung (bei dieser Funktion bin ich mir noch nicht sicher, ob sie so funktioniert)

NL.UD <- function(x){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=runif(ncol(x)-1,min=-0.4,max=0.4)
  }
  x
}

# Korrelationen aus Gleichverteilung

Phi.UD <- function(x){
  for (i in 1:nrow(x)) {
    x[i,-i]=runif(ncol(x)-1,min=-0.3,max=0.3)
  }
  x
}

nobs <- 48
# diese Objekte werden für die Simulation gebraucht:


syst.var  <-  function(x) {
  result=matrix(NA,nrow=nrow(x),ncol=(ncol(x)-1))
  for (i in 1:nrow(x)) {
    result[i,] <- x[i,-which.max(x[i,])]^2
  }
  rowSums(abs(result))
}


syst.var  <-  function(x) {
  result=matrix(NA,nrow=nrow(x),ncol=(ncol(x)-1))
  for (i in 1:nrow(x)) {
    result[i,] <- x[i,-which.max(x[i,])]^2
  }
  rowSums(abs(result))
}


NL.equal <- function(x){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=sqrt(syst.var(x)[i]/4)
  }
  x
}

NL.one <- function(x){
  y <- x
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])]=0
  }
  
  set.seed(1000)
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])][sample(dim(x)[2]-1,1,replace=T)]=sqrt(syst.var(x)[i])
  }
  y
}


NL.two <- function(x){
  y <- x
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])]=0
  }
  
  set.seed(1000)
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])][sample(dim(x)[2]-1,2,replace=T)]=sqrt(syst.var(x)[i]/2)
  }
  y
}

# Ladungsmatrix erstellen; m ist Wert, der für NL eingesetzt wird 
NL.fixed <- function(x,m){
  y <- x
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])]=0
  }
  
  set.seed(1000)
  for (i in 1:nrow(y)) {
    y[i,-which.max(y[i,])][sample(1:4,m,replace=T)]=sqrt(syst.var(x)[i]/2)
  }
  y
}

# Ladungsmatrix erstellen; m ist Wert, der für NL eingesetzt wird 
NL <- function(x,m){
  for (i in 1:nrow(x)) {
    x[i,-which.max(x[i,])]=m
  }
  x
}

# Korrelationsmatrix erstellen; m ist Wert, der für Phi eingesetzt wird 
phi <- function(x,m){
  for (i in 1:nrow(x)) {
    x[i,-(1:i)]=m
  }
  x
}


#NL
#1 : nur Hauptladungen
#2: 0.01 Nebenladung
#3: 0.05 Nebenladung
#4: 0.1 Nebenladung
#5: NV-verteilte Nebenladung mit sd=0.1
#6  NV-verteilte Nebenladung mit sd=0.241
#7  NV-verteilte Nebenladung mit sd=0.3
#phi
#1: keine Korrelationen
#2: 0.01 Korrelationen
#3: 0.05 Korrelationen
#4: NV-verteilte Korrelationen
setCorrelationMatrix <- function(NLmean, NLsd, phimean,phisd, constant.NL = F, addError=F) {
  options1 <-""
  options2 <- ""
  
  original.fx <- NL.fixed(fa.ges$loadings, 0.05)
  
  Phi <- Phi.0(fa.ges$loadings)
  if(NLsd == 0) {
    fx <- NL.fixed(fa.ges$loadings,NLmean)
    options1 <- paste(NLmean, " Nebenladung")
  } else {
    fx <- NL.NV(fa.ges$loadings,NLmean, NLsd)
    options1 <- paste("Nebenladung nv mit MW ", NLmean, " und sd ", NLsd)
  }
  
  if(phisd == 0) {
    Phi <- Phi.fixed(fa.ges$Phi, phimean)
    options2 <- paste(phimean, " Korrelation")
  } else {
    Phi <- Phi.NV(fa.ges$Phi, phimean, phisd)
    options1 <- paste("Korrelation nv mit MW ", phimean, " und sd ", phisd)
  }
  
  if(!constant.NL) {
    for (i in 1:nrow(fx)) {
      cat("original: ", fx[i,-which.max(fx[i,])][sample(1:(ncol(fx)-1),1,replace=T)], " now  ", sqrt(syst.var(original.fx)[i]))
    #fx[i,-which.max(fx[i,])][sample(1:(ncol(fx)-1),1,replace=T)]<-sqrt(syst.var(original.fx)[i])
    }
  }
  
  #fx <- fa.ges$loadings
  
  optionstext <<-paste(options1,options2)
  
#  print(fx)
  
  #print(Phi)
  # mit dem sim.structure-Befehl werden dann die Korrelationsmatrizen simuliert
  #sp <- sim.structure(fx=fx,Phi=Phi,n=nobs,uniq=fa.ges$uniquenesses)  
  unique <- fa.ges$residual
  #sp <<- sim.structure.stella(fx=fx,Phi=Phi,uniq=unique,n=0)  
  if(addError) {
  sp <<- sim.structure.stella(fx=fx,Phi=Phi,uniq=unique,n=0) 
  corM <- cov2cor(sp$model)
  #  sp <<- sim.structure.stella(fx=fx,Phi=Phi,n=0)  
  #  corM <- cor(sp$model)
  # corM <- fa.ges$residual
  } else {
  sp <<- sim.structure(fx=fx,Phi=Phi,n=0)  
  corM <- sp$model
  }

  
  
  #bei Korrelation:
  
  #sp <<- sim.structure.stella(fx=fx,Phi=Phi,n=0)  
  #corM <- cor(sp$model)
  corM
}
