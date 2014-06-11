#nur einmal ausführen
source("clusterfunctions.R")

#install.packages("plsgenomics")
#install.packages("ppcor")
#install.packages("corpcor")


library(ClustOfVar)
library(plsgenomics)
library(ppcor)
library(corpcor)
library(psych)


getNewVarsFromLoadings <- function(x,loadings) {
  partialLoads <- factor.scores(x, loadings, method="Thurstone")$scores
  partialLoads
}



getVar <- function(method,x,cor, k) {

  
  if(method=="kmeansmds") {
    loads.kmeansmds <- cmdsolve.loading(cor,k=k) 
    vars.kmeansmds <- getNewVarsFromLoadings(x, loads.kmeansmds)
  } else if(method=="kmeanscor") {
    loads.kmeanscor <- kmeansCor.loading(cor, k=k)
    vars.kmeanscor <- getNewVarsFromLoadings(x, loads.kmeanscor)
  } else if(method=="varClust") {
    loads.varClust <- varClust.loadings(cor, k=k)
    vars.varClust <- getNewVarsFromLoadings(x, loads.varClust)
  } else if(method=="pca") {
    loads.pca <- pca.loadings(cor, k= k)
    vars.pca <- getNewVarsFromLoadings(x, loads.pca)
  }
}




colon <- data(Colon)

x <- Colon$X

y <- Colon$Y

cor <- cor(x)


nvars <- c(5,25,50,100,500)
result.matrix <- matrix(nrow=5, ncol=length(nvars), 0)

colnames(result.matrix) <- nvars
rownames(result.matrix) <- c("full", "pca", "kmeansmds", "kmeanscor", "varClust")
  

vars.kmeansmds <- getVar("kmeansmds",x,cor,k=nvar) 
vars.kmeanscor <- getVar("kmeanscor",x,cor,k=nvar) 
vars.varClust <- getVar("varClust",x,cor,k=nvar) 
vars.pca <- getVar("pca",x,cor,k=nvar) 

data.all <- as.data.frame(cbind(y = Colon$Y, Colon$X))
data.kmeansmds <- as.data.frame(cbind(y,vars.kmeansmds))


ld1 <- lda(y ~ ., data=data.kmeansmds)
ld2 <- lda(y ~ ., data=data.all)

#wichtige Befehle zum Simulieren:
#ld <- lda(train.y ~ ., data=train.data)

#pred <- predict( ld, test.data)$class

#dazu sinnvoll Kreuzvalidierung programmieren und gute GRafiken erstellen
