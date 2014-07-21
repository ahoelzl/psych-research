source("Daten_einlesen.R")

source("cmdsolve.R")

cor.sp <- corM

k <- 5

faclustering <- fclustering(cor.sp,k)

completecor <- completeCor(cor.sp,k)

averagecor <-  averageCor(cor.sp,k)


kmeansmds <- cmdsolve(cor.sp,k)

kmeansCor <- kmeansCor(cor.sp,k)

  varClust1 <- varClust(facs,k)

varClust2  <-  varClust2(facs,k)

vc1 <- varClust1
vc2 <- varClust2
listOfClusters <- list("faclustering" = faclustering, "completecor" = completecor, "averagecor"=averagecor,
                       "kmeansmds"=kmeansmds, "kmeansCor"=kmeansCor, "varClust1"=vc1, "varClust2"=vc2)

resultMatrix <- matrix(0,nrow=7,ncol=7)

for(i in 1:length(listOfClusters)) {
  for(j in 1:length(listOfClusters)) {

result <- vergleich(listOfClusters[[i]],listOfClusters[[j]], 1)[1] +  vergleich(listOfClusters[[i]],listOfClusters[[j]], 1)[2]
resultMatrix[i,j] <- result
}
}

rownames(resultMatrix) <-  names(listOfClusters)
colnames(resultMatrix) <-  names(listOfClusters)

resultMatrix


faclustering


kmeansmds

leukaemie <- read.table("http://www.statistik.lmu.de/institut/lehrstuhl/semsto/Lehre/Multivariate2014/Tutorium/leukaemie.txt",header=T)



