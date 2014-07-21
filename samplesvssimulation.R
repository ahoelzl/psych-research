
ntry <- 50
nobs <- 100

results.samples <-c()
results.simulations <- c()

for(i in 1:ntry) {
loads <- NL.one(fa.ges$loadings)


Phi <- fa.ges$Phi

# for(z in 1:nrep) {
set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )

  sim <-  sim.structure(fx=loads,Phi=Phi, uniq=fa.ges$uniquenesses, n=nobs,  items=T, cat=5)

cor.simulation <- sim$r

samps <- sample((x=1:nrow(data)), size=nobs, replace=T)

daten.sp <- data[samps,]

cor.samples <- cor(as.matrix(daten.sp), use="pairwise.complete.obs", method="pearson")


dim.samples <- dim(cor.samples)[1] - 1
dist.samples <- getDist(cor.samples, F)
fit.samples <- cmdscale(d=dist.samples,eig=TRUE, k=dim) # k is the number of dim
points.samples <- fit.samples$points

dim.simulation <- dim(cor.simulation)[1] - 1
dist.simulation <- getDist(cor.simulation, F)
fit.simulation <- cmdscale(d=dist.simulation,eig=TRUE, k=dim) # k is the number of dim
points.simulation <- fit.simulation$points



result.samples <- clValid(obj=points.samples, nClust=2:minv,clMethods="kmeans", validation=c(validationtype), verbose=F)
number.samples <-  as.numeric(as.character(optimalScores(result.samples)[3,3]))

result.simulation <- clValid(obj=points.simulation, nClust=2:minv,clMethods="kmeans", validation=c(validationtype), verbose=F)
number.simulation<- as.numeric(as.character(optimalScores(result.simulation)[3,3]))

results.samples  <- c(results.samples , number.samples)
results.simulations <- c(results.simulations, number.simulation)


}

mean(results.samples)
mean(results.simulations)
