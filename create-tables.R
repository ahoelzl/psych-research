source("Daten_einlesen.R")
source("GesamtDatenAnalysen.R")

source("number-clusters-Samples.R")
source("number-clusters-Simulation.R")
source("SampleClusterSimiliaritySimulation.R")
source("EFASimulationClusterSimiliaritySimulation.R")
source("faktorensimulation.R")
source("CFA.R")
library(psych)
library(gridExtra)
addError <- F

paintTable <- function(table, title,footnote) {

  title <- gsub(" ", "", title)
  filename <- paste0("currentresults2/",title,".jpeg")
  print(filename)
  jpeg(filename, width=1000, height=1000)
  
  ###rounds if possible
  
  for(i in 1:dim(table)[1]) {
    for(j in 1:dim(table)[2]) {
      if(is.numeric(table[i,j])) {
      table[i,j] <- round(table[i,j], 2)
      }
    }
  }
  
  table <- tableGrob(table)
  
  grid.newpage()
  h <- grobHeight(table)
  w <- grobWidth(table)
  title2 <- textGrob(title, y=unit(0.5,"npc") + 0.5*h, 
                    vjust=0, gp=gpar(fontsize=20))
  footnote <- textGrob(footnote, 
                       x=unit(0.5,"npc") - 0.5*w,
                       y=unit(0.5,"npc") - 0.5*h, 
                       vjust=1, hjust=0,gp=gpar( fontface="italic"))
  gt <- gTree(children=gList(table, title2, footnote))
 

 

  grid.draw(gt)
  dev.off()
  
}




toSimulate <- c("faclust", "averagecor", "completecor", "averagecorcor", "completecorcor", "kmeansmds")
method.names.EFA <<- c("MAP", "Paralell-mcomp", "Paralell-nfact", "AIC")
method.names.normal <<- c("APN" ,"AD" ,"ADM" ,"FOM","Connectivity", "Dunn" ,"Silhouette")


allnobs <-c(100,200, 500,1000)


test1 <- getClusterSimiliarity.simulation.methods(methods=c(1,2,3),zuordnung.ges, toSimulate, fa.ges) 


getClusterSimiliarity.samples(nrep=200,numbercluster=5)

getClusterNumberBias.simulation.methods(types= c("kmeans", "average", "complete", "faclust"), 
                                        methods=c(1,2,3), fa.ges)


#getClusterNumberBias.simulation.methods.original(method=1, fa.ges)
#getClusterNumberBias.simulation.methods.original(method=2, fa.ges)
#getClusterNumberBias.simulation.methods.original(method=3, fa.ges)


r1 <- getClusterNumberBiasVariance.samples(200,types= c("kmeans", "average", "complete", "faclust"))


runCFR(nrep=200, nobs=1000, efa=T)

runCFR(nrep=200, nobs=1000, efa=F)

