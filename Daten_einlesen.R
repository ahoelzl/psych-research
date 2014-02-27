#######im Moment von hier einlesen!!

library(psych)
### Daten öffnen kurz  ###

## Datensatz besteht aus Ziffern, die durch kein Zeichen getrennt sind

raw <- readLines("NEOPIR_Ost.TXT")

# Variablen werden einzeln extrahiert:

klinisch <- as.factor(substring(raw,275,275))

NEO <- substring(raw, 23, 262)
NEO <- lapply(strsplit(NEO, ""), as.numeric)
NEO <- do.call(rbind, NEO)

## Werte des NEO gehen von 1 bis 5 

is.na(NEO) <- NEO>5
summary(NEO)

# Anzahl fehlender Werte
sum(is.na(NEO))

klinisch <- factor( klinisch, levels= c(1,2,3), labels = c("stationär","normal","ambulant"))
age <- as.numeric(substring(raw,272,273))
sex <- as.factor(substring(raw,271,271))

# hier werden die Items umkodiert, die andersherum gepolt sind

NEOk <- NEO
NEOk[,61] <- 6-NEO[,61]
NEOk[,1] <- 6-NEO[,1]
NEOk[,121] <- 6-NEO[,121]
NEOk[,181] <- 6-NEO[,181]
NEOk[,36] <- 6-NEO[,36]
NEOk[,96] <- 6-NEO[,96]
NEOk[,156] <- 6-NEO[,156]
NEOk[,11] <- 6-NEO[,11]
NEOk[,71] <- 6-NEO[,71]
NEOk[,46] <- 6-NEO[,46]
NEOk[,106] <- 6-NEO[,106]
NEOk[,166] <- 6-NEO[,166]
NEOk[,21] <- 6-NEO[,21]
NEOk[,81] <- 6-NEO[,81]
NEOk[,231] <- 6-NEO[,231]
NEOk[,141] <- 6-NEO[,141]
NEOk[,56] <- 6-NEO[,56]
NEOk[,116] <- 6-NEO[,116]
NEOk[,176] <- 6-NEO[,176]
NEOk[,206] <- 6-NEO[,206]
NEOk[,236] <- 6-NEO[,236]
NEOk[,32] <- 6-NEO[,32]
NEOk[,92] <- 6-NEO[,92]
NEOk[,7] <- 6-NEO[,7]
NEOk[,67] <- 6-NEO[,67]
NEOk[,127] <- 6-NEO[,127]
NEOk[,187] <- 6-NEO[,187]
NEOk[,42] <- 6-NEO[,42]
NEOk[,102] <- 6-NEO[,102]
NEOk[,162] <- 6-NEO[,162]
NEOk[,222] <- 6-NEO[,222]
NEOk[,17] <- 6-NEO[,17]
NEOk[,77] <- 6-NEO[,77]
NEOk[,137] <- 6-NEO[,137]
NEOk[,52] <- 6-NEO[,52]
NEOk[,112] <- 6-NEO[,112]
NEOk[,27] <- 6-NEO[,27]
NEOk[,87] <- 6-NEO[,87]
NEOk[,147] <- 6-NEO[,147]
NEOk[,207] <- 6-NEO[,207]
NEOk[,33] <- 6-NEO[,33]
NEOk[,93] <- 6-NEO[,93]
NEOk[,153] <- 6-NEO[,153]
NEOk[,183] <- 6-NEO[,183]
NEOk[,213] <- 6-NEO[,213]
NEOk[,8] <- 6-NEO[,8]
NEOk[,68] <- 6-NEO[,68]
NEOk[,128] <- 6-NEO[,128]
NEOk[,43] <- 6-NEO[,43]
NEOk[,103] <- 6-NEO[,103]
NEOk[,163] <- 6-NEO[,163]
NEOk[,18] <- 6-NEO[,18]
NEOk[,78] <- 6-NEO[,78]
NEOk[,138] <- 6-NEO[,138]
NEOk[,198] <- 6-NEO[,198]
NEOk[,228] <- 6-NEO[,228]
NEOk[,53] <- 6-NEO[,53]
NEOk[,113] <- 6-NEO[,113]
NEOk[,173] <- 6-NEO[,173]
NEOk[,28] <- 6-NEO[,28]
NEOk[,88] <- 6-NEO[,88]
NEOk[,148] <- 6-NEO[,148]
NEOk[,208] <- 6-NEO[,208]
NEOk[,238] <- 6-NEO[,238]
NEOk[,4] <- 6-NEO[,4]
NEOk[,64] <- 6-NEO[,64]
NEOk[,124] <- 6-NEO[,124]
NEOk[,39] <- 6-NEO[,39]
NEOk[,99] <- 6-NEO[,99]
NEOk[,159] <- 6-NEO[,159]
NEOk[,189] <- 6-NEO[,189]
NEOk[,219] <- 6-NEO[,219]
NEOk[,14] <- 6-NEO[,14]
NEOk[,74] <- 6-NEO[,74]
NEOk[,134] <- 6-NEO[,134]
NEOk[,49] <- 6-NEO[,49]
NEOk[,109] <- 6-NEO[,109]
NEOk[,169] <- 6-NEO[,169]
NEOk[,199] <- 6-NEO[,199]
NEOk[,229] <- 6-NEO[,229]
NEOk[,24] <- 6-NEO[,24]
NEOk[,84] <- 6-NEO[,84]
NEOk[,144] <- 6-NEO[,144]
NEOk[,234] <- 6-NEO[,234]
NEOk[,59] <- 6-NEO[,59]
NEOk[,119] <- 6-NEO[,119]
NEOk[,35] <- 6-NEO[,35]
NEOk[,95] <- 6-NEO[,95]
NEOk[,155] <- 6-NEO[,155]
NEOk[,10] <- 6-NEO[,10]
NEOk[,70] <- 6-NEO[,70]
NEOk[,130] <- 6-NEO[,130]
NEOk[,190] <- 6-NEO[,190]
NEOk[,220] <- 6-NEO[,220]
NEOk[,45] <- 6-NEO[,45]
NEOk[,105] <- 6-NEO[,105]
NEOk[,20] <- 6-NEO[,20]
NEOk[,80] <- 6-NEO[,80]
NEOk[,140] <- 6-NEO[,140]
NEOk[,55] <- 6-NEO[,55]
NEOk[,115] <- 6-NEO[,115]
NEOk[,175] <- 6-NEO[,175]
NEOk[,205] <- 6-NEO[,205]
NEOk[,30] <- 6-NEO[,30]
NEOk[,90] <- 6-NEO[,90]
NEOk[,150] <- 6-NEO[,150]

## alle Rekodieren

NEOkod <- NEOk-1

NEOdf <- as.data.frame(NEOkod)

NEOdfa <- cbind(NEOdf,age)
# NEOdfs <- cbind(NEOdf,sex)

# nur noch den Teil des Datensatzes nehmen ohne die klinischen Patienten
NEOdf <- NEOdf[klinisch=="normal",]
NEOdfa <- NEOdfa[klinisch=="normal",]
# NEOdfs <- NEOdfs[klinisch=="normal",]


## Summerwerte für Facetten bilden

N1<-rowSums(NEOdf[,c(1,31,61,91,121,151,181,211)])
E1<-rowSums(NEOdf[,c(2,32,62,92,122,152,182,212)])
O1<-rowSums(NEOdf[,c(3,33,63,93,123,153,183,213)])
A1<-rowSums(NEOdf[,c(4,34,64,94,124,154,184,214)])
C1<-rowSums(NEOdf[,c(5,35,65,95,125,155,185,215)])
N2<-rowSums(NEOdf[,c(6,36,66,96,126,156,186,216)])
E2<-rowSums(NEOdf[,c(7,37,67,97,127,157,187,217)])
O2<-rowSums(NEOdf[,c(8,38,68,98,128,158,188,218)])
A2<-rowSums(NEOdf[,c(9,39,69,99,129,159,189,219)])
C2<-rowSums(NEOdf[,c(10,40,70,100,130,160,190,220)])
N3<-rowSums(NEOdf[,c(11,41,71,101,131,161,191,221)])
E3<-rowSums(NEOdf[,c(12,42,72,102,132,162,192,222)])
O3<-rowSums(NEOdf[,c(13,43,73,103,133,163,193,223)])
A3<-rowSums(NEOdf[,c(14,44,74,104,134,164,194,224)])
C3<-rowSums(NEOdf[,c(15,45,75,105,135,165,195,225)])
N4<-rowSums(NEOdf[,c(16,46,76,106,136,166,196,226)])
E4<-rowSums(NEOdf[,c(17,47,77,107,137,167,197,227)])
O4<-rowSums(NEOdf[,c(18,48,78,108,138,168,198,228)])
A4<-rowSums(NEOdf[,c(19,49,79,109,139,169,199,229)])
C4<-rowSums(NEOdf[,c(20,50,80,110,140,170,200,230)])
N5<-rowSums(NEOdf[,c(21,51,81,111,141,171,201,231)])
E5<-rowSums(NEOdf[,c(22,52,82,112,142,172,202,232)])
O5<-rowSums(NEOdf[,c(23,53,83,113,143,173,203,233)])
A5<-rowSums(NEOdf[,c(24,54,84,114,144,174,204,234)])
C5<-rowSums(NEOdf[,c(25,55,85,115,145,175,205,235)])
N6<-rowSums(NEOdf[,c(26,56,86,116,146,176,206,236)])
E6<-rowSums(NEOdf[,c(27,57,87,117,147,177,207,237)])
O6<-rowSums(NEOdf[,c(28,58,88,118,148,178,208,238)])
A6<-rowSums(NEOdf[,c(29,59,89,119,149,179,209,239)])
C6<-rowSums(NEOdf[,c(30,60,90,120,150,180,210,240)])

set.seed(089)
#N1, E2, O3, A4, C5

cols <- c(1,31,61,91,121,151,181,211, 7,37,67,97,127,157,187,217, 13,43,73,103,133,163,193,223,
                    19,49,79,109,139,169,199,229, 25,55,85,115,145,175,205,235)

####das wird normalerweise untersucht
#facs <- as.matrix((NEOdf[,cols]), ncols=48)
#facs <- as.data.frame(cbind(N1,N2,N3,N4,N5,N6,E1,E2,E3,E4,E5,E6,O1,O2,O3,O4,O5,O6,A1,A2,A3,A4,A5,A6,C1,C2,C3,C4,C5,C6))

###C1, C2, C3, C4, C5
#cols.c.only <- c(5,35,65,95,125,155,185,215,10,40,70,100,130,160,190,220,15,45,75,105,135,165,195,225
#                 ,20,50,80,110,140,170,200,230, 25,55,85,115,145,175,205,235,30,60,90,120,150,180,210,240)

#N2, E3, O4, A5, C1
#cols.different.variation <- c(6,36,66,96,126,156,186,216, 12,42,72,102,132,162,192,222,
#                              18,48,78,108,138,168,198,228, 24,54,84,114,144,174,204,234, 
#                              5,35,65,95,125,155,185,215)

#facs <- as.matrix((NEOdf[,cols]), ncols=48)
facs <- as.matrix((NEOdf[,cols]), ncols=48)


## Analyse fehlender Werte

table(rowSums(is.na(NEOdf))) # seltsamerweise keiner mit mehr als 24 fehlenden Werten; 10 haben nur fehlende Werte, die fallen sowieso immer raus
hist(rowSums(is.na(NEOdf)))
sum(rowSums(is.na(NEOdf))>24) # wieviele Fälle mit über 10% fehlender Werte
rowSums(is.na(NEOdf))[rowSums(is.na(NEOdf))>24] #  wieviele fehlende Werte haben die (alle 240)
which(rowSums(is.na(NEOdf))>24) # wer ist das?
is.na(NEOdf[which(rowSums(is.na(NEOdf))>24),]) # wo haben sie die fehlenden Werte (überall)

## Datensatz ohne unter 18jährigen





corM <- cor(facs, use="pairwise.complete.obs", method="pearson")
corcorM <- cor(corM, use="pairwise.complete.obs", method="pearson")


fa.ges <- fa(facs, nfactors=5, max.iter=100, fm="ml", rotate="promax", method="pearson")
comparing <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))

zuordnung.ges <- apply(fa.ges$loadings,1,function(x) which.max(abs(x)))


