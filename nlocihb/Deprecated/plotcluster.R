library(data.table)
library(ggraptR)
setwd("~/Documents/GitHub/EpistasisSim/nlocihb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "positive") fitfun <- c(fitfun, "B. Positive Epistasis")
  if(jaccmatrix[1,c] == "negative") fitfun <- c(fitfun, "C. Negative Epistasis")
  if(jaccmatrix[1,c] == "multiplicative") fitfun <- c(fitfun, "A. Multiplicative")
  if(jaccmatrix[1,c] == "directional") fitfun <- c(fitfun, "D. Directional QT")
  if(jaccmatrix[1,c] == "diminishingReturns") fitfun <- c(fitfun, "E. Truncating QT")
  if(jaccmatrix[1,c] == "stabilizing") fitfun <- c(fitfun, "F. Stabilizing QT")
}

fitfun <- c(fitfun, rep('G. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac), empjaccs)

data = data.frame(fitfun, Generation, nloci, meanjac)

ggplot(data, aes(y=meanjac, x=nloci)) + 
  geom_point(aes(shape=as.factor(Generation), colour=fitfun), 
             stat="identity", 
             position="identity", 
             alpha=0.5, 
             size=4) + 
  geom_line(aes(colour=fitfun, shape = as.factor(Generation)), 
            stat="identity", 
            position="identity", 
            alpha=0.5,
            size = 1.25) + 
  theme_grey() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=15, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  guides(shape=guide_legend(title="Generation")) + 
  guides(colour=guide_legend(title="Fitness Function")) + 
  xlab("Number of loci") + 
  ylab("Mean Jaccard Score")

