################################
## Main Text: # loci analysis ##
################################

library(data.table)
library(ggraptR)
library(gridExtra)
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
  if(jaccmatrix[1,c] == "diminishingReturns") fitfun <- c(fitfun, ". Truncating QT")
  if(jaccmatrix[1,c] == "stabilizing") fitfun <- c(fitfun, "E. Stabilizing QT")
}

fitfun <- c(fitfun, rep('D. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac), empjaccs)

data <- data.frame(fitfun, Generation, nloci, meanjac)

data <- data[-which(data$fitfun == ". Truncating QT" | 
                      data$fitfun == "D. Directional QT" | 
                      data$fitfun == "E. Stabilizing QT"),]

dev.off()
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

##################################
## Main Text: Epistasis Analysis ##
##################################

setwd("~/Documents/GitHub/EpistasisSim/EmpT0hb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000),
               rep("D. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)

data <- data[-which(data$treatment == "E. Truncating QT" | 
                      data$treatment == "D. Directional QT" | 
                      data$treatment == "F. Stabilizing QT"),]
# ggraptR(data)
ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + # ggtitle("121 Hap Blocks. Empirical T0") + 
  xlab("Fitness Function") + 
  ylab("Jaccard Score") + 
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#################################
## Supplement: # loci analysis ##
#################################

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

plot1 <- ggplot(data, aes(y=meanjac, x=nloci)) + 
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
  guides(shape=F) + 
  guides(colour=F) + 
  ylab("Mean Jaccard Score") + 
  xlab("Number of Loci") + 
  ylim(c(0,1)) + 
  ggtitle("121 Haplotype Blocks")

#################################

setwd("~/Documents/GitHub/EpistasisSim/nlocisnps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
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
nloci <- c(as.numeric(jaccmatrix[3,]), 4977,4977)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac, na.rm = T), empjaccs)

data = data.frame(fitfun, Generation, nloci, meanjac)

plot2 <- ggplot(data, aes(y=meanjac, x=nloci)) + 
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
  guides(colour = guide_legend(title = "Fitness Function")) + 
  xlab("Numer of loci") + 
  theme(axis.text.y.left = element_blank()) +
  ylab("") + 
  ylim(c(0,1)) + 
  ggtitle("4977 SNPs on 121 Haplotype Blocks")

#################################

grid.arrange(plot1, plot2, ncol=2, widths = c(1,1.35))

####################################
## Supplement: Mut Drift analysis ##
####################################

setwd("~/Documents/GitHub/EpistasisSim/EmpT0hb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000),
               rep("G. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)
# ggraptR(data)
plot1 <- ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=F) + 
  ggtitle("121 Hap Blocks. Empirical T0") + 
  theme(axis.text.x.bottom = element_blank()) +
  xlab("") + 
  ylab("Jaccard Score") + 
  ylim(c(0,1)) 

####################################

setwd("~/Documents/GitHub/EpistasisSim/MutDrifthb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000),
               rep("G. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)
# ggraptR(data)
plot2 <- ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + 
  ggtitle("121 Hap Blocks. Neutral T0") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank()) +
  ylim(c(0,1)) 

####################################

setwd("~/Documents/GitHub/EpistasisSim/EmpT0snps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = '../results/EmpT0/sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000),
               rep("G. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)
# ggraptR(data)
plot3 <- ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=F) + 
  ggtitle("4977 snps on 121 hap blocks. Empirical T0") + 
  xlab("Fitness Function") + 
  ylab("Jaccard Score") + 
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

####################################

setwd("~/Documents/GitHub/EpistasisSim/MutDriftsnps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000),
               rep("G. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)
# ggraptR(data)
plot4 <- ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + ggtitle("4977 snps on 121 hap blocks. Neutral T0") + 
  xlab("Fitness Function") + 
  ylab("") + 
  ylim(c(0,1)) +
  theme(axis.text.y.left = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

####################################

grid.arrange(plot1, plot2, plot3, plot4, ncol=2,
             heights = c(1,1.35), widths = c(1,1.15))

####################################