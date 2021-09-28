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

fitfun <- c(fitfun, rep('E. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac), empjaccs)

data <- data.frame(fitfun, Generation, nloci, meanjac)

data <- data[-which(data$fitfun == ". Truncating QT" | 
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
  ylab("Mean Jaccard Score") +
  scale_color_viridis(discrete = T, option = "turbo") +
  scale_fill_viridis(discrete = T, option = "turbo")

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
               rep("E. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)

data <- data[-which(data$treatment == "E. Truncating QT" | 
                      data$treatment == "F. Stabilizing QT"),]
# ggraptR(data)
ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.3) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + # ggtitle("121 Hap Blocks. Empirical T0") + 
  xlab("Fitness Function") + 
  ylab("Jaccard Score") + 
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T)
  


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
  ggtitle("121 Haplotype Blocks") +
  # scale_color_brewer(palette = "Pastel1")
  scale_color_viridis(discrete = T, option = "turbo") +
  scale_fill_viridis(discrete = T, option = "turbo")

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
  xlab("Number of Loci") + 
  theme(axis.text.y.left = element_blank()) +
  ylab("") + 
  ylim(c(0,1)) + 
  ggtitle("4977 SNPs on 121 Haplotype Blocks") + 
  # scale_color_brewer(palette = "RdYlBu")
  scale_color_viridis(discrete = T, option = "turbo") +
  scale_fill_viridis(discrete = T, option = "turbo")
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
  ylim(c(0,1)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T)
totaldata <- array()


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
title <- rep("121 Hap Blocks. Empirical T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)


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
title <- rep("121 Hap Blocks. Neutral T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata,data)



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
title <- rep("4977 snps on 121 hap blocks. Empirical T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)



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
title <- rep("4977 snps on 121 hap blocks. Neutral T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)


totaldata <- totaldata[-1,]

ggplot(totaldata, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), 
               stat="boxplot", position="dodge", alpha=0.5, width=0.3) + 
  facet_wrap(~ title) + 
  theme_grey() + 
  theme(text=element_text(family="sans", 
                          face="plain", color="#000000", 
                          size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + 
  xlab("Fitness function") + 
  ylab("Jaccard Score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T)


###################################
## Supplement: Fitness Functions ##
###################################

par(mar=c(4.0, 4.0, 1.5, 1.5), mfrow = c(2,2))
###################################

nmuts <- 1:121
x <- nmuts / 121


positive <- 1.13^(nmuts) * exp(8 * x^2)
multiplicative <- 1.13^nmuts
negative <- 1.13^(nmuts) * exp(-8 * x^2)
lb <- min(positive, negative, multiplicative)
ub <- max(positive, negative, multiplicative)

positive <- positive / max(positive)
multiplicative <- multiplicative / max(multiplicative)
negative <- negative / max(negative)


plot(x = x, y = log(positive), col = turbo(7)[2], type = "l" , xlab = "Phenotype",
     ylab = "log(Relative Fitness)", lwd = 2, ylim = c(-23,0), cex.lab = 1.25, 
     cex.axis = 1.5)
lines(x = x, y = log(multiplicative), col = turbo(7)[1], lwd = 2)
lines(x = x, y = log(negative), col = turbo(7)[3], lwd = 2)
title("121 Haplotype Blocks. Multiplicative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("A. Multiplicative", "B. Positive Epistasis", 
#                                  "C. Negative Epistasis", "Initial Distribution"), 
#        col = c("orange", "cyan", "pink", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))
abline(v = 0.223786, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)
abline(v = 0.325956, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)

ci <- seq(from = 0.380902, to = 0.495603, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[2])

ci <- seq(from = 0.331542, to = 0.445071, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[1])

ci <- seq(from = 0.29318, to = 0.402245, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[3])

###################################
x <- seq(from = 0, to = 1, length.out = 200)

fmin = 0
fmax = 1
b <- -0.3
s <- 0.1
r <- -15
DirectionalEpistasis = fmin + ((fmax - fmin) / ((1 + s * exp(r *(x + b))) ^(1/s)));


a <- 10
b <- -0.25
tmp = fmin + (fmax - fmin) * (1 - 1 / exp(a * (x + b)));
tmp[tmp < 0.0] = 1e-30;
TruncatingEpistasis = tmp;

mu <- 0.4
std <- 0.07
StabilizingEpistasis = exp(-0.5 * ((x - mu)^2 / std ^ 2));


plot(x = x, y = log(DirectionalEpistasis), type = "l", col = turbo(7)[4], 
     xlab = "Phenotype", ylab = "log(Relative Fitness)", lwd = 2,
     ylim = c(-23,0), cex.lab = 1.25, cex.axis = 1.5)
lines(x = x, y = log(TruncatingEpistasis), col = turbo(7)[5], lwd = 2)
lines(x = x, y = log(StabilizingEpistasis), col = turbo(7)[6], lwd = 2)
title("121 Haplotype Blocks. Quantative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("D. Stabiliing Epistasis",
#                                  "E. Directional Epistasis", 
#                                  "F. Truncating Epistasis", 
#                                  "Initial Distribution"), 
#        col = c("red", "green", "Blue", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))
abline(v = 0.223786, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)
abline(v = 0.325956, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)

ci <- seq(from = 0.300193, to = 0.406439, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[4])

ci <- seq(from = 0.297201, to = 0.402461, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[6])

ci <- seq(from = 0.314602, to = 0.418423, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[5])

###################################

nmuts <- 1:121
x <- nmuts / 121

a <- 8
positive <- 1.13^(nmuts) * exp(8 * x^2)
multiplicative <- 1.13^nmuts
negative <- 1.13^(nmuts) * exp(-8 * x^2)
lb <- min(positive, negative, multiplicative)
ub <- max(positive, negative, multiplicative)

positive <- positive / max(positive[1:121])
multiplicative <- multiplicative / max(multiplicative[1:121])
negative <- negative / max(negative[1:121])


plot(x = x, y =  log(positive), col = turbo(7)[2], type = "l" , xlab = "Phenotype",
     ylab = " (Relative Fitness)", lwd = 2, 
     xlim = c(0, 1), cex.lab = 1.25, cex.axis = 1.5, ylim = c(-23,0))
lines(x = x, y =  log(multiplicative), col = turbo(7)[1], lwd = 2)
lines(x = x, y =  log(negative), col = turbo(7)[3], lwd = 2)
title("4977 SNPs on 121 Haplotype Blocks. Multiplicative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("A. Multiplicative", "B. Positive Epistasis", "
#                                  C. Negative Epistasis", "Initial Distribution"),
#        col = c("orange", "cyan", "pink", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.392263, to = 0.403815, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))
abline(v = 0.392263, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)
abline(v = 0.403815, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)

ci <- seq(from = 0.433199, to = 0.443856, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[2])

ci <- seq(from = 0.432522, to = 0.444504, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[1])

ci <- seq(from = 0.430782,to = 0.441271, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[3])

###################################

x <- seq(from = 0, to = 1, length.out = 200)

fmin = 0
fmax = 1
b <- -0.41
s <- 0.1
r <- -75
DirectionalEpistasis = fmin + ((fmax - fmin) / ((1 + s * exp(r *(x + b))) ^(1/s)));


a <- 40
b <- -0.395
tmp = fmin + (fmax - fmin) * (1 - 1 / exp(a * (x + b)));
tmp[tmp < 0.0] = 1e-30;
TruncatingEpistasis = tmp;

mu <- 0.435
std <- 0.0175
StabilizingEpistasis = exp(-0.5 * ((x - mu)^2 / std ^ 2));


plot(x = x, y =  log(DirectionalEpistasis), type = "l", col = turbo(7)[4], 
     xlab = "Phenotype", ylab = "log(Relative Fitness)", lwd = 2,
     ylim = c(-23,0), cex.lab = 1.25, cex.axis = 1.5, xlim = c(0, 1))
lines(x = x, y =  log(TruncatingEpistasis), col = turbo(7)[5], lwd = 2)
lines(x = x, y =  log(StabilizingEpistasis), col = turbo(7)[6], lwd = 2)
title("4977 SNPs on 121 Haplotype Blocks. Quantative Fitness Functions", adj = 0)
legend("bottomright", legend = c("A. Multiplicative", 
                                 "B. Positive Epistasis", 
                                 "C. Negative Epistasis", 
                                 "D. Directional QT", 
                                 "E. Truncating QT", 
                                 "F. Stabilizing QT",
                                 "Initial Distribution"),
       col = c(turbo(7)[1:6], rgb(0,0,0,0.15)), 
       lwd = c(2,2,2,2,2,2,8),
       bg = "white", cex = 1)

ci <- seq(from = 0.392263, to = 0.403815, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))
abline(v = 0.392263, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)
abline(v = 0.403815, lty = 2, col = rgb(0,0,0,0.15), lwd = 0.5)

ci <- seq(from = 0.401406, to = 0.412646, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[4])

ci <- seq(from = 0.401266, to = 0.412542, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[6])

ci <- seq(from = 0.401696, to = 0.412732, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(7, alpha = 0.15)[5])

###################################


