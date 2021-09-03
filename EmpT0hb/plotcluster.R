library(data.table)
library(ggraptR)
setwd("~/Documents/GitHub/EpistasisSim/Empt0snps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = './results/EmpT0/sim.results.csv'))
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
ggraptR(data)
