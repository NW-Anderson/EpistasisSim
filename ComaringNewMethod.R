########################
## Internal Functions ##
########################
GetJaccards <- function(rawout, cutoff, results){
  for(gen in 1:2){
    signsnps <- vector("list",10)
    names(signsnps) <- paste('pop', 1:10)
    for(pop in 1:10){
      popdat <- rawout[which(rawout[,2] == pop),]
      if(gen == 1){
        signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoff[1:nloci])
      }else if(gen == 2){
        signsnps[[pop]] <- which((popdat[3,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoff[1:nloci])
      }
    }
    jaccmat <- array(dim = c(10,10))
    for(j in 1:9){
      for(k in (j+1):10){
        jaccmat[j,k] <- length(intersect(signsnps[[j]], signsnps[[k]])) / length(union(signsnps[[j]], signsnps[[k]]))
        if(length(union(signsnps[[j]], signsnps[[k]])) == 0) jaccmat[j,k] <- 0
        cat("Intersect: ", length(intersect(signsnps[[j]], signsnps[[k]])))
        cat("    Union: ", length(union(signsnps[[j]], signsnps[[k]])), "\n")
      }
    }
    jaccmat <- na.omit(as.vector(jaccmat))
    # results[(5 * i + 1):(5 *  i + 5)] <- quantile(jaccmat)
    if(gen == 1) results[1,] <- mean(jaccmat)
    if(gen == 2) results[2,] <- mean(jaccmat)
    
  }
  return(results)
}
################
## Parameters ##
################
library(doMC)
library(stringr)
library(foreach)
library(ggraptR)
library(data.table)
opts <- list(preschedule = FALSE)
registerDoMC(5)
setwd("~/Documents/GitHub/EpistasisSim")
haplotypedata <- fread(file = "sortedsnpdata.csv")
hap_blocks.jaccard.sim10.RDS <- readRDS(file = 'hap_block_snps.jaccard.sim10.RDS')

iter <- 2

seeds <- sample(1:2^15, 5 * iter)
npops = 10
nloci = 500
popsize = 1000
fmin = 0
fmax = 1
ahat <- 8
bhat <- 0
scaleT0 <- 0
scales <- 0

########################
## Positive Epistasis ##
########################
ExpJaccards <- foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'exponential', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d a=", ahat,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " New.slim | tail -n +14 > output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                               cutoff = haplotypedata$Gen10_AFC10,
                               results = results)
  system(paste("rm ./output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
  return(results)
}
rownames(ExpJaccards) <- c("Gen6", "Gen10")
########################
## Negative Epistasis ##
########################
ahat <- -ahat
NegJaccards <- foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'exponential', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d a=", ahat,
               " -d b=", bhat,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " New.slim | tail -n +14 > output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                         cutoff = haplotypedata$Gen10_AFC10,
                         results = results)
  system(paste("rm ./output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
  return(results)
}
rownames(NegJaccards) <- c("Gen6", "Gen10")
#############################
## Multiplicative Fitness  ##
#############################
SSJaccards <- foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[2 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " sim.slim | tail -n +14 > output/polygenic.csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/polygenic.csv', sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                         cutoff = haplotypedata$Gen10_AFC10,
                         results = results)
  system("rm ./output/polygenic.csv")
  return(results)
}
rownames(SSJaccards) <- c("Gen6", "Gen10")
###########################
## Directional Epistasis ##
###########################
# load("~/Documents/GitHub/EpistasisSim/DirectionalEpistasis.RData")
# 
# shat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,1])
# rhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,2])
# bhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,3]) * 100
shat <- 0.1
rhat <- -75
bhat <- -0.41

DirectionalJaccards <- foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'directional', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d s=", shat,
               " -d r=", rhat,
               " -d b=", bhat,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " New.slim | tail -n +14 > output/directional",
               "_seed=", seeds[i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/directional',
                                            "_seed=", seeds[i],
                                            "_s=", shat,
                                            "_r=", rhat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                         cutoff = haplotypedata$Gen10_AFC10,
                         results = results)
  system(paste("rm ./output/directional",
               "_seed=", seeds[i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ''))
  return(results)
}
rownames(DirectionalJaccards) <- c("Gen6", "Gen10")
###################################
## Diminishing Returns Epistasis ##
###################################

ahat <- 40
bhat <- -0.395

DimRetJaccards <- foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[2 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'diminishingReturns', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d b=", bhat,
               " -d a=", ahat,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " New.slim | tail -n +14 > output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/diminishingReturns',
                                            "_seed=", seeds[2 * iter + i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                         cutoff = haplotypedata$Gen10_AFC10,
                         results = results)
  system(paste("rm ./output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
  return(results)
}
rownames(DimRetJaccards) <- c("Gen6", "Gen10")
###########################
## Stabilizing Epistasis ##
###########################
mu <- 0.435
std <- 0.0175

StabJaccards <-  foreach(i = 1:iter, .options.multicore=opts, .combine = 'cbind') %dopar%{
  results <- array(dim = c(2,1))
  system(paste("slim -d seed=", seeds[2 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'stabilizing', "'", '"',
               " -d mu=", mu,
               " -d std=", std,
               " -d scaleT0=", scaleT0,
               " -d scales=", scales,
               " New.slim | tail -n +14 > output/stabilizing",
               "_seed=", seeds[2 * iter + i],
               "_mu=", mu,
               "_std=", std, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/stabilizing',
                                            "_seed=", seeds[2 * iter + i],
                                            "_mu=", mu,
                                            "_std=", std, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  results <- GetJaccards(rawout = rawout,
                         cutoff = haplotypedata$Gen10_AFC10,
                         results = results)
  system(paste("rm ./output/stabilizing",
               "_seed=", seeds[2 * iter + i],
               "_mu=", mu,
               "_std=", std, ".csv", sep = ''))
  return(results)
}
rownames(StabJaccards) <- c("Gen6", "Gen10")
################################################
if(scales == 1 | scaleT0 == 1){
  DimRetJaccards$Gen6 <- rep(0.0, times = 45)
  DimRetJaccards$Gen10 <- rep(0.0, times = 45)
}
treatment <- c(rep("G. Empirical", times = 73),
               rep("A. Multitplicative", times = 90),
               rep("B. Positive Epistasis", times = 90),
               rep("C. Negative Epistasis", times = 90),
               rep("D. Directional Epistasis", times = 90),
               rep("E. Truncating Epistasis", times = 90),
               rep("F. Stabilizing Epistasis", times = 90))

generation <- c(rep(6, each = 45), rep(10, times = 28),
                rep(c(6,10), times = 6, each = 45))

Jaccards <- c(hap_blocks.jaccard.sim10.RDS[[1]],
              hap_blocks.jaccard.sim10.RDS[[2]],
              SSJaccards$Gen6,
              SSJaccards$Gen10,
              ExpJaccards$Gen6,
              ExpJaccards$Gen10,
              NegJaccards$Gen6,
              NegJaccards$Gen10,
              DirectionalJaccards$Gen6,
              DirectionalJaccards$Gen10,
              DimRetJaccards$Gen6,
              DimRetJaccards$Gen10,
              StabJaccards$Gen6,
              StabJaccards$Gen10)
data <- data.frame(treatment, generation, Jaccards)
save.image(file = "fourkrun.RData")
# png(filename = "EmpericalBoth.png", width = 1800, height = 1000)
# ggplot(data, aes(y=Jaccards, x=as.factor(treatment))) + 
#   geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
#   theme_grey() + 
#   theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5)) + 
#   guides(fill=guide_legend(title="generation")) + 
#   ggtitle("Pairwise Similarity in Allele Frequency Shifts") + 
#   xlab("Fitness Function") + 
#   ylab("Jaccard Score") + 
#   ylim(0,1) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# dev.off()

# ggraptR(data)
