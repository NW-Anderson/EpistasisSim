GetJaccards <- function(rawout, cutoffs){
  nloci <- length(cutoffs)
  results <- array(dim = c(1,2))
  for(gen in 1:2){
    signsnps <- vector("list",10)
    names(signsnps) <- paste('pop', 1:10)
    for(pop in 1:10){
      popdat <- rawout[which(rawout[,2] == pop),]
      if(gen == 1){
        signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoffs[1:nloci])
      }else if(gen == 2){
        signsnps[[pop]] <- which((popdat[3,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoffs[1:nloci])
      }
    }
    jaccmat <- array(dim = c(10,10))
    for(j in 1:9){
      for(k in (j+1):10){
        jaccmat[j,k] <- length(intersect(signsnps[[j]], signsnps[[k]])) / length(union(signsnps[[j]], signsnps[[k]]))
        if(length(union(signsnps[[j]], signsnps[[k]])) == 0) jaccmat[j,k] <- 0
        # cat("Intersect: ", length(intersect(signsnps[[j]], signsnps[[k]])))
        # cat("    Union: ", length(union(signsnps[[j]], signsnps[[k]])), "\n")
      }
    }
    jaccmat <- na.omit(as.vector(jaccmat))
    # results[(5 * i + 1):(5 *  i + 5)] <- quantile(jaccmat)
    if(gen == 1) results[,1] <- mean(jaccmat)
    if(gen == 2) results[,2] <- mean(jaccmat)
    
  }
  return(results)
}
GetCutoffs <- function(rawout){
  loci <- colnames(rawout)[3:ncol(rawout)]
  loci <- substr(loci,5, nchar(loci))
  loci <- as.numeric(loci) + 1
  cutoffs <- empdata$Gen10_neutAFC999[loci]
  names(cutoffs) <- loci
  return(cutoffs)
}
library(data.table)
library(doMC)
library(stringr)
library(foreach)
library(ggraptR)
opts <- list(preschedule = FALSE)
registerDoMC(7)
setwd("/media/lee/HDD_Array/nwanderson/EpistasisSim/jacccalc/alphaComparison/")
# setwd("~/Documents/GitHub/EpistasisSim/alphaComparison/")
empdata <- fread(file = "sortedhbdata.csv")
# setwd("~/Documents/GitHub/EpistasisSim/ABC")
load(file = "ABCoutput.RData")
abcEst <- sum(out$param * out$weights)
# setwd("/media/nathan/T7/EpistasisSim/alphaComparison/")
# setwd("/Volumes/T7/EpistasisSim/alphaComparison")
alphas <- list.files(path = './SLiMouts/')
sim.results <- array(dim = c(1002, 2 * length(alphas)))
sim.results[1,] <- rep(paste("alpha=", c(8,-8,abcEst,0), sep = ""), 
                       each = 2, times = 1)
sim.results[2,] <- rep(c(6,10), each = 1, times = 4)
for(a in alphas){
  files <- list.files(path = paste("./SLiMouts/", a, "/", sep = ''))
  meanjaccs <- foreach(sim = files, 
                       .options.multicore=opts, 
                       .combine = 'rbind') %dopar%{
                         path <- paste("./SLiMouts/", a, "/", sim, sep = '')
                         rawout <- read.csv(file = path)
                         rawout <- rawout[,-ncol(rawout)]
                         cutoffs <- GetCutoffs(rawout = rawout)
                         return(GetJaccards(rawout = rawout, 
                                            cutoffs = cutoffs))
                       }
  sim.results[3:1002, which(sim.results[1,] == a)] <- meanjaccs
}
write.csv(sim.results, file = "sim.results.csv")
