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
setwd("/media/lee/HDD_Array/nwanderson/EpistasisSim/jacccalc/nlocihb/")
# setwd("~/Documents/GitHub/EpistasisSim/nlocihb/")
empdata <- fread(file = "sortedhbdata.csv")
# setwd("/media/nathan/T7/EpistasisSim/nlocihb/")
# setwd("/Volumes/T7/EpistasisSim/nlocihb")
FFs <- list.files(path = './SLiMouts/')
sim.results <- array(dim = c(103, 2 * length(FFs) * 13))
sim.results[1,] <- rep(c("positive",
                         "negative", # 1
                         "multiplicative", # 2
                         "directional", # 3
                         "diminishingReturns", # 4
                         "stabilizing"), # 5
                       each = 2, times = 13)
sim.results[2,] <- rep(c(6,10), each = 1, times = 6*13)
sim.results[3,] <- rep(c(seq(from = 10, to = 120, by = 10), 121), each = 12, times = 1)
for(ff in FFs){
  files <- list.files(path = paste("./SLiMouts/", ff, "/", sep = ''))
  splt <- do.call(rbind, strsplit(files, "_"))
  for(nloci in c(seq(from = 10, to = 120, by = 10), 121)){
    nlocisims <- which(splt[,9] == paste("nloci=", nloci, sep = ""))
    meanjaccs <- foreach(sim = files[nlocisims], .options.multicore=opts, .combine = 'rbind') %dopar%{
      path <- paste("./SLiMouts/", ff, "/", sim, sep = '')
      rawout <- read.csv(file = path)
      rawout <- rawout[,-ncol(rawout)]
      cutoffs <- GetCutoffs(rawout = rawout)
      return(GetJaccards(rawout = rawout, cutoffs = cutoffs))
    }
    sim.results[4:103, which(sim.results[1,] == ff & sim.results[3,] == nloci)] <- meanjaccs
  }
}
write.csv(sim.results, file = "sim.results.csv")
