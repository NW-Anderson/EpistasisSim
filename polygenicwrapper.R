library(doMC)
library(stringr)
library(foreach)
library(EasyABC)
# setwd(("/media/lee/HDD_Array/nwanderson/EpistasisSim"))
setwd("~/Documents/GitHub/EpistasisSim")
system(paste("slim Sim.slim | tail -n +14 > output/polygenic.csv", sep = ""))
rawout <- as.matrix(read.csv(file = paste('./output/polygenic.csv', sep = "")))[,-1159]
results <- vector(length = 2)
for(i in 0:1){
  signsnps <- vector("list",10)
  names(signsnps) <- paste('pop', 1:10)
  for(pop in 1:10){
    popdat <- rawout[which(rawout[,2] == pop),]
    if(i == 0){
      signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
    }else if(i == 1){
      signsnps[[pop]] <- which((popdat[3,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
    }
  }
  jaccmat <- array(dim = c(10,10))
  for(j in 1:9){
    for(k in (j+1):10){
      jaccmat[j,k] <- length(intersect(signsnps[[j]], signsnps[[k]])) / length(union(signsnps[[j]], signsnps[[k]]))
    }
  }
  jaccmat <- na.omit(as.vector(jaccmat))
  # results[(5 * i + 1):(5 *  i + 5)] <- quantile(jaccmat)
  results[i + 1] <- c(mean(jaccmat))
}
# system(paste("rm ./output/polygenic.csv", sep = ''))
results
