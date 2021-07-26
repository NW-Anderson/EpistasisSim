setwd("~/Documents/GitHub/EpistasisSim")

iter <- 100

seeds <- sample(1:2^15, 2 * iter)
npops = 10
nloci = 1156
RR = 0.5
popsize = 100
fmin = 0
fmax = 1


###########################
## Directional Epistasis ##
###########################
load("~/Documents/GitHub/EpistasisSim/DirectionalEpistasis.RData")

shat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,1])
rhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,2])
bhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,3]) * 100

DirectionalJaccards <- list("Gen6" = c(),
                            "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'directional', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d s=", shat,
               " -d r=", rhat,
               " -d b=", bhat,
               " EpistaticSim.slim | tail -n +14 > output/directional",
               "_seed=", seeds[i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/directional',
                                            "_seed=", seeds[i],
                                            "_s=", shat,
                                            "_r=", rhat,
                                            "_b=", bhat, ".csv", sep = "")))[,-1159]
  for(gen in 1:2){
    signsnps <- vector("list",10)
    names(signsnps) <- paste('pop', 1:10)
    for(pop in 1:10){
      popdat <- rawout[which(rawout[,2] == pop),]
      if(gen == 1){
        signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
      }else if(gen == 2){
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
    DirectionalJaccards[[gen]] <- c(DirectionalJaccards[[gen]], jaccmat)
  }
  system(paste("rm ./output/directional",
               "_seed=", seeds[i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ''))
}

###########################
## Exponential Epistasis ##
###########################
load("~/Documents/GitHub/EpistasisSim/ExponentialEpistasis.RData")

ahat <- sum(Exponential$weights * Exponential$param[,1])
bhat <- sum(Exponential$weights * Exponential$param[,2]) * 100

ExponentialJaccards <- list("Gen6" = c(),
                            "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'exponential', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d a=", ahat,
               " -d b=", bhat,
               " EpistaticSim.slim | tail -n +14 > output/exponential",
               "_seed=", seeds[iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[iter+i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))[,-1159]
  for(gen in 1:2){
    signsnps <- vector("list",10)
    names(signsnps) <- paste('pop', 1:10)
    for(pop in 1:10){
      popdat <- rawout[which(rawout[,2] == pop),]
      if(gen == 1){
        signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
      }else if(gen == 2){
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
    ExponentialJaccards[[gen]] <- c(ExponentialJaccards[[gen]], jaccmat)
  }
  system(paste("rm ./output/exponential",
               "_seed=", seeds[iter+i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}
