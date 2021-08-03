################
## Parameters ##
################
setwd("~/Documents/GitHub/EpistasisSim")

iter <- 1

seeds <- sample(1:2^15, 5 * iter)
npops = 10
nloci = 1156
RR = 0.5
popsize = 2000
fmin = 0
fmax = 1
ahat <- 5
bhat <- -1
#####################
## New Exponential ##
#####################
NewExpJaccards <- list("Gen6" = c(),
                       "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'exponential', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d a=", ahat,
               " -d b=", bhat,
               " New.slim | tail -n +14 > output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[i],
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
    NewExpJaccards[[gen]] <- c(NewExpJaccards[[gen]], jaccmat)
  }
  system(paste("rm ./output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}

#####################
## Old Exponential ##
#####################
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
               "_seed=",seeds[iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[iter + i],
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
               "_seed=", seeds[iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
  
}

#####################
## New Directional ##
#####################

shat <- 1
rhat <- -10
bhat <- -0.4

NewDirectionalJaccards <- list("Gen6" = c(),
                            "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[2 * iter + i],
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
               " New.slim | tail -n +14 > output/directional",
               "_seed=", seeds[2 * iter + i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/directional',
                                            "_seed=", seeds[2 * iter + i],
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
    NewDirectionalJaccards[[gen]] <- c(NewDirectionalJaccards[[gen]], jaccmat)
  }
  system(paste("rm ./output/directional",
               "_seed=", seeds[2 * iter + i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ''))
}

#################
## Directional ##
#################
shat <- 1
rhat <- -10
bhat <- -0.4

DirectionalJaccards <- list("Gen6" = c(),
                            "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[3 * iter + i],
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
               "_seed=", seeds[3 * iter + i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/directional',
                                            "_seed=", seeds[3 * iter + i],
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
               "_seed=", seeds[3 * iter + i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ''))
}

#####################
## Selective Sweep ##
#####################
SSJaccards <- list("Gen6" = c(),
                   "Gen10" = c())

for(i in 1:iter){
  system(paste("slim -d seed=", seeds[4 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " Sim.slim | tail -n +14 > output/polygenic.csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/polygenic.csv', sep = "")))[,-1159]
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
    SSJaccards[[gen]] <- c(SSJaccards[[gen]], jaccmat)
  }
}
################################################

treatment <- c(rep("New", times = 90), 
               rep("Old", times = 90))

generation <- c(rep(c(6,10), times = 2, each = 45))
Jaccards <- c(NewJaccards$Gen6,
              NewJaccards$Gen10,
              ExponentialJaccards$Gen6,
              ExponentialJaccards$Gen10)
data <- data.frame(treatment, generation, Jaccards)

ggraptR(data)
