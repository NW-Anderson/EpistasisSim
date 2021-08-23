setwd("~/Documents/GitHub/EpistasisSim")

iter <- 1

seeds <- sample(1:2^15, 4 * iter)
npops = 10
nloci = 100
RR = 0.25
popsize = 2000
fmin = 0
fmax = 1


###########################
## Directional Epistasis ##
###########################
# load("~/Documents/GitHub/EpistasisSim/DirectionalEpistasis.RData")
# 
# shat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,1])
# rhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,2])
# bhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,3]) * 100
shat <- 1
rhat <- -10
bhat <- -0.4

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
# load("~/Documents/GitHub/EpistasisSim/ExponentialEpistasis.RData")
# 
# ahat <- sum(Exponential$weights * Exponential$param[,1])
# bhat <- sum(Exponential$weights * Exponential$param[,2]) * 100

ahat <- 5
bhat <- -1

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

###################################
## Diminishing Returns Epistasis ##
###################################

ahat <- 10
bhat <- -0.4

DimRetJaccards <- list("Gen6" = c(),
                            "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[2 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'diminishingReturns', "'", '"',
               " -d fmin=", fmin,
               " -d fmax=", fmax,
               " -d b=", bhat,
               " -d a=", ahat,
               " EpistaticSim.slim | tail -n +14 > output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/diminishingReturns',
                                            "_seed=", seeds[2 * iter + i],
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
    DimRetJaccards[[gen]] <- c(DimRetJaccards[[gen]], jaccmat)
  }
  system(paste("rm ./output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}

#####################
## Selective Sweep ##
#####################
SSJaccards <- list("Gen6" = c(),
                   "Gen10" = c())

for(i in 1:iter){
  system(paste("slim -d seed=", seeds[3 * iter + i],
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

##############
## Plotting ##
##############
library(vioplot)
treatment <- c(rep("Diminishing Returns", times = 90), 
               rep("Directional", times = 90), 
               rep("Exponential", times = 90), 
               rep("Selective Sweep", times = 90))

generation <- c(rep(c(6,10), times = 4, each = 45))
Jaccards <- c(DimRetJaccards$Gen6, 
              DimRetJaccards$Gen10,
              DirectionalJaccards$Gen6,
              DirectionalJaccards$Gen10,
              ExponentialJaccards$Gen6,
              ExponentialJaccards$Gen10,
              SSJaccards$Gen6,
              SSJaccards$Gen10)
data <- data.frame(treatment, generation, Jaccards)

ggraptR(data)

####################
## Finished Plots ##
####################

ggplot(data, aes(y=Jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), 
               stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  geom_violin(aes(fill=as.factor(generation)), stat="ydensity", 
              position="dodge", alpha=0.5, trim=TRUE, scale="area") + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", 
                          size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="generation")) + 
  ggtitle("1156 Loci: No Linkage") + 
  xlab("Treatment") + 
  ylab("Jaccard Score")




