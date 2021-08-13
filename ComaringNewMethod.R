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
        signsnps[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoff)
      }else if(gen == 2){
        signsnps[[pop]] <- which((popdat[3,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoff)
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
    results[[gen]] <- c(results[[gen]], jaccmat)
  }
  return(results)
}
################
## Parameters ##
################
library(ggraptR)
setwd("~/Documents/GitHub/EpistasisSim")
haplotypedata <- readRDS(file = 'hap_blocks.RDS')
hap_blocks.jaccard.sim10.RDS <- readRDS(file = 'hap_blocks.jaccard.sim10.RDS')

iter <- 1

seeds <- sample(1:2^15, 5 * iter)
npops = 10
nloci = 121
RR = 0.5
popsize = 1750
fmin = 0
fmax = 1
ahat <- 3
bhat <- 0

########################
## Positive Epistasis ##
########################
ExpJaccards <- list("Gen6" = c(),
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
               " New.slim | tail -n +14 > output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/exponential',
                                            "_seed=", seeds[i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  ExpJaccards <- GetJaccards(rawout = rawout,
                               cutoff = haplotypedata$Gen10_AFC10,
                               results = ExpJaccards)
  system(paste("rm ./output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}
########################
## Negative Epistasis ##
########################
ahat <- -ahat
NegJaccards <- list("Gen6" = c(),
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
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  NegJaccards <- GetJaccards(rawout = rawout,
                             cutoff = haplotypedata$Gen10_AFC10,
                             results = NegJaccards)
  system(paste("rm ./output/exponential",
               "_seed=", seeds[i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}
#############################
## Multiplicative Fitness  ##
#############################
SSJaccards <- list("Gen6" = c(),
                     "Gen10" = c())

for(i in 1:iter){
  system(paste("slim -d seed=", seeds[3 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " sim.slim | tail -n +14 > output/polygenic.csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/polygenic.csv', sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  SSJaccards <- GetJaccards(rawout = rawout,
                              cutoff = haplotypedata$Gen10_AFC10,
                              results = SSJaccards)
  system("rm ./output/polygenic.csv")
}

###########################
## Directional Epistasis ##
###########################
# load("~/Documents/GitHub/EpistasisSim/DirectionalEpistasis.RData")
# 
# shat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,1])
# rhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,2])
# bhat <- sum(ABC_SLiM$weights * ABC_SLiM$param[,3]) * 100
shat <- 0.1
rhat <- -15
bhat <- -0.3

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
  DirectionalJaccards <- GetJaccards(rawout = rawout,
                            cutoff = haplotypedata$Gen10_AFC10,
                            results = DirectionalJaccards)
  
  system(paste("rm ./output/directional",
               "_seed=", seeds[i],
               "_s=", shat,
               "_r=", rhat,
               "_b=", bhat, ".csv", sep = ''))
}
###################################
## Diminishing Returns Epistasis ##
###################################

ahat <- 10
bhat <- -0.25

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
               " New.slim | tail -n +14 > output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/diminishingReturns',
                                            "_seed=", seeds[2 * iter + i],
                                            "_a=", ahat,
                                            "_b=", bhat, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  DimRetJaccards <- GetJaccards(rawout = rawout,
                                     cutoff = haplotypedata$Gen10_AFC10,
                                     results = DimRetJaccards)
  system(paste("rm ./output/diminishingReturns",
               "_seed=", seeds[2 * iter + i],
               "_a=", ahat,
               "_b=", bhat, ".csv", sep = ''))
}


###########################
## Stabilizing Epistasis ##
###########################

mu <- 0.4
std <- 0.07

StabJaccards <- list("Gen6" = c(),
                       "Gen10" = c())
for(i in 1:iter){
  system(paste("slim -d seed=", seeds[2 * iter + i],
               " -d npops=", npops,
               " -d nloci=", nloci,
               " -d RR=",RR,
               " -d popsize=", popsize,
               " -d " ,'"', 'fitnessFunction=', "'", 'stabilizing', "'", '"',
               " -d mu=", mu,
               " -d std=", std,
               " New.slim | tail -n +14 > output/stabilizing",
               "_seed=", seeds[2 * iter + i],
               "_mu=", mu,
               "_std=", std, ".csv", sep = ""))
  rawout <- as.matrix(read.csv(file = paste('./output/stabilizing',
                                            "_seed=", seeds[2 * iter + i],
                                            "_mu=", mu,
                                            "_std=", std, ".csv", sep = "")))
  rawout <- rawout[,-ncol(rawout)]
  StabJaccards <- GetJaccards(rawout = rawout,
                                cutoff = haplotypedata$Gen10_AFC10,
                                results = StabJaccards)
  system(paste("rm ./output/stabilizing",
               "_seed=", seeds[2 * iter + i],
               "_mu=", mu,
               "_std=", std, ".csv", sep = ''))
}


################################################

treatment <- c(rep("A Empirical", times = 73),
               rep("B Multitplicative", times = 90),
               rep("C Positive Epistasis", times = 90),
               rep("D Negative Epistasis", times = 90),
               rep("E Directional Epistasis", times = 90),
               rep("F Truncating Epistasis", times = 90),
               rep("G Stabilizing Epistasis", times = 90))

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

ggraptR(data)
