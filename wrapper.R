library(doMC)
library(stringr)
library(foreach)
library(EasyABC)
setwd(("/media/lee/HDD_Array/nwanderson/EpistasisSim"))
# setwd("~/Documents/GitHub/EpistasisSim")

# opts <- list(preschedule = FALSE)
# registerDoMC(5)
#################
## ABC wrapper ##
#################
## test params ##
#################
# seed = 1
# npops = 10
# nloci = 1156
# RR = 0.5
# popsize = 100
# fmin = 0
# fmax = 1
# s = 1
# r = -0.05
# a = 0.2
# b = -1.3
# par <- c(seed, npops, nloci, RR, popsize, fmin, fmax, s, r, a, b)
# rm(list=ls()[-7])

####################
## model function ##
####################
model <- function(par){
  seed <- par[1]
  npops <- par[2]
  nloci <- par[3]
  RR <- par[4]
  popsize <- par[5]
  fmin <- par[6]
  fmax <- par[7]
  s <- par[8]
  r <- par[9]
  a <- par[10]
  b <- par[11] * 100
  
  # seed <- par[1]
  # s <- par[2]
  # r <- par[3]
  # b <- par[4]
  # 
  # npops <- 10 
  # nloci <- 1156
  # RR <- 0.5
  # popsize <- 100
  # fmin <- 0
  # fmax <- 1
  # a <- 1
  
  tryCatch( {
    system(paste("slim -d seed=", seed,
                 " -d npops=", npops,
                 " -d nloci=", nloci,
                 " -d RR=",RR,
                 " -d popsize=", popsize,
                 " -d " ,'"', 'fitnessFunction=', "'", 'directional', "'", '"',
                 " -d fmin=", fmin,
                 " -d fmax=", fmax,
                 " -d s=", s,
                 " -d r=", r,
                 " -d a=", a,
                 " -d b=", b,
                 " EpistaticSim.slim | tail -n +14 > output/directional",
                 "_seed=", seed,
                 "_s=", s,
                 "_r=", r,
                 "_a=", a,
                 "_b=", b, ".csv", sep = ""))
    rawout <- as.matrix(read.csv(file = paste('./output/directional',
                                              "_seed=", seed,
                                              "_s=", s,
                                              "_r=", r,
                                              "_a=", a,
                                              "_b=", b, ".csv", sep = "")))[,-1159]
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
    system(paste("rm ./output/directional",
                 "_seed=", seed,
                 "_s=", s,
                 "_r=", r,
                 "_a=", a,
                 "_b=", b, ".csv", sep = ''))
    return(results)
  },error = function(e) {return(c(rep(0, 2)))} )
}
###########
## prior ##
###########
#par <- c(seed, npops, nloci, RR, popsize, fmin, fmax, s, r, a, b)

prior <- list(c("unif",10,10), # npops
              c("unif",1156,1156), # nloci
              c("unif",0.5,0.5), # RR
              c("unif",2000,2000), # popsize
              c("unif",0,0), # fmin
              c("unif",1,1), # fmax
              c("unif",0,1), # s
              c("unif",-5,0), # r
              c("unif",1,1), # a
              c("unif",-1.50,-1.20)) # b

# prior <- list(c("unif",0,1), # s
#               c("unif",-5,0), # r
#               c("unif",-150,-120)) # b)

##############
## Observed ##
##############
rawout <- readRDS("jaccard.empirical.AFC_01.RDS")
observed <- vector(length = 2)
# for(i in 0:1){
#   observed[(5 * i + 1):(5 * i + 5)] <- quantile(rawout[[i+1]])
# }
observed[1] <- mean(rawout[[1]])
observed[2] <- mean(rawout[[2]])
rm(rawout) # , i)
##################
## Call easyABC ##
##################
ABC_SLiM <- ABC_sequential(method="Lenormand", use_seed=T,
                           model=model, prior=prior, summary_stat_target=observed,
                           nb_simul=5, n_cluster = 1) 

save(ABC_SLiM, file = "ABCoutput.RData")


