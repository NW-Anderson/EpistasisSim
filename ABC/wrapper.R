library(doMC)
library(stringr)
library(foreach)
library(EasyABC)
setwd("/media/lee/HDD_Array/nwanderson/EpistasisSim/ABC")
# setwd("~/Documents/GitHub/EpistasisSim/ABC")

# opts <- list(preschedule = FALSE)
# registerDoMC(5)
#################
## ABC wrapper ##
#################
## test params ##
#################
seed = 1
npops = 10
nloci = 121
popsize = 1750
scaleT0 = 0
scales = 0
fitnessFunction = 1
fmin = 0
fmax = 1
s = 1
r = 1
a = 8
b = 1
mu = 1
std = 1 
par <- c(seed, npops, nloci, popsize, scaleT0, scales,
         fitnessFunction, fmin, fmax, s, r, a, b, mu, std)
rm(list=ls()[-9])

####################
## model function ##
####################
model <- function(par){
  seed <- par[1]
  npops <- par[2]
  nloci <- par[3]
  popsize <- par[4]
  scaleT0 <- par[5]
  scales <- par[6]
  fitnessFunction <- par[7]
  fmin <- par[8]
  fmax <- par[9]
  s <- par[10]
  r <- par[11]
  a <- par[12]
  b <- par[13] * 100
  mu <- par[14]
  std <- par[15]

  
  tryCatch( {
    system(paste("slim -d fmin=", fmin,
                 " -d fmax=", fmax,
                 " -d npops=", npops,
                 " -d nloci=", nloci,
                 " -d popsize=", popsize,
                 " -d scaleT0=", scaleT0,
                 " -d scales=0", scales,
                 " -d seed=", seed,
                 " -d fitnessFunction=", fitnessFunction,
                 " -d a=", a,
                 " -d s=", s,
                 " -d r=", r, 
                 " -d b=", b,
                 " -d mu=", mu,
                 " -d std=", std,
                 " New.slim | tail -n +14 > output/ff=", fitnessFunction,
                 "_seed=", seed,
                 "_a=", a,
                 "_s=", s,
                 "_r=", r,
                 "_b=", b,
                 "_mu=", mu,
                 "_std=", std,
                 ".csv", sep = ""))
    rawout <- as.matrix(read.csv(file = paste("output/ff=", fitnessFunction,
                 "_seed=", seed,
                 "_a=", a,
                 "_s=", s,
                 "_r=", r,
                 "_b=", b,
                 "_mu=", mu,
                 "_std=", std,
                 ".csv", sep = "")))[,-1159]
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
    system(paste("rm output/ff=", fitnessFunction,
                 "_seed=", seed,
                 "_a=", a,
                 "_s=", s,
                 "_r=", r,
                 "_b=", b,
                 "_mu=", mu,
                 "_std=", std,
                 ".csv", sep = ""))
    return(results)
  },error = function(e) {return(c(rep(0, 2)))} )
}
###########
## prior ##
###########
#par <- c(seed, npops, nloci, RR, popsize, fmin, fmax, s, r, a, b)

prior <- list(c("unif",10,10), # npops
              c("unif",121,121), # nloci
              c("unif",1750,1750), # popsize
              c("unif", 0, 0), # scale T0
              c("unif", 0, 0), # scales
              c("unif", 1,1), # fitnessFunction
              c("unif",0,0), # fmin
              c("unif",1,1), # fmax
              c("unif",1,1), # s
              c("unif",1,1), # r
              c("unif",0,15), # a
              c("unif",1,1), # b
              c("unif", 1,1), # mu
              c("unif", 1,1)) # std

##############
## Observed ##
##############
rawout <- readRDS("hap_blocks.jaccard.neutral999.RDS")
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
out <- ABC_sequential(method="Lenormand", use_seed=T,
                              model=model, prior=prior, summary_stat_target=observed,
                              nb_simul=1000, n_cluster = 5) 

save(out, file = "ABCoutput.RData")


