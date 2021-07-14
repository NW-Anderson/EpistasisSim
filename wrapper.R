library(doMC)
library(stringr)
library(foreach)
opts <- list(preschedule = FALSE)
registerDoMC(5)
#################
## ABC wrapper ##
#################
## test params ##
#################
seed = 1
npops = 10
nloci = 1156
RR = 0.5
popsize = 100
fmin = 0
fmax = 1
s = 1
r = -0.05
a = 0.2
b = -130
par <- c(seed, npops, nloci, RR, popsize, fmin, fmax, s, r, a, b)
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
  b <- par[11]
  
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
               "_b=", b, ".csv", sep = ""), intern = T)
  
  rawout <- as.matrix(read.csv(file = paste('./output/', list.files(path = './output')[1], sep = '')))[,-1159]
  sigsnps6 <- vector("list",10)
  names(sigsnps6) <- paste('pop', 1:10)
  for(pop in 1:10){
    popdat <- rawout[which(rawout[,2] == pop),]
    sigsnps6[[pop]] <- which((popdat[2,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= 0.01)
  }
  jaccmat <- array(dim = c(10,10))
  for(j in 1:9){
    for(k in (j+1):10){
      jaccmat[j,k] <- length(intersect(sigsnps6[[j]], sigsnps6[[k]])) / length(union(sigsnps6[[j]], sigsnps6[[k]]))
    }
  }
  jaccmat <- na.omit(as.vector(jaccmat))
  system(paste("rm ./output/", list.files(path = './output')[1], sep = ''))
  return(jaccmat)
}
###########
## prior ##
###########
#par <- c(seed, npops, nloci, RR, popsize, fmin, fmax, s, r, a, b)

prior <- list(c("unif",10,11), # npops 
              c("unif",1156,1156), # nloci
              c("unif",0.5,0.5), # RR
              c("unif",2000,2000), # popsize
              c("unif",0,0), # fmin
              c("unif",1,1), # fmax
              c("unif",0,1), # s
              c("unif",-5,0), # r
              c("unif",1,1), # a
              c("unif",120,150)) # b

##############
## Observed ##
##############
observed <- readRDS("jaccard.empirical.AFC_01.RDS")







################
## Deprecated ##
################