library(doMC)
library(stringr)
library(foreach)
library(EasyABC)
library(data.table)
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
  
  cutoffs <- c(0.03628571, 0.06996693, 0.05728657, 0.06885771, 0.08735743, 
               0.07732164, 0.08610714, 0.05750029, 0.08378586, 0.08664343, 
               0.03192971, 0.06521457, 0.06039286, 0.06771457, 0.08325029, 
               0.05935714, 0.07535714, 0.08007157, 0.04628629, 0.07892957, 
               0.08807186, 0.07421429, 0.05407229, 0.07221486, 0.07200143, 
               0.07589421, 0.07810743, 0.07525036, 0.06264329, 0.07800029, 
               0.08628571, 0.06121429, 0.05046429, 0.07867857, 0.05942943, 
               0.08096521, 0.05050086, 0.07657314, 0.08582179, 0.06650057, 
               0.08685771, 0.08793, 0.06182314, 0.08428571, 0.08614343, 
               0.07621486, 0.08246479, 0.04250057, 0.055501, 0.07707229, 
               0.06982143, 0.06571457, 0.08121514, 0.066608, 0.08817971, 
               0.04132171, 0.08789314, 0.07900029, 0.03942886, 0.06392914, 
               0.08618021, 0.05592943, 0.04435771, 0.06071429, 0.07010714, 
               0.08385857, 0.06817964, 0.06935714, 0.09207257, 0.06785857, 
               0.05721514, 0.04814343, 0.08742914, 0.06828571, 0.08032221, 
               0.04539314, 0.073572, 0.07407186, 0.04114286, 0.04942886, 
               0.053286, 0.07150114, 0.082036, 0.06667886, 0.08542943, 
               0.07910807, 0.08221429, 0.06764343, 0.08235886, 0.08578571, 
               0.06017886, 0.07600057, 0.062, 0.06843, 0.08875, 0.06650014, 
               0.06021514, 0.06417914, 0.06657229, 0.08571486, 0.06364329, 
               0.04621514, 0.08017857, 0.08735771, 0.03721486, 0.07642914, 
               0.06357143, 0.08678814, 0.055536, 0.05742886, 0.070072, 
               0.08528629, 0.080072, 0.07728571, 0.08592929, 0.07432171, 
               0.0855, 0.08403679, 0.08028629, 0.04435771, 0.08339343)

  
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
                 ".csv", sep = "")))
    rawout <- rawout[,-ncol(rawout)]
    
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


