############################################################
##' parammaker.R example.
##' this R script has two purposes in this analysis. 
##' 1st, it creates a matrix containing all of the parameter values fed into slim by the cluster.
##' These include seed values, number of loci, and parameter values
##' 2nd, it creates a dataframe of randomly selected loci to be included in each simulation
############################################################

library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/ExampleAnalysis")

## deprecated code 
# mult <- 6
# iter = 1000
# fitnessFunction <- rep(c("exponential", # 1
#                          "exponential", # 1
#                          "multiplicative", # 2
#                          "directional", # 3
#                          "diminishingReturns", # 4
#                          "stabilizing" # 5
#                          ), each = 1000)

## randomly sampling seed values for each of the simulations performed.
seed <- sample(1:2^34, 7800)
##' the fitness functions multiplicative, directional (positive and negative), shifted optimum, truncating and directional are labeled
##' 1-5 
fitnessFunction <- rep(c(1,1,2,3,4,5), each = 1300)
## creating parameter values for each model, parameters are left as 0 in the models where they are not used.
a <- rep(c(8,-8,0, 0, 40, 0), each = 1300)
s <- rep(c(0,0,0,0.1,0,0), each = 1300)
r <- rep(c(0,0,0,-75,0,0), each = 1300)
b <- rep(c(0,0,0,-0.41,-0.395,0), each = 1300)
mu <- rep(c(0,0,0,0,0,0.435), each = 1300)
std <- rep(c(0,0,0,0,0,0.0175), each = 1300)
## selecting the number of loci that will be used in each simulation.
nloci <- rep(c(round(seq(from = 10, to = 4970, length.out = 12)), 4977), each = 100, times = 6)
## assigning a simulation "number" to each set of seeds and parameter values, this helps with making sure the correct row of some data
## frames are read in by slim within the cluster
sim <- 0:(7800-1)


## deprecated code block
# a <- rep(c(8,-8,0, 0, 40, 0), each = 1000)
# s <- rep(c(0,0,0,0.1,0,0), each = 1000)
# r <- rep(c(0,0,0,-75,0,0), each = 1000)
# b <- rep(c(0,0,0,-0.41,-0.395,0), each = 1000)
# mu <- rep(c(0,0,0,0,0,0.435), each = 1000)
# std <- rep(c(0,0,0,0,0,0.0175), each = 1000)

## saving the parameter data frame for later
dat <- cbind(seed, fitnessFunction, a, s, r, b, mu, std, nloci, sim)
fwrite(dat, file = 'params.txt', col.names = F)

## Here we are randomly selecting loci to be in each simulation
## make the data frame have as many rows as simulations we are running and as many columns as the maximum number of loci 
## include in any analysis
sampledloci <- array(dim = c(7800, 4977))
for(n in 1:length(nloci)){
  ## filling in the data frame with snps which are selected, the rest are left as NAs that are delt with in slim
  sampledloci[n,1:nloci[n]] <- sample(0:4976,nloci[n])
}
# saving this dataframe
fwrite(sampledloci, file = "sampledloci.csv", col.names = F)




## deprecated below this point.
fmin <- 0
fmax <- 1
npops <- rep(10, times = mult * iter)
popsize <- 1750
scaleT0 <- 0
scales <- 0


###############
## Test data ##
###############

mult <- 3
iter = 2
seed <- sample(1:2^15, mult * iter)
# fitnessFunction <- rep(c("exponential",
#                          "exponential",
#                          "multiplicative"), each = 2)
fitnessFunction <- rep(c(1,1,2), each = 2)
a <- rep(c(8,-8,0), each = 2)
s <- rep(c(0,0,0), each = 2)
r <- rep(c(0,0,0), each = 2)
b <- rep(c(0,0,0), each = 2)
mu <- rep(c(0,0,0), each = 2)
std <- rep(c(0,0,0), each = 2)
dat <- cbind(seed, fitnessFunction, a,s,r,b,mu,std)
# fwrite(dat, file = 'testparams.txt', col.names = F)

# # read.csv(file = "figure1parms.txt")
# 
# 
# 
# "slim -d seed=", seeds[i],
# " -d npops=", npops,
# " -d nloci=", nloci,
# " -d popsize=", popsize,
# " -d " ,'"', 'fitnessFunction=', "'", 'exponential', "'", '"',
# " -d fmin=", fmin,
# " -d fmax=", fmax,
# " -d a=", ahat,
# " -d scaleT0=", scaleT0,
# " -d scales=", scales,
# " New.slim | tail -n +14 > output/exponential",
# "_seed=", seeds[i],
# "_a=", ahat,
# "_b=", bhat, ".csv", sep = ""))