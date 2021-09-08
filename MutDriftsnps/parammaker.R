library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/MutDriftsnps")
mult <- 6
iter = 1000
seed <- sample(1:2^15, mult * iter)
# fitnessFunction <- rep(c("exponential", # 1
#                          "exponential", # 1
#                          "multiplicative", # 2
#                          "directional", # 3
#                          "diminishingReturns", # 4
#                          "stabilizing" # 5
#                          ), each = 1000)
fitnessFunction <- rep(c(1,1,2,3,4,5), each = 1000)
a <- rep(c(8,-8,0, 0, 40, 0), each = 1000)
s <- rep(c(0,0,0,0.1,0,0), each = 1000)
r <- rep(c(0,0,0,-75,0,0), each = 1000)
b <- rep(c(0,0,0,-0.41,-0.395,0), each = 1000)
mu <- rep(c(0,0,0,0,0,0.435), each = 1000)
std <- rep(c(0,0,0,0,0,0.0175), each = 1000)

dat <- cbind(seed, fitnessFunction, a,s,r,b,mu,std)
fwrite(dat, file = 'params.txt', col.names = F)


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