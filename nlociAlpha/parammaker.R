library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/nlociAlpha")
mult <- 6
iter = 1000
# fitnessFunction <- rep(c("exponential", # 1
#                          "exponential", # 1
#                          "multiplicative", # 2
#                          "directional", # 3
#                          "diminishingReturns", # 4
#                          "stabilizing" # 5
#                          ), each = 1000)
seed <- sample(1:2^34, 9100)
fitnessFunction <- rep(c(1,1,1,1,1,1,2), each = 1300)
a <- rep(c(8,-8,1,-1,15,-15,0), each = 1300)
s <- rep(0, times = 9100)
r <- rep(0, times = 9100)
b <- rep(0, times = 9100)
mu <- rep(0, times = 9100)
std <- rep(0, times = 9100)
nloci <- rep(c(seq(from = 10, to = 120, by = 10), 121), each = 100, times = 7)
sim <- 0:(9100-1)

dat <- cbind(seed, fitnessFunction, a, s, r, b, mu, std, nloci, sim)
fwrite(dat, file = 'params.txt', col.names = F)


sampledloci <- array(dim = c(9100, 121))
for(n in 1:length(nloci)){
  sampledloci[n,1:nloci[n]] <- sample(0:120,nloci[n])
}
fwrite(sampledloci, file = "sampledloci.csv", col.names = F)





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