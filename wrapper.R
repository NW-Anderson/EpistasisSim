library(doMC)
library(stringr)
library(foreach)
opts <- list(preschedule = FALSE)
registerDoMC(5)

setwd("~/Documents/GitHub/EpistasisSim")
npops = 10
nloci = 1156
RR = 0.5
popsize = 2000
fmin = 0
fmax = 1
s = 1
r = -0.05
a = 0.2
b = -130

for(s in seq(from = 0.5, to = 1.5, length.out = 5)){
  for(r in seq(from = -.15, to = 0.05, length.out = 5)){
    # for(b in seq(from = 120, to = 140, length.out = 10)) {
      foreach(b = seq(from = 122, to = 140, length.out = 10), .options.multicore=opts) %dopar%{
        system(paste("slim -d npops=", npops,
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
                     " EpistaticSim.slim > output/directional",
                     "_s=", s,
                     "_r=", r,
                     "_a=", a,
                     "_b=", b, ".txt", sep = ""))
        
      }
    # }
  }
}
