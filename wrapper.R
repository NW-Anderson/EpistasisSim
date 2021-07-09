library(doMC)
library(stringr)
library(foreach)
opts <- list(preschedule = FALSE)
registerDoMC(5)

setwd("~/Documents/GitHub/EpistasisSim")

foreach(i = 1:1000, .options.multicore=opts) %dopar%{
  slim -d
}