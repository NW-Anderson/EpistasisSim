GetNumberofLines <- function(rawout, cutoffs){
  nloci <- 121
  signsnps <- vector("list",10)
  names(signsnps) <- paste('pop', 1:10)
  for(pop in 1:8){
    popdat <- rawout[which(rawout[,2] == pop),]
    signsnps[[pop]] <- which((popdat[3,3:ncol(popdat)] - popdat[1,3:ncol(popdat)]) >= cutoffs[1:nloci])
  }
  tmparray <- array()
  for(i in 1:8){
    tmpvec <- c()
    for(j in 1:121){
      tmpvec <- c(tmpvec, j %in% signsnps[[i]])
    }
    tmparray <- rbind(tmparray, tmpvec)
  }
  tmparray <- tmparray[-1,]
  results <- colSums(tmparray)
  return(results)
}

GetCutoffs <- function(rawout){
  loci <- colnames(rawout)[3:ncol(rawout)]
  loci <- substr(loci,5, nchar(loci))
  loci <- as.numeric(loci) + 1
  cutoffs <- empdata$Gen10_neutAFC999[loci]
  names(cutoffs) <- loci
  return(cutoffs)
}
library(data.table)
library(doMC)
library(stringr)
library(foreach)
library(ggraptR)
opts <- list(preschedule = FALSE)
registerDoMC(7)
# setwd("/media/lee/HDD_Array/nwanderson/EpistasisSim/jacccalc/alphaComparison/")
setwd("~/Documents/GitHub/EpistasisSim/alphaComparison/")
empdata <- fread(file = "sortedhbdata.csv")
setwd("~/Documents/GitHub/EpistasisSim/ABC")
load(file = "ABCoutput.RData")
abcEst <- sum(out$param * out$weights)
setwd("~/Documents/GitHub/EpistasisSim/SFS")
hap_blocks.RFS <- readRDS("~/Documents/GitHub/EpistasisSim/SFS/hap_blocks.RFS.rds")
# setwd("/media/nathan/T7/EpistasisSim/alphaComparison/")
setwd("/Volumes/T7/EpistasisSim/alphaComparison")
alphas <- list.files(path = './SLiMouts/')
alphas <- alphas[2:3] 
sim.results <- array(dim = c(6,121))
rownames(sim.results) <- paste(rep(c("alpha = 36.5", 
                                     "alpha = 0", 
                                     "Directional QT"), 
                                   each = 2), 
                               c("mean", "SEM"))
colnames(sim.results) <- paste("Loci", 1:121, sep = "")

sim.results <- list()
for(a in alphas){
  files <- list.files(path = paste("./SLiMouts/", a, "/", sep = ''))
  tmp <- foreach(sim = files, 
                 .options.multicore=opts, 
                 .combine = 'rbind') %dopar%{
                   path <- paste("./SLiMouts/", a, "/", sim, sep = '')
                   rawout <- read.csv(file = path)
                   rawout <- rawout[,-ncol(rawout)]
                   cutoffs <- GetCutoffs(rawout = rawout)
                   return(GetNumberofLines(rawout = rawout,
                                    cutoffs = cutoffs) / 8)
                 }
  sim.results[[a]] <- tmp
}


setwd("/Volumes/T7/EpistasisSim/EmpT0hb")
FFs <- list.files(path = './SLiMouts/')
ff <- FFs[2]
files <- list.files(path = paste("./SLiMouts/", ff, "/", sep = ''))
tmp <- foreach(sim = files, 
               .options.multicore=opts, 
               .combine = 'rbind') %dopar%{
                 path <- paste("./SLiMouts/", ff, "/", sim, sep = '')
                 rawout <- read.csv(file = path)
                 rawout <- rawout[,-ncol(rawout)]
                 cutoffs <- GetCutoffs(rawout = rawout)
                 return(GetNumberofLines(rawout = rawout,
                                         cutoffs = cutoffs) / 8)
               }
sim.results[[ff]] <- tmp
setwd("~/Documents/GitHub/EpistasisSim/SFS")
# write.csv(sim.results, file = "sim.results.csv")
sim.results[["empirical"]] <- hap_blocks.RFS[[2]]
save(sim.results, file = "sim.results.RData")

