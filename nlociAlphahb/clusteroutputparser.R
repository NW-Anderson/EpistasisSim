setwd("~/Documents/GitHub/EpistasisSim/ABC")
load(file = "ABCoutput.RData")
abcEst <- sum(out$param * out$weights)
# setwd("/media/nathan/T7/EpistasisSim/nlociAlpha")
setwd("/Volumes/T7/EpistasisSim/nlociAlphahb")

filenames <- list.files(path = './SLiMouts/')
parsefilenames(filenames = filenames, abcEst = abcEst)
########################
## Internal Functions ##
########################
parsefilenames <- function(filenames, abcEst){
  alphas <- paste("alpha=", c(8,-8,1,-1,abcEst,-abcEst,0), sep = "")
  for(a in alphas){
    system(paste("mkdir ./SLiMouts/", a, sep = ""))
  }
  for(i in 1:length(filenames)){
    split <- strsplit(filenames[i], "_")[[1]]
    split[length(split)] <- substr(split[length(split)], 1, nchar(split[length(split)]) - 4)
    if(split[3] == "a=8") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=8/", sep = ""))
    if(split[3] == "a=-8") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=-8/", sep = ""))
    if(split[3] == "a=1") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=1/", sep = ""))
    if(split[3] == "a=-1") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=-1/", sep = ""))
    if(split[3] == paste("a=",abcEst, sep = "")) system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=", abcEst, "/", sep = ""))
    if(split[3] == paste("a=",-abcEst, sep = "")) system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=", -abcEst, "/", sep = ""))
    if(split[3] == "a=0") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=0/", sep = ""))
  }
}
