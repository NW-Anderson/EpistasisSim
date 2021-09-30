# setwd("/media/nathan/T7/EpistasisSim/nlociAlpha")
setwd("/Volumes/T7/EpistasisSim/nlociAlpha")

filenames <- list.files(path = './SLiMouts/')
parsefilenames(filenames = filenames)
########################
## Internal Functions ##
########################
parsefilenames <- function(filenames){
  alphas <- paste("alpha=", c(8,-8,1,-1,15,-15,0), sep = "")
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
    if(split[3] == "a=15") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=15/", sep = ""))
    if(split[3] == "a=-15") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=-15/", sep = ""))
    if(split[3] == "a=0") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/alpha=0/", sep = ""))
  }
}
