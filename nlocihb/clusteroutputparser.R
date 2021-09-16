# setwd("/media/nathan/T7/EpistasisSim/nlocihb")
setwd("/Volumes/T7/EpistasisSim/nlocihb")

filenames <- list.files(path = './SLiMouts/')
parsefilenames(filenames = filenames)
########################
## Internal Functions ##
########################
parsefilenames <- function(filenames){
  FFs <- c("positive",
          "negative", # 1
          "multiplicative", # 2
          "directional", # 3
          "diminishingReturns", # 4
          "stabilizing") # 5
  for(ff in FFs){
    system(paste("mkdir ./SLiMOuts/", ff, sep = ""))
  }
  for(i in 1:length(filenames)){
    split <- strsplit(filenames[i], "_")[[1]]
    split[length(split)] <- substr(split[length(split)], 1, nchar(split[length(split)]) - 4)
    if(split[1] == "ff=1" & split[3] == "a=-8") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/negative/", sep = ""))
    if(split[1] == "ff=1" & split[3] == "a=8") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/positive/", sep = ""))
    if(split[1] == "ff=2") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/multiplicative/", sep = ""))
    if(split[1] == "ff=3") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/directional/", sep = ""))
    if(split[1] == "ff=4") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/diminishingReturns/", sep = ""))
    if(split[1] == "ff=5") system(paste("mv ./SLiMouts/", filenames[i], " ./SLiMouts/stabilizing/", sep = ""))
  }
}
