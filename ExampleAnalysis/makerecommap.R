############################################################
##' makerecommap.R example.
##' this R script has several purposes in this analysis. 
##' 1st, it takes in the RDS data structures containing the empirical results from david 
##' and creates a csv containing all of the data needed by slim
##' 2nd, it takes in the position loci chosen to be included in each simulation (see parammaker.R)
##' and creates a vector of recombination fractions b/w adjacent SNPs included in each of the models
############################################################

##  first let's read in the empirical data from david
library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/ExampleAnalysis/")
## this is the main data frame that contains most of the data we need. However, it includes antiquated AFC cuttoffs 
## since this data frame was created, we used the AFC > 99.9% of neutral simulations to identify candidate snps
## but the rest of the info contained should be correct
x <- readRDS(file = "haplotype_blocks.snp_res.RDS")
## the newer cutoff values are in this other dataframe. Along with some information to identify the snps
y <- readRDS("hap_block_snps.neutral_AFC_cutoffs.RDS")

## Here I am taking the useful information from these two dataframes and putting together.
tmpdata <- cbind(x$chr, x$pos, x$haplotype_block, y$T0, x$selCoef, x$cov, y$Gen10_neutAFC99, y$Gen10_neutAFC999)
tmpdata <- data.frame(tmpdata)
colnames(tmpdata) <- c("chr", 'pos', 'haplotype_block', 'T0', 'selCoef', 'cov', 'Gen10_neutAFC99', 'Gen10_neutAFC999')
## fixing some of the classes in the new dataframe.
tmpdata$pos <- as.integer(tmpdata$pos)
tmpdata$T0 <- as.numeric(tmpdata$T0)
tmpdata$selCoef <- as.numeric(tmpdata$selCoef)
tmpdata$cov <- as.integer(tmpdata$cov)
tmpdata$Gen10_neutAFC99 <- as.numeric(tmpdata$Gen10_neutAFC99)
tmpdata$Gen10_neutAFC999 <- as.numeric(tmpdata$Gen10_neutAFC999)
## snps were ordered based on thier position on thr chomosome, not position on the haplotype blocks,
## this can cause issues for how the data is read into slim, so we reorder the snps first based on which
## hplotype block they belong to, then by their position on these haplotypeblocks
sorteddata <- c()
for(haps in unique(tmpdata$haplotype_block)){
  hapdata <- tmpdata[which(tmpdata$haplotype_block == haps),]
  sorteddata <- rbind(sorteddata,hapdata[order(as.integer(hapdata$pos)),])
}
## this concludes our making the new dataframe. It is saved as sortedsnpdata.csv below.

## Next we're gonna make the recombination map for all of the simulations
## T californiculs recombination rate per base pair
RF.bp <- 0.016 / 1e6
## this is a matrix containing all of the loci selected for each simulation in rows
sampledloci <- fread(file = "sampledloci.csv")

## this data frame is formatted the same way as the sampledloci.csv one
RecomFractions <- array(dim = c(7800,4977))
## looping through each simulation (rows)
for(sim in 1:nrow(sampledloci)){
  ## first lets sort the snps we sampled
  simloci <- sort(na.omit(unlist(sampledloci[sim,])))+1
  simloci <- sorteddata[simloci,]
  simrecommap <- c()
  ## then we will sort through the haplotype blocks the snps are possibly located on
  for(i in unique(simloci$haplotype_block)){
    ## finding the position of all of the snps which are located on the current haplotype block
    pos <- as.numeric(simloci$pos[simloci$haplotype_block == i])
    ## this checks that they are indeed in order
    if(min(diff(pos)) < 0) stop()
    ## calculating recombination fraction between snps which share a haplotype block
    RF <- RF.bp * diff(pos)
    ## adding a recombination fraction of 0.5 at the end of the haplotype block to unlink them
    simrecommap <- c(simrecommap, RF, 0.5)
  }
  ## getting rid of the final 0.5 because it is unneeded
  simrecommap <- simrecommap[-length(simrecommap)]
  ## saving data into the larger dataframe
  RecomFractions[sim,1:length(simrecommap)] <- simrecommap
}
## saving data
# write.csv(sorteddata, file = 'sortedsnpdata.csv')
# fwrite(RecomFractions, file = "RecomMap.csv", col.names = F)
# write.table(RecomFractions, file = "RecomMap.csv", sep = ",", col.names = F)
