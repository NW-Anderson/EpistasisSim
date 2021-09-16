library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/nlocisnps/")
x <- readRDS(file = "haplotype_blocks.snp_res.RDS")
y <- readRDS("hap_block_snps.neutral_AFC_cutoffs.RDS")

tmpdata <- cbind(x$chr, x$pos, x$haplotype_block, y$T0, x$selCoef, x$cov, y$Gen10_neutAFC99, y$Gen10_neutAFC999)
tmpdata <- data.frame(tmpdata)
colnames(tmpdata) <- c("chr", 'pos', 'haplotype_block', 'T0', 'selCoef', 'cov', 'Gen10_neutAFC99', 'Gen10_neutAFC999')
tmpdata$pos <- as.integer(tmpdata$pos)
tmpdata$T0 <- as.numeric(tmpdata$T0)
tmpdata$selCoef <- as.numeric(tmpdata$selCoef)
tmpdata$cov <- as.integer(tmpdata$cov)
tmpdata$Gen10_neutAFC99 <- as.numeric(tmpdata$Gen10_neutAFC99)
tmpdata$Gen10_neutAFC999 <- as.numeric(tmpdata$Gen10_neutAFC999)
sorteddata <- c()
for(haps in unique(tmpdata$haplotype_block)){
  hapdata <- tmpdata[which(tmpdata$haplotype_block == haps),]
  sorteddata <- rbind(sorteddata,hapdata[order(as.integer(hapdata$pos)),])
}



RF.bp <- 0.016 / 1e6
sampledloci <- fread(file = "sampledloci.csv")

RecomFractions <- array(dim = c(7800,4977))
for(sim in 1:nrow(sampledloci)){
  simloci <- sort(na.omit(unlist(sampledloci[sim,])))
  simloci <- sorteddata[simloci,]
  simrecommap <- c()
  for(i in unique(simloci$haplotype_block)){
    pos <- as.numeric(simloci$pos[simloci$haplotype_block == i])
    if(min(diff(pos)) < 0) stop()
    RF <- RF.bp * diff(pos)
    simrecommap <- c(simrecommap, RF, 0.5)
  }
  simrecommap <- simrecommap[-length(simrecommap)]
  RecomFractions[sim,1:length(simrecommap)] <- simrecommap
}

# write.csv(sorteddata, file = 'sortedsnpdata.csv')
# fwrite(RecomFractions, file = "RecomMap.csv", col.names = F)
write.table(RecomFractions, file = "RecomMap.csv", sep = ",", col.names = F)
