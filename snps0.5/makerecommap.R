library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/snps0.5/")
x <- readRDS(file = "haplotype_blocks.snp_res.RDS")
y <- readRDS("hap_block_snps.neutral_AFC_cutoffs.RDS")
z <- fread(file = "mutdriftinfo.csv")

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

sorteddata$T0 <- rep(0.5, times = 4977)
sorteddata$Gen10_neutAFC99 <- rep(mean(z$Gen10_AFC99[which(z$AF == 0.5)]), times = 4977)
sorteddata$Gen10_neutAFC999 <- rep(mean(z$Gen10_AFC999[which(z$AF == 0.5)]), times = 4977)

RF.bp <- 0.016 / 1e6

RecomFractions <- c()
for(i in unique(sorteddata$haplotype_block)){
  pos <- as.numeric(sorteddata$pos[sorteddata$haplotype_block == i])
  if(min(diff(pos)) < 0) stop()
  RF <- RF.bp * diff(pos)
  RecomFractions <- c(RecomFractions, RF, 0.5)
}
RecomFractions <- RecomFractions[-length(RecomFractions)]
# write.csv(sorteddata, file = 'sortedsnpdata.csv')
# write.csv(RecomFractions, file = 'RecomMap.csv')
