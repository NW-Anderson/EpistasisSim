library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/EmpT0hbLinkage/")
x <- readRDS(file = "hap_blocks.res.RDS")
y <- readRDS("hap_blocks.neutral_AFC_cutoffs.RDS")

tmpdata <- cbind(x$chr, x$start, x$stop, x$tag, y$T0, x$selCoef, x$cov, y$Gen10_neutAFC99, y$Gen10_neutAFC999)
tmpdata <- data.frame(tmpdata)
colnames(tmpdata) <- c("chr", 'start', 'stop', 'haplotype_block', 'T0', 'selCoef', 'cov', 'Gen10_neutAFC99', 'Gen10_neutAFC999')
tmpdata$start <- as.integer(tmpdata$start)
tmpdata$stop <- as.integer(tmpdata$stop)
tmpdata$T0 <- as.numeric(tmpdata$T0)
tmpdata$selCoef <- as.numeric(tmpdata$selCoef)
tmpdata$cov <- as.integer(tmpdata$cov)
tmpdata$Gen10_neutAFC99 <- as.numeric(tmpdata$Gen10_neutAFC99)
tmpdata$Gen10_neutAFC999 <- as.numeric(tmpdata$Gen10_neutAFC999)
tmpdata$midpoint <- (tmpdata$stop + tmpdata$start) / 2
tmpdata$chr <- gsub("\\_0..*", "", tmpdata$haplotype_block)
sorteddata <- c()
for(chrom in unique(tmpdata$chr)){
  chromdata <- tmpdata[which(tmpdata$chr == chrom),]
  sorteddata <- rbind(sorteddata, chromdata[order(chromdata$midpoint),])
}



write.csv(sorteddata, file = "sortedhbdata.csv")



RF.bp <- 0.016 / 1e6 /4 

RecomFractions <- c()
for(chrom in unique(sorteddata$chr)){
  chromdata <- sorteddata[which(sorteddata$chr == chrom),]
  pos <- chromdata$midpoint
  if(min(diff(pos)) < 0) stop()
  RF <- RF.bp * diff(pos)
  RecomFractions <- c(RecomFractions, RF, 0.5)
}

RecomFractions <- RecomFractions[-length(RecomFractions)]
write.csv(RecomFractions, file = 'RecomMap4xred.csv')
