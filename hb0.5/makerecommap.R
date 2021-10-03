library(data.table)
setwd("~/Documents/GitHub/EpistasisSim/hb0.5/")
x <- readRDS(file = "hap_blocks.res.RDS")
y <- readRDS("hap_blocks.neutral_AFC_cutoffs.RDS")
z <- fread(file = "mutdriftinfo.csv")

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
sorteddata <- c()
for(haps in unique(tmpdata$haplotype_block)){
  hapdata <- tmpdata[which(tmpdata$haplotype_block == haps),]
  sorteddata <- rbind(sorteddata,hapdata[order(as.integer(hapdata$start)),])
}

sorteddata$T0 <- rep(0.5, times = 121)
sorteddata$Gen10_neutAFC99 <- rep(mean(z$Gen10_AFC99[which(z$AF == 0.5)]), times = 121)
sorteddata$Gen10_neutAFC999 <- rep(mean(z$Gen10_AFC999[which(z$AF == 0.5)]), times = 121)



write.csv(sorteddata, file = "sortedhbdata.csv")



