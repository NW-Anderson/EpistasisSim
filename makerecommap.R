library(data.table)
setwd("~/Documents/GitHub/EpistasisSim")
x <- fread(file = "haplotype_blocks.snp_res.csv")

sorteddata <- data.frame()
for(haps in unique(x$haplotype_block)){
  hapdata <- x[x$haplotype_block == haps]
  sorteddata <- rbind(sorteddata,hapdata[order(hapdata$pos),])
}



RF.bp <- 0.016 / 1e6

RecomFractions <- c()
for(i in unique(sorteddata$haplotype_block)){
  pos <- sorteddata$pos[sorteddata$haplotype_block == i]
  if(min(diff(pos)) < 0) stop()
  RF <- RF.bp * diff(pos)
  RecomFractions <- c(RecomFractions, RF, 0.5)
}
RecomFractions <- RecomFractions[-length(RecomFractions)]
# write.csv(sorteddata, file = 'sortedsnpdata.csv')
# write.csv(RecomFractions, file = 'RecomMap.csv')
