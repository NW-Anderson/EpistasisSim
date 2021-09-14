library(data.table)
library(ggraptR)
setwd("~/Documents/GitHub/EpistasisSim/MutDriftsnps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. Positive Epistasis",
                   "C. Negative Epistasis",
                   "A. Multiplicative",
                   "D. Directional QT",
                   "E. Truncating QT",
                   "F. Stabilizing QT"), each = 2000),
               rep("G. Empirical", each = 2))

generation <- c(rep(c(6,10), each = 1000, times = 6), 
                c(6,10))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards, empjaccs))
data <- data.frame(treatment, generation, jaccards)
# ggraptR(data)
ggplot(data, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_boxplot(aes(fill=as.factor(generation)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="generation")) + ggtitle("4977 snps on 121 hap blocks. Neutral T0") + 
  xlab("Fitness Function") + 
  ylab("Jaccard Score") + 
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

