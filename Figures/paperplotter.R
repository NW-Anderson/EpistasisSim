library(data.table)
library(ggraptR)
library(gridExtra)
library(viridis)
library(matrixStats)
library(tidyr)
library(patchwork)
dev.off()

################################
## Main Text: # loci Analysis ##
################################
totaldata <- array()
setwd("~/Documents/GitHub/EpistasisSim/nlociAlphahb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "alpha=8") fitfun <- c(fitfun, "E. alpha = 8")
  if(jaccmatrix[1,c] == "alpha=-8") fitfun <- c(fitfun, "B. alpha = - 8")
  if(jaccmatrix[1,c] == "alpha=1") fitfun <- c(fitfun, "D. alpha = 1")
  if(jaccmatrix[1,c] == "alpha=-1") fitfun <- c(fitfun, "B. alpha = - 1")
  if(jaccmatrix[1,c] == "alpha=36.5629892630851") fitfun <- c(fitfun, "C. alpha = 36.5")
  if(jaccmatrix[1,c] == "alpha=-36.5629892630851") fitfun <- c(fitfun, "A. alpha = -36.5")
  if(jaccmatrix[1,c] == "alpha=0") fitfun <- c(fitfun, "B. alpha = 0\n    (Multiplicative)")
}

fitfun <- c(fitfun, rep('E. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 182)
meanjac <- c(colMeans(meanjac), empjaccs)

sdjac <- jaccmatrix[4:103,]
sdjac <- mapply(sdjac, FUN=as.numeric)
sdjac <- matrix(data = sdjac, nrow = 100, ncol = 182)
sdjac <- c(colSds(sdjac), empjaccs)

U95 <- c()
L95 <- c()
tmp <- jaccmatrix[4:103,]
tmp <- mapply(tmp, FUN=as.numeric)
tmp <- matrix(data = tmp, nrow = 100, ncol = 182)
for(c in 1:ncol(tmp)){
  U95 <- c(U95, quantile(tmp[,c], probs = c(0.25,0.75))[2])
  L95 <- c(L95, quantile(tmp[,c], probs = c(0.25,0.75))[1])
}
rm(tmp)
U95 <- c(U95, empjaccs)
L95 <- c(L95, empjaccs)
data <- data.frame(fitfun, Generation, nloci, meanjac, sdjac, U95, L95)
data <- data[-which(data$fitfun == "B. alpha = - 8" | data$fitfun == "E. alpha = 8" |
                      data$fitfun == "D. alpha = 1" | data$fitfun == "B. alpha = - 1"),]

totaldata <- rbind(totaldata, data)


setwd("~/Documents/GitHub/EpistasisSim/nlocihb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "positive") fitfun <- c(fitfun, "B. Positive Epistasis")
  if(jaccmatrix[1,c] == "negative") fitfun <- c(fitfun, "C. Negative Epistasis")
  if(jaccmatrix[1,c] == "multiplicative") fitfun <- c(fitfun, "A. Multiplicative")
  if(jaccmatrix[1,c] == "directional") fitfun <- c(fitfun, "D. Directional QT")
  if(jaccmatrix[1,c] == "diminishingReturns") fitfun <- c(fitfun, "E. Truncating QT")
  if(jaccmatrix[1,c] == "stabilizing") fitfun <- c(fitfun, "F. Stabilizing QT")
}

fitfun <- c(fitfun, rep('G. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac), empjaccs)

sdjac <- jaccmatrix[4:103,]
sdjac <- mapply(sdjac, FUN=as.numeric)
sdjac <- matrix(data = sdjac, nrow = 100, ncol = 156)
sdjac <- c(colSds(sdjac), empjaccs)

U95 <- c()
L95 <- c()
tmp <- jaccmatrix[4:103,]
tmp <- mapply(tmp, FUN=as.numeric)
tmp <- matrix(data = tmp, nrow = 100, ncol = 156)
for(c in 1:ncol(tmp)){
  U95 <- c(U95, quantile(tmp[,c], probs = c(0.25,0.75))[2])
  L95 <- c(L95, quantile(tmp[,c], probs = c(0.25,0.75))[1])
}
rm(tmp)
U95 <- c(U95, empjaccs)
L95 <- c(L95, empjaccs)
data = data.frame(fitfun, Generation, nloci, meanjac, sdjac, U95, L95)
data = data[which(data$fitfun == "D. Directional QT"),]

totaldata <- rbind(totaldata, data)
totaldata <- totaldata[-which(totaldata$nloci == 120),]
totaldata <- totaldata[-which(totaldata$fitfun == "E. Empirical"),]
totaldata <- totaldata[-1,]

totaldata$fitfun <- factor(totaldata$fitfun, 
                           levels = sort(unique(totaldata$fitfun)), 
                           labels = c("alpha = -36.5", 
                                      "alpha = 0\n    (Multiplicative)",
                                      "alpha = 36.5", "Directional QT"))

################################

### \u03b1
plot1 <- ggplot(totaldata, aes(y=meanjac, x=nloci)) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], size = 0.75) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) +
  geom_point(aes(shape=as.factor(Generation), 
                 colour=fitfun), 
             stat="identity", 
             position="identity", 
             alpha=1, 
             size=2.5) + 
  geom_line(aes(colour=fitfun,
                shape = as.factor(Generation)),
            stat="identity",
            position="identity",
            alpha=0.6,
            size = 1) +
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=20, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  guides(shape=guide_legend(title = "Generation",
                            override.aes = list(alpha=1))) + 
  theme(legend.position = "top") +
  guides(colour = "none") + 
  xlab("Number of Loci") + 
  ylab("Mean Jaccard Score") +
  labs(tag = "b",
       title = "Genetic Parallelism\n ~ Number of Loci") +
  scale_color_manual(values = turbo(11)[c(1,6,9,4,11)]) +
  scale_x_continuous(breaks = c(seq(from = 20, 
                                    to = 100,
                                    by = 20),121)) +
  # scale_y_continuous(sec.axis = sec_axis(~ ., breaks = empjaccs, 
  #                                        labels = c("Empirical\nGen 6",
  #                                                   "Empirical\nGen 10"))) +
  geom_errorbar(aes(colour=fitfun, 
                    ymin=L95, 
                    ymax=U95), 
                width=0, alpha=1, size=1.2) + 
  geom_ribbon(aes(colour=fitfun, 
                  shape = as.factor(Generation),
                  fill = fitfun,
                  ymin=L95, 
                  ymax=U95), 
              alpha = 0.2, colour = NA,
              show.legend = F) + 
  coord_cartesian(xlim = c(10,121), clip = 'off') +
  scale_fill_manual(values = turbo(11)[c(1,6,9,4,11)]) +
  theme(axis.text.y.right = element_text(size = 12, 
                                         face = "bold", 
                                         color = turbo(11)[11])) +
  theme(axis.text.y.left = element_blank()) +
  theme(axis.title.y.left = element_blank()) +
  theme(panel.grid = element_blank()) + 
  theme(plot.tag = element_text(face = "bold"))





###################################
## Main Text: Epistasis Analysis ##
##################################

totaldata <- array()

setwd("~/Documents/GitHub/EpistasisSim/alphaComparison")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]

treatment <- c(rep(c("B. alpha = 8",
                     "C. alpha = - 8",
                     "D. alpha = 36.5",
                     "A. alpha = 0\n    (Multiplicative)"), each = 2000))

generation <- c(rep(c(6,10), each = 1000, times = 4))

jaccards <- c()
for(i in 1:8){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(c(jaccards))
data <- data.frame(treatment, generation, jaccards)

totaldata <- rbind(totaldata, data)

setwd("~/Documents/GitHub/EpistasisSim/EmpT0hb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("121 Hap Blocks. Empirical T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards)

totaldata <- rbind(totaldata, data)
totaldata <- totaldata[-1,]
totaldata <- totaldata[which(totaldata$treatment == "D. alpha = 36.5" |
                               totaldata$treatment == "A. alpha = 0\n    (Multiplicative)" |
                               totaldata$treatment == "D. Directional QT"),]

totaldata$treatment <-  factor(totaldata$treatment, 
                               levels = sort(unique(totaldata$treatment)), 
                               labels = c("Multiplicative\n(\u03b1 = 0)", 
                                          "Positive Epistasis\n(\u03b1 = 36.5)", "Directional QT"))

##################################

plot2 <-  ggplot(totaldata, aes(y=jaccards, 
                               x=as.factor(treatment))) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], 
             size = .75) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], 
             size  = .75) +
  geom_boxplot(aes(fill=as.factor(generation)), 
               stat="boxplot", 
               position="dodge", 
               alpha=1, width=0.3) + 
  theme_bw() + 
  labs(tag = "a",
       title = "Genetic Parallelism\n ~ Fitness Function") +
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=20, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) +
  xlab("Fitness Function") + 
  ylab("Mean Jaccard Score") + 
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = empjaccs, 
                                         labels = c("Empirical\nGen 6",
                                                    "Empirical\nGen 10")),
                     limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust=1)) +
  scale_color_manual(values = viridis(12)[c(3,12)]) + 
  scale_fill_manual(values = viridis(12)[c(3,12)]) +
  theme(axis.text.y.right = element_text(size = 12, 
                                         face = "bold", 
                                         color = turbo(11)[11],
                                         hjust = 0.5)) +
  theme(panel.grid = element_blank()) + 
  theme(legend.position = "top") +
  theme(plot.tag = element_text(face = "bold"))






###################################
## Main: Replicate Freq Spectrum ##
###################################

setwd("~/Documents/GitHub/EpistasisSim/SFS")
load(file = "sim.results.RData")

plotdat <- data.frame("Group" = c(rep("B. Multiplicative", length(sim.results$`alpha=0`)),
                                  rep("C. Positive Epistasis", length(sim.results$`alpha=36.5629892630851`)),
                                  rep("A. Directional QT", length(sim.results$directional)),
                                  rep("D. Empirical", length(sim.results$empirical))),
                      "Proportion" = c(sim.results$`alpha=0`,
                                       sim.results$`alpha=36.5629892630851`,
                                       sim.results$directional,
                                       sim.results$empirical))
plotdat$Group <- factor(plotdat$Group, 
                        levels = sort(unique(plotdat$Group)),
                        labels = c("Directional QT", 
                                   "Multiplicative (\u03b1 = 0)",
                                   "Positive Epistasis (\u03b1 = 36.5)",
                                   "Empirical"))

###################################

# ggplot(plotdat, aes(x=Proportion, fill=Group)) +
#   geom_histogram(aes(y=..density..),
#                  stat="bin",
#                  position=position_dodge2(reverse = FALSE,preserve ="single"),
#                  alpha=1,
#                  breaks = sort(c(unique(plotdat$Proportion) + 0.045,
#                                  unique(plotdat$Proportion) - 0.045)),
#                  center = 1) +
#   theme_bw() +
#   theme(text=element_text(family="sans",
#                           face="plain",
#                           color="#000000",
#                           size=15, hjust=0.5, vjust=0.5)) +
#   guides(fill=guide_legend(title="Fitness Function")) +
#   xlab("Proportion") +
#   ylab("Density") +
#   scale_fill_manual("Group",values=turbo(11)[c(4,6,9,11)]) +
#   scale_x_continuous(breaks = sort(unique(plotdat$Proportion))) +
#   theme(panel.grid = element_blank())

plotdat2 <- plotdat %>% 
  dplyr::count(Proportion, Group) %>% 
  complete(Group, Proportion, fill = list(n=0)) %>%
  group_by(Group) %>% 
  mutate(Prop = n/sum(n))

plotdat2$Generation <- rep(c(6,10), times = 18)
plotdat2$Gen <- rep(10, times = 36)
# plotdat2 <- rbind(plotdat2,
#                   data.frame(Group = "Directional QT",
#                              Proportion = 0,
#                              n = 0, 
#                              Prop = 0,
#                              Generation = 6,
#                              Gen = 6))

plot3 <-  ggplot(plotdat2, aes(x=(8 * Proportion), fill=Group)) +
  geom_bar(aes(y=Prop),
           stat = "identity",
           position=position_dodge2(reverse = FALSE,preserve ="single"),
           width = .85) +
  # geom_histogram(aes(y=..density..),
  #                stat="bin",
  #                position=position_dodge2(reverse = FALSE,preserve ="single"),
  #                alpha=1,
  #                breaks = sort(c(unique(plotdat$Proportion) + 0.045,
  #                                unique(plotdat$Proportion) - 0.045)),
  #                center = 1) +
  theme_bw() +
  labs(tag = "c",
       title = "Replicate\n Frequency Spectra") +
  theme(text=element_text(family="sans",
                          face="plain",
                          color="#000000",
                          size=20, hjust=0.5, vjust=0.5)) +
  guides(color=guide_legend(title="Fitness Function")) +
  xlab("Number of Lines") +
  ylab("Proportion of Hap Blocks") +
  # scale_x_continuous(breaks = sort(unique(plotdat$Proportion))) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(3,"lines"),
        legend.text = element_text(vjust = 0.5)) +
  scale_fill_manual("Group",values=turbo(11)[c(1,4,6,9,11)],
                    limits = c("Negative Epistasis (\u03b1 = - 36.5)", "Directional QT",
                               "Multiplicative (\u03b1 = 0)",
                               "Positive Epistasis (\u03b1 = 36.5)",
                               "Empirical")) +
  scale_color_manual("Group",values=turbo(11)[c(1,4,6,9,11)],
                     limits = c("Negative Epistasis (\u03b1 = - 36.5)", "Directional QT",
                                "Multiplicative (\u03b1 = 0)",
                                "Positive Epistasis (\u03b1 = 36.5)",
                                "Empirical")) +
  # geom_point(aes(shape=as.factor(Generation), 
  #                colour=Group,
  #                y=Prop), 
  #            stat="identity", 
  #            position="identity", 
  #            alpha=0, 
  #            size=4) +
  # geom_point(aes(shape=as.factor(Gen),
  #                color=Group,
  #                y=Prop), 
  #            stat="identity", 
#            position=position_dodge2(reverse = FALSE,preserve ="single", width = 0.7), 
#            alpha=c(rep(1, times = 36)), 
#            size=4) +
guides(shape="none") + 
  # geom_line(aes(y=Prop, x=(8 * Proportion), color=Group),
  #                       stat="identity",
  #                       position=position_dodge2(reverse = FALSE,preserve ="single", width = 0.7),
  #                       alpha=0.6,
  #                       size = 0.75,
  #           linetype = "dashed") +
  guides(fill=guide_legend(title="Fitness Function")) +
  theme(plot.tag = element_text(face = "bold"))




# plot3 <- ggplot(plotdat2, aes(x=(8 * Proportion), color=Group)) +
#   geom_point(aes(y=Prop),
#              stat="identity",
#              position="identity",
#              alpha=0.6,
#              size=4) +
#   geom_line(aes(y=Prop),
#             stat="identity",
#             position="identity",
#             alpha=0.6,
#             size = 1.25) +
#   theme_bw() +
#   theme(text=element_text(family="sans",
#                           face="plain",
#                           color="#000000",
#                           size=15, hjust=0.5, vjust=0.5)) +
#   guides(fill=guide_legend(title="Fitness Function")) +
#   xlab("Number of Lines") +
#   ylab("Proportion of SNPs") +
#   scale_color_manual("Group",values=turbo(11)[c(1,4,6,9,11)],
#                      limits = c("alpha = - 36.5", "Directional QT",
#                                 "Multiplicative\n(alpha = 0)",
#                                 "Positive Epistasis\n(alpha = 36.5)",
#                                 "Empirical")) +
#   scale_x_continuous(breaks = 8 * sort(unique(plotdat$Proportion))) +
#   theme(panel.grid = element_blank())
# 



###################################

plot2 + plot1 + plot3
################################
## Supplement: alpha analysis ##
################################

setwd("~/Documents/GitHub/EpistasisSim/nlociAlphahb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "alpha=8") fitfun <- c(fitfun, "F. alpha = 8")
  if(jaccmatrix[1,c] == "alpha=-8") fitfun <- c(fitfun, "B. alpha = - 8")
  if(jaccmatrix[1,c] == "alpha=1") fitfun <- c(fitfun, "E. alpha = 1")
  if(jaccmatrix[1,c] == "alpha=-1") fitfun <- c(fitfun, "C. alpha = - 1")
  if(jaccmatrix[1,c] == "alpha=36.5629892630851") fitfun <- c(fitfun, "G. alpha = 36.5")
  if(jaccmatrix[1,c] == "alpha=-36.5629892630851") fitfun <- c(fitfun, "A. alpha = - 36.5")
  if(jaccmatrix[1,c] == "alpha=0") fitfun <- c(fitfun, "D. alpha = 0\n    (Multiplicative)")
}

fitfun <- c(fitfun, rep('H. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 182)
meanjac <- c(colMeans(meanjac), empjaccs)

data <- data.frame(fitfun, Generation, nloci, meanjac)

data <- data[-which(data$fitfun == "H. Empirical"),]
data <- data[-which(data$nloci == 120),]
data$fitfun <-  factor(data$fitfun, 
                       levels = sort(unique(data$fitfun)), 
                       labels = substr(sort(unique(data$fitfun)), 
                                       4, nchar(sort(unique(data$fitfun)))))

################################

ggplot(data, aes(y=meanjac, x=nloci)) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], size = 0.75) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) +
  geom_point(aes(shape=as.factor(Generation), colour=fitfun), 
             stat="identity", 
             position="identity", 
             alpha=1, 
             size=4) + 
  geom_line(aes(colour=fitfun, shape = as.factor(Generation)), 
            stat="identity", 
            position="identity", 
            alpha=0.6,
            size = 1.25) + 
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=20, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  guides(shape=guide_legend(title="Generation")) + 
  guides(colour=guide_legend(title="Fitness Function")) + 
  xlab("Number of Loci") + 
  ylab("Mean Jaccard Score") +
  scale_color_manual(values = turbo(11)[c(1,2,3,6,7,8,9,11)]) +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = empjaccs, 
                                         labels = c("Empirical\nGen 6",
                                                    "Empirical\nGen 10")),
                     limits = c(0,1)) +
  scale_x_continuous(breaks = c(seq(from = 20, 
                                    to = 100,
                                    by = 20),121)) +
  theme(axis.text.y.right = element_text(size = 12, 
                                         face = "bold", 
                                         color = turbo(11)[11])) +
  theme(panel.grid = element_blank())






#################################
## Supplement: # loci analysis ##
#################################

setwd("~/Documents/GitHub/EpistasisSim/nlocihb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "positive") fitfun <- c(fitfun, "C. Positive Epistasis\n    (a = 8)")
  if(jaccmatrix[1,c] == "negative") fitfun <- c(fitfun, "A. Negative Epistasis\n    (a = - 8)")
  if(jaccmatrix[1,c] == "multiplicative") fitfun <- c(fitfun, "B. Multiplicative\n    (a = 0)")
  if(jaccmatrix[1,c] == "directional") fitfun <- c(fitfun, "F. Directional QT")
  if(jaccmatrix[1,c] == "diminishingReturns") fitfun <- c(fitfun, "E. Truncating QT")
  if(jaccmatrix[1,c] == "stabilizing") fitfun <- c(fitfun, "D. Stabilizing QT")
}

fitfun <- c(fitfun, rep('G. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 121,121)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac), empjaccs)

data = data.frame(fitfun, Generation, nloci, meanjac)
data <- data[-which(data$fitfun == "G. Empirical"),]
data <- data[-which(data$nloci == 120),]

plot1 <- ggplot(data, aes(y=meanjac, x=nloci)) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], size = 0.75) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], size  = 0.75) +
  geom_point(aes(shape=as.factor(Generation), colour=fitfun), 
             stat="identity", 
             position="identity", 
             alpha=1, 
             size=4) + 
  geom_line(aes(colour=fitfun, shape = as.factor(Generation)), 
            stat="identity", 
            position="identity", 
            alpha=0.6,
            size = 1.25) + 
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=20, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  guides(shape=F) + 
  guides(colour=F) + 
  ylab("Mean Jaccard Score") + 
  xlab("Number of Loci") + 
  ylim(c(0,1)) + 
  ggtitle("121 Haplotype Blocks") +
  scale_color_manual(values = turbo(11)[c(2,6,8,10,5,4,11)]) +
  annotate("text", x = 110, y = 0.6275, 
           label = "Empirical Gen 6" , color=turbo(11)[11], 
           size=4 , fontface="bold") +
  annotate("text", x = 110, y = .775, 
           label = "Empirical Gen 10" , color=turbo(11)[11], 
           size=4 , fontface="bold")+
  scale_x_continuous(breaks = c(seq(from = 20, 
                                    to = 100,
                                    by = 20),121)) +
  theme(panel.grid = element_blank())






#################################

setwd("~/Documents/GitHub/EpistasisSim/nlocisnps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


fitfun <- c()
for(c in 1:ncol(jaccmatrix)){
  if(jaccmatrix[1,c] == "positive") fitfun <- c(fitfun, "C. Positive Epistasis\n    (a = 8)")
  if(jaccmatrix[1,c] == "negative") fitfun <- c(fitfun, "A. Negative Epistasis\n    (a = - 8)")
  if(jaccmatrix[1,c] == "multiplicative") fitfun <- c(fitfun, "B. Multiplicative\n    (a = 0)")
  if(jaccmatrix[1,c] == "directional") fitfun <- c(fitfun, "F. Directional QT")
  if(jaccmatrix[1,c] == "diminishingReturns") fitfun <- c(fitfun, "E. Truncating QT")
  if(jaccmatrix[1,c] == "stabilizing") fitfun <- c(fitfun, "D. Stabilizing QT")
}

fitfun <- c(fitfun, rep('G. Empirical', times = 2))
Generation <- c(as.numeric(jaccmatrix[2,]), 6,10)
nloci <- c(as.numeric(jaccmatrix[3,]), 4977,4977)
meanjac <- jaccmatrix[4:103,]
meanjac <- mapply(meanjac, FUN=as.numeric)
meanjac <- matrix(data = meanjac, nrow = 100, ncol = 156)
meanjac <- c(colMeans(meanjac, na.rm = T), empjaccs)

data = data.frame(fitfun, Generation, nloci, meanjac)

data <- data[-which(data$fitfun == "G. Empirical"),]
data <- data[-which(data$nloci == 4970),]
data$fitfun <-  factor(data$fitfun, 
                       levels = sort(unique(data$fitfun)), 
                       labels = substr(sort(unique(data$fitfun)), 
                                       4, nchar(sort(unique(data$fitfun)))))


plot2 <- ggplot(data, aes(y=meanjac, x=nloci)) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], size = 1.) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], size  = 1.) +
  geom_point(aes(shape=as.factor(Generation), colour=fitfun), 
             stat="identity", 
             position="identity", 
             alpha=1, 
             size=4) + 
  geom_line(aes(colour=fitfun, shape = as.factor(Generation)), 
            stat="identity", 
            position="identity", 
            alpha=0.6,
            size = 1.25) + 
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=20, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  scale_size(range=c(1, 3)) + 
  guides(shape=guide_legend(title="Generation")) + 
  guides(colour = guide_legend(title = "Fitness Function")) + 
  xlab("Number of Loci") + 
  theme(axis.text.y.left = element_blank()) +
  ylab("") + 
  ylim(c(0,1)) + 
  ggtitle("4977 SNPs on 121 Haplotype Blocks") + 
  scale_color_manual(values = turbo(11)[c(2,6,8,10,5,4,11)]) +
  annotate("text", x = 4500, y = 0.555, 
           label = "Empirical Gen 6" , color=turbo(11)[11], 
           size=4 , fontface="bold") +
  annotate("text", x = 4500, y = .595, 
           label = "Empirical Gen 10" , color=turbo(11)[11], 
           size=4 , fontface="bold") +
  scale_x_continuous(breaks = sort(unique(data$nloci))[which(1:12 %% 2 == 0)]) +
  theme(panel.grid = element_blank())







#################################

grid.arrange(plot1, plot2, ncol=2, widths = c(1,1.35))

#####################
## Supplement: RFS ##
#####################
setwd("~/Documents/GitHub/EpistasisSim/SFS")
load(file = "sim.results.RData")

plotdat <- data.frame("Group" = c(rep("B. Multiplicative", length(sim.results$`alpha=0`)),
                                  rep("C. Positive Epistasis", length(sim.results$`alpha=36.5629892630851`)),
                                  rep("A. Directional QT", length(sim.results$directional)),
                                  rep("D. Empirical", length(sim.results$empirical))),
                      "Proportion" = c(sim.results$`alpha=0`,
                                       sim.results$`alpha=36.5629892630851`,
                                       sim.results$directional,
                                       sim.results$empirical))
plotdat$Group <- factor(plotdat$Group, 
                        levels = sort(unique(plotdat$Group)),
                        labels = c("Directional QT", 
                                   "Multiplicative\n(alpha = 0)",
                                   "Positive Epistasis\n(alpha = 36.5)",
                                   "Empirical"))
plotdat <- plotdat %>% 
  dplyr::count(Proportion, Group) %>% 
  complete(Group, Proportion, fill = list(n=0)) %>%
  group_by(Group) %>% 
  mutate(Prop = n/sum(n))
plotdat$Gen <- rep(10, times = 36)
plotdat$Generation <- rep("Generation 10", times = 36)
totaldata <- plotdat

load(file = "sim.results.gen6.RData")

plotdat <- data.frame("Group" = c(rep("B. Multiplicative", length(sim.results$`alpha=0`)),
                                  rep("C. Positive Epistasis", length(sim.results$`alpha=36.5629892630851`)),
                                  rep("A. Directional QT", length(sim.results$directional)),
                                  rep("D. Empirical", length(sim.results$empirical))),
                      "Proportion" = c(sim.results$`alpha=0`,
                                       sim.results$`alpha=36.5629892630851`,
                                       sim.results$directional,
                                       sim.results$empirical))
plotdat$Group <- factor(plotdat$Group, 
                        levels = sort(unique(plotdat$Group)),
                        labels = c("Directional QT", 
                                   "Multiplicative\n(alpha = 0)",
                                   "Positive Epistasis\n(alpha = 36.5)",
                                   "Empirical"))
plotdat <- plotdat %>% 
  dplyr::count(Proportion, Group) %>% 
  complete(Group, Proportion, fill = list(n=0)) %>%
  group_by(Group) %>% 
  mutate(Prop = n/sum(n))


plotdat$Gen <- rep(6, times = 44)
plotdat$Generation <- rep("Generation 6", times = 44)

totaldata <- rbind(totaldata, plotdat)

totaldata$Proportion[which(totaldata$Generation == "Generation 6")] <- 
  totaldata$Proportion[which(totaldata$Generation == "Generation 6")] * 10 
totaldata$Proportion[totaldata$Generation == "Generation 10"] <- 
  totaldata$Proportion[totaldata$Generation == "Generation 10"] * 8

totaldata$Generation <- factor(totaldata$Generation,
                               levels = c("Generation 6",
                                          "Generation 10"),
                               labels = c("Generation 6",
                                          "Generation 10"))
#####################


ggplot(totaldata, aes(x=(Proportion), color=Group)) +
  geom_point(aes(y=Prop,
                 shape=as.factor(Gen)),
             stat="identity",
             position="identity",
             alpha=1,
             size=4) +
  geom_line(aes(y=Prop,
                shape=as.factor(Gen)),
            stat="identity",
            position="identity",
            alpha=0.6,
            size = 1.25) +
  theme_bw() +
  theme(text=element_text(family="sans",
                          face="plain",
                          color="#000000",
                          size=15, hjust=0.5, vjust=0.5)) +
  guides(fill=guide_legend(title="Fitness Function")) +
  guides(shape=guide_legend(title = "Generation")) + 
  xlab("Number of Lines") +
  ylab("Proportion of SNPs") +
  scale_color_manual("Group",values=turbo(11)[c(4,6,9,11)],
                     limits = c("Directional QT",
                                "Multiplicative\n(alpha = 0)",
                                "Positive Epistasis\n(alpha = 36.5)",
                                "Empirical")) +
  scale_x_continuous(breaks = seq(from = 0, to = 10, by  = 2)) +
  theme(panel.grid = element_blank()) + 
  facet_wrap(~ as.factor(Generation), 
             scales = "free") +
  ylim(c(0,.4))

#########################
## Supplement: Linkage ##
#########################
setwd("~/Documents/GitHub/EpistasisSim/EmpT0hb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("Freely Recombining", times = length(jaccards))
totaldata <- data.frame(treatment, generation, jaccards, title)

setwd("~/Documents/GitHub/EpistasisSim/EmpT0hbLinkage/2xRed")

jaccmatrix <- as.matrix(fread(file = 'sim.results.csv'))
jaccmatrix <- jaccmatrix[-1,-1]


treatment <- c(rep(c("B. Positive Epistasis",
                     "C. Negative Epistasis",
                     "A. Multiplicative",
                     "D. Directional QT",
                     "E. Truncating QT",
                     "F. Stabilizing QT"), each = 2000))

generation <- c(rep(c(6,10), each = 1000, times = 6))

jaccards <- c()
for(i in 1:12){
  jaccards <- c(jaccards, jaccmatrix[3:1002,i])
}
jaccards <- as.numeric(jaccards)
title <- rep("Linked 2x Reduced RR", times = length(jaccards))
totaldata <- rbind(totaldata, data.frame(treatment, generation, jaccards, title))



setwd("~/Documents/GitHub/EpistasisSim/EmpT0hbLinkage")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("Linked", times = length(jaccards))
totaldata <- rbind(totaldata, data.frame(treatment, generation, jaccards, title))






totaldata <- totaldata[-which(totaldata$treatment == "G. Empirical"),]
totaldata$treatment <- factor(totaldata$treatment, 
                              levels = sort(unique(totaldata$treatment)), 
                              labels = c("Multiplicative",
                                         "Positive Epistasis", 
                                         "Negative Epistasis", 
                                         "Directional QT",    
                                         "Truncating QT",
                                         "Stabilizing QT" ))

#########################

ggplot(totaldata, aes(y=jaccards, x=as.factor(treatment))) + 
  geom_hline(yintercept=empjaccs[1], 
             linetype="dotdash", 
             color = turbo(11)[11], 
             size = .75) +
  geom_hline(yintercept=empjaccs[2], 
             linetype="dashed", 
             color = turbo(11)[11], 
             size  = .75) +
  geom_boxplot(aes(fill=as.factor(generation)), 
               stat="boxplot", position="dodge", 
               alpha=1, width=0.3) + 
  facet_wrap(~ title) + 
  theme_bw() +  
  theme(text=element_text(family="sans", face="plain", 
                          color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="Generation")) + 
  xlab("Fitness Function") + 
  ylab("Mean Jaccard Score") + 
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = empjaccs, 
                                         labels = c("Empirical\nGen 6",
                                                    "Empirical\nGen 10")),
                     limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust=1)) +
  scale_color_manual(values = viridis(12)[c(3,12)]) + 
  scale_fill_manual(values = viridis(12)[c(3,12)]) +
  theme(axis.text.y.right = element_text(size = 12, 
                                         face = "bold", 
                                         color = turbo(11)[11])) +
  theme(panel.grid.minor = element_blank())

####################################
## Supplement: Mut Drift analysis ##
####################################

totaldata <- array()


setwd("~/Documents/GitHub/EpistasisSim/EmpT0hb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("121 Hap Blocks. Empirical T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)


setwd("~/Documents/GitHub/EpistasisSim/MutDrifthb")
empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("121 Hap Blocks. Neutral T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata,data)

setwd("~/Documents/GitHub/EpistasisSim/hb0.5")

empjaccs <- readRDS(file = 'hap_blocks.jaccard.neutral999.RDS')
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
title <- rep("121 Hap Blocks. T0=0.5", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata,data)

setwd("~/Documents/GitHub/EpistasisSim/EmpT0snps")
empjaccs <- readRDS(file = 'hap_block_snps.jaccard.neutral999.RDS')
tmp <- c(mean(empjaccs[[1]]), mean(empjaccs[[2]]))
empjaccs <- tmp
rm(tmp)
jaccmatrix <- as.matrix(fread(file = '../results/EmpT0/sim.results.csv'))
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
title <- rep("4977 SNPs. Empirical T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)



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
title <- rep("4977 SNPs. Neutral T0", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)



setwd("~/Documents/GitHub/EpistasisSim/snps0.5")
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
title <- rep("4977 SNPs. T0=0.5", times = length(jaccards))
data <- data.frame(treatment, generation, jaccards, title)
totaldata <- rbind(totaldata, data)

totaldata <- totaldata[-1,]
####################################
hlinedf <- totaldata[which(totaldata$treatment == "G. Empirical"),]
hlinegen10df <- hlinedf[which(hlinedf$generation == 10),]
hlinegen6df <- hlinedf[which(hlinedf$generation == 6),]
rm(hlinedf)

totaldata <- totaldata[-which(totaldata$treatment == "G. Empirical"),]
totaldata$treatment <- factor(totaldata$treatment, 
                              levels = sort(unique(totaldata$treatment)), 
                              labels = c("Multiplicative",
                                         "Positive Epistasis", 
                                         "Negative Epistasis", 
                                         "Directional QT",    
                                         "Truncating QT",
                                         "Stabilizing QT" ))

ann_text1 <- data.frame(treatment = factor("Truncating QT", 
                                           levels = c("Multiplicative",
                                                      "Positive Epistasis", 
                                                      "Negative Epistasis", 
                                                      "Directional QT",    
                                                      "Truncating QT",
                                                      "Stabilizing QT" )),
                        jaccards = 0.62, 
                        title = "121 Hap Blocks. T0=0.5",
                        lab = "Empirical Gen 6")

ann_text2 <- data.frame(treatment = factor("Truncating QT", 
                                           levels = c("Multiplicative",
                                                      "Positive Epistasis", 
                                                      "Negative Epistasis", 
                                                      "Directional QT",    
                                                      "Truncating QT",
                                                      "Stabilizing QT" )),
                        jaccards = 0.76, 
                        title = "121 Hap Blocks. T0=0.5",
                        lab = "Empirical Gen 10")

ann_text3 <- data.frame(treatment = factor("Truncating QT", 
                                           levels = c("Multiplicative",
                                                      "Positive Epistasis", 
                                                      "Negative Epistasis", 
                                                      "Directional QT",    
                                                      "Truncating QT",
                                                      "Stabilizing QT" )),
                        jaccards = 0.54, 
                        title = "4977 SNPs. T0=0.5",
                        lab = "Empirical Gen 6")

ann_text4 <- data.frame(treatment = factor("Truncating QT", 
                                           levels = c("Multiplicative",
                                                      "Positive Epistasis", 
                                                      "Negative Epistasis", 
                                                      "Directional QT",    
                                                      "Truncating QT",
                                                      "Stabilizing QT" )),
                        jaccards = 0.65, 
                        title = "4977 SNPs. T0=0.5",
                        lab = "Empirical Gen 10")
totaldata <- totaldata[-which(is.na(totaldata$jaccards)),]
pointrangedata <- totaldata %>% group_by(treatment, generation, title) %>% 
  summarize(jaccard.mean = mean(jaccards),
            jaccard.max = max(jaccards),
            jaccard.min = min(jaccards))

####################################

ggplot(pointrangedata, aes(x=as.factor(treatment))) +
  geom_hline(data = hlinegen10df, aes(yintercept = jaccards), 
             linetype="dashed",
             color = turbo(11)[11], 
             size = .75) + 
  geom_hline(data = hlinegen6df, aes(yintercept = jaccards),
             linetype = "dotdash",
             color = turbo(11)[11], 
             size = .75) +
  geom_pointrange(aes(y=jaccard.mean, 
                      ymin=jaccard.min, 
                      ymax=jaccard.max,
                      color=as.factor(generation)), 
                  position=position_dodge(width = 0.3), 
                  size = 1,
                  fatten = 2,
                  alpha = 0.6)+
  facet_wrap(~ title) + 
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", color="#000000", 
                          size=20, hjust=0.5, vjust=0.5)) + 
  guides(color=guide_legend(title="Generation")) + 
  xlab("Fitness function") + 
  ylab("Jaccard Score") +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_manual(values = viridis(12)[c(3,7)]) + 
  scale_fill_manual(values = viridis(12)[c(3,7)]) +
  geom_text(data = ann_text1,
            aes(y=jaccards),
            label = "Empirical Gen 6", 
            color=turbo(11)[11], 
            size=4 , fontface="bold") +
  geom_text(data = ann_text2,
            aes(y=jaccards),
            label = "Empirical Gen 10", 
            color=turbo(11)[11], 
            size=4 , fontface="bold") + 
  geom_text(data = ann_text3,
            aes(y=jaccards),
            label = "Empirical Gen 6", 
            color=turbo(11)[11], 
            size=4 , fontface="bold") +
  geom_text(data = ann_text4,
            aes(y=jaccards),
            label = "Empirical Gen 10", 
            color=turbo(11)[11], 
            size=4 , fontface="bold") + 
  theme(panel.grid.minor  = element_blank())



###################################
## Supplement: Fitness Functions ##
###################################

nmuts <- 1:121
x <- nmuts / 121


positive <- 1.13^(nmuts) * exp(8 * x^2)
multiplicative <- 1.13^nmuts
negative <- 1.13^(nmuts) * exp(-8 * x^2)
lb <- min(positive, negative, multiplicative)
ub <- max(positive, negative, multiplicative)

positive <- positive / max(positive)
multiplicative <- multiplicative / max(multiplicative)
negative <- negative / max(negative)
###################################
linesdat = data.frame(title = rep(c("121 Haplotype Blocks. Multiplicative Fitness Functions"), times = (3 * 121)),
                      relfit = c(positive, multiplicative, negative),
                      pheno = rep(x, times = 3),
                      fitfun = rep(c("Positive Epistasis", 
                                     "Multiplicative",
                                     "Negative"), each = 121))

###################################
plot(x = x, y = log10(positive), col = turbo(11)[8], type = "l" , xlab = "Phenotype",
     ylab = "log(Relative Fitness)", lwd = 2, ylim = c(-10,0), cex.lab = 1.25, 
     cex.axis = 1.5)
lines(x = x, y = log10(multiplicative), col = turbo(11)[6], lwd = 2)
lines(x = x, y = log10(negative), col = turbo(11)[2], lwd = 2)
title("121 Haplotype Blocks. Multiplicative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("A. Multiplicative", "B. Positive Epistasis", 
#                                  "C. Negative Epistasis", "Initial Distribution"), 
#        col = c("orange", "cyan", "pink", rgb(0,0,0,0.25)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-10, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.25))
abline(v = 0.223786, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)
abline(v = 0.325956, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)

# positive
ci <- seq(from = 0.380902, to = 0.495603, length.out = 3)
y <- rep(-8.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[8])

# multiplicative
ci <- seq(from = 0.331542, to = 0.445071, length.out = 3)
y <- rep(-9, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[6])

# negative
ci <- seq(from = 0.29318, to = 0.402245, length.out = 3)
y <- rep(-9.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[2])
###################################
rangesdat <- data.frame(title = rep(c("121 Haplotype Blocks. Multiplicative Fitness Functions"), times = 4),
                        xmin = c(0.223786, 0.380902, 0.331542, 0.29318),
                        xmax = c(0.325956, 0.495603, 0.445071, 0.402245),
                        y = c(-10, -8.5, -9, -9.5),
                        fitfun = c("Initial", 
                                   "Positive Epistasis", 
                                   "Multiplicative",
                                   "Negative Epistasis"))
vlinesdat <- data.frame(title = rep(c("121 Haplotype Blocks. Multiplicative Fitness Functions"), times = 2),
                        xint = c(0.223786, 0.325956),
                        fitfun = rep("Initial", times = 2))

rm(list = ls()[-c(3,8,10)])

###################################
x <- seq(from = 0, to = 1, length.out = 200)

fmin = 0
fmax = 1
b <- -0.3
s <- 0.1
r <- -15
DirectionalEpistasis = fmin + ((fmax - fmin) / ((1 + s * exp(r *(x + b))) ^(1/s)));


a <- 10
b <- -0.25
tmp = fmin + (fmax - fmin) * (1 - 1 / exp(a * (x + b)));
tmp[tmp < 0.0] = 1e-30;
TruncatingEpistasis = tmp;

mu <- 0.4
std <- 0.07
StabilizingEpistasis = exp(-0.5 * ((x - mu)^2 / std ^ 2));

###################################



linesdatadd <- data.frame(title = rep("121 Haplotype Blocks. Quantative Fitness Functions", times = (3 * 200)),
                          relfit = c(DirectionalEpistasis, TruncatingEpistasis, StabilizingEpistasis),
                          pheno = rep(x, times = 3),
                          fitfun = rep(c("Directional QT",
                                         "Truncating QT",
                                         "Stabilizing QT"), each = 200))

linesdat <- rbind(linesdat, linesdatadd)

###################################

plot(x = x, y = log10(DirectionalEpistasis), type = "l", col = turbo(11)[4], 
     xlab = "Phenotype", ylab = "log(Relative Fitness)", lwd = 2,
     ylim = c(-10,0), cex.lab = 1.25, cex.axis = 1.5)
lines(x = x, y = log10(TruncatingEpistasis), col = turbo(11)[5], lwd = 2)
lines(x = x, y = log10(StabilizingEpistasis), col = turbo(11)[10], lwd = 2)
title("121 Haplotype Blocks. Quantative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("D. Stabiliing Epistasis",
#                                  "E. Directional Epistasis", 
#                                  "F. Truncating Epistasis", 
#                                  "Initial Distribution"), 
#        col = c("red", "green", "Blue", rgb(0,0,0,0.25)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-10, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.25))
abline(v = 0.223786, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)
abline(v = 0.325956, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)

# directional
ci <- seq(from = 0.300193, to = 0.406439, length.out = 3)
y <- rep(-9.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[4])

# stabilizing
ci <- seq(from = 0.297201, to = 0.402461, length.out = 3)
y <- rep(-9, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[10])

# truncating
ci <- seq(from = 0.314602, to = 0.418423, length.out = 3)
y <- rep(-8.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[5])


###################################

rangesdatadd <- data.frame(title = rep(c("121 Haplotype Blocks. Quantative Fitness Functions"), times = 4),
                           xmin = c(0.223786, 0.300193, 0.297201, 0.314602),
                           xmax = c(0.325956, 0.406439, 0.402461, 0.418423),
                           y = c(-10, -9.5, -9, -8.5),
                           fitfun = c("Initial",
                                      "Directional QT",
                                      "Truncating QT",
                                      "Stabilizing QT"))
rangesdat <- rbind(rangesdat, rangesdatadd)


vlinesdatadd <- data.frame(title = rep("121 Haplotype Blocks. Quantative Fitness Functions", times = 2),
                           xint = c(0.223786, 0.325956),
                           fitfun = rep("Initial", times = 2))
vlinesdat <- rbind(vlinesdat, vlinesdatadd)

rm(list = ls()[-c(7,11,18)])

###################################

nmuts <- 1:121
x <- nmuts / 121

a <- 8
positive <- 1.13^(nmuts) * exp(8 * x^2)
multiplicative <- 1.13^nmuts
negative <- 1.13^(nmuts) * exp(-8 * x^2)
lb <- min(positive, negative, multiplicative)
ub <- max(positive, negative, multiplicative)

positive <- positive / max(positive[1:121])
multiplicative <- multiplicative / max(multiplicative[1:121])
negative <- negative / max(negative[1:121])

###################################


linesdatadd <- data.frame(title = rep("4977 SNPs on 121 Haplotype Blocks. Multiplicative Fitness Functions", times = 121),
                          relfit = c(positive, multiplicative, negative),
                          pheno = rep(x, times = 3),
                          fitfun = rep(c("Positive Epistasis", 
                                         "Multiplicative",
                                         "Negative"), each = 121))
linesdat <- rbind(linesdat, linesdatadd)

###################################


plot(x = x, y =  log10(positive), col = turbo(11)[8], type = "l" , xlab = "Phenotype",
     ylab = "log(Relative Fitness)", lwd = 2, 
     xlim = c(0, 1), cex.lab = 1.25, cex.axis = 1.5, ylim = c(-10,0))
lines(x = x, y =  log10(multiplicative), col = turbo(11)[6], lwd = 2)
lines(x = x, y =  log10(negative), col = turbo(11)[2], lwd = 2)
title("4977 SNPs on 121 Haplotype Blocks. Multiplicative Fitness Functions", adj = 0)
# legend("bottomright", legend = c("A. Multiplicative", "B. Positive Epistasis", "
#                                  C. Negative Epistasis", "Initial Distribution"),
#        col = c("orange", "cyan", "pink", rgb(0,0,0,0.25)), lwd = c(2,2,2,8),
#        bg = "white", cex = 1)

ci <- seq(from = 0.392263, to = 0.403815, length.out = 3)
y <- rep(-10, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.25))
abline(v = 0.392263, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)
abline(v = 0.403815, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)

# positive
ci <- seq(from = 0.433199, to = 0.443856, length.out = 3)
y <- rep(-8.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[8])

# mult
ci <- seq(from = 0.432522, to = 0.444504, length.out = 3)
y <- rep(-9, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[6])

# negative
ci <- seq(from = 0.430782,to = 0.441271, length.out = 3)
y <- rep(-9.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[2])

###################################

rangesdatadd <- data.frame(title = rep("4977 SNPs on 121 Haplotype Blocks. Multiplicative Fitness Functions", times = 4),
                           xmin = c(0.392263, 0.433199, 0.432522, 0.430782),
                           xmax = c(0.403815, 0.443856, 0.444504, 0.441271),
                           y = c(-10, -8.5, -9, -9.5),
                           fitfun = c("Initial", 
                                      "Positive Epistasis", 
                                      "Multiplicative",
                                      "Negative Epistasis"))
rangesdat <- rbind(rangesdat, rangesdatadd)


vlinesdatadd <- data.frame(title = rep("4977 SNPs on 121 Haplotype Blocks. Multiplicative Fitness Functions", times = 2),
                           xint = c(0.392263, 0.403815),
                           fitfun = rep("Initial", times = 2))
vlinesdat <- rbind(vlinesdat, vlinesdatadd)

rm(list = ls()[-c(4,10,13)])
###################################

x <- seq(from = 0, to = 1, length.out = 200)

fmin = 0
fmax = 1
b <- -0.41
s <- 0.1
r <- -75
DirectionalEpistasis = fmin + ((fmax - fmin) / ((1 + s * exp(r *(x + b))) ^(1/s)));


a <- 40
b <- -0.395
tmp = fmin + (fmax - fmin) * (1 - 1 / exp(a * (x + b)));
tmp[tmp < 0.0] = 1e-30;
TruncatingEpistasis = tmp;

mu <- 0.435
std <- 0.0175
StabilizingEpistasis = exp(-0.5 * ((x - mu)^2 / std ^ 2));

###################################

linesdatadd <- data.frame(title = rep("4977 SNPs on 121 Haplotype Blocks. Quantative Fitness Functions", times = (3 * 200)),
                          relfit = c(DirectionalEpistasis, TruncatingEpistasis, StabilizingEpistasis),
                          pheno = rep(x, times = 3),
                          fitfun = rep(c("Directional QT",
                                         "Truncating QT",
                                         "Stabilizing QT"), each = 200))

linesdat <- rbind(linesdat, linesdatadd)


###################################

plot(x = x, y =  log10(DirectionalEpistasis), type = "l", col = turbo(11)[4], 
     xlab = "Phenotype", ylab = "log(Relative Fitness)", lwd = 2,
     ylim = c(-10,0), cex.lab = 1.25, cex.axis = 1.5, xlim = c(0, 1))
lines(x = x, y =  log10(TruncatingEpistasis), col = turbo(11)[5], lwd = 2)
lines(x = x, y =  log10(StabilizingEpistasis), col = turbo(11)[10], lwd = 2)
title("4977 SNPs on 121 Haplotype Blocks. Quantative Fitness Functions", adj = 0)
legend("bottomright", legend = c("Multiplicative", 
                                 "Positive Epistasis", 
                                 "Negative Epistasis", 
                                 "Directional QT", 
                                 "Truncating QT", 
                                 "Stabilizing QT",
                                 "Initial Distribution"),
       col = c(turbo(11)[c(6,8,2,4,5,10)], rgb(0,0,0,0.25)), 
       lwd = c(2,2,2,2,2,2,8),
       bg = "white", cex = 1)

ci <- seq(from = 0.392263, to = 0.403815, length.out = 3)
y <- rep(-10, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.25))
abline(v = 0.392263, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)
abline(v = 0.403815, lty = 2, col = rgb(0,0,0,0.25), lwd = 0.5)

# dir
ci <- seq(from = 0.401406, to = 0.412646, length.out = 3)
y <- rep(-9.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[4])

# stab
ci <- seq(from = 0.401266, to = 0.412542, length.out = 3)
y <- rep(-9, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[10])

# trunc
ci <- seq(from = 0.401696, to = 0.412732, length.out = 3)
y <- rep(-8.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = turbo(11, alpha = 0.25)[5])

###################################

rangesdatadd <- data.frame(title = rep(c("4977 SNPs on 121 Haplotype Blocks. Quantative Fitness Functions"), times = 4),
                           xmin = c(0.392263, 0.401406, 0.401266, 0.401696),
                           xmax = c(0.403815, 0.412646, 0.412542, 0.412732),
                           y = c(-10, -9.5, -9, -8.5),
                           fitfun = c("Initial",
                                      "Directional QT",
                                      "Truncating QT",
                                      "Stabilizing QT"))
rangesdat <- rbind(rangesdat, rangesdatadd)


vlinesdatadd <- data.frame(title = rep("4977 SNPs on 121 Haplotype Blocks. Quantative Fitness Functions", times = 2),
                           xint = c(0.392263, 0.403815),
                           fitfun = rep("Initial", times = 2))
vlinesdat <- rbind(vlinesdat, vlinesdatadd)

rm(list = ls()[-c(7,11,18)])
###################################

linesdat$fitfun <- factor(linesdat$fitfun, 
                          levels = c("Multiplicative",
                                     "Positive Epistasis", 
                                     "Negative", 
                                     "Directional QT",    
                                     "Truncating QT",
                                     "Stabilizing QT"), 
                          labels = c("Multiplicative",
                                     "Positive Epistasis", 
                                     "Negative Epistasis", 
                                     "Directional QT",    
                                     "Truncating QT",
                                     "Stabilizing QT"))
rangesdat$fitfun <- factor(rangesdat$fitfun, 
                           levels = c("Multiplicative",
                                      "Positive Epistasis", 
                                      "Negative Epistasis", 
                                      "Directional QT",    
                                      "Truncating QT",
                                      "Stabilizing QT"), 
                           labels = c("Multiplicative",
                                      "Positive Epistasis", 
                                      "Negative Epistasis", 
                                      "Directional QT",    
                                      "Truncating QT",
                                      "Stabilizing QT"))

vlinesdat$fitfun <- factor(vlinesdat$fitfun, 
                           levels = c("Multiplicative",
                                      "Positive Epistasis", 
                                      "Negative Epistasis", 
                                      "Directional QT",    
                                      "Truncating QT",
                                      "Stabilizing QT"), 
                           labels = c("Multiplicative",
                                      "Positive Epistasis", 
                                      "Negative Epistasis", 
                                      "Directional QT",    
                                      "Truncating QT",
                                      "Stabilizing QT"))

###################################

ggplot(linesdat, aes(y=log10(relfit), x=pheno)) + 
  geom_line(aes(colour=fitfun), 
            stat="identity", 
            position="identity", 
            alpha=0.6,
            size = 1.25) + 
  facet_wrap(~ title) + 
  theme_bw() + 
  theme(text=element_text(family="sans", 
                          face="plain", 
                          color="#000000", 
                          size=15, 
                          hjust=0.5, 
                          vjust=0.5)) + 
  xlab("Phenotype") + 
  ylab("log(Relative Fitness)") +
  theme(panel.grid = element_blank()) + 
  guides(color=guide_legend(title="Fitness Function")) +
  coord_cartesian(ylim = c(-10,0), clip = 'on') +
  geom_vline(data = vlinesdat, aes(color = fitfun,
                                   xintercept=xint),
             linetype="dashed") + 
  scale_color_manual("Group",values=c(turbo(11)[c(6,8,2,4,5,10)], rgb(0,0,0,0.15)),
                     limits = c(as.character(sort(unique(linesdat$fitfun))), 
                                "Initial")) + 
  geom_segment(data = rangesdat, aes(color = fitfun,
                                     x = xmin,
                                     xend = xmax,
                                     y = y,
                                     yend = y),
               size = 3,
               lineend = "round")
## unique adaptive genomic architecture

