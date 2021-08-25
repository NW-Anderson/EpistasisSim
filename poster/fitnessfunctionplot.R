################
## Quantative ##
################
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


plot(x = x, y = log(DirectionalEpistasis), type = "l", col = "green", 
     xlab = "Phenotype", ylab = "log(Relative Fitness)", lwd = 2,
     ylim = c(-23,0), cex.lab = 1.25, cex.axis = 1.5)
lines(x = x, y = log(TruncatingEpistasis), col = "blue", lwd = 2)
lines(x = x, y = log(StabilizingEpistasis), col = "red", lwd = 2)
title("Quantative Fitness Functions", adj = 0)
legend("bottomright", legend = c("D. Stabiliing Epistasis","E. Directional Epistasis", "F. Truncating Epistasis", "Initial Distribution"), 
       col = c("red", "green", "Blue", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
       bg = "white", cex = 1.5)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))

ci <- seq(from = 0.300193, to = 0.406439, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(0,1,0,0.15))

ci <- seq(from = 0.297201, to = 0.402461, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(1,0,0,0.15))

ci <- seq(from = 0.314602, to = 0.418423, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(0,0,1,0.15))

####################
## Multiplicative ##
####################
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


plot(x = x, y = log(positive), col = "cyan", type = "l" , xlab = "Phenotype",
     ylab = "log(Relative Fitness)", lwd = 2, ylim = c(-23,0), cex.lab = 1.25, cex.axis = 1.5)
lines(x = x, y = log(multiplicative), col = 'orange', lwd = 2)
lines(x = x, y = log(negative), col = 'pink', lwd = 2)
title("Multiplicative Fitness Functions", adj = 0)
legend("bottomright", legend = c("A. Multiplicative", "B. Positive Epistasis", "C. Negative Epistasis", "Initial Distribution"), 
       col = c("orange", "cyan", "pink", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
       bg = "white", cex = 1.5)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(-23, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))

ci <- seq(from = 0.380902, to = 0.495603, length.out = 3)
y <- rep(-21.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,1,1,0.15))

ci <- seq(from = 0.331542, to = 0.445071, length.out = 3)
y <- rep(-22, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(1, .647, 0,0.15))

ci <- seq(from = 0.29318, to = 0.402245, length.out = 3)
y <- rep(-22.5, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(1, .753, .796,0.15))

