x <- seq(from = 0, to = 1, length.out = 200)

fmin = 0
fmax = 1
b <- -0.3
s <- 0.1
r <- -15
DirectionalEpistasis = fmin + ((fmax - fmin) / ((1 + s * exp(r *	(x + b))) ^(1/s)));


a <- 10
b <- -0.25
tmp = fmin + (fmax - fmin) * (1 - 1 / exp(a * (x + b)));
tmp[tmp < 0.0] = 0.0;
TruncatingEpistasis = tmp;

mu <- 0.4
std <- 0.07
StabilizingEpistasis = exp(-0.5 * ((x - mu)^2 / std ^ 2));

plot(x = x, y = DirectionalEpistasis, type = "l", col = "red", xlab = "Phenotype", ylab = "Relative Fitness", lwd = 2)
title("Quantative Fitness Functions", adj = 0)
legend("topleft", legend = c("Directional", "Truncating", "Stabiliing", "Initial 95% CI"), 
       col = c("red", "green", "Blue", rgb(0,0,0,0.15)), lwd = c(2,2,2,8),
       bg = "white")
lines(x = x, y = TruncatingEpistasis, col = "green", lwd = 2)
lines(x = x, y = StabilizingEpistasis, col = "Blue", lwd = 2)

ci <- seq(from = 0.223786, to = 0.325956, length.out = 3)
y <- rep(0.0, times = length(ci))
lines(x=ci, y = y, lwd = 8, col = rgb(0,0,0,0.15))

ci <- seq(from = 0.300193, to = 0.406439, length.out = 3)
y <- rep(0.04, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(0,1,0,0.15))

ci <- seq(from = 0.297201, to = 0.402461, length.out = 3)
y <- rep(0.02, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(1,0,0,0.15))

ci <- seq(from = 0.314602, to = 0.418423, length.out = 3)
y <- rep(0.06, times = length(ci))
lines(x = ci, y = y, lwd = 8, col = rgb(0,0,1,0.15))
