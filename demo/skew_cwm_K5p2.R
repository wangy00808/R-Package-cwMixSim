
#####################################################
# 		DEMO (skew-cwm)
# Simulating a bivariate five-component transformation
# cluster-weighted model
#
# Overlap characteristics: Maximum Omega = 0.05
# Average skewness: 0.50
# Mixing proportions [equal]: 0.5 and 0.5
#####################################################

# Simulate a mixture model with specified characteristics
set.seed(11)
A <- MixSim.skewcwm(BarGamma = 0.50, MaxOmega = 0.05, K = 5, p = 2)

# Construct a contour plot for the simulated mixture   
contour2d(A, method = "skewcwm", mar = c(1, 1, 1, 1), xlim = c(-7, 2), ylim = c(-15, 80), nlevels = 250, lwd = 1.5)

# Generate 200 realizations from the simulated mixture model
data.skewcwm <- simdataset.skewcwm(n = 200, Pi = A$Pi, Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y, lambda = A$lambda)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(6, 21, 17, 19, 15)
points(data.skewcwm$X, data.skewcwm$Y, pch = plotchar[data.skewcwm$id], cex = 1.8, lwd = 2)

