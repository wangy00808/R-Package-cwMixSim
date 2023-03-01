
#####################################################
# 		DEMO 1 (skew-cwm)
# Simulating a bivariate two-component transformation
# cluster-weighted model
#
# Overlap characteristics: Omega = 0.10
# Maximum skewness: 1.00
# Mixing proportions [equal]: 0.5 and 0.5
#####################################################

# Simulate a mixture model with specified characteristics
set.seed(10)
A <- MixSim.skewcwm(MaxGamma = 1.00, BarOmega = 0.10, K = 2, p = 2)

# Construct a contour plot for the simulated mixture 
contour2d(A, method = "skewcwm", mar = c(1, 1, 1, 1), xlim = c(-4, 2), ylim = c(-3.5, 2), nlevels = 55, lwd = 1.5, xaxt = "n", yaxt = "n")

# Generate 100 realizations from the simulated mixture model
data.skewcwm <- simdataset.skewcwm(n = 100, Pi = A$Pi, Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y, lambda = A$lambda)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(17, 21)
points(data.skewcwm$X, data.skewcwm$Y, pch = plotchar[data.skewcwm$id], cex = 1.8, lwd = 2)

