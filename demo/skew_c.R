
#####################################################
# 		DEMO (skew)
# Simulating a bivariate two-component transformation
# mixture model
#
# Overlap characteristics: Omega = 0.05
# Average Skewness: 2.0
# Mixing proportions [equal]: 0.5 and 0.5
#####################################################

# Simulate a mixture model with specified characteristics
set.seed(1)
A <- MixSim.skew(BarGamma = 2.0, BarOmega = 0.05, K = 2, p = 2)

# Construct a contour plot for the simulated mixture   
contour2d(A, method = "skew", mar = c(1, 1, 1, 1), xlim = c(0, 0.35), ylim = c(0, 0.38), nlevels = 50, lwd = 1.5, xaxt = "n", yaxt = "n")

# Generate 100 realizations from the simulated mixture model
data.skew <- simdataset.skew(n = 100, Pi = A$Pi, Mu = A$Mu, S = A$S, lambda = A$lambda)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(17, 21)
points(data.skew$X, pch = plotchar[data.skew$id], cex = 1.8, lwd = 2)
