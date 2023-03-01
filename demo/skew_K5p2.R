
#####################################################
# 		DEMO (skew)
# Simulating a bivariate five-component transformation
# mixture model
#
# Overlap characteristics:
#   - Maximum Omega = 0.20
#   - Average Omega = 0.05
# Maximum Skewness: 5.0
# Mixing proportions [equal]: 0.5 and 0.5
#####################################################

# Simulate a mixture model with specified characteristics
set.seed(2)
A <- MixSim.skew(MaxGamma = 5.00, MaxOmega = 0.20, BarOmega = 0.05, K = 5, p = 2)

# Construct a contour plot for the simulated mixture   
contour2d(A, method = "skew", mar = c(1, 1, 1, 1), xlim = c(-0.15, 0.35), ylim = c(0.165, 0.28), nlevels = 100, lwd = 1.5, xaxt = "n", yaxt = "n")

# Generate 200 realizations from the simulated mixture model
data.skew <- simdataset.skew(n = 200, Pi = A$Pi, Mu = A$Mu, S = A$S, lambda = A$lambda)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(6, 19, 21, 17, 15)
points(data.skew$X, pch = plotchar[data.skew$id], cex = 1.8, lwd = 2)

# Show the overlap map and the pair of components with the highest overlap
A$OmegaMap
A$rcMax
A$OmegaMap[A$rcMax[1], A$rcMax[2]] + A$OmegaMap[A$rcMax[2], A$rcMax[1]]
