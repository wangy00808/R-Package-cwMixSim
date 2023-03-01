
######################################################
# 		DEMO (cwm)
# Simulating a bivariate five-component
# cluster-weighted model
#
# Overlap characteristics:
#   - Average Overlap = 0.01
#   - Maximum Overlap = 0.04
# Mixing proportions [equal]: 0.25, 0.25, 0.25, and 0.25
#######################################################

# Simulate a mixture model with specified characteristics
set.seed(1)
A <- MixSim.cwm(BarOmega = 0.01, MaxOmega = 0.04, K = 5, p = 2)

# Construct a contour plot for the simulated mixture   
contour2d(A, method = "cwm", mar = c(1, 1, 1, 1), xlim = c(-5.57, 4.5), ylim = c(-26, 25), nlevels = 100, lwd = 1.5, xaxt = "n", yaxt = "n")

# Generate 200 realizations from the simulated mixture model
data.cwm <- simdataset.cwm(n = 200, Pi = A$Pi, Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(17, 19, 21, 6, 15)
points(data.cwm$X, data.cwm$Y, pch = plotchar[data.cwm$id], cex = 1.8, lwd = 2)
