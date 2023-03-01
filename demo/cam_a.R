
###############################################
# 		DEMO (cwm)
# Simulating a bivariate two-component
# cluster-weighted model
#
# Overlap characteristics [high]: Omega = 0.01
# Mixing proportions [unequal]: 0.5 and 0.5
###############################################

# Simulate a mixture model with specified characteristics
set.seed(1)
A <- MixSim.cwm(BarOmega = 0.01, K = 2, p = 2)

# Construct a contour plot for the simulated mixture   
contour2d(A, method = "cwm", mar = c(1, 1, 1, 1), nlevels = 50, lwd = 1.5, xaxt = "n", yaxt = "n")

# Generate 100 realizations from the simulated mixture model
data.cwm <- simdataset.cwm(n = 100, Pi = A$Pi, Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y)

# Add the generated data points to the contour
# plotting characters represent the true id of the points
plotchar <- c(17, 21)
points(data.cwm$X, data.cwm$Y, pch = plotchar[data.cwm$id], cex = 1.8, lwd = 2)

