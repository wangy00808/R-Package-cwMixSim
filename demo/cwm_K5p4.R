
######################################################
# 		DEMO (cwm)
# Simulating a four-variate five-component
# cluster-weighted model
#
# Overlap characteristics:
#   - Maximum Overlap = 0.01
#   - Average Overlap is not specified
# Mixing proportions [equal]: 0.20, 0.20, 0.20, 0.20, and 0.20
# At the data simulation step, mixture proportions are changed
#######################################################

# Simulate a mixture model with specified characteristics
set.seed(1)
A <- MixSim.cwm(MaxOmega = 0.01, N = 1000, K = 5, p = 4)
# Present the overlap map and the pair of components with highest overlap
A$OmegaMap
A$rcMax
A$OmegaMap[A$rcMax[1], A$rcMax[2]] + A$OmegaMap[A$rcMax[2], A$rcMax[1]]

# Generate 200 realizations from the simulated mixture model
# with mixing proportions changed to be 0.1, 0.1, 0.1, 0.1, and 0.6
data.cwm <- simdataset.cwm(n = 200, Pi = c(0.1, 0.1, 0.1, 0.1, 0.6),
Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y)


