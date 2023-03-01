
#####################################################
# 		DEMO (skew-cwm)
# Simulating a four-variate five-component transformation
# cluster-weighted model
#
# Overlap characteristics: Maximum Omega = 0.05
# Average skewness: 0.50
# Mixing proportions [equal]: 0.5 and 0.5
#####################################################

# Simulate a mixture model with specified characteristics
set.seed(1)
A <- MixSim.skewcwm(BarGamma = 0.50, BarOmega = 0.05, K = 5, p = 4)

# Show the overlap map and the pair of components with the highest overlap
A$OmegaMap
A$rcMax
A$OmegaMap[A$rcMax[1], A$rcMax[2]] + A$OmegaMap[A$rcMax[2], A$rcMax[1]]
