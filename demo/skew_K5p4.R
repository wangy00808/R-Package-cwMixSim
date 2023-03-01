
#####################################################
# 		DEMO (skew)
# Simulating a four-variate five-component transformation
# mixture model
#
# Overlap characteristics: Maximum Omega = 0.20
# Maximum skewness: 3.0
# Mixing proportions [unequal]: 0.5 and 0.5
#####################################################

# Simulate a transformation mixture model with equal spherical covariance matrices
# and mixing proportions greater or equal than 0.1
set.seed(1)
A <- MixSim.skew(MaxGamma = 3.00, MaxOmega = 0.20, K = 5, p = 4, hom = TRUE, sph = TRUE, PiLow = 0.1)

# Show the overlap map and the pair of components with the highest overlap
A$OmegaMap
A$rcMax
A$OmegaMap[A$rcMax[1], A$rcMax[2]] + A$OmegaMap[A$rcMax[2], A$rcMax[1]]

