\name{cwMixSim-package}
\alias{cwMixSim-package}
\docType{package}
\title{
Simulating skewed mixture models and cluster-weighted models to study performance of clustering algorithms
}
\description{
Simulates transformation mixture models with pre-specified skewness and overlap. Simulates  cluster-weighted models with pre-specified skewness and overlap. Simulates data from generated models. Provides visualization for 2-dimensional cases in the form of a contour plot.}

\details{
Function MixSim.skew() - simulates transformation mixture models with pre-specified skewness and overlap.
Function MixSim.cwm() - simulates cluster-weighted models with pre-specified overlap.
Function MixSim.skewcwm() - simulates skewed cluster-weighted models with pre-specified skewness and overlap.
Function simdataset.skew() - generates data from transformation mixture models simulated by MixSim.skew().
Function simdataset.cwm() - generates data from cluster-weighted models simulated by MixSim.cwm().
Function simdataset.skewcwm() - generates data from skewed cluster-weighted models simulated by MixSim.skewcwm().
Function contour2d() - constructs a contour plot for bivariate models simulated from MixSim.skew(), MixSim.cwm(), or MixSim.skewcwm().
}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}

