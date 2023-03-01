\name{simdataset.skewcwm}
\alias{simdataset.skewcwm}
\title{ Simulating data from cluster-weighted models}
\description{
Generates a dataset of size n given the parameters of a skewed cluster-weighted model simulated by MixSim.skewcwm().}

\usage{simdataset.skewcwm(n, Pi, Mu.x, S2.x, Beta, S2.y, lambda)}


\arguments{
  \item{n }{ sample size}
  \item{Pi }{ vector of mixing propotions [length K]}
  \item{Mu.x }{ maxtrix consisting of components' mean vectors [K x (p-1)]}
  \item{S2.x }{ array consisting of components' covariance matrices for X [(p - 1) x (p-1) x K]}
  \item{Beta }{ matrix of regression coefficients [K x p] (first column represents intercepts)}
  \item{S2.y }{ vector of components' variances for y [length K]}
  \item{lambda }{ value of transformation (skewness) parameter}
 }

\value{
  \item{X }{ simulated data set for X}
  \item{Y }{ simulated data set for y}
  \item{id }{ classification vector}
}

\examples{
set.seed(1)
# controls levels of average skewness and average overlap
A <- MixSim.skewcwm(BarGamma = 0.5, BarOmega = 0.1, K = 3, p = 3)

# similates the data from obtained model
W <- simdataset.skewcwm(n = 100, Pi = A$Pi, Mu.x = A$Mu.x, S2.x = A$S2.x, Beta = A$Beta, S2.y = A$S2.y, lambda = A$lambda)
head(W$X)
head(W$Y)
table(W$id)
}

\keyword{MixSim, skewness, cluster-weighted model}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}

