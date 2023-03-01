\name{simdataset.skew}
\alias{simdataset.skew}
\title{ Simulating data from transformation mixtures}
\description{
Generates a dataset of size n given the parameters of a transformation model simulated by MixSim.skew().}

\usage{simdataset.skew(n, Pi, Mu, S2, lambda)}


\arguments{
  \item{n }{ sample size}
  \item{Pi }{ vector of mixing proportions (length K)}
  \item{Mu }{ matrix consisting of components' mean vectors (K x p)}
  \item{S2 }{ set of components' covariance matrices (p x p x K)}
  \item{lambda }{ value of transformation (skewness) parameter}
}


\value{
 \item{X }{ simulated data set}
 \item{id }{ classification vector}
}


\examples{

set.seed(1)
# controls levels of average skewness and overlap
A <- MixSim.skew(BarGamma = 0.5, BarOmega = 0.1, K = 3, p = 3)

# simulates the data from the obtained model
W <- simdataset.skew(n = 100, Pi = A$Pi, Mu = A$Mu, S2 = A$S2, lambda = A$lambda)
head(W$X)
table(W$id)

}

\keyword{simulation, skewness}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
