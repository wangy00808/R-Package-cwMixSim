\name{MixSim.skewcwm}
\alias{MixSim.skewcwm}
\title{Simulating skewed cluster weighted model}
\description{
Generates skewed cluster-weighted models based on the pre-specified levels of maximum or average skewness as well as maximum and/or average overlap }

\usage{
MixSim.skewcwm(BarOmega = NULL, MaxOmega = NULL, BarGamma = NULL,
      MaxGamma = NULL, K, p, Pi = NULL, Mu.x_lim = c(-5, 5),
      Beta_lim = c(-5, 5), S2.y_lim = c(0, 1), S2.x_lim = c(0, 1),
      N = 500, nM = 200, resN = 100, tol = 1e-06)
}

\arguments{
  \item{BarOmega }{ value of desired average overlap}
  \item{MaxOmega }{ value of desired maximum overlap}
  \item{BarGamma }{ value of desired average skewness among components}
  \item{MaxGamma }{ value of desired maximum skewness among components}
  \item{K }{ number of mixture components}
  \item{p }{ number of dimensions}
  \item{Pi }{ mixing proportions; equal if not specified}
  \item{Mu.x_lim }{ limits of a Uniform distribuion for simulating means of X}
  \item{Beta_lim }{ limits of a Uniform distribuion for simulating coefficients Beta}
  \item{S2.y_lim }{ limits of a Uniform distribuion for simulating initial variance of y}
  \item{S2.x_lim }{ limits of a Uniform distribuion for simulating initial eigenvalues for the variance of X}
  \item{N }{ sample size for estimating the expected value of overlap}
  \item{nM }{ sample size for simulated data to estimate skewness per Mardia's definition}
  \item{resN }{ maximum number of mixture resimulations if desired charasteristics are not met}
  \item{tol }{ desired accuracy (convergence tolerance)}
}

\details{
If 'BarOmega' and 'MaxOmega' are both specified, a mixture will be simulated based on both characteristics (if achievable). If only 'BarOmega' or 'MaxOmega' is specified, a mixture is simulated solely based on the specified characteristic. The user can specify either desired average skewness 'BarGamma' or maximum skewness 'MaxGamma'.
}

\value{
  \item{lambda }{ value of skewness parameter}
  \item{Pi }{ vector of mixing propotions [length K]}
  \item{Mu.x }{ maxtrix consisting of components' mean vectors [K x (p-1)]}
  \item{S2.x }{ array consisting of components' covariance matrices for X [(p - 1) x (p-1) x K]}
  \item{Beta }{ matrix of regression coefficients [K x p] (first column represents intercepts)}
  \item{S2.y }{ vector of components' variances for y [length K]}
  \item{OmegaMap }{ matrix of misclassification probabilities}
  \item{BarOmega }{ value of average overlap}
  \item{MaxOmega }{ value of maximum overlap}
  \item{rcMax }{ row and column numbers for the pair of components producing maximum overlap 'MaxOmega'}
  \item{Gamma }{ skewness value for K components}
  \item{BarGamma }{ value of average skewness}
  \item{MaxGamma }{ value of maximum skewness}
}


\examples{
set.seed(1)
# controls levels of average skewness and average overlap
A <- MixSim.skewcwm(BarGamma = 0.5, BarOmega = 0.1, K = 3, p = 3)

}

\keyword{MixSim, skewness, cluster-weighted model}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}

