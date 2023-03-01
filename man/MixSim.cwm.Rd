\name{MixSim.cwm}
\alias{MixSim.cwm}
\title{Simulating cluster weighted model}
\description{
Generates cluster-weighted models based on pre-specified levels of maximum and/or average overlap }

\usage{
MixSim.cwm(BarOmega = NULL, MaxOmega = NULL, K, p, Pi = NULL,
        Mu.x_lim = c(-5, 5), Beta_lim = c(-5, 5), S2.y_lim = c(0, 1),
        S2.x_lim = c(0, 1), N = 500, resN = 10, tol = 1e-06)
}

\arguments{
  \item{BarOmega }{ value of desired average overlap}
  \item{MaxOmega }{ value of desired maximum overlap}
  \item{K }{ number of mixture components}
  \item{p }{ number of dimensions in the joint distribution; i.e., X is (p-1)-dimensional, y is univariate}
  \item{Pi }{ mixing proportions; equal if not specified}
  \item{Mu.x_lim }{ limits of a Uniform distribuion for simulating means of X}
  \item{Beta_lim }{ limits of a Uniform distribuion for simulating coefficients Beta}
  \item{S2.y_lim }{ limits of a Uniform distribuion for simulating initial variance of y}
  \item{S2.x_lim }{ limits of a Uniform distribuion for simulating initial eigenvalues for the variance of X}
  \item{N }{ sample size for estimating the expected value of overlap}
  \item{resN }{ maximum number of mixture resimulations if desired charasteristics are not met}
 \item{tol }{ desired accuracy (convergence tolerance)}
}

\details{
If 'BarOmega' and 'MaxOmega' are both specified, a mixture will be simulated based on both characteristics (if achievable). If only 'BarOmega' or 'MaxOmega' is specified, a mixture is simulated solely based on the specified characteristic.
}

\value{
  \item{Pi }{ vector of mixing propotions [length K]}
  \item{Mu.x }{ maxtrix consisting of components' mean vectors [K x (p-1)]}
  \item{S2.x }{ array consisting of components' covariance matrices for X [(p - 1) x (p-1) x K]}
  \item{Beta }{ matrix of regression coefficients [K x p] (first column represents intercepts)}
  \item{S2.y }{ vector of components' variances for y [length K]}
  \item{OmegaMap }{ matrix of misclassification probabilities}
  \item{BarOmega }{ value of average overlap}
  \item{MaxOmega }{ value of maximum overlap}
  \item{rcMax }{ row and column numbers for the pair of components producing maximum overlap 'MaxOmega'}
}

\examples{
set.seed(1)

# controls the level of average overlap
A <- MixSim.cwm(BarOmega = 0.01, K = 3, p = 3)

# controls the level of maximum overlap
B <- MixSim.cwm(MaxOmega = 0.1, K = 3, p = 3)

}

\keyword{MixSim, cluster-weighted model}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
