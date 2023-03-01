\name{MixSim.skew}
\alias{MixSim.skew}
\title{Simulating transformation mixture}
\description{
Generates transformation mixture models (with skewness) based on the pre-specified levels of maximum or average skewness as well as maximum and/or average overlap }

\usage{
MixSim.skew(BarOmega = NULL, MaxOmega = NULL, BarGamma = NULL,
            MaxGamma = NULL, K, p, sph = FALSE, hom = FALSE,
            ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), nM = 1000,
            resN = 100, tol = 1e-06, eps = 1e-06, lim = 1e06)}

\arguments{
  \item{BarOmega }{ value of desired average overlap}
  \item{MaxOmega }{ value of desired maximum overlap}
  \item{BarGamma }{ value of desired average skewness among components}
  \item{MaxGamma }{ value of desired maximum skewness among components}
  \item{K }{ number of mixture components}
  \item{p }{ number of dimensions}
  \item{sph }{ covariance matrix structure (FALSE - non-spherical [default], TRUE - spherical)}
  \item{hom }{ heteroscedastic or homoscedastic components (FALSE - heteroscedastic [default], TRUE - homoscedastic)}
  \item{ecc }{ maximum eccentricity allowed (controls elongation of components through the ratio of maximum and minimum covariance matrix eigenvalues)}
  \item{PiLow }{ value of the lowwer bound for mixing proportions (if PiLow is not reachable with respect to K, equal proportions are taken; PiLow = 1.0 implies equal proportions [default]}
  \item{int }{ mean vectors are simulated uniformly on a hypercube with sides specified by int = (lower.bound, upper.bound) [(0, 1) by default]}
  \item{nM }{ sample size for simulated data to estimate skewness per Mardia's definition}
  \item{resN }{ maximum number of mixture resimulations if desired charasteristics are not met}
  \item{tol }{ desired accuracy (convergence tolerance)}
  \item{eps }{ error bound for overlap computations}
  \item{lim }{ maximum number of integration terms (per Davies (1980))}
}

\details{
If 'BarOmega' and 'MaxOmega' are both specified, a mixture will be simulated based on both characteristics (if achievable). If only 'BarOmega' or 'MaxOmega' is specified, a mixture is simulated solely based on the specified characteristic. For details on parameters 'PiLow', 'ecc', 'hom', 'sph', etc. we refer the reader to the JSS paper by Melnykov et al. (2012). The user can specify either desired average skewness 'BarGamma' or maximum skewness 'MaxGamma'.
}

\value{
  \item{lambda }{ value of skewness parameter}
  \item{Pi }{ vector of mixing propotions}
  \item{Mu }{ maxtrix consisting of components' mean vectors K*p}
  \item{S2 }{ array consisting of components' covariance matrices p*p*K}
  \item{OmegaMap }{ matrix of misclassification probabilities}
  \item{BarOmega }{ value of average overlap}
  \item{MaxOmega }{ value of maximum overlap}
  \item{rcMax }{ row and column numbers for the pair of components producing maximum overlap 'MaxOmega'}
  \item{Gamma }{ skewness value for K components}
  \item{BarGamma }{ value of average skewness}
  \item{MaxGamma }{ value of maximum skewness}
}


\examples{
set.seed(123)

# controls levels of average skewness and average overlap
A <- MixSim.skew(BarGamma = 0.5, BarOmega = 0.1, K = 3, p = 3)

# controls levels of maximum skewness and maximum overlap
B <- MixSim.skew(MaxGamma = 1, MaxOmega = 0.2, K = 3, p = 3)

# controls level of average skewness as well as average and maximum overlap
C <- MixSim.skew(BarGamma = 2, BarOmega = 0.2, MaxOmega = 0.4, K = 3, p = 3)

}

\keyword{MixSim, skewness}

\references{
Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', Journal of Statistical Software, 51:12, 1-25.
}
