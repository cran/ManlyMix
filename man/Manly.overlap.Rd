\name{Manly.overlap}
\alias{Manly.overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimates the overlap for a Manly mixture}
\description{
Estimates the pairwise overlap matrix for a Manly mixture by simulating samples based on user-specified parameters. 
}
\usage{
Manly.overlap(tau, Mu, S, la, N = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{la }{ matrix of transformation parameters (K x p)}
  \item{tau }{ vector of mixing proportions (length K)}
  \item{Mu }{ matrix of mean vectors (K x p)}
  \item{S }{ array of covariance matrices (p x p x K)}
  \item{N }{ number of samples simulated}
}

\details{
Estimates the pairwise overlap matrix for a Manly mixture. Overlap is defined as sum of two misclassification probabilities.}
\value{
  \item{OmegaMap}{matrix of misclassification probabilities (K x K); OmegaMap[i,j] is the probability that X coming from the i-th component is classified to the j-th component.}
  \item{BarOmega}{value of average overlap.}
  \item{MaxOmega}{value of maximum overlap.}
}
\examples{
set.seed(123)
#sets the number of components, dimensionality and sample size
K <- 3
p <- 2

#sets the mixture parameters
tau <- c(0.25, 0.3, 0.45)
Mu <- matrix(c(4.5,4,5,7,8,5.5),3)
la <- matrix(c(0.2,0.5,0.3,0.25,0.35,0.4),3)
S <- array(NA, dim = c(p,p,K))
S[,,1] <- matrix(c(0.4,0,0,0.4),2)
S[,,2] <- matrix(c(1,-0.2,-0.2,0.6),2)
S[,,3] <- matrix(c(2,-1,-1,2),2)

#computes the overlap
A <- Manly.overlap(tau, Mu, S, la)
print(A)
}
\references{
Maitra, R. and Melnykov, V. (2010) ``Simulating data to study performance of finite mixture modeling and clustering algorithms'', \emph{Journal of Computational and Graphical Statistics}, 2:19, 354-376.

Melnykov, V., Chen, W.-C., and Maitra, R. (2012) ``MixSim: An R Package for Simulating Data to Study Performance of Clustering Algorithms'', \emph{Journal of Statistical Software}, 51:12, 1-25.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ Overlap }
