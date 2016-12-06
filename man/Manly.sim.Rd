\name{Manly.sim}
\alias{Manly.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulates Manly mixture dataset}
\description{
Simulates Manly mixture dataset given the mixture parameters and sample size.
}
\usage{
Manly.sim(n, la, tau, Mu, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n }{ sample size}
  \item{la }{ matrix of transformation parameters (K x p)}
  \item{tau }{ vector of mixing proportions (length K)}
  \item{Mu }{ matrix of mean vectors (K x p)}
  \item{S }{ array of covariance matrices (p x p x K)}
}
\details{
Simulates a Manly mixture dataset. Manly mixture data points are computed from back-transforming Gaussian distributed data points using user-specified transformation parameters 'la'. 
}
\value{
 \item{X }{ the simulated Manly mixture dataset}
 \item{id }{ the simulated membership of the data}
}


\examples{
set.seed(123)

#sets the number of components, dimensionality and sample size
K <- 3
p <- 2
n <- 1000

#sets the parameters to simulate data from 
tau <- c(0.25, 0.3, 0.45)
Mu <- matrix(c(12,4,4,12,4,10),3)
la <- matrix(c(1.2,0.5,1,0.5,0.5,0.7),3)
S <- array(NA, dim = c(p,p,K))
S[,,1] <- matrix(c(4,0,0,4),2)
S[,,2] <- matrix(c(5,-1,-1,3),2)
S[,,3] <- matrix(c(2,-1,-1,2),2)

#use function Manly.sim to simulate dataset with membership
A <- Manly.sim(n, la, tau, Mu, S)

#plot the data
plot(A$X, col = A$id)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ Data simulation }
