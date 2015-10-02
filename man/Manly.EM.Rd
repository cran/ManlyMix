\name{Manly.EM}
\alias{Manly.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for Manly mixture model }
\description{
Runs the EM algorithm for a Manly mixture model with specified initial membership and transformation parameters.
}
\usage{
Manly.EM(X, id = NULL, la = NULL, tau = NULL, Mu = NULL, S = NULL, 
tol = 1e-5, max.iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{id }{ initial membership vector (length n)}
  \item{la }{ initial transformation parameters (K x p)}
  \item{tau }{ initial vector of mixing proportions (length K)}
  \item{Mu }{ initial matrix of mean vectors (K x p)}
  \item{S }{ initial array of covariance matrices (p x p x K)}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
}
\details{
Runs the EM algorith for a Manly mixture model for a provided dataset. Manly mixture model assumes that a multivariate Manly transformation applied to each component allows to reach near-normality. 
A user has a choice to specify either initial id vector 'id' and transformation parameters 'la' or initial mode parameters 'la', 'tau', 'Mu', and 'S'. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a Gaussian mixtuire model. If some transformation parameters have to be excluded from the consideration, in the corresponding positions of matrix 'la', the user has to specify value 0.
Notation: n - sample size, p - dimensionality of the dataset X, K - number of mixture components.
}
\value{
  \item{la }{ matrix of the estimated transformation parameters (K x p)}
  \item{tau }{ vector of mixing proportions (length K)}
  \item{Mu }{ matrix of the estimated mean vectors (K x p)}
  \item{S }{ array of the estimated covariance matrices (p x p x K)}
  \item{gamma }{ matrix of posterior probabilities (n x K)}
  \item{id }{ estimated membership vector (length n)}
  \item{ll }{ log likelihood value}
  \item{bic }{ Bayesian Information Criterion}
 \item{iter }{ number of EM iterations run}
 \item{flag }{ convergence flag (0 - success, 1 - failure)}
}

\seealso{Manly.select}
\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the EM algorithm for a Gaussian mixture model based on K-means solution
A <- Manly.EM(X, id.km)
id.Gauss <- A$id

table(id.true, id.Gauss)

# Run the EM algorithm for a Manly mixture model based on Gaussian mixture solution
la <- matrix(0.1, K, p)
B <- Manly.EM(X, id.Gauss, la)
id.Manly <- B$id

table(id.true, id.Manly)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ EM algorithm }
\keyword{ Forward and backward selection }% __ONLY ONE__ keyword per line
