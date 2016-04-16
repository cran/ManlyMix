\name{Manly.select}
\alias{Manly.select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Manly transformation selection }
\description{
Runs forward or backward model selection procedures for finding the optimal model in terms of BIC.
}
\usage{
Manly.select(X, model, method, tol = 1e-5, max.iter = 1000, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{model }{ list containing parameters of the initial model}
  \item{method }{ model selection method (options 'forward' and 'backward')}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{silent }{ output control}
}
\details{
Runs Manly forward and backward model selection procedures for a provided dataset. Forward and backward selection can be started from any ManlyMix object provided in 'model'. Manly transformation parameters are provided in matrix 'model$la'. If some transformations are not needed for specific components, zeros have to be specified in corresponding poisition. When all transformation parameters are set to zero, Manly mixture model degenerates to a Gaussian mixture model.  
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

\seealso{Manly.EM}
\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the EM algorithm for a Gaussian mixture model based on K-means solution
G <- Manly.EM(X, id = id.km)
id.G <- G$id

# Run FORWARD SELECTION ('silent' is on)
F <- Manly.select(X, model = G, method = "forward", silent = TRUE)

# Run the EM algorithm for a full Manly mixture model based on Gaussian mixture solution
la <- matrix(0.1, K, p)
M <- Manly.EM(X, id = id.G, la = la)

# Run BACKWARD SELECTION ('silent' is off)
B <- Manly.select(X, model = M, method = "backward")

BICs <- c(G$bic, M$bic, F$bic, B$bic)
names(BICs) <- c("Gaussian", "Manly", "Forward", "Backward")
BICs

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ EM algorithm }
\keyword{ Forward and backward selection }% __ONLY ONE__ keyword per line
