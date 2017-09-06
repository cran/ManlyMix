\name{Manly.Kmeans}
\alias{Manly.Kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ k-means algorithm with Manly transformation }
\description{
Runs the CEM algorithm for k-means clustering with specified initial membership and transformation parameters.
}
\usage{
Manly.Kmeans(X, id = NULL, la = NULL, Mu = NULL, S = NULL, 
initial = "k-means", K = NULL, nstart = 100, method = "ward.D", 
tol = 1e-5, max.iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{id }{ initial membership vector (length n)}
  \item{la }{ initial transformation parameters (K x p)}
  \item{Mu }{ initial matrix of mean vectors (K x p)}
  \item{S }{ initial vector of variances (K)}
   \item{initial }{ initialization strategy of the EM algorithm ("k-means" - partition obtained by k-means clustering, "hierarchical" - partition obtained by hierarchical clustering)}
  \item{K }{ number of clusters for the k-means initialization}
  \item{nstart }{ number of random starts for the k-means initialization}
  \item{method }{ linkage method for the hierarchical initialization}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
}
\details{
Runs the CEM algorithm for k-means clustering with Manly transformation for a provided dataset. The model assumes that a multivariate Manly transformation applied to each component allows to reach near-normality. 
A user has a choice to specify either initial id vector 'id' and transformation parameters 'la' or initial mode parameters 'la', 'Mu', and 'S'. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a k-means model. If some transformation parameters have to be excluded from the consideration, in the corresponding positions of matrix 'la', the user has to specify value 0.
Notation: n - sample size, p - dimensionality of the dataset X, K - number of mixture components.
}
\value{
  \item{la }{ matrix of the estimated transformation parameters (K x p)}
  \item{Mu }{ matrix of the estimated mean vectors (K x p)}
  \item{S }{ array of the estimated covariance matrices (K)}
  \item{id }{ estimated membership vector (length n)}
 \item{iter }{ number of EM iterations run}
 \item{flag }{ convergence flag (0 - success, 1 - failure)}
}

\seealso{Manly.EM}
\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

# Obtain initial memberships based on the traditional K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the CEM algorithm for k-means with Manly transformation based on traditional k-means solution
la <- matrix(0.1, K, p)
B <- Manly.Kmeans(X, id.km, la)
id.Manly <- B$id

ClassAgree(id.Manly, id.true)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ k-means }
\keyword{ Manly transformation }
\keyword{ CEM algorithm }% __ONLY ONE__ keyword per line
