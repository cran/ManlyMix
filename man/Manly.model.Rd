\name{Manly.model}
\alias{Manly.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Manly mixture model }
\description{
Runs all the functionality of a Manly mixture model. 
}
\usage{
Manly.model(X, K = 1:5, Gaussian = FALSE, initial = "k-means", 
nstart = 100, method = "ward.D",  short.iter = 5, 
select = "none", silent = TRUE, plot = FALSE, var1 = NULL, 
var2 = NULL, VarAssess = FALSE, overlap = FALSE, N = 1000, 
tol = 1e-5, max.iter = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{K }{ number of components tested}
  \item{Gaussian }{ whether Gaussian mixture models are run or not}
  \item{initial }{ initialization strategy of the EM algorithm ("k-means" - partition obtained by k-means clustering, "hierarchical" - partition obtained by hierarchical clustering, "emEM" - parameters estimated by the emEM algorithm)}
  \item{nstart }{ number of random starts for the k-means or the emEM initialization}
  \item{method }{ linkage method for the hierarchical initialization}
  \item{short.iter }{ number of short emEM iterations to run}
  \item{select }{ control to run Manly.select or not ("none" - do not run Manly.select , "forward" - run forward selection, "backward" - run backward selection) }
  \item{silent }{ control the output from Manly.select}
  \item{plot }{ control to construct the density or contour plot or not}
  \item{var1 }{ x-axis variable for contour plot or variable for density plot}
  \item{var2 }{ y-axis variable for contour plot}
  \item{VarAssess }{ run the variability assessment of the Manly mixture model or not}
  \item{overlap }{ estimate the overlap of Manly mixture components or not}
  \item{N }{ number of Monte Carlo simulations to run in the Manly.overlap function}  
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{...}{further arguments related to \code{\link{Manly.plot}}}
}
\details{
Wrapper function that incorporates all functionality associated with Manly mixture modeling.
}
\value{
  \item{Model }{ best mixture model obtained}
  \item{VarAssess }{ estimated variance-covariance matrix for model parameter estimates}
  \item{Overlap }{ estimated overlap of Manly mixture components}
}

\seealso{Manly.EM}

\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

Obj <- Manly.model(X, K = 1:5, initial = "emEM", nstart = 1, short.iter = 5)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
