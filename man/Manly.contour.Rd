\name{Manly.contour}
\alias{Manly.contour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Contour plot for Manly mixture model }
\description{
Provides a contour plot for the fitted data with Manly mixture model.
}
\usage{
Manly.contour(X, var1 = 1, var2 = 2, model = NULL, x.slice = 100, 
y.slice = 100, x.mar = 1, y.mar = 1, col = "lightgrey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{var1 }{ x-axis variable}
  \item{var2 }{ y-axis variable}
  \item{model }{ fitted Manly mixture model}
  \item{x.slice }{ number of slices in the first variable sequence in the contour}
  \item{y.slice }{ number of slices in the second variable sequence in the contour}
  \item{x.mar }{ value to be subtracted/added to the smallest/largest observation in the x-axis}
  \item{y.mar }{ value to be subtracted/added to the smallest/largest observation in the y-axis}
  \item{col }{ color of the contour lines}
  \item{...}{further arguments related to \code{\link{contour}}}
}
\details{
Provides the contour plot for the fitted data by Manly mixture model. This plot is designed for two-dimensional data. 
}

\seealso{Manly.EM}

\examples{

set.seed(123)

K <- 2; p <- 2
X <- as.matrix(faithful)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the EM algorithm for a Manly mixture model based on K-means solution
la <- matrix(0.1, K, p)
B <- Manly.EM(X, id.km, la)

Manly.contour(X, model = B, x.mar = 1, y.mar = 2,
xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", 
ylab = "", nlevels = 10, drawlabels = FALSE, 
lwd = 3.2, col = "lightgrey", pch = 19)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Contour plot }
\keyword{ Manly mixture }% __ONLY ONE__ keyword per line
