\name{Manly.var}
\alias{Manly.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variability assessment of Manly mixture model }
\description{
Runs the variability assessment for a Manly mixture model.
}
\usage{
Manly.var(X, model = NULL, conf.CI = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{ dataset matrix (n x p)}
  \item{model }{ Manly mixture model}
  \item{conf.CI }{ confidence level, say 95 percent confidence}
}
\details{
Returns the estimated variance-covariance matrix and confidence intervals for model parameter estimates.
}
\value{
  \item{V}{variance-covariance matrix.}
  \item{CI}{confidence intervals for each parameter.}
}

\seealso{Manly.EM}

\examples{

set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result 
#all skewness parameters set to be 0.1 as the initialization of the EM algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0.1, K, p)

#Run the EM algorithm with Manly mixture model
M.EM <- Manly.EM(X, id.km, la)
     
# Run the variability assessment
Manly.var(X, M.EM, conf.CI = 0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ Variability assessment }% __ONLY ONE__ keyword per line
