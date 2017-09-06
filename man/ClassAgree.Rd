\name{ClassAgree}
\alias{ClassAgree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the confusion matrix and number of misclassifications}
\description{
Calculates the confusion matrix and number of misclassifications.
}
\usage{
ClassAgree(est.id, trueid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est.id }{ estimated membership vector}
  \item{trueid }{ true membership vector}
}
\value{
  \item{ClassificationTable }{ confusion table between true and estimated partitions}
  \item{MisclassificationNum }{ number of misclassifications}
}

\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

ClassAgree(id.km, id.true)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manly mixture }
\keyword{ EM algorithm }
\keyword{ Forward and backward selection }% __ONLY ONE__ keyword per line
