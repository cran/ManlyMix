\name{ManlyMix-package}
\alias{ManlyMix-package}
\docType{package}
\title{
Finite mixture modeling and model-based clustering based on Manly mixture models.
}
\description{
The utility of this package includes finite mixture modeling and model-based clustering based on Manly mixtures as well as forward and backward model selection procedures.
}
\details{
\tabular{ll}{
Package: \tab ManlyMix\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2015-10-03\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}
Function 'Manly.sim' simulates Manly mixture components.

Function 'Manly.EM' runs the EM algorithm for Manly mixture models.

Function 'Manly.select' runs forward and backward model selection procedures.

Function 'Manly.CEM' runs k-means model with Manly transformation.

Function 'Manly.var' produces the variance-covariance matrix of the parameter estimates from Manly mixture model.

Function 'Manly.contour' produces the contour plot of Manly mixture.
}
\author{ Xuwen Zhu and Volodymyr Melnykov.

Maintainer: Xuwen Zhu <xzhu20@crimson.ua.edu>
}

\examples{

set.seed(123)

K <- 3; p <- 4
X <- as.matrix(iris[,-5])
id.true <- rep(1:K, each = 50)

# Obtain initial memberships based on the K-means algorithm
id.km <- kmeans(X, K)$cluster

# Run the CEM algorithm for Manly K-means model
la <- matrix(0.1, K, p)
C <- Manly.Kmeans(X, id = id.km, la = la)

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
