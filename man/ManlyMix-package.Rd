\name{ManlyMix-package}
\alias{ManlyMix-package}
\alias{ManlyMix-package}
\docType{package}
\title{
Model-based clustering with Manly mixture models 
}
\description{
The package runs finite mixture modeling and model-based clustering using Manly mixture, Manly forward and Manly backward models
}
\details{
\tabular{ll}{
Package: \tab ManlyMix\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}
Function 'Manly.EM' runs the EM algorithm for finite mixture models with Manly mixture components.

}
\author{ Xuwen Zhu, Volodymyr Melnykov

Maintainer: Xuwen Zhu <xzhu20@crimson.ua.edu>
}
\references{
Zhu, X. and Melnykov, V. Manly transformation in finite mixture modeling, under revision in Computational Statistics & Data Analysis.

}
\examples{

require(ManlyMix)
set.seed(777)

X <- as.matrix(faithful)

#get initial id from K-means
id <- kmeans(X,2)$cluster

#run Manly.EM
la <- matrix(0.1, 2, 2)
B <- Manly.EM(X, id, la)

}
