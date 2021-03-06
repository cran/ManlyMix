set.seed(123)

#Application to dataset bankruptcy
data("bankruptcy")
X <- as.matrix(bankruptcy[,2:3])
id <- as.numeric(bankruptcy[,1]) + 1


n <- dim(X)[1]
p <- dim(X)[2]
K <- max(id) 

#run the traditional K-means algorithm
M.K <- kmeans(X, K)
id.km <- M.K$cluster
ClassAgree(id.km, id)

#run the Manly K-means algorithm
M.MK <- Manly.Kmeans(X, id = id.km, la = matrix(0.1, K, p))
ClassAgree(M.MK$id, id)

#run Gaussian mixture model
M.Gauss <- Manly.EM(X, id = id.km, la = matrix(0, K, p))
ClassAgree(M.Gauss$id, id)

#run the EM algorithm
M.EM <- Manly.EM(X, id = id.km, la = matrix(0.1, K, p))
ClassAgree(M.EM$id, id)

#run the forward selection
M.F <- Manly.select(X, M.Gauss, method = "forward", silent = TRUE)
ClassAgree(M.F$id, id)

#run the backward algorithm
M.B <- Manly.select(X, M.EM, method = "backward", silent = TRUE)
ClassAgree(M.B$id, id)


#plot the results
M.K$id <- id.km
M.K$tau <- rep(1/K, K)
M.K$Mu <- M.K$centers
M.K$la <- matrix(0, K, p)
M.K$S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	diag(M.K$S[,,k]) <- M.K$tot.withinss / n / p
}
Manly.plot(X, var1 = 1, var2 = 2, model = M.K, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 4, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

M.MK$tau <- rep(1/K, K)
S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	diag(S[,,k]) <- M.MK$S[k]
}
M.MK$S <- S
Manly.plot(X, var1 = 1, var2 = 2, model = M.MK, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.plot(X, var1 = 1, var2 = 2, model = M.Gauss, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 100, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.plot(X, var1 = 1, var2 = 2, model = M.EM, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.plot(X, var1 = 1, var2 = 2, model = M.F, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)

Manly.plot(X, var1 = 1, var2 = 2, model = M.B, x.mar = 60, y.mar = 50, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab = "", nlevels = 30, drawlabels = FALSE, lwd = 3.2, col = "lightgrey", pch = 19)
