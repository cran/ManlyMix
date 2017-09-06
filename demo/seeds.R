set.seed(14)

#Application to dataset seeds
data("seeds")
X <- as.matrix(seeds[,1:7])
id <- as.numeric(seeds[,8])


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

