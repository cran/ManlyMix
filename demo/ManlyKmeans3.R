set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result and all skewness parameters set to be 0 as the initialization of the Sigma K-means algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0, K, p)

#Run the Sigma K-means algorithm
M.MK3 <- Manly.Kmeans(X, id.km, la)
print(M.MK3)