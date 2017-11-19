set.seed(123)

#Use iris dataset
K <- 3; p <- 4
X <- as.matrix(iris[,-5])

#Use k-means clustering result and all skewness parameters set to be 0.1 as the initialization of the EM algorithm  
id.km <- kmeans(X, K)$cluster
la <- matrix(0.1, K, p)

#Run the EM algorithm with Manly mixture model
M.EM <- Manly.EM(X, id.km, la)
     
# Run the variability assessment
result <- Manly.var(X, M.EM, conf.CI = 0.95)
print(result$V)
print(result$CI)