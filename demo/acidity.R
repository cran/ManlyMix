set.seed(123)

#Application to dataset AIS
data("acidity")

K <- 2
p <- 1
X <- acidity
n <- length(X)
#run the traditional K-means algorithm
M.K <- kmeans(X, K)
id.km <- M.K$cluster


#run the Manly K-means algorithm
M.MK <- Manly.Kmeans(X, id = id.km, la = matrix(0.1, K, p))


#run Gaussian mixture model
M.Gauss <- Manly.EM(X, id = id.km, la = matrix(0, K, p))


#run the EM algorithm
M.EM <- Manly.EM(X, id = id.km, la = matrix(0.1, K, p))

#run the forward selection
M.F <- Manly.select(X, M.Gauss, method = "forward", silent = TRUE)

#run the backward algorithm
M.B <- Manly.select(X, M.EM, method = "backward", silent = TRUE)

x11(height = 4, width = 4)
par(mar = rep(0.1, 4))
#plot the results
M.K$id <- id.km
M.K$tau <- rep(1/K, K)
M.K$Mu <- M.K$centers
M.K$la <- matrix(0, K, p)
M.K$S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	M.K$S[,,k] <- M.K$tot.withinss / n / p
}

Manly.plot(X= acidity, model = M.K, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")

dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Kmeans.pdf")
dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Kmeans.eps")


M.MK$tau <- rep(1/K, K)
S <- array(0, dim = c(p,p,K))
for(k in 1:K){
	S[,,k] <- M.MK$S[k]
}
M.MK$S <- S

Manly.plot(X= acidity, model = M.MK, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")

dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyKmeans.pdf")
dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyKmeans.eps")

Manly.plot(X= acidity, model = M.Gauss, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")

dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Gaussian.pdf")
dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Gaussian.eps")


Manly.plot(X= acidity, model = M.EM, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")
dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Manly.pdf")
dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_Manly.eps")


Manly.plot(X= acidity, model = M.F, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")

dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyF.pdf")

dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyF.eps")

Manly.plot(X= acidity, model = M.B, var1 = 1, main = "", ylim = c(0, 0.75), xlab = "", xaxt = "n", ylab = "", yaxt = "n", x.slice = 200, col = "red")
dev.copy2pdf(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyB.pdf")
dev.copy2eps(file = "/Users/zhu/Desktop/Xuwen/R_journal3/Figures/acidity_ManlyB.eps")

