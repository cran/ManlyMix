options(warn=-1)
Manly.model <- function(X, K = 1:5, Gaussian = FALSE, initial = "k-means", nstart = 100, method = "ward.D", short.iter = 5, select = "none", silent = TRUE, plot = FALSE, var1 = NULL, var2 = NULL, VarAssess = FALSE, conf.CI = NULL, overlap = FALSE, N = 1000, tol = 1e-5, max.iter = 1000, ...){
	if(!is.matrix(X)){
		if(is.vector(X)){
			n <- length(X)
			X <- matrix(X, n, 1)
		}
	}
	VAR <- NULL
	n <- dim(X)[1]
	p <- dim(X)[2]
	if(select == "forward"){Gaussian <- TRUE}
	if(select == "backward"){Gaussian <- FALSE}

	if(initial == "k-means"){
		bic <- Inf
		M <- NULL
		for(k1 in K){
			id.km <- kmeans(X, centers = k1, nstart = nstart)$cluster
			if(Gaussian){
				res1 <- try(temp <- Manly.EM(X, id = id.km, tol = tol, max.iter = max.iter))
				if(!inherits(res1, "try-error")){
					if(temp$bic < bic){
						M <- temp
						bic <- M$bic
					}
				}
			}
			else{
				res1 <- try(temp <- Manly.EM(X, id = id.km, la = matrix(0.1, k1, p), tol = tol, max.iter = max.iter))
				if(!inherits(res1, "try-error")){
					if(temp$bic < bic){
						M <- temp
						bic <- M$bic
					}
				}


			}
		}
	}
	else if(initial == "hierarchical"){

		bic <- Inf
		M <- NULL
		for(k1 in K){
			H <- hclust(dist(X), method = method)
			id.hier <- cutree(H, k1)
			if(Gaussian){
				res1 <- try(temp <- Manly.EM(X, id = id.hier, tol = tol, max.iter = max.iter))
				if(!inherits(res1, "try-error")){
					if(temp$bic < bic){
						M <- temp
						bic <- M$bic
					}
				}
			}
			else{
				res1 <- try(temp <- Manly.EM(X, id = id.hier, la = matrix(0.1, k1, p), tol = tol, max.iter = max.iter))
				if(!inherits(res1, "try-error")){
					if(temp$bic < bic){
						M <- temp
						bic <- M$bic
					}
				}
			}
		}


	}
	else if(initial == "emEM"){
		options(warn=-1)
		bic <- Inf
		M <- NULL
		for(k1 in K){
			iter <- 0
			ll <- -Inf
			M1 <- NULL
			repeat{
				iter <- iter + 1
				rs <- sample(1:n, k1)
				
				id.km <- kmeans(X, centers = matrix(X[rs,], k1, p), iter.max = 1)$cluster
				if(Gaussian){
					res1 <- try(temp <- Manly.EM(X, id = id.km, tol = tol, max.iter = short.iter))
				}
				else{

					res1 <- try(temp <- Manly.EM(X, id = id.km, la = matrix(0.1, k1, p), tol = tol, max.iter = max.iter))

				}
				if(!inherits(res1, "try-error")){
					if(temp$ll > ll){
						M1 <- temp
						ll <- M1$ll
					}
				}	
	
				if(iter == nstart){break}
			}
			res1 <- try(temp <- Manly.EM(X, tau = M1$tau, Mu = M1$Mu, S = M1$S, la = M1$la, tol = tol, max.iter = max.iter))
			if(!inherits(res1, "try-error")){
				if(temp$bic < bic){
					M <- temp
					bic <- M$bic
					
				}
			}
			
		}

	}
	if(!Gaussian){
		if(VarAssess){
			VAR <- Manly.var(X, model = M, conf.CI = conf.CI)
		}
	}

	if(select == "forward"){

		M <- Manly.select(X, model = M, method = "forward", tol = tol, max.iter = max.iter, silent = silent)
	}
	else if(select == "backward"){

		M <- Manly.select(X, model = M, method = "backward", tol = tol, max.iter = max.iter, silent = silent)
	}
	O <- NULL
	if(overlap){

		O <- Manly.overlap(M$tau, M$Mu, M$S, M$la, N = N)

	}
	if(plot){

		Manly.plot(X, var1 = var1, var2 = var2, model = M, ...)

	}
	return(list(Model = M, VarAssess = VAR, Overlap = O))

}
Manly.EM <- function(X, id = NULL, la = NULL, tau = NULL, Mu = NULL, S = NULL, tol = 1e-5, max.iter = 1000){
	if(!is.matrix(X)){
		if(is.vector(X)){
			n <- length(X)
			X <- matrix(X, n, 1)
		}
	}


	if (tol <= 0) stop("Wrong value of tol...\n")
	if (max.iter < 1) stop("Wrong number of iterations iter...\n")


	xdims <- dim(X)
	n <- xdims[1]
	p <- xdims[2]


	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")



	if(is.null(id) && (is.null(Mu) || is.null(tau) || is.null(S))) stop("Must provide one initialization method...\n")
	if(!is.null(id)){


		K <- max(id)	

		if(K < 1) stop("Wrong number of mixture components K...\n")
		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}
		if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
		if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")


		x1 <- as.vector(t(X))
		la1 <- as.vector(t(la))
		gamma1 <- rep(0, n*K)
		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)
		tau <- rep(0, K)
		Mu1 <- rep(0, K*p)
		S1 <- rep(0, K*p*p)


		result <- .C("run_Manly", x1 = as.double(x1), id = as.integer(id), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), tau = as.double(tau), Mu1 = as.double(Mu1), S1 = as.double(S1), gamma1 = as.double(gamma1), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "ManlyMix")
	}
	else{
		K <- length(tau)
		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}

		equal.K <- c(dim(la)[1], dim(Mu)[1], dim(S)[3])
		equal.p <- c(dim(la)[2], dim(Mu)[2], dim(S)[1], dim(S)[2])

		if (K < 1) stop("Wrong number of mixture components K...\n")
		if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3])) stop("Inconsistent number of mixture components K...\n")
		if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")


		x1 <- as.vector(t(X))
		la1 <- as.vector(t(la))
		gamma1 <- rep(0, n*K)
		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)
		id <- rep(0, n)
		Mu1 <- as.vector(t(Mu))
		S1 <- as.vector(S)


		result <- .C("run_Manly2", x1 = as.double(x1), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), tau = as.double(tau), Mu1 = as.double(Mu1), S1 = as.double(S1), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv),PACKAGE = "ManlyMix")

	}



	ind <- result$la1==0
	
	M <- K - 1 + 2 * K * p + K * p * (p + 1) / 2 - sum(ind)
	BIC <- Manly.bic(result$ll[1], n, M)

	
	ret <- list(la = matrix(result$la1, nrow=K, byrow=TRUE), tau = result$tau, Mu = matrix(result$Mu1, nrow = K, byrow = TRUE), S = array(result$S1, dim = c(p, p, K)), gamma = matrix(result$gamma1, nrow = n, byrow = TRUE), id = result$id, ll = result$ll[1], bic = BIC, iter = result$conv[1], flag = result$conv[2])

	class(ret) <- "ManlyMix"
	if(result$conv[2] == 1) {
		warning("The EM algorithm does not converge...\n")
	}	
	return(ret)



}


Manly.bic <- function(logl, n, M){
	return(-2 * logl + M * log(n))
}





