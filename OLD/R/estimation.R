# R HYPERPARAMETER ESTIMATION EXPERIMENTS
library(MASS)

print('Maximum likelihood estimation of hyperparameters...')

# squared exponential covariance function handle
ker <- function(x1,x2,p) {
	distsq <- sum((x1-x2)^2)
	(p[2]^2) * exp(-distsq/(2*(p[1]^2)))
}

# handle for generating a covariance matrix
kermat <- function(xs,ker,p) {
	size <- length(xs)
	mat <- matrix(0, nrow = size, ncol = size)
	for (i in 1:size) {
		for (j in 1:size) {
			mat[i,j] <- ker(xs[i],xs[j],p)
		}
	}
	return(mat)
}

# negative log likelihood to minimize
# parameters are always >= 0!
nll <- function(p,x,y) {
	if (length(p) != 3) 
		vec <- y - p[5] * x - p[4] 
	else
		vec <- y
	mat <- kermat(x,ker,p[1:2]) + (p[3]^2) * diag(n)
	rhks <- vec %*% solve(mat,vec)
	cplex <- 1/2 * log(abs(det(mat)))
	nll <- rhks + cplex
	return(nll)
}

## TODO: make this work!
dnll <- function(p,x,y) {
	
	K <- kermat(x,ker,p[1:2])
	dKdsf <- 2 * K/ p[2]
	dKdsn <- 2 * p[3] * diag(n)
	mat1 <- matrix(rep.int(x,n), nrow = n)
	mat2 <- matrix(rep.int(x,n), nrow = n, byrow = TRUE)
	mat <- (mat1 - mat2)^2
	dKdl  <- -mat * K /(p[1]^3)
	Kn <- K + (p[3]^2) * diag(n)
	if (length(p) == 3) {
		vec <- solve(Kn, y)
		dnll1 <- sum(diag(solve(Kn, dKdl)))/2 - vec %*% dKdl %*% vec
		dnll2 <- sum(diag(solve(Kn, dKdsf)))/2 - vec %*% dKdsf %*% vec
		dnll3 <- sum(diag(solve(Kn, dKdsn)))/2 - vec %*% dKdsn %*% vec
		dnll <- c(dnll1,dnll2,dnll3)
	}
	else {
		vec <- solve(Kn, y - p[5] * x - p[4])
		dnll1 <- sum(diag(solve(Kn, dKdl)))/2 - vec %*% dKdl %*% vec
		dnll2 <- sum(diag(solve(Kn, dKdsf)))/2 - vec %*% dKdsf %*% vec
		dnll3 <- sum(diag(solve(Kn, dKdsn)))/2 - vec %*% dKdsn %*% vec
		dnll4 <- -2 * rep.int(1,n) %*% vec
		dnll5 <- -2 * x %*% vec
		dnll <- c(dnll1,dnll2,dnll3,dnll4,dnll5)
	}
	return(dnll)
}

# load data saved from matlab into a text file
# in windows
# file <- "C:/Users/okankoc/Learning-and-Adaptation/hyper/estimation1D.txt"
file <- "/home/okankoc/Learning-and-Adaptation/hyper/estimation1D.txt"
values <- read.table(file)
x <- values$V1
Y <- values[,2:6]
# values with mean
Y2 <- values[,7:11] 
n <- length(x)

# minimization
for (i in 1:5) {
	print(sprintf('Dataset %d', i))
	print('Using conjugate gradient method ...')
	out1 <- optim(par = c(0.1,0.1,0.1), nll, gr = dnll, x, Y[,i], method = "CG")
	print('Estimated zero-mean GP hyperparameters:')
	print(sprintf('ell = %f, sf = %f, sn = %f', abs(out1$par[1]), abs(out1$par[2]), abs(out1$par[3])))
	
	#print('Using Newtons method ...')
	#out2 <- nlm(nll, p = c(0.5,0.5,0.5), x, Y[,i])
	#print('Estimated zero-mean GP hyperparameters:')
	#print(sprintf('ell = %f, sf = %f, sn = %f', abs(out2$estimate[1]), abs(out2$estimate[2]), abs(out2$estimate[3])))
}

# REPEAT BY INCLUDING LINEAR MEAN

# minimization
for (i in 1:5) {
	print(sprintf('Dataset %d', i))
	print('Using conjugate gradient method ...')
	out1 <- optim(par = c(0.5,0.5,0.5,0.5,0.5), nll, gr = dnll, x, Y2[,i], method = "CG")
	print('Covariance hyperparameters:')
	print(sprintf('ell = %f, sf = %f, sn = %f', abs(out1$par[1]), abs(out1$par[2]), abs(out1$par[3])))
	print('Mean hyperparameters:')
	print(sprintf('slope = %f, intercept = %f', out1$par[5], out1$par[4]))

	#print('Using Newtons method ...')
	#out2 <- nlm(nll, p = c(0.5,0.5,0.5,0.5,0.5), x, Y2[,i])
	#print('Covariance hyperparameters:')
	#print(sprintf('ell = %f, sf = %f, sn = %f', abs(out2$estimate[1]), abs(out2$estimate[2]), abs(out2$estimate[3])))
	#print('Mean hyperparameters:')
	#print(sprintf('slope = %f, intercept = %f', out2$estimate[5], out2$estimate[4]))
}
