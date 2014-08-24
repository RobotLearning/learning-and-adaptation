# R hyperparameter estimation experiments
library(MASS)
library(geoR)

print('Maximum likelihood estimation of hyperparameters...')

# load data saved from matlab into a text file
# in windows
# file <- "C:/Users/okankoc/Learning-and-Adaptation/hyper/estimation2D.txt"
file <- "/home/okankoc/Learning-and-Adaptation/hyper/estimation2D.txt"
values <- read.table(file)
X <- values[,1:2]
Y <- values[,3:7]
# values with mean
Y2 <- values[,8:12] 
n <- length(X)/2

## ML LIKELIHOOD FITTING 
# using function 'likfit' from package 'geoR'
sn0 <- seq(from = 0.01, to = 0.2, length = 10)
L0 <- seq(from = 0.01, to = 1.1, length = 5) 
sf0 <- seq(from = 0.01, to = 2.2, length = 5)
psiR0 <- c(0,1,2,3)
# a matrix can be used to provide several initial values
par0 <- matrix(c(L0, sf0), ncol = 2)

for (i in 1:5) {
  ml <- likfit(coords = X, data = Y[,i], trend = "cte", 
  	 	fix.nugget = FALSE, nugget = 0.04,
  	 	ini.cov.pars = c(0.3, 0.5),
  	 	cov.model = "gaussian", 
      fix.psiA = TRUE, fix.psiR = FALSE, psiR = 2,
  	 	lik.method = "ML")
  # output the results
  print(summary(ml))
  # covariance length scale parameters must be divided by sqrt(2) to get 
  # comparable value with matlab
}

## ML LIKELIHOOD FITTING FOR ADDED MEAN CASE
# using function 'likfit' from package 'geoR'
for (i in 1:5) {
  ml <- likfit(coords = X, data = Y2[,i], trend = "cte", 
  	 	fix.nugget = FALSE, nugget = 0.04,
  	 	ini.cov.pars = c(0.3, 0.5),
  	 	cov.model = "gaussian", 
      fix.psiA = TRUE, fix.psiR = FALSE, psiR = 2,
  	 	lik.method = "ML")
  # output the results
  print(summary(ml))
}
## REML LIKELIHOOD FITTING FOR ADDED MEAN CASE
for (i in 1:5) {
  reml <- likfit(coords = X, data = Y2[,i], trend = "1st", 
               fix.nugget = FALSE, nugget = 0.04,
               ini.cov.pars = c(0.3, 0.5),
               cov.model = "gaussian", 
               fix.psiA = TRUE, fix.psiR = FALSE, psiR = 2,
               lik.method = "REML")
  # output the results
  print(summary(reml))
}
