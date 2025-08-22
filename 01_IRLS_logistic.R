setwd("~/dev/comp_stat")

#########################################################################
### EXAMPLE 2.4 NEWTON'S METHOD (BIVARIATE)
#########################################################################

#########################################################################
# x = initial value
# itr = number of iterations to run
# x.values = contains values of x for each iteration
# g = objective function
# g.prime = first derivative of objective function
# g.2prime = second derivative of objective function
#########################################################################
# Read in the data
data <- read.table("datasets/facerecognition.dat", 
                   header = TRUE,  # Set to TRUE if your file has a header row
                   sep = " ",       # Specify the separator (e.g., " ", ",", "\t")
                   na.strings = c("NA"))  # Specify strings that represent missing values
# Assign y and X
y = data$match
X = data[, c("eyediff", "nosecheekdiff", "variabilityratio")]
# X = data[, c("eyediff")]
# Add a new column of ones
X = cbind(1, X)
X = as.matrix(X)

# Dimension of the problem
n = nrow(X)
k = ncol(X)

## INITIAL VALUES
beta_init = runif(k, -0.1, 0.1)
itr = 40
beta_vals = matrix(0,k,itr+1)
beta_vals[, 1] = beta_init # Adds initial values to the matrix

# Define b_theta.
get_b_theta = function(X, beta_val){
  b_theta = matrix(0,itr,1)
  for (i in 1:n) {
    b_theta[i,1] = log(1 + exp( X[i,] %*% beta_val ))
  }
  return(b_theta)
}

# Define W and pi, the probability that y = 1 given the covariates X.
get_W = function(X, betas){
  n = nrow(X)
  c = X %*% betas
  pi = exp(c) / ( exp(c) + 1 )
  W = matrix(0, n, n)
  for (i in 1:n) {
    W[i, i] = pi[i]*(1-pi[i])
  }
  return(list(W=W, pi=pi))
}

## OBJECTIVE FUNCTION AND DERIVATIVES
log_lik = function(y, X, betas) {
  # y := n x 1 vector
  # X := n x k matrix
  # beta_vec := k x1 vector
  b_theta = get_b_theta(X, betas)
  out = t(y) %*% X %*% betas - t(b_theta) %*% rep(1, nrow(b_theta))
  return(out)
}

score_func = function(y, X, pi) {
  out = t(X) %*% (y - pi)
  return(out)
}
hessian = function(X, W) {
  out = -t(X) %*% W %*% X
  return(out)
}

## MAIN
for(i in 2:itr){
  get_W_result = get_W(X, beta_vals[, i-1])
  W = get_W_result$W
  pi = get_W_result$pi
  beta_vals[, i] = beta_vals[, i-1] - solve(hessian(X, W))%*%score_func(y, X, pi)
}

# compare the industry solver with the one we just made:
glm(match ~ eyediff + nosecheekdiff + variabilityratio, family=binomial, data=data)$coef
beta_vals[, itr]

# > # compare the industry solver with the one we just made:
#   > glm(match ~ eyediff + nosecheekdiff + variabilityratio, family=binomial, data=data)$coef
# (Intercept)          eyediff    nosecheekdiff variabilityratio 
# 1.077480        -9.244228       -13.406903         1.543547 
# > beta_vals[, itr]
# [1]   1.077480  -9.244228 -13.406903   1.543547

# Wow!! A perfect match to the precision stated. I'm happy with that.