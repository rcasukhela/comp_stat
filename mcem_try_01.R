# My goal is to create a general EM solver
# That can accept the complete log-likelihood (that is simulated partially
# by sampling the missing values, so Monte Carlo E-Step),
# Use numerical solving techniques for the M-Step,
# And bootstrapping for the variance estimation step.

# After some thought, I decided that my goal for Givens's and Hoeting's book
# should be to make generic solvers and leave plenty of room for code optimization
# and efficiency because I don't have strong use cases for these algorithms yet.
# This approach is the simplest way to demonstrate basic understanding of the
# fundamentals of the algorithms and not bog myself down.

neg_log_lik = function(n, lambda, y_vec){
  out = -1*(n*log(lambda) - lambda*sum(y_vec))
  return(out)
}

# Simulation Hyperparameters (user-defined)
n = 100
lambda = 2
i_itr = 50
j_itr = 100

# Data generating function.
complete_data = rexp(n, rate=lambda) # complete data
c = runif(n, 0.5, 0.6) # censoring times
delta = ifelse(complete_data <= c, 1, 0) # was data censored?
obs = pmin(complete_data, c)
hist(obs)
full_data = cbind(complete_data, c, obs, delta)

# Generate observed data (what we would see in the experiment)
obs_data = full_data[, c(3, 4)]










lambdas = c()

# Make a good first guess for lambda.
lambda_init = 1 / mean(obs_data[delta == 1, c(1)])

lambdas = c(lambda_init)

# Grab the observed values where we are censored.
censored_values = obs_data[delta == 0, c(1)]

num_cen = length(censored_values)

# Simulate the censored values. We know from properties of the exponential
# distribution that the missing values will have an Exp(lambda_guess) + c_i
# distribution, owing to the memory-less property of that distribution.
# So we can forget about the past, simulate using the guess for the rate that
# we have, and just add the censoring time back to get what the value should've
# been.

lambda_itr = lambda_init
for (i in 1:i_itr){
  ys = matrix(
    nrow=n,
    ncol=j_itr
  )
  for (j in 1:j_itr){
    z = rexp(num_cen, lambda_itr) + censored_values
    y = c(z, obs_data[delta == 1, c(1)])
    ys[,j] = y
  }
  
  Q_hat = function(n, lambda, j_itr, ys){
    out = -1*(n*log(lambda) - lambda * sum(colSums(ys)))/j_itr
    return(out)
  }
  
  lambda_itr = optim(lambda_itr, Q_hat, n=n, j_itr=j_itr, ys=ys, method="BFGS")$par
  
  lambdas = c(lambdas, lambda_itr)
}

lambda_final = lambdas[length(lambdas)]

# Now, bootstrap the variance estimates. Start lambda around 2.5.
var_lambdas = c()
var_lambda_init = 2.5
var_lambdas = c(var_lambdas, var_lambda_init)











