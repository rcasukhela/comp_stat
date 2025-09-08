# https://www.r-bloggers.com/2020/08/generating-data-from-a-truncated-distribution/
rtruncpois <- function(n, lambda, a, b) {
  u <- runif(n, ppois(a - 1, lambda), ppois(b, lambda))
  x <- qpois(u, lambda)
  return(x)
}

# read in data
dat = read.table("datasets/hivrisk.dat", header=TRUE, sep = " ")


data = dat

# Number of data points in study
n = sum(data[,2])

# Group Probabilities
alphas = c(1/3)
betas = c(1/3)

# Draw parameters
group_z_p = 0
group_t_p = 1.7
group_p_p = 6.2


# Test Run ----------------------------------------------------------------
# First, assign an individual to a group.
data$freq_z = 0
data$freq_t = 0 
data$freq_p = 0

groups = c()
counts = c()

for (i in 1:n){
  group_assignment = sample(c("Group Z", "Group T", "Group P"), size=1, prob=c(alphas[1], betas[1], (1-alphas[1]-betas[1])))
  groups = c(groups, group_assignment)
  count = -1
  if (group_assignment == "Group Z"){
    count = 0
    data[data$encounters == count, ]$freq_z = data[data$encounters == count, ]$freq_z + 1
  }
  else if (group_assignment == "Group T"){
    # while (count < 0 | count > 16){
    #   count = rpois(1, group_t_p)
    # }
    count = rtruncpois(1, group_t_p, 0, 16)
    data[data$encounters == count, ]$freq_t = data[data$encounters == count, ]$freq_t + 1
  }
  else {
    # while (count < 0 | count > 16){
    #   count = rpois(1, group_p_p)
    # }
    count = rtruncpois(1, group_p_p, 0, 16)
    data[data$encounters == count, ]$freq_p = data[data$encounters == count, ]$freq_p + 1
  }
  counts = c(counts, count)
}
data


# Clean Up Simulated Data -------------------------------------------------
# First, create the promiscuous sample.
data$freq_p = data$freq_p + pmax(data$frequency - data$freq_z - data$freq_t, 0)
data
# Rescale frequencies such that they add up to the original observed frequency.
new_freq_z = round((data$freq_z / (data$freq_z + data$freq_t + data$freq_p)) * data$frequency)
new_freq_t = round((data$freq_t / (data$freq_z + data$freq_t + data$freq_p)) * data$frequency)
new_freq_p = round((data$freq_p / (data$freq_z + data$freq_t + data$freq_p)) * data$frequency)

data$freq_z = new_freq_z
data$freq_t = new_freq_t
data$freq_p = new_freq_p

# The previous step can introduce NaNs, so fill those with 0.
data[is.na(data)] = 0
data

# Optimize Log Likelihood Given Simulated Data ----------------------------
neg_log_lik <- function(alpha, beta, mu, lambda, data) {
  i <- 0:16
  n_i <- data
  log_prob <- log(alpha * (i == 0) + beta * dpois(i, mu) + (1 - alpha - beta) * dpois(i, lambda))
  ll <- sum(n_i * log_prob)
  return(-ll)
}

# Define the function to optimize
nll_opt <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  mu <- params[3]
  lambda <- params[4]
  
  # Calculate the log-likelihood or other objective function
  nll <- neg_log_lik(alpha, beta, mu, lambda, data)
  
  # Return the negative log-likelihood (for minimization)
  return(nll)
}

# Initial values for the parameters
params_init <- c(alpha = alphas[1], beta = betas[1], mu = group_t_p, lambda = group_p_p)

# Run the optimization
fit <- optim(params_init, nll_opt, method = "L-BFGS-B",
             lower=c(0, 0, 0, 0), upper=c(1, 1, Inf, Inf))
fit$par
