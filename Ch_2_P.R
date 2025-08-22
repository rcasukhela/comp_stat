# Problem 2.1

# The following data are an i.i.d. sample from a Cauchy(\theta, 1) distribution.
# `c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40,
#    4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)`

# A note for myself, gamma (the scale parameter) is set to 1.


## 2.1a

### Graph the log-likelihood function.

# The log-likelihood function for an i.i.d. Cauchy random sample is:

log_lik = function(x, theta) {
  # x := nx1 vector
  # theta := scalar
  n = length(x)
  out = -n*log(pi) - sum(log(1 + (x - theta)^2))
  return(out)
}

plot_results = c()
plot_theta = seq(-10, 10, length.out=100)
plot_x = runif(10,-10,10)
for (i in 1:length(plot_theta)) {
  plot_results = c(plot_results, log_lik(plot_x, plot_theta[i]))
}

plot(plot_theta, plot_results)
