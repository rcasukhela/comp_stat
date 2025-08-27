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
plot_x = rcauchy(1000, -1, 1)
plot_x = c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40,
          4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
for (i in 1:length(plot_theta)) {
  plot_results = c(plot_results, log_lik(plot_x, plot_theta[i]))
}

plot(plot_theta, plot_results)

# Observations:
# At small n, the log likelihood is very noisy and can have multiple optima
# / stationary points.
# At large n, the log likelihood stabilizes and there is a clear peak at
# approximately the true location parameter.
# For our specific data, we have a highly non-convex function!!! This can lead
# to very serious issues with estimation.

score_func = function(x, theta) {
  out = 2 * sum( (x - theta) / (1 + (x - theta)^2) )
  return(out)
}

hessian   = function(x, theta) {
  out = sum( 2*((x - theta)^2 -1) / (1 + (x - theta)^2)^2 )
  return(out)
}

# Now, we can find the MLE of theta using the Newton-Raphson method.
# Initial Values
# Try the following starting points for theta: -11, -1, 0, 1.5, 4, 4.7, 7, 8, 38
# Is the mean of the data a good starting point?
x = c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40,
      4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
theta = mean(x)
itr = 40

# Main
for(i in 1:itr){theta = theta - score_func(x, theta)/hessian(x, theta)}
theta

# Observations:
# The MLE approximation changes drastically by starting point. Clearly, intelligent
# starting points are crucial for the algorithm to properly converge.
# I think -1 and 0 worked the best, but this is also the region where the convexity
# of the function was guaranteed. Beyond this, the mean, and other starting points,
# sucked.


# 2.1b:
#########################################################################
### EXAMPLE 2.1 BISECTION
#########################################################################
univariate_bisection = function(a, b, data, x, itr, g, g.prime){
    # a = initial left endpoint
    # b = initial right endpoint
    # x = initial value
    # itr = number of iterations to run
    # g = objective function
    # g.prime = first derivative of objective function
    #########################################################################
    
    ## MAIN
    for (i in 1:itr){
      if (g.prime(data, a)*g.prime(data, x) < 0) {b = x}
      else {a = x}
      x = a+(b-a)/2
    }
    
    ## OUTPUT
    return(list(
    theta=x		# FINAL ESTIMATE
    ,g.theta=g(data, x)		# OBJECTIVE FUNCTION AT ESTIMATE
    ,g.prime.theta=g.prime(data, x) 	# GRADIENT AT ESTIMATE
    ))
    
    
}

## INITIAL VALUES
a = -1
b = 1
theta = a+(b-a)/2
itr = 40

univariate_bisection(a, b, x, a+(b-a)/2, itr, log_lik, score_func)

# $theta
# [1] -0.1922866
# 
# $g.theta
# [1] -72.91582
# 
# $g.prime.theta
# [1] -1.871461e-12

# So, the bisection method got us the same answer as well.

# Use additional runs to illustrate manners in which the bisection method
# may fail to find the global maximum.

# A simple case is choosing an interval where the optimum cannot be.
a = -10
b = -5
univariate_bisection(a, b, x, a+(b-a)/2, itr, log_lik, score_func)









# 2.1c: Apply fixed-point iterations as in 2.29, starting from -1,
# with scaling choices of alpha = 1, 0.64, 0.25. Investigate other choices
# of starting values and scaling factors.

## INITIAL VALUES
alpha = 0.25
theta = 3.5
itr = 100000

## MAIN
for(i in 1:itr){theta = alpha*score_func(x, theta) + theta}

## OUTPUT
theta		# FINAL ESTIMATE
log_lik(x, theta) 		# OBJECTIVE FUNCTION AT ESTIMATE
score_func(x, theta) 	# GRADIENT AT ESTIMATE

# Reducing the scaling factor really helped the algorithm to converge to the
# right value. Different starting points lead to different optima.






# 2.1d: using (-2, -1) as starting points, apply the secant method to estimate theta.
# Then try (-3, 3) as starting points. And other ones.
## INITIAL VALUES
thetas = c(-3, 3)
itr = 100

## MAIN
for(t in 2:itr){
  new_theta = ( thetas[t] - score_func(x, thetas[t])
    * ( 
        (thetas[t] - thetas[t-1])
        / (score_func(x, thetas[t]) - score_func(x, thetas[t-1])) ))
  thetas[t+1] = new_theta
  if (abs(thetas[t+1] - thetas[t]) < 1E-6) break
  }

# Added in a simple stopping condition just so I don't get NaNs.
# Yeah, about the same conclusions as the other algorithms.



# 2.1e: For an iid Cauchy sample, it's tough to say what's better or worse.
# I mean all algorithms seem to perform well, but you have to be careful with
# hyperparameters. That's always true.

# In the case of an iid Normal sample, I'd expect Newtwon-Raphson to blow the
# other algorithms out of the water because the likelihood function is more
# strictly convex and the N-R method exploits more information about the objective
# function, which also happens to be more globally accurate.

