# https://www.r-bloggers.com/2020/08/generating-data-from-a-truncated-distribution/
rtruncpois <- function(n, lambda, a, b) {
  u <- runif(n, ppois(a - 1, lambda), ppois(b, lambda))
  x <- qpois(u, lambda)
  return(x)
}

# read in data
data = read.table("datasets/hivrisk.dat", header=TRUE, sep = " ")

# Number of data points in study
n = sum(data[,2])

# Group Probabilities
alphas = c(1/3)
betas = c(1/3)

# Draw parameters
group_z_p = 0
group_t_p = 5
group_p_p = 12

# Test run:
# First, assign an individual to a group.
groups = c()
counts = c()
for (i in 1:n){
  group_assignment = sample(1:3, size=1, prob=c(alphas[1], betas[1], (1-alphas[1]-betas[1])))
  groups = c(groups, group_assignment)
  count = -1
  while ( count < 0 | count > 16 ) {
    count = switch(group_assignment,
               rpois(1, group_z_p),
               rpois(1, group_t_p),
               rpois(1, group_p_p))
  }
  counts = c(counts, count)
}
g_c = as.data.frame(t(table(cbind(groups, counts))))[, c(2,3)]
sim_df = as.data.frame(t(table(counts)))[,c(2,3)]
sum(g_c$Freq)

# Given the observed data, sample the 
data[2,2]