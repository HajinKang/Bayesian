# Understanding the Metropolis-Hastings Algorithm
# by Siddhartha Chib and Edward Greenberg



# ### 7.1
mu <- c(1,2)
cov_mat <- matrix(c(1, 0.9, 0.9, 1), 2, 2, TRUE)
set.seed(2019121102)
sample_chol <- rmvnorm(3000, mu, cov_mat)

# 7.1.1 Method 1 : MH Random walk with uniform increments
alpha <- function(x, y) {
  num <- exp(-0.5 * ((y-mu)%*% solve(cov_mat)%*%(y-mu)))
  den <- exp(-0.5 * ((x-mu)%*% solve(cov_mat)%*%(x-mu)))
  min(num/den, 1)
}
sample_val <- function(x) {
  delta_1 <- 0.75
  delta_2 <- 1
  return( x + c(runif(1, -delta_1, delta_1), runif(1, -delta_2, delta_2)))
}

N <- 10000
x_start <- 0
results_1 <- matrix(NA, nrow = N, ncol = 2)
results_1[1,] <- x_start
i <- 2
for( i in 2:N ) {
  new_val <- sample_val(results_1[i - 1, ])
  alpha_x_y <- alpha(results_1[i - 1], new_val)
  if(runif(1) < alpha_x_y){
    results_1[i, ] <- new_val
  }else{
    results_1[i,] <- results_1[i-1, ]
  }
}
results_1 <- results_1[round(0.1*N):N,]

# 7.1.2 Method 2 : MH Random walk with normal increments
sample_val <- function(x) {
  delta_1 <- 0.6
  delta_2 <- 0.4
  return(x + c(rnorm(1, 0, sqrt(delta_1)), rnorm(1, 0, sqrt(delta_2))))
}
results_2 <- matrix(NA, nrow = N, ncol = 2)
results_2[1,] <- x_start
i <- 2
for( i in 2:N ) {
  new_val <- sample_val(results_2[i - 1, ])
  alpha_x_y <- alpha(results_2[i - 1], new_val)
  if(runif(1) < alpha_x_y){
    results_2[i, ] <- new_val
  }else{
    results_2[i,] <- results_2[i-1, ]
  }
}
results_2 <- results_2[round(0.1*N):N,]

# 7.1.3 Method 3 : Psuedo Dominating Density method
D <- diag(c(2, 2))
k <- 0.9
f <- function(z) 1/(2*pi)*1/sqrt(det(cov_mat)) * exp(-0.5*(z-mu)%*%solve(cov_mat)%*%(z-mu))
h <- function(z) dmvnorm(z, c(0, 0), D)
results_3 <- matrix(NA, nrow = N, ncol = 2)
results_3[1,] <- x_start
i <- 2
for(i in 2:N){
  cur_val <- results_3[i - 1, ]
  accept <- FALSE
  while(accept == FALSE){
    z <- c(rmvnorm(1, c(0, 0), D))
    if(runif(1) <= f(z)/(k*(h(z)))) {
      new_val <- z
      accept <- TRUE
    }
  }
  C1 <- f(cur_val) < k*h(cur_val)
  C2 <- f(new_val) < k*h(new_val)
  if(C1 == 1) alpha_mult <- 1
  if(C1 == 0 & C2 ==1) alpha_mult <- k*h(cur_val)/f(cur_val)
  if(C1 == 1 & C2 ==0) alpha_mult <- min(f(new_val)*h(cur_val)/(f(cur_val)*h(new_val)),1)
  if(runif(1) <= alpha_mult){
    results_3[i, ] <- new_val
  }else{
    results_3[i,] <- results_3[i-1, ]
  }
}
results_3 <- results_3[round(0.1*N):N,]

# 7.1.4 Method 4 : MH Random walk with autoregressive chains
sample_val <- function(x) {
  delta_1 <- 1
  delta_2 <- 1
  return(2*mu - x + c(runif(1, -delta_1, delta_1), runif(1, -delta_2, delta_2)))
  }
results_4 <- matrix(NA, nrow = N, ncol = 2)
results_4[1,] <- x_start
i <- 2
for( i in 2:N ) {
  new_val <- sample_val(results_4[i - 1, ])
  alpha_x_y <- alpha(results_4[i - 1], new_val)
  if(runif(1) < alpha_x_y){
    results_4[i, ] <- new_val
  }else{
    results_4[i,] <- results_4[i-1, ]
  }
}
results_4 <- results_4[round(0.1*N):N,]

# 7.1 Summary: Autocorrelation of the chain across four methods

corr_chain <- function(temp){
  n <- dim(temp)[1]
  cor(temp[1:(n-1)], temp[2:n])
}
cor_results <- data.frame( "Method 1" = corr_chain(results_1),
                           "Method 2" = corr_chain(results_2),
                           "Method 3" = corr_chain(results_3),
                           "Method 4" = corr_chain(results_4))
cor_results <- t(cor_results)
colnames(cor_results) <- c("rho")


windows()
par(mfrow=c(2,2))
plot(results_1, ylim=c(-6,6), xlim=c(-6,6), main="Method 1")+points(sample_chol, col="orange")
plot(results_2, ylim=c(-6,6), xlim=c(-6,6), main="Method 2")+points(sample_chol, col="orange")
plot(results_3, ylim=c(-6,6), xlim=c(-6,6), main="Method 3")+points(sample_chol, col="orange")
plot(results_4, ylim=c(-6,6), xlim=c(-6,6), main="Method 4")+points(sample_chol, col="orange")



# ### 7.2
library(mvtnorm)
library(ordinalLBM)

set.seed(2019121102)
n <- 100
ar2_data <- arima.sim( n = n, list( order=c(2,0,0), ar = c(1,-0.5)), sd = 1)

Y <- ar2_data[1:2]
N <- 5000
params <- matrix(NA, nrow = N, ncol = 3)
y <- cbind(ar2_data[3:n])
wt <- cbind(ar2_data[2:(n-1)], ar2_data[1:(n-2)])
V_inv_mat <- function(phi){
  matrix( c( 1-phi[2]^2, -phi[1]*(1+phi[2]),
             -phi[1]*(1+phi[2]),1-phi[2]^2) , 2, 2, T)
}
Psi <- function(phi, sig2) {
  V_inv <- V_inv_mat(phi)
  1/sig2 *1 /sqrt(abs(det(V_inv))) * exp(-1/(2*sig2) * t(Y)%*%V_inv%*%Y)
}
params[1,] <- c(2, 0.5, 0.2)
i <- 2
for(i in 2:N) {
  #simulate gamma rv
  cur_phi <- params[(i-1), 2:3]
  V_inv <- V_inv_mat(cur_phi)
  params[i,1] <- 1/rgamma(1, n/2, 0.5*(Y%*%V_inv%*%Y + sum((y - wt%*%cur_phi)^2)))
  
  G <- t(wt)%*%wt
  phi_mu <- solve(G) %*%(t(wt)%*%y)
  phi_sig <- params[i,1]*solve(G)
  new_phi <- c(rmvnorm(1, phi_mu, phi_sig))
  alpha_mult <- min(Psi(new_phi, params[i, 1])/ Psi(cur_phi, params[i, 1]), 1)
  if(runif(1) < alpha_mult){
    params[i, c(2,3)] <- new_phi
  }else{
    params[i, c(2,3)] <- cur_phi
  }
}
params <- params[round(0.1*N):N,]
params <- cbind(params, sqrt(params[,1]))
results <- apply(params,2, function(z){
  c(mean(z), sd(z), median(z), quantile(z, prob=c(0.025, 0.975)),
    cor(z[2:4501],z[1:4500]))
})
results <- t(results)
colnames(results) <- c("mean","sd", "median", "lower", "upper","correlation")
results <- as.data.frame(results)
# results$"o.b.m sd" <- sqrt(c(apply(params, 2,function(z) olbm(z, 50))))
results <- results[-1,]
rownames(results) <- c("phi1","phi2","sigma^2")
# results <- results[,c("mean", "o.b.m sd", "sd",
#                      "median", "lower", "upper","correlation")]


