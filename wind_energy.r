library(numDeriv)

setwd("~/Desktop/statistical_modelling/Assignments")

# 1. read the data tuno.txt into R
data <- read.table("~/Desktop/statistical_modelling/Assignments/tuno.txt", header = TRUE)
head(data)
data$pow.obs_norm = data$pow.obs / 5000 # max capacity 

#2. Make a graphical presentation of data or parts of the data, and present some summary statistics.

qqnorm(data$pow.obs_norm)
qqline(data$pow.obs_norm)
summary(data$pow.obs_norm)

library(MASS)
boxcox(lm(pow.obs ~ 1, data = data))
#lambda is between 0-0.5


# box-cox transformation
bc <- function(lambda, y) {
  y.i <- (y ^(lambda) - 1)/lambda
  if(lambda==0){y.i <- log(y)}
  return(y.i)
}

ll.lambda <- function(lambda,y){
  n <- length(y) ## number of obs
  y.lambda <- bc(lambda, y) ## Transform data

  yl.bar <- mean(y.lambda) 
  sigmasq <- mean((y.lambda-yl.bar)^2) 
 
  ## profile log-likelihood
  - n/2 * log(sigmasq) - n/2 + (lambda-1)*sum(log(y))
}

## profile likelihood of lambda

# lambda <- seq(0,0.5,by=0.01)
# ll <- sapply(lambda, ll.lambda, y=data$pow.obs_norm)
# plot(lambda, ll-max(ll),type="l",ylim=c(-5,0))
# lines(range(lambda), -qchisq(0.95,df=1)/2 * c(1,1),lty=2)
# mle.lambda <- optimize(ll.lambda, c(0, 2), y=data$pow.obs, maximum=TRUE)
# mle.lambda
#0.3467421

lambda <- mle.lambda$maximum
# qqnorm(bc(lambda, y))
# qqline(bc(lambda, y))


# Wind Energy equation 1 transformation
eq1 <- function(lambda, y) {
  1/lambda * log(y ^ lambda / (1 - y ^ lambda))
}

eq1.ll.lambda <- function(lambda,y){
  n <- length(y) ## number of obs
  y.lambda <- eq1(lambda, y) ## Transform data
  
  yl.bar <- mean(y.lambda) 
  sigmasq <- mean((y.lambda-yl.bar)^2) 
  
  ## profile log-likelihood
  - n/2 * log(sigmasq) - n/2 + sum(log(1 / (y * (-y ^ lambda + 1))))
}

## profile likelihood of lambda

lambda <- seq(0.001, 0.5, by=0.01)
ll <- sapply(lambda, eq1.ll.lambda, y=data$pow.obs_norm)

plot(lambda, ll-max(ll), type="l")
lines(range(lambda), -qchisq(0.95,df=1)/2 * c(1,1),lty=2)

mle.eq1.lambda <- optimize(eq1.ll.lambda, c(0.1, 0.5),y=data$pow.obs_norm, maximum=TRUE)
mle.eq1.lambda
#0.3467421

lambda <- mle.lambda$maximum
qqnorm(eq1(lambda, data$pow.obs_norm))
qqline(eq1(lambda, data$pow.obs_norm))


# Wind Energy equation 2 transformation


eq2 <- function(lambda, y) {
  2 * log(y^lambda / (1-y)^(1-lambda))
}

eq2.ll.lambda <- function(lambda,y){
  n <- length(y) ## number of obs
  y.lambda <- eq2(lambda, y) ## Transform data
  yl.bar <- mean(y.lambda) 
  sigmasq <- mean((y.lambda-yl.bar)^2) 
  ## profile log-likelihood
  - n/2 * log(sigmasq) - n/2 + sum(log(2 * lambda + 2 * y - (4 * y * lambda)) - log(y * (-y + 1)))
}

## profile likelihood of lambda

lambda <- seq(0, 1, by=0.01)
ll <- sapply(lambda, eq2.ll.lambda, y=data$pow.obs_norm)

plot(lambda, ll-max(ll), type="l")
lines(range(lambda), -qchisq(0.95,df=1)/2 * c(1,1),lty=2)

mle.lambda <- optimize(eq2.ll.lambda, c(0.1, 0.5),y=data$pow.obs_norm,maximum=TRUE)
mle.lambda
#0.2523292

lambda <- mle.lambda$maximum
qqnorm(eq2(lambda, data$pow.obs_norm))
qqline(eq2(lambda, data$pow.obs_norm))


#profile likelihood with respect to mean with the first transformation
nll <- function(theta, lambda, y) {
  mu <- theta[1]
  sigma2 <- theta[2]
  y_lambda <- eq1(lambda, y)
  - (sum(dnorm(y_lambda, mu, sd=sqrt(sigma2), log=TRUE)) + sum(log(1 / (y * (-y ^ lambda + 1)))))
}

opt <- nlminb(c(0, 1), nll, lambda=mle.eq1.lambda$maximum, y=data$pow.obs_norm, lower=c(-Inf, 0), upper=c(Inf, Inf))
opt

H <- hessian(nll, opt$par, lambda=mle.eq1.lambda$maximum, y=data$pow.obs_norm)

# Mean confidence interval
opt$par[1] + c(-1,1) * 4 * sqrt(solve(H)[1,1])
# [1] 1.427681 3.642367

# Sigma confidence interval
opt$par[2] + c(-1,1) * 4 * sqrt(solve(H)[2,2])
# [1] 14.71449 29.42899

par(mfrow=c(1,1))
mu <- seq(1, 4, by=0.1)
sigma2 <- seq(15, 30, by=0.1)
ll <- matrix(nrow = length(mu), ncol = length(sigma2))

## Find the values of the log likelihood
for(i in 1:length(mu)){
  for(j in 1:length(sigma2)){
      ll[i,j] <- nll(c(mu[i], sigma2[j]), lambda=mle.eq1.lambda$maximum, y=data$pow.obs_norm)
  }
}

Nll <- -ll + opt$objective

## Define the levels
alpha <- c(0.5, 0.25, 0.1, 0.05, 0.01)
contour(x = mu, y = sigma2, z = Nll,
        levels = -qchisq(1 - alpha, df = 2) / 2,
        labels = alpha, xlab="mu", ylab="sigma2")

#Wind speed = Weibull distribution

# TODO try Weibull 

lognormal_l <- function(theta, y){
  log_y <- log(y)
  sum(dnorm(log_y, mean=theta, sd=sd(log_y), log=TRUE)) + sum(log(1/y))
}



normal_l <- function(theta, y){
  sum(dnorm(y, mean=theta, sd=sd(y), log=TRUE))
}

par(mfrow=c(2, 2))
theta <- seq(-100, 100)
normal_l_at_theta <- sapply(theta, normal_l, y = data$ws30)
plot(theta, normal_l_at_theta / min(normal_l_at_theta), type="l", main="Normal")

mle_normal.mu <- optimize(normal_l, c(0, 100), y=data$ws30, maximum=TRUE)
mle_normal.mu

abline(v=mle_normal.mu)

theta <- seq(-100, 100)
ln <- sapply(theta, lognormal_l, y = data$ws30)
plot(theta, ln / max(ln), type="l", main="Lognormal")

mle_lognormal.mu <- optimize(lognormal_l, c(0, 100), y=data$ws30, maximum=TRUE)
mle_lognormal.mu


abline(v=mle_lognormal.mu)

par(mfrow=c(1, 1))
qqnorm(data$ws30)
qqline(data$ws30)
fit <- qqnorm(log(data$ws30))
qqline(log(data$ws30))

# Wind direction
par(mfrow=c(1, 1))
plot(data$wd30)
hist(data$wd30)

qqnorm(data$wd30)

# Fitting normal distribution

theta <- seq(-10, 10)
l <- sapply(theta, normal_l, y = data$wd30)
plot(theta, l / max(l), type="l", main="normal")

mle_normal.mu <- optimize(normal_l, c(0, 100), y=data$wd30, maximum=TRUE)
mle_normal.mu
# 3.60239

# Fitting Von Mises distribution
library(circular)

von_mises_nll <- function(theta, y) {
  mu <- theta[1]
  kappa <- theta[2]
  -sum(dvonmises(y, mu, kappa, log=TRUE))
}

opt <- nlminb(c(3, 1), von_mises_nll, y=data$wd30, lower=c(-Inf, 0), upper=c(Inf, Inf))
opt
library(numDeriv)
# Confidence intervals
H <- hessian(von_mises_nll, opt$par, y=data$wd30) 
opt$par[1] + c(-1,1) * 4 * sqrt(solve(H)[1,1])
opt$par[2] + c(-1,1) * 4 * sqrt(solve(H)[2,2])

mu <- seq(4, 5.5, by=0.1)
kappa <- seq(0.2, 1, by=0.1)
ll <- matrix(nrow = length(mu), ncol = length(kappa))

## Find the values of the log likelihood
for(i in 1:length(mu)){
  for(j in 1:length(kappa)){
    ll[i,j] <- von_mises_nll(c(mu[i], kappa[j]), y=data$wd30)
  }
}

Nll <- -ll + opt$objective

## Define the levels
alpha <- c(0.5, 0.25, 0.1, 0.05, 0.01)
contour(x = mu, y = kappa, z = Nll,
        levels = -qchisq(1 - alpha, df = 2) / 2,
        labels = alpha, xlab="mu", ylab="kappa")

# TODO: use qqplot with arbitrary theoretical distribution comparison. Might need to install a new package.

#PART TWO

# TODO

fit <- glm(data$pow.obs_norm ~ data$wd30 + data$ws30, data = data, family = gaussian(link = "identity"))
summary(fit)

fit1 <- glm(data$pow.obs_norm ~ data$ws30, data = data, family = gaussian(link = "identity"))
summary(fit1)



#PART THREE AUTORERESSION

setwd("~/Desktop/statistical_modelling/Assignments")

# 1. read the data tuno.txt into R
data <- read.table("~/Desktop/statistical_modelling/Assignments/tuno.txt", header = TRUE)
head(data)
data$pow.obs_norm = data$pow.obs / 5000 # max capacity 


# 1. Fit simple linear model Y(0.2) = β0 +β1ws +β2ws2 +ε; ε ∼ N(0,σ^2)
lambda <- 0.2

# m1 (without log):
pow <- eq1(lambda, data$pow.obs_norm)

nll <- function(theta, pow, ws) {
  b0 <- theta[1]
  b1 <- theta[2]
  b2 <- theta[3]
  sigma2 <- theta[4]
  mu <- b0 + b1 * ws + b2 * (ws ^ 2)
  -(sum(dnorm(pow, mu, sd=sqrt(sigma2), log=TRUE)))
}

opt1 <- nlminb(c(0, 0, 0, 1), nll, pow=pow, ws=data$ws30, lower=c(-Inf, -Inf, -Inf, 0), upper=c(Inf, Inf, Inf, Inf))
opt1
# b0=-7.14619103  b1=1.70479162 b2=-0.03210044 sigma=13.90059881
#Likelihood 787.6524

m1 <- glm(pow ~ data$ws30 + I(data$ws30 ^ 2))
summary(m1)
# -7.1461 1.7048  -0.0321 
confint(m1)


# m2 (with ws log):
log_ws <- log(data$ws30)
m2 <- glm(pow ~ log_ws + I(log_ws ^ 2), family=gaussian("identity"))
m2

# m3 (ws^3 term):
m3 <- glm(pow ~ data$ws30 + I(data$ws30 ^ 2) + I(data$ws30 ^ 3), family=gaussian("identity"))
m3

# m4 (without ws^2):
m4 <- glm(pow ~ data$ws30, family=gaussian("identity"))
m4

# m5 (with wd)
m5 <- glm(pow ~ data$ws30 + I(data$ws30 ^ 2) + data$wd30, family=gaussian("identity"))
m5

# m6 (with wd using von misses) - Not sure if this is how to combine them :(
von_misses_wd <- dvonmises(data$wd30, opt$par[1], opt$par[2])
m6 <- glm(pow ~ data$ws30 + I(data$ws30 ^ 2) + von_misses_wd, family=gaussian("identity"))
m6

AIC(m1, m2, m3, m4, m5, m6)
# m1 is the best!

residuals <- resid(m1)
qqnorm(residuals)
shapiro.test(residuals)
matrix <- cbind(residuals[-1], residuals[-length(residuals)])
plot(matrix)

library(mclust)
mvn('XII', matrix)
#sigmasq = 13.935
#rho = 0
#loklik -1570.563

# 2. From the residuals of the model in 1, fit a multivariate gaussian 
#    with mean 0. This will tell us something about the relationship between
#    residuals. Show CIs and standard error analysis, profile likelihoods etc..
# 3. Algebraic fisher information vs numerical fisher information that we got
#    above

# 4. For the residual model in 2. and 3. find better models using 
#    transformations on the parameters. Redo analysis 3 and 4 on the best model.

library(mvtnorm)
library(numDeriv)

## negative log likelihood
nll <- function(theta, y){
  ## set up parameters
  sigma2 <- theta[1]
  rho <- theta[2]
  Sigma <- matrix(1, ncol = 2, nrow = 2)
  Sigma[1, 2] <- Sigma[2, 1] <- rho
  Sigma <- sigma2 * Sigma
  ## negative log-likelihood
  - sum(dmvnorm(y, sigma = Sigma, log=TRUE))
}

## Initialize parameters
theta0 <- c(1,0)
## Find optimal parameters and the Hessian
opt <- nlminb(theta0, nll,
              lower = c(0, -1),
              upper = c(Inf, 1), y = matrix)
opt

H <- hessian(nll, opt$par, y=matrix)

n <- dim(matrix)[1]
sigma2.hat <- opt$par[1]
# [1] 13.93593
rho.hat <- opt$par[2]
# [1] 0.3230793
H
# [,1]       [,2]
# [1,]  1.477781  -7.429019
# [2,] -7.429019 395.142075

n / sigma2.hat ^ 2
# [1] 1.477781
n * (1 + rho.hat ^ 2) / (1 - rho.hat ^ 2) ^ 2
# [1] 395.142
- n * rho.hat / (sigma2.hat * (1 - rho.hat^2))
# [1] -7.429018

## Plot of likelihood
## Define the parameters where the likelihood
## should be calculated
sigma2 <- seq(10, 20, by = 0.01)
rho <- seq(0.2, 0.5, by = 0.01)
ll <- matrix(nrow = length(sigma2), ncol = length(rho))

## Find the values of the log likelihood
for(i in 1:length(sigma2)){
  for(j in 1:length(rho)){
    ll[i,j] <- nll(c(sigma2[i], rho[j]), y=matrix)
  }
}
## normalized log-likelihood
Nll <- -ll + opt$objective
## Define the levels
alpha <- c(0.5,0.25,0.1,0.05,0.01)
contour(x = sigma2, y = rho,
        z = Nll, ylim=c(0,1),
        levels = -qchisq(1 - alpha, df = 2) / 2,
        labels = alpha, xlab="sigma2",ylab="rho")

# TODO add confidence interval and profile likelihood of parameters
# TODO make profile likelihood plot of rho and compare with quadratic approximation
# TODO use fisher's z transform to improve model and report hessian

# ----
# 5. Estimate AR on residuals of the model in 1. compare it with model from 4.

x <- residuals
  
plot(x, type="l")
acf(x)

## Moment estimates
(phi.hat <- acf(x, plot=FALSE)$acf[2]) ## ACF in lag 1
cor(x[-1], x[-length(x)])
(sigma <- var(x) * (1-phi.hat^2))


## Likelihood estimation
nll <- function(theta, y){
  n <- length(y) - 1
  n/2 * log(theta[1]) + 1/(2*theta[1]) * sum((y[-1]-theta[2]*y[-(n+1)])^2)
}

## MLE

# numerical
(opt <- nlminb(c(1, 1/2), nll, y=x, lower=c(0.001 ,-0.999),upper=c(Inf, 0.999)))

# analytical
(phi <- sum(x[-1] * x[-length(x)]) / sum(x[-length(x)] ^ 2))

## standard errors
library(numDeriv)
V <- solve(hessian(nll, opt$par, y=x))
(se <- sqrt(diag(V)))

## Profile likelihood for phi
llp.phi <- function(phi, x, x0){
  n <- length(x) - 1
  -n/2 * log(sum((x[-1]-phi*x[-(n+1)])^2))
}

## Plot profile likelihood
phi <- seq(max(opt$par[2]-3*se[2],-1),
           min(opt$par[2]+3*se[2],1),
           length=200)

llp <- sapply(phi, llp.phi, x=x[-1], x0)

plot(phi, exp(llp - max(llp)),type="l", main = "Profile likelihood for phi")
lines(range(phi), exp(-c(1,1) * qchisq(0.95,df=1)/2), col=2,lty=2)

## Directly in R
arima(x, order=c(1,0,0))
opt$par ## rather close. 

## Did it help?
e <- x[-1] - opt$par[2] * x[-length(x)]
acf(e)

# ----
# 6. Combine both model 4. and model 5.

ar_nll <- function(theta, y){
  sigma <- theta[1]
  phi <- theta[2]
  n <- length(y) - 1
  n/2 * log(sigma) + 1/(2*sigma) * sum((y[-1] - phi*y[-(n+1)])^2)
}

nll <- function(theta, pow, ws) {
  b0 <- theta[1]
  b1 <- theta[2]
  b2 <- theta[3]
  # sigma2_norm <- theta[4]
  
  sigma2_ar <- theta[4]
  phi <- theta[5]
  
  mu <- b0 + b1 * ws + b2 * (ws ^ 2)
  # linear_nll <- -(sum(dnorm(pow, mu, sd=sqrt(sigma2_norm), log=TRUE)))

  residuals <- pow - mu
  ar_theta <- c(sigma2_ar, phi)
  ar_nll <- ar_nll(ar_theta, residuals)
  
  ar_nll
}

opt <- nlminb(c(0, 0, 0, 1, 1/2), nll, pow=pow, ws=data$ws30, lower=c(-Inf, -Inf, -Inf, 0.001 ,-0.999), upper=c(Inf, Inf, Inf, Inf, 0.999))
opt
# [1] 505.5884

# 7. Decode short and long term predictions from model in 6.

## AR(1) model

# x <- data$pow.obs_norm # TODO change this to residuals from linear model
# 
# plot(x, type="l")
# acf(x)
# 
# ## Moment estimates
# (phi.hat <- acf(x, plot=FALSE)$acf[2]) ## ACF in lag 1
# cor(x[-1], x[-length(x)])
# (sigma <- var(x) * (1-phi.hat^2))
# 
# 
# ## Likelihood estimation
# nll <- function(theta, y){
#   n <- length(y) - 1
#   n/2 * log(theta[1]) + 1/(2*theta[1]) * sum((y[-1]-theta[2]*y[-(n+1)])^2)
# }
# 
# ## MLE
# 
# # numerical
# (opt <- nlminb(c(1, 1/2), nll, y=x, lower=c(0.001 ,-0.999),upper=c(Inf, 0.999)))
# 
# # analytical
# (phi <- sum(x[-1] * x[-length(x)]) / sum(x[-length(x)] ^ 2))
# 
# ## standard errors
# library(numDeriv)
# V <- solve(hessian(nll, opt$par, y=x))
# (se <- sqrt(diag(V)))
# 
# ## Profile likelihood for phi
# llp.phi <- function(phi, x, x0){
#   n <- length(x) - 1
#   -n/2 * log(sum((x[-1]-phi*x[-(n+1)])^2))
# }
# 
# ## Plot profile likelihood
# phi <- seq(max(opt$par[2]-3*se[2],-1),
#            min(opt$par[2]+3*se[2],1),
#            length=200)
# 
# llp <- sapply(phi, llp.phi, x=x[-1], x0)
# 
# plot(phi, exp(llp - max(llp)),type="l")
# lines(range(phi), exp(-c(1,1) * qchisq(0.95,df=1)/2), col=2,lty=2)
# 
# ## Directly in R
# arima(x, order=c(1,0,0))
# opt$par ## rather close. 
# 
# ## Did it help?
# e <- x[-1] - opt$par[2] * x[-length(x)]
# acf(e)
# 
# 
# 
