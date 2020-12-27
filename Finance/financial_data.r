library(tidyverse)
library(fitdistrplus)
library(car)

#PROJECT 3: FINANCIAL DATA


# Setting working directory 
setwd("~/Desktop/statistical_modelling/Assignments/")

finance_data <- read.table('finance_data.csv', header = TRUE, sep =';')

head(finance_data)
hist(finance_data$SLV)
qqnorm(finance_data$SLV)
qqline(finance_data$SLV)

finance_data[which.min(finance_data$SLV), ]
finance_data[which.max(finance_data$SLV), ]

head(finance_data[order(finance_data$SLV), ])
tail(finance_data[order(finance_data$SLV), ])

normal_ll <- function(mu, sigma, y) {
  sum(dnorm(y, mean=mu, sd=sigma, log=TRUE))
}

mu_theta <- seq(-0.5, 0.5, by=0.01)
normal_ll_at_theta <- sapply(mu_theta, normal_ll, sigma=sd(finance_data$SLV), y=finance_data$SLV)
plot(mu_theta, exp(normal_ll_at_theta / max(normal_ll_at_theta)), type="l")

#The Cauchy model is useful as a model for data with heavy tails, characterized by the presence of outliers
cauchy_model <- function(x, location, scale){
  dcauchy(x, location = mle.location, scale = mle.scale, log = FALSE)
}

mle_location <- optimize(cauchy_model, c(-3, 3), x = finance_data$SLV, location = 1, scale = 0, maximum=TRUE)

x<- finance_data$SLV
fitdist(x, "t", start= list(df=length(x), ncp=mean(x)))
plot(fit1)

df <- 1
t <- function(x, df ){
prod(dt(x, df, log = FALSE))
}
 t <- optimise(t, c(-3,3), x= finance_data$SLV)$maximum
 


theta <- seq(-0.5, 0.5, by=0.01)
cm <- sapply(theta, cauchy_model, x=finance_data$SLV)
plot(theta, cm, type="l")
length(theta)
length(cm)

mle_slv.mu <- optimize(normal_ll, c(-3, 3), sigma=sd(finance_data$SLV), y=finance_data$SLV, maximum=TRUE)
mle_slv.mu

abline(v=mle_slv.mu)

sd(finance_data$SLV)

par(mfrow=c(2, 1))
qqPlot(finance_data$SLV, distribution = "norm")
qqPlot(finance_data$SLV, distribution = "t", df=2)

# TODO try modelling as binomial with either positive or negative return as a proportion


### Summary statistics ###
mean(finance_data$SLV)
#0.00146765
sd(finance_data$SLV)
#0.04830741

summary(finance_data$SLV)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.238893 -0.026350  0.002226  0.001468  0.033122  0.267308 

ggplot(data = finance_data, mapping = aes(x =SLV)) +
  geom_histogram(fill = 'steelblue', color ="black", bins = 100)




#normMLE <- fitdist(finance_data$SLV, "norm", method = "mle")
#summary(normMLE)


finance_data$SLV_normalized <- (finance_data$SLV / max(finance_data$SLV))
qqplot <- qqnorm(finance_data$SLV_normalized)
plot(qqplot)

slv_len <- count(finance_data)

y <- finance_data$SLV_normalized
n <- slv_len$n
x <- seq(1, n)

eps = 1e-5
mu = seq(-min(y) + eps, max(y), length = n)
lambda = seq(0, 5, length = n)
lik = matrix(0, n, n)
for (i in 1:n) lik[, i] = boxcox((y + mu[i]) ~ x, lambda = lambda, plotit = FALSE)$y
image(lik, xlab = "mu", ylab = "lambda", main = "likelihood")


#There are a couple outliers above and below so it is not normal

#The probability that a random variable deviates by a given amount from its 
#expectation is referred to as a tail probability.

#plot(normMLE)



  





#PART 2
# Setting working directory 
setwd("~/Desktop/statistical_modelling/Assignments/")

finance_data <- read.table('finance_data.csv', header = TRUE, sep =';')
par(mfrow=c(1,1))
plot(hist(finance_data$SLV, breaks = 454))

# Mixture model
finance_data$date <- strptime(finance_data$time, format='%Y-%m-%d')

#a) With the same data as in Assignment 1, fit normal mixture models with 2 and 3 components, 
#argue which one is better and compare with the best model from assignment 1.
plot(finance_data$date, finance_data$SLV, type='l')

## Normal mixture: transform
## natural to working parameters
normal.mix.pn2pw <- function(n.dist, mu, sd, delta){
  if(n.dist == 1){
    return (list(mu=mu, sd=sd, tau=1))
  }
  if(sum(delta) >= 1){
    print("sum(delta) should be < 1")
    return()
  }
  if(length(mu) != n.dist){
    print("length(mu) should be n.dist")
    return()
  }
  if(length(sd) != n.dist){
    print("length(sd) should be n.dist")
    return()
  }
  if(length(delta) != (n.dist-1)){
    print("length(delta) should be n.dist-1")
    return()
  }
  tau <- log(delta/(1-sum(delta)))
  return(list(mu=mu, sd=sd, tau=tau))
}

## Normal mixture: transform
## working to natural parameters
normal.mix.pw2pn <- function(n.dist, mu, sd, tau){
  if(n.dist == 1){
    return (list(mu=mu, sd=sd, delta=1))
  }
  if(length(mu) != n.dist){
    print("length(mu) should be n.dist")
    return()
  }
  if(length(sd) != n.dist){
    print("length(sd) should be n.dist")
    return()
  }
  if(length(tau) != (n.dist-1)){
    print("length(delta) should be n.dist-1")
    return()
  }
  delta <- exp(tau)/(1+sum(exp(tau)))
  delta <- c(1-sum(delta),delta)
  return(list(mu=mu, sd=sd, delta=delta))
}

## negative log likelihood function
nll <- function(theta, n.dist, y){
  mu <- theta[1:n.dist]
  sd <- theta[(n.dist+1):(2*n.dist)]
  tau <- theta[(2*n.dist + 1):(3*n.dist-1)]
  n.pars <- normal.mix.pw2pn(n.dist, mu, sd, tau)
  n <- length(y)
  nll <- 0
  for(i in 1:n){
    nll <- nll - log(sum(n.pars$delta * dnorm(y[i], mean=n.pars$mu, sd=n.pars$sd)))
  }
  return(nll)
}

# Estimation with one distribution
n.dist <- 1
mu <- c(mean(finance_data$SLV))
sd <- c(sd(finance_data$SLV))
delta <- c(1.0)
wpars <- normal.mix.pn2pw(n.dist, mu, sd, delta)
theta <- c(wpars$mu, wpars$sd, wpars$tau)
(opt1 <- nlminb(theta, nll, n.dist=n.dist, y=finance_data$SLV))
(npars1 <- normal.mix.pw2pn(n.dist, opt1$par[1:n.dist], opt1$par[(n.dist+1):(2*n.dist)], opt1$par[(2*n.dist+1):(3*n.dist-1)]))
# $mu 0.001467675, $sd 0.04825421 $delta 1

## Estimation with two distributions
n.dist <- 2
mu <- mean(finance_data$SLV) * c(1/2, 3/2)
sd <- sd(finance_data$SLV) * c(1/2, 3/2)
delta <- c(0.5)
wpars <- normal.mix.pn2pw(n.dist, mu, sd, delta)
theta <- c(wpars$mu, wpars$sd, wpars$tau)
(opt2 <- nlminb(theta, nll, n.dist=n.dist, y=finance_data$SLV))
(npars2 <- normal.mix.pw2pn(n.dist, opt2$par[1:n.dist], opt2$par[(n.dist+1):(2*n.dist)], opt2$par[(2*n.dist+1):(3*n.dist-1)]))
# $mu 0.00394719 -0.02514932 $sd 0.04046551 0.09471918 $delta 0.9147822 0.0852178

## Estimation with three distributions
n.dist <- 3
mu <- mean(finance_data$SLV) * c(1/2, 1, 3/2)
sd <- sd(finance_data$SLV) * c(1/2, 1, 3/2)
delta <- c(1, 1)/3
wpars <- normal.mix.pn2pw(n.dist, mu, sd, delta)
theta <- c(wpars$mu, wpars$sd, wpars$tau)
(opt3 <- nlminb(theta, nll, n.dist=n.dist, y=finance_data$SLV))
(npars3 <- normal.mix.pw2pn(n.dist, opt3$par[1:n.dist], opt3$par[(n.dist+1):(2*n.dist)], opt3$par[(2*n.dist+1):(3*n.dist-1)]))
#$mu 0.0485287351  0.0004392401 -0.0197941796 $sd 0.02136200 0.03866547 0.09172425 $delta 0.06457801 0.83276478 0.10265721

c(opt1$objective,opt2$objective,opt3$objective)
#[1] -731.9998 -749.8222 -750.1281

(AIC <- 2 * c(opt1$objective, opt2$objective,
              opt3$objective) +
    2 * c(length(opt1$par), length(opt2$par),
          length(opt3$par)))

## Small function to calculate mixture dist
mix.dist <- function(x, npars, idx) {
  sum(npars$delta[idx] * dnorm(x, mean=npars$mu[idx], sd=npars$sd[idx]))
}

par(mfrow=c(1,1))
hist(finance_data$SLV, breaks=100)
# abline(v=npars2$mu[1], col=2)
# abline(v=npars2$mu[2], col=3)

x <- seq(-0.5, 0.5, by=0.001)
gaussian_mix_1 <- sapply(x, mix.dist, npars=npars2, idx=1)
lines(x, gaussian_mix_1 * 3, type="l", col=2, pch=2)
gaussian_mix_2 <- sapply(x, mix.dist, npars=npars2, idx=2)
lines(x, gaussian_mix_2 * 3, type="l", col=3, pch=2)

hist(finance_data$SLV, breaks=100)
# abline(v=npars3$mu[1], col=2)
# abline(v=npars3$mu[2], col=3)
# abline(v=npars3$mu[3], col=4)

x <- seq(-0.5, 0.5, by=0.001)
gaussian_mix_1 <- sapply(x, mix.dist, npars=npars3, idx=1)
lines(x, gaussian_mix_1 * 3, type="l", col=2, pch=2)
gaussian_mix_2 <- sapply(x, mix.dist, npars=npars3, idx=2)
lines(x, gaussian_mix_2 * 3, type="l", col=3, pch=2)
gaussian_mix_3 <- sapply(x, mix.dist, npars=npars3, idx=3)
lines(x, gaussian_mix_3 * 3, type="l", col=4, pch=2)

#b) For the best mixture model, report confidence interval for the parameters, and give an interpretation of these intervals.
#c) For the two component model make a profile likelihood plot of one of the variance parameters.
#d) In the previous question you should see multiple maxima, reprametrize the model such that you only see one maximum.
#e) Discuss the interpretaion of the models.


#Hidden Markov Models
#a) Fit two and three state normal Hidden Markov Models to the data and conclude on the best choice
norm.HMM.pn2pw <- function(m, mu, sd, gamma) {                                              
  tgamma  <- NULL
  
  # Transform sd to unconstrained parameter as sd > 0
  sd <- log(sd)
  if (m > 1) {
    
    # Transform gamma(mxm) to unconstrained parameter 
    # This is just softmax inverse
    # gamma:
    # [,1] [,2]
    # [1,] 0.95 0.05
    # [2,] 0.05 0.95
    gamma_div_diag <- gamma / diag(gamma)
    # [,1]       [,2]
    # [1,] 1.00000000 0.05263158
    # [2,] 0.05263158 1.00000000
    foo <- log(gamma_div_diag) 
    # [,1]      [,2]
    # [1,]  0.000000 -2.944439
    # [2,] -2.944439  0.000000
    tgamma<- as.vector(foo[!diag(m)]) # Selects the non-diag elements 
    # [1] -2.944439 -2.944439
  }                                             
  parvect <- c(mu, sd, tgamma)                    
  parvect                                         
}  

norm.HMM.pw2pn <- function(m, parvect) { 
  mu <- parvect[1:m]
  
  # Transform sd to constrained param, sd > 0
  sd <- exp(parvect[(m+1):(2*m)])
  gamma <- diag(m) # Identity matrix with dim m
  if(m > 1) {  
    
    # Transform gamma(mxm) to constrained parameter (i.e. sum(gamma_ij)) == 1)
    # This is just a softmax.
    gamma[!gamma] <- exp(parvect[(2*m+1):(m*m)])
    #            [,1]       [,2]
    # [1,] 1.00000000 0.05263158
    # [2,] 0.05263158 1.00000000
    
    gamma <- gamma / apply(gamma, 1, sum)
    # [,1] [,2]
    # [1,] 0.95 0.05
    # [2,] 0.05 0.95
  }                                                
  
  # Delta is stationary distribution. Solve for δ(Im − Γ + U) = 1
  delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  list(mu=mu, sd=sd, gamma=gamma, delta=delta)           
}  

norm.HMM.mllk <- function(parvect, x, m) {
  pn <- norm.HMM.pw2pn(m, parvect)
  if (m==1) return (-sum(dnorm(x, mean=pn$mu, sd=pn$sd, log=TRUE))) 

  n <- length(x)
  phi <- pn$delta * dnorm(x[1], mean=pn$mu, sd=pn$sd)
  phi_sum <- sum(phi)
  lscale <- log(phi_sum)
  phi <- phi / phi_sum
  for (i in 2:n) {
    phi <- phi %*% pn$gamma * dnorm(x[i], mean=pn$mu, sd=pn$sd)
    phi_sum <- sum(phi)
    lscale <- lscale + log(sum(phi_sum))
    phi <- phi / phi_sum
  }
  mllk <- -lscale
  mllk
}

norm.HMM.mle <- function(x, m, mu, sd, gamma) {                                                      
  parvect0  <- norm.HMM.pn2pw(m, mu, sd, gamma) 
  mod <- nlm(norm.HMM.mllk, parvect0, x=x, m=m, stepmax=1, iterlim=300)       
  pn <- norm.HMM.pw2pn(m, mod$estimate)            
  mllk <- mod$minimum                              
  np <- length(parvect0)                          
  AIC <- 2 * (mllk + np)                              
  n <- sum(!is.na(x))                           
  BIC <- 2 * mllk + np * log(n)                         
  list(mu=pn$mu, sd=pn$sd, gamma=pn$gamma, delta=pn$delta,   
       code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC)   
}  

y <- finance_data$SLV
m <- 1
mu0 <- mean(y)
sd0 <- sd(y)
gamma0 <- 0

## optimize
fit1 <- norm.HMM.mle(x=y, m, mu0, sd0, gamma0)
fit1
# $mu [1] 0.001467152
# $sd [1] 0.04825411
# $mllk [1] -731.9998
# $AIC [1] -1460

m <- 2
mu0 <- quantile(y,c(0.25,0.75)) # TODO sweep multiple init values
sd0 <- c(0.8, 1.2)
gamma0 <- matrix(0.45, ncol=m, nrow=m)
diag(gamma0) <- 1 - (m-1) * gamma0[1,1]

## optimize
fit2 <- norm.HMM.mle(y, m, mu0, sd0, gamma0)
fit2

m <- 3
mu0 <- quantile(y,c(0.2, 0.5, 0.8))
sd0 <- c(0.8, 1.2)
gamma0 <- matrix(0.45, ncol=m, nrow=m)
diag(gamma0) <- 1 - (m-1) * gamma0[1,1]

## optimize
fit3 <- norm.HMM.mle(y, m, mu0, sd0, gamma0)
fit3

par(mfrow=c(1, 1))
plot(finance_data$date, finance_data$SLV, type='l')

# fit1
abline(h=fit1$mu, col=2)

# fit2
abline(h=fit2$mu[1], col=3)
abline(h=fit2$mu[2], col=3)

# fit3
abline(h=fit3$mu[1], col=4)
abline(h=fit3$mu[2], col=4)
abline(h=fit3$mu[3], col=4)

(AIC <- c(fit1$AIC, fit2$AIC, fit3$AIC))
# [1] -1459.9989 -1499.3636  -646.7287
# [1] -1460.000 -1499.364 -1415.361

(ll <- -c(fit1$mllk, fit2$mllk, fit3$mllk))

m <- c(1, 2, 3)
df <- m + m + (m ^ 2 - m) ## mu + sigma + gamma
cbind(df, ll, AIC)

##################################################
## Working with the 2-state mmodel

## Finding the standard errors
m <- 2
parvect  <- norm.HMM.pn2pw(m, fit2$mu, fit2$sd, fit2$gamma)
mod <- nlm(norm.HMM.mllk, parvect, x=y, m=m, hessian=TRUE)  
mod

## Organize the result
parvect <- mod$estimate
names(parvect) <- c("mu1","mu2","sd1", "sd2", "tau21","tau12")

# Calculating se from hessian causes nans. This might because
# some of the parameters (e.g. gamma) is on the edge of the parameter
# space.
# se <- sqrt(diag(solve(mod$hessian)))
## Working pars + standard error
# round(cbind(parvect, se), digits=2) ## note se of tau31
# fit2$gamma

norm.HMM.generate_sample <- function(n, m, mu, sd, gamma, delta=NULL) {                                                           
    if(is.null(delta)) delta<-solve(t(diag(m) - gamma + 1), rep(1,m))    
    mvect <- 1:m                                               
    state <- numeric(n)                                        
    state[1] <- sample(mvect, 1, prob=delta)                     
    for (i in 2:n)                                             
      state[i] <- sample(mvect, 1, prob=gamma[state[i - 1],])
    x <- rnorm(n, mean=mu[state], sd=sd[state])                      
    x                                                          
}       

n <- length(y)
k <- 100
m <- 2
mu0 <- quantile(y, c(0.25, 0.75))
sd0 <- c(0.8, 1.2)
gamma0 <- matrix(0.45, ncol=m, nrow=m)
diag(gamma0) <- 1 - (m-1) * gamma0[1,1]

## Matrices to stores bootstap res.
GAMMA <- matrix(ncol=m * m, nrow=k)
Mu <- matrix(ncol=m, nrow=k)
Sd <- matrix(ncol=m, nrow=k)
Delta <- matrix(ncol=m, nrow=k)
Code <- numeric(k)

## Boot strap with two different optimizers
for(i in 1:k){
  set.seed(i)
  y.sim <- norm.HMM.generate_sample(n, m, fit2$mu, fit2$sd, fit2$gamma)
  ## fit model to sample
  fit2.tmp <- norm.HMM.mle(y.sim, m, mu0, sd0, gamma0)
  ## Store result
  GAMMA[i, ] <- c(fit2.tmp$gamma[1, ], fit2.tmp$gamma[2, ])
  Mu[i, ] <- fit2.tmp$mu
  Sd[i, ] <- fit2.tmp$sd
  Delta[i, ] <- fit2.tmp$delta
  Code[i] <- fit2.tmp$code
  print(c(i, Code[i]))
}
sum(Code!=1)

## Plot the results (i.e. statistical distribution of
## estimates
par(mfrow=c(1,2))
hist(Mu[ ,1])
rug(fit2$mu,lwd=2,col=2)
hist(Mu[ ,2])
rug(fit2$mu,lwd=2,col=2)

par(mfrow=c(1,2))
hist(Sd[ ,1])
rug(fit2$sd, lwd=2,col=2)
hist(Sd[ ,2])
rug(fit2$sd, lwd=2,col=2)

par(mfrow=c(2,2))
hist(GAMMA[ ,1])
rug(fit2$gamma,lwd=2,col=2)
hist(GAMMA[ ,2])
rug(fit2$gamma,lwd=2,col=2)
hist(GAMMA[ ,3])
rug(fit2$gamma,lwd=2,col=2)
hist(GAMMA[ ,4])
rug(fit2$gamma,lwd=2,col=2)

## Confidence intervals (90%)
## 95% CI for Mu and Sd
apply(Mu, 2, quantile,prob=c(0.025,0.975))
apply(Sd, 2, quantile,prob=c(0.025,0.975))

## 95% CI for gamma
t(round(apply(GAMMA, 2, quantile,prob=c(0.025,0.975)), digits=3))

## 95% CI for delta
round(t(apply(Delta, 2, quantile, prob=c(0.025,0.975))),digits=3)

## Correlation of Mu and Sd
round(cov2cor(cov(Mu)), digits=2)
round(cov2cor(cov(Sd)), digits=2)

## Plot the interdepence
plot(data.frame(Mu))
plot(data.frame(Sd))
plot(data.frame(Delta))


#b) For the chosen model report 95% confidence intervals for the working parameters.


#c) Report the natural parameters and interpret the result.
#d) Compare the following distributions (by plots)
#The long term distribution of the return.
#The 1-step ahead distribution given that you known the current state.
fit2$delta * fit2$gamma
plot(fit2$delta * fit2$gamma)


#e) Report 95% confidence intervals for some (or all) natural parameters (note that the natural paramters include the stationary distribution).
#Some options for finding these CI’s are
#– Formula (3.2) in the (HMM) textbook.
#– The bootstrap method in Section 3.6.2.
#– Profile likelihood
#f) Discuss what would be needed in order to make short term (say 1-2 weeks) predictions.
par(mfrow=c(1, 1))
## depmix
library(depmixS4)
hmm <- depmix(returns ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=finance_data$SLV))
hmmfit <- fit(hmm, verbose = FALSE)

post_probs <- posterior(hmmfit)

plot(post_probs$state, type='s', main='True Regimes', xlab='', ylab='Regime')

matplot(post_probs[,-1], type='l', main='State Probability', ylab='Probability')

legend(x='topright', c('Bull','Bear'), fill=1:2, bty='n')




