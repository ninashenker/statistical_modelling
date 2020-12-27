#PROJECT 2: SURVIVAL DATA
library(LaplacesDemon)
library(survival)

# Setting working directory 
setwd("~/Desktop/statistical_modelling/Assignments/")

#Analysis of the binary data

# Read the data Logistic.txt into R.

data <- read.table("Logistic.txt", header = TRUE)
head(data)

#  AZT AIDS_yes   n
#  Yes   25      170
#  No    44      168

# Fit the Binomial distribution to the data (i.e. consider all data as coming from the same population).

ss <- 338
probs <- dbinom(x = 0:ss, size = ss, prob = (25+44)/ss)
barplot(probs,names.arg = 0:ss, xlab = "Patients with AIDS (X)", ylab = "Probability")


# Fit the Binomial separately to the two distributions and test if there is a difference between the groups.

y1 <- 25
n1 <- 170
y2 <- 44
n2 <- 168

prob1 <- dbinom(x = 0:n1, size = n1, y1/n1)
prob2 <- dbinom(x = 0:n2, size = n2, y2/n2)
par(mfrow=c(1,2))
barplot(prob1/max(prob1),names.arg = 0:n1)
barplot(prob2/max(prob2),names.arg = 0:n2)

# Test if there is a difference between groups
prop.test(c(y1,y2), c(n1, n2))

#OUTPUT:
#2-sample test for equality of proportions with continuity correction
#data:  c(y1, y2) out of c(n1, n2)
#X-squared = 6.171, df = 1, p-value = 0.01299
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.20593715 -0.02375473
#sample estimates:
#  prop 1    prop 2 
#0.1470588 0.2619048 

#p-value = 0.01299 <0.05 so we can say they are significantly different

## Log odds-ratio
theta <- log(y1/(n1-y1)*(n2-y2)/y2)
theta
#-0.721766

#Estimate parameters in the model (p0 probability of AIDS in control group, p1 probability of AIDS in treatment group) 
#and report a confidence interval for the parameter describing the difference, compare with the result above.


## Control group (AZTno)
x0 <- 44
n0 <- 168
p0 <- 44/168

## Treatment group (AZTyes)
x1 <- 25
n1 <- 170
p1 <- 25/170

## parameters
theta <- log (p1/(1-p1))
eta <- log (p0/(1-p0))

b0 = log(p0/(1-p0))

b1 = log(p1/(1-p1)) - log(p0/(1-p0))
# -0.7217
# Therefore there is a negative correlation when given the treatment i.e. when given AZT your log odds of having Aids are -0.721

p0 = exp(b0) / (1 + exp(b0))
p1 = exp(b0 + b1) / (1 + exp(b0 + b1))
# Projecting back into probability scale, p0 = 0.261, p1 = 0.1470, therefore, when given AZT you have a lower probability of having Aids 

loglik.fun <- function(b0, b1, y, n, t) {
  p = exp(b0 + b1 * t) / (1 + exp(b0 + b1 * t))
  prod(dbinom(y, n, p))
}

## profile likelihood
lp.theta <- function(theta, y, n, t){
  fun.tmp <- function(b0, b1, y, n, t){
    loglik.fun(b0, b1, y, n, t)
  }
  ## interval from plot
  opt <- optimise(fun.tmp, c(-20, 20), b1=theta, y=y, n=n, t=t, maximum=TRUE)
  opt$objective
}

t <- c(0, 1)
y <- c(x0, x1)
n <- c(n0, n1)
theta_b1 <- seq(-2, 2, 0.001)
theta_b1_profile_l <- sapply(theta_b1, lp.theta, y=y, n=n, t=t)
plot(theta_b1, theta_b1_profile_l / max(theta_b1_profile_l), type='l')
max <- max(theta_b1_profile_l)

chisq_l <- exp(-qchisq(0.95, df=1) / 2)
abline(h=chisq_l, col=2)
abline(v=0, col=3)

fun <- function(theta, y, n, t){
  lp.theta(theta, y, n, t) / max - chisq_l
}

theta_b1_shifted_l <- sapply(theta_b1, fun, y=y, n=n, t=t)
plot(theta_b1, theta_b1_shifted_l, type='l')

CI <- c(uniroot(f = fun, interval = c(-2, -1), y=y, n=n, t=t)$root,  
        uniroot(f = fun, interval = c(-0.5, 1), y=y, n=n, t=t)$root)
par(mfrow=c(1,1))
plot(theta_b1, theta_b1_profile_l / max(theta_b1_profile_l), type='l', xlab = "AZT yes coefficient", ylab = "Likelihood", main = "Profile likelihood for AZT-yes coefficient")
max <- max(theta_b1_profile_l)

text(CI[1] + 0.2, chisq_l + 0.05, round(CI[1], 2))
text(CI[2] + 0.2, chisq_l + 0.05, round(CI[2], 2))
text(b1 + 0.3, 1, round(b1, 2))
abline(h=chisq_l,col=2)
points(CI[1], chisq_l, col=2, type = 'p', pch = 16)
points(CI[2], chisq_l, col=2, type = 'p', pch = 16)
points(b1, 1, type = 'p', pch = 16)


# PART 2

#Fit a logistic regression model for the binary outcome AIDS=”yes” versus AIDS=”no” 
#with the explanatory variable treatment with AZT (Yes, NO). 

logit <- glm(cbind(data$AIDS_yes, data$n - data$AIDS_yes) ~ data$AZT, data = data, family = "binomial")
summary(logit)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.0361     0.1755  -5.904 3.54e-09 ***
# data$AZTYes  -0.7218     0.2787  -2.590  0.00961 ** 

coef(logit)
exp(-1.036092)
exp(-0.721766)

#Present the odds ratio for the effect of AZT on AIDS with 95% confidence interval and interpret the result in words.

predict(logit, type="response")

# Traditional CI (not using likelihood)
exp(-1.0361)
exp(-0.7218)
# -------------
exp(summary(logit)$coefficients[2,1] + 1.96 * summary(logit)$coefficients[2,2]) # upper
exp(summary(logit)$coefficients[2,1] - 1.96 * summary(logit)$coefficients[2,2]) # lower
exp(summary(logit)$coefficients[1,1] + 1.96 * summary(logit)$coefficients[1,2]) # upper
exp(summary(logit)$coefficients[1,1] - 1.96 * summary(logit)$coefficients[1,2]) # lower
#               2.5 %    97.5 %
#(Intercept) 0.2515723 0.5004942
#data$AZTYes 0.2813743 0.8390689

# Traditional CI with confint
exp(confint.default(logit))
#               2.5 %    97.5 %
#(Intercept) 0.2515723 0.5004942
#data$AZTYes 0.2813743 0.8390689

# Wald CI for profile likelihood (using quadratic apporximation of likelihood)
# -----------------------------
exp(confint(logit)) 
#               2.5 %    97.5 %
#(Intercept) 0.2489862 0.4962491
#data$AZTYes 0.2782713 0.8330139

# LR Test (raio of null model and test model likelihoods * 2)
-------------------------------
(x <- lrtest(logit))
LRTest <- 2 * (x$LogLik[1] - x$LogLik[2])
LRTest # 6.929332
anova(logit, test = "LRT") # LR = 6.929332, Pr(>chi) = 0.008479
(lr_pval <- 1 - pchisq(LRTest, 1)) # Convert LR Test to P value
LRTest < qchisq(p = 0.95, df = 1) # LRTest < Pr(chi)

# TODO to get confidence intervals we need to plot 
# likelihood with qchisq line and find roots. Alternatively 
# we can create a for loop and just search over an interval
# until we reach the chisquare bound.

# Score Test
anova(logit, test="Rao")
-------------------------------
# score test = 6.8597, Pr(>Chi) = 0.008816

#It is the estimated amount by which the log odds of AIDS would decrease if AZT is given. The log odds of leaves.presence when AIDs is 0 is just above in the first row.

#likelihood function

nll <- function(params, xi, y, n) {
  b_0 <- params[1]
  b_1 <- params[2]
  -sum((b_0 + b_1*xi) * y - log(1 + exp(b_0 + (b_1 * xi))))
}

p<- (data$AIDS_yes / data$n) 
opt <- nlminb(c(0.1,0.99), nll, xi=c(1,0), y = p, n=data$n)
opt

library(numDeriv)
H <- hessian(nll, opt$par, n=data$n, y=p, xi=c(1,0)) ## Negative Hessian

# From stackexchange (we think it's a Wald CI)

# Estimate se from hessian (positive hessian because 
# we're minimizing negative loglikelihood), then find CI 
# using se and parameters.
variance <- solve(H)
se <- sqrt(diag(variance))

# Summary stats (params, stderror, Wald z statisti)
tab <- cbind(opt$par, se, z.stat=opt$par/se)
rownames(tab) <- c("intercept","AZT_coef")
tab

# Wald CI for both parameters
upper <- exp(opt$par + 1.96 * se)
lower <- exp(opt$par - 1.96 * se)
# c(lower[1], opt$par[1], upper[1])
c(lower[2], opt$par[2], upper[2])


# Profile likelihood of B1
pll<- function(b_1, xi, y, n){
  fun <- function(b_0, b_1, xi, y, n){
    nll(params=c(b_0, b_1), xi=xi, y=y, n=n)
  }
  optimise(fun, c(-10, 10), b_1=b_1, xi=xi, y=y, n=n)$objective
}

theta <- seq(-10, 8, by=0.01)
profile <- sapply(theta, pll, xi=c(1, 0), y = p, n=data$n)
pp <- plot(theta, exp(-(profile - min(profile))), type="l")


#Test the hypothesis of no effect of AZT on AIDS using:

# 1. The Wald test
## Wald test:
H <- hessian(pll, opt$par[2], n=data$n, y=p, xi=c(1,0)) ## Negative Hessian
sd.theta <- as.numeric(sqrt(1/H)) # stderr
W <- opt$par[2] / sd.theta
1 - pchisq(W^2, df = 1)


# 2. The likelihood ratio test

(chisq <- 2 * (-opt$objective - (-pll(theta, xi=c(1, 0), y = p, n=data$n))))
pchisq(chisq, df=1, lower.tail=FALSE)

# 3. The score test 
z.exp <- score(lambda0,y)/sqrt(ExpInfo(lambda0,y))
(p.score.exp <- 1-pchisq(z.exp^2,df=1)) ## p-value


## Score test (using observed information)
#library(numDeriv)
#(I <- -hessian(ll,theta0,y=y))
#ObsInfo(theta0,y)
#(z.obs2 <- score(theta0,y)/sqrt(I))
#(p.score.obs2 <- 1-pchisq(z.obs2^2,df=1))## P-value







# Setting working directory 
setwd("~/Desktop/statistical_modelling/Assignments/")


# Analysis of the survival time data

# Read the data actg320.txt into R. "Import Dataset" button.

actg320 <- read.table("actg320.txt", header = TRUE)

# How many patients got AIDS or died in the two treatment groups?
event_count <- count(actg320, vars = 'event')
event_count
#1151

frequencies <- table(actg320['event'])
# 0      1 
# 1055   96

#What is the proportion of patients that got AIDS or died in the two group? 

x2 <- frequencies[2]
y2 <- unname(x2)
print(y2)

x1 <- frequencies[1]
y1 <- unname(x1)
print(y1)

proportion_with_AIDS = y2/(y2+y1)
print(proportion_with_AIDS)
#0.08340573

#Other relevant number that could be calculated?
#How long was the total follow-up time in the two groups?
max(actg320$time)
#364 so about a year
  
# Fit an exponential distribution, using numerical methods, to the time of event (time) in the data set, 
# remember to take into account that some of the data is censored (i.e. we only know that the time to the event is longer that the reported time).
#1: Using all data (i.e. ignore the treatment effect)

#FIT EXPONENTIAL MODEL
mod1 <- survreg(Surv(actg320$time) ~ event, data = actg320,
                dist="exponential")
summary(mod1)
confint(mod1)

#2: Separately for the two treatments
mod2 <- survreg(Surv(actg320$time) ~ tx + event, data = actg320,
                dist="exponential")
summary(mod2)
confint(mod2)

# Compare the likelihood for the above models and conclude

## Exponential regression model
## likelihood function
nll.exp <- function(theta,t,d,X){
  mu <- exp(X %*% theta)
  - sum(d * dexp(t, rate = 1/mu, log=TRUE) + (1-d) *
          pexp(t, rate = 1/mu, lower.tail = FALSE, log.p=TRUE))
}

# Design matrix
X <- cbind(1,actg320$tx==1)
head(X)
t <- actg320$time
d <- actg320$event

## Find MLE
opt.exp <- nlminb(c(0,0), nll.exp, t=t, d = d, X=X)

par(mfrow=c(1,1))

mu <- seq(7, 8.5, by=0.01)
kappa <- seq(-0.5, 1.5, by=0.01)
ll <- matrix(nrow = length(mu), ncol = length(kappa))

## Find the values of the log likelihood
for(i in 1:length(mu)){
  for(j in 1:length(kappa)){
    ll[i,j] <- nll.exp(c(mu[i], kappa[j]), t=t, d = d, X=X)
  }
}

Nll <- -ll + opt.exp$objective
## Define the levels
alpha <- c(0.5, 0.25, 0.1, 0.05, 0.01)
contour(x = mu, y = kappa, z = Nll,
        levels = -qchisq(1 - alpha, df = 2) / 2,
        labels = alpha, xlab="mu", ylab="kappa")



# Formulate a model where one parameter indicate the treatment effect and find the MLE and compare with the result above.

#survival <- glm(family = Gamma(link="log")))

# Find the Wald confidence interval for the treatment parameter in the model above.
#Derive the theoretical results for the models above, including the standard error estimates, use this to formulate and implement the profile likelihood function for the treatment parameter



#PART TWO
#Descriptive statistics
setwd("~/Desktop/statistical_modelling/Assignments/")
actg320 <- read.table("actg320.txt", header = TRUE)
#Plot the survival functions in the two treatment groups, which group seems to be doing best?
## Non parametric analysis (Kaplan Meier)
library(survival)
dat <- data.frame(t=actg320$time, treatment = actg320$tx, event = actg320$event, cd4 = actg320$cd4)
Surv <- survfit(Surv(t, event) ~ treatment, conf.type = "log-log",
                   data = dat)




## Kaplan-Meier plot
par(mfrow=c(1,1))
plot(Surv, conf.int = FALSE, las = 1, xlab = "Survival time in Days", 
     ylab = "Estimated Survival Prob.", main="Survival Distributions by Treatment", col=3:4, lwd = 2, mark.time = TRUE)
legend("bottomleft", col = 3:4, c("2 drug treatment","3 drug treatment"), lwd = 2)





#Plot the cumulative incidence functions for the two groups, which plot would you prefer?


#PLOT THE CUMULATIVE INCIDENCE
plot(Surv, fun=function(x) { 1- x }, conf.int = F,
     las = 1, xlab = "Days",
     ylab = "Estimated Failure Prob.", main = "Cumulative Incidence Functions", col=3:4, lwd = 2)
legend("topleft", col=3:4, c("2-drug treatment","3-drug treatment"), lwd=2)


#Compare the survival in the two treatment groups using a log-rank test.
#LOG RANK TESTING THE SURVIVAL CURVES ARE THE SAME
survdiff(Surv(t, event) ~ treatment, data = dat)

#survdiff(formula = Surv(t, event) ~ treatment, data = dat)

#N Observed Expected (O-E)^2/E (O-E)^2/V
#treatment=0 577       63     47.1      5.37      10.5
#treatment=1 574       33     48.9      5.17      10.5

#Chisq= 10.5  on 1 degrees of freedom, p= 0.001 

#SIgnificant so that the curves are indeed NOT the same meaning that the drug makes a significant difference

#Parametric survival models
#Fit parametric survival models containing treatment (tx) and CD4 count (cd4) as explanatory variables.
#Try using the exponential, Weibull and log-logistic models, which one gave the best fit (and why)?

#FIT EXPONENTIAL MODEL
mod2 <- survreg(Surv(t, event) ~ treatment + cd4, data = dat,
                dist="exponential")
summary(mod2)
confint(mod2)

#FIT WEIBULL MODEL
mod3 <- survreg(Surv(t, event) ~ treatment + cd4, data = dat,
                dist="weibull")
summary(mod3)
confint(mod3)

#FIT LOG-LOGISTIC
mod4 <- survreg(Surv(t, event) ~ treatment + cd4, data = dat,
                dist = "loglogistic")
summary(mod4)
confint(mod4)

#compare AIC
AIC(mod2)
#1645.838
AIC(mod3)
#1640.671
AIC(mod4)
#1639.655

#AIC for loglogistic is the lowest. 

#Using the survival model you chose, make a table of estimates and their 95% confidence intervals.



#Using your model compute the time ratio for the treatment effect. Similarly, compute the time ratio for the effect of increasing the CD4 count
#with 50. In both cases uncertainty evaluation (e.g. confidence intervals) should be included. Interpret the results in words.

exp(cbind(coef(mod4),confint(mod4)))

#                           2.5 %      97.5 %
#(Intercept) 921.350304 559.460759 1517.329622
#treatment     2.323208   1.316460    4.099857
#cd4           1.021019   1.013541    1.028553


#The median time to survival time for new treatment is 2.32 times the median for old treatment for the same load 95% CI (1.32- 4.10)

exp(c(0.02080,0.01344965, 0.02815334)*50)
#[1] 2.829217 1.959095 4.086411
#The median time to failure for increasing CD4 with 50 is 2.82 times the median for the same bearing with 95% CI 1.959095-4.086411.

#log odds
exp(coef(mod4)[3]/mod4$scale)
#The odds of cd4 new are 1.017 times the odds for cd4 old with the same load.

#Assess the goodness of fit of this model using a plot based on the Cox snell residuals.
dat$z <- (log(dat$t)-mod4$linear.predictors) / mod4$scale
dat$CS4 <- log(1+exp(dat$z))
surv4 <- survfit(Surv(dat$CS4, dat$event==1)~1 , data = dat)
plot(surv4$time, -log(surv4$surv))
abline(a=0, b=1, lty=2)


#Give a graphical presentation of your model

