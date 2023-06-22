############################################################################
## PoE - Chapter 12 - Regression with Time-Series Data: non-stationary Variables
## Author: Robson Tigre (T.A)
## Created: Nov 26 2017
## Last modified: Dec 15 2017
## R version 3.4.0 (2017-04-21) and RStudio 1.0.143 on Windows 10
############################################################################

# Before sending an email or posting a question to Piazza regarding the content of this script, please:
  ## Try to run line by line to get what is going on in a context; 
  ## Pay attention to the comments. Answers to your question may be a couple of lines below or above your current line;
  ## Check the help documentation for the function you are trying to understand.
# If the question persists after you've followed those steps, please don't hesitate to discuss it.

# Main sections in this script are:
	# 12.1 STATIONARY AND NON-STATIONARY VARIABLES
	# 12.2 SPURIOUS REGRESSIONS 
  # 12.3 UNIT ROOT TESTS FOR STATIONARITY
	# 12.4 COINTEGRATION


# Loading the data
rm(list=ls())
if (!require(devtools)) {         # if R can't load and attach the devtools package (i.e. !require)...
  install.packages("devtools")    # ...then R must install this package...
  library(devtools)               # ...and then try to load and attach it again. Same for other packages below.
}
if (!require(PoEdata)) {
  devtools::install_git("https://github.com/ccolonescu/PoEdata")
  library(PoEdata)
}

data("usa")
str(usa) # gdp "real US gross domestic product"; inf "annual inflation rate"; f "federal funds rate"; b "3-year Bond rate"
is.ts(usa) # FALSE
usa.ts <- ts(usa)

# Instead of referring to variables as usa.ts[,"gdp"], I prefer to store them as individual objects. Do as you prefer.
gdp <- usa.ts[,"gdp"]           # gdp
inf <- usa.ts[,"inf"]           # inf
f <- usa.ts[,"f"]               # f
b <- usa.ts[,"b"]               # b
D1_gdp <- diff(usa.ts[,"gdp"])  # first difference of gdp
D1_inf <- diff(usa.ts[,"inf"])  # first difference of inf
D1_f <- diff(usa.ts[,"f"])      # first difference of f
D1_b <- diff(usa.ts[,"b"])      # first difference of b


#######################################################################################################
# 12.1 STATIONARY AND NON-STATIONARY VARIABLES
#######################################################################################################

# Formally, a time series y_t is stationary if its MEAN AND VARIANCE are CONSTANT OVER TIME, and if the COVARIANCE BETWEEN TWO VALUES FROM THE SERIES DEPENDS ONLY ON THE LENGTH OF TIME SEPARATING THE TWO VALUES, and not on the actual times at which the variables are observed:
  # E(y_t) = \mu (constant mean - doesn't depend on t) - Non-stationary series with non-constant means are often described as not having the property of mean reversion. That is, stationary series have the property of mean reversion.
  # var(y_t) = \sigma^2 (constant variance - doesn't depend on t)
  # cov(y_t, y_{t-s}) = cov(y_t, y_{t+s}) = \gamma_s (covariance between two values depends only on the length of time, s, separating them - doesn't depend on t itself)

# inflation and GDP, both their levels and their changes display characteristics of non-stationarity.
plot(gdp, main = "Real gross domestic product (GDP)")
plot(D1_gdp, main="Change in GDP")

plot(inf, main="Inflation rate")
plot(D1_inf, main="Change in the inflation rate")

# federal funds rate and the bond rate display characteristics of non-stationarity, their changes display characteristics of stationarity
plot(f, main="Federal funds rate")
plot(D1_f, main="Change in the federal funds rate")

plot(b, main= "Three-year bond rate")
plot(D1_b, main="Change in the bond rate")


#################################################
## 12.1.1 the first-order autoregressive model
#################################################

# the AR(1) model, is a useful univariate time series model for explaining the difference between stationary and non-stationary series: y_t = \rho y_{t-1} + v_t with |\rho| < 1 -> |\rho| < 1 implies that y_t is stationary. The AR(1) process shows that each realization of the random variable y_t contains a proportion \rho of last period's value y_{t-1} plus an error v_t drawn from a distribution with mean zero and variance \sigma^2_v.

# Setting parameters to create artificial data and generate the AR models below
N <- 500
a <- 1
l <- 0.01
rho <- 0.7
set.seed(1)
v <- ts(rnorm(N,0,1))

# AR(1) model with zero mean
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- rho*y[t-1]+v[t] # (12.2a) y_t = \rho y_{t-1} + v_t with |\rho|<1 -> E(y_t)=0
}
plot(y,type='l', ylab="rho*y[t-1]+v[t]",  ylim = c(-6, 6))
abline(h= mean(rho*y+v, na.rm=TRUE), lty=2, col="red")


# AR(1) model with non-zero mean
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+rho*y[t-1]+v[t] # (12.2b) y_t = \alpha + \rho y_{t-1} + v_t with |\rho|<1 -> E(y_t) = \alpha/(1-\rho)
}
plot(y,type='l', ylab="a+rho*y[t-1]+v[t]",  ylim = c(-2, 10))
abline(h= mean(a+rho*y+v, na.rm=TRUE), lty=2, col="blue") 


# AR(1) model fluctuating around a linear trend
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+l*time(y)[t]+rho*y[t-1]+v[t] # (12.2c) y_t = \alpha + \rho y_{t-1} + ∆ t + v_t -> E(y_t) = \mu + ∆*t (depends on t)
}
plot(y,type='l', ylab="a+l*time(y)[t]+rho*y[t-1]+v[t]")


#################################################
## 12.1.2 random walk models 
#################################################

# Random walk
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- y[t-1]+v[t] # (12.3a) y_t = y_{t-1} + v_t (i.e. \rho=1) -> E(y_t)=y_0 **BUT** var(y_t) = t*\sigma^2_v (depends on t)
}
plot(y,type='l', ylab="y[t-1]+v[t]")
abline(h=mean(y+v), lty=2)

# Notice that recursively we can depart from y_1 = y_0 + v_1 and rewrite y_t as  y_t = y_{t-1} + v_t = t*\alpha + y_0 + \sum{v_s}. This latter component, \sum{v_s}, is often called the stochastic trend because a stochastic component v_t is added for each time t, and because it causes the time series to trend in unpredictable directions. If the variable y_t is subjected to a sequence of positive shocks, v_t > 0, followed by a sequence of negative shocks, v_t < 0, it will have the appearance of wandering upward, then downward. PoE p.481


# Random walk with drift -> stochastic trend + deterministic trend
a <- 0.1 # each realization of y_t now contains an intercept (i.e., the drift component) \alpha
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+y[t-1]+v[t] # (12.3b)  y_t = \alpha + y_{t-1} + v_t -> E(y_t)= t*\alpha + y_0 **AND** var(y_t) = t*\sigma^2_v (both depend on t)
}
plot(y,type='l', ylab="a+y[t-1]+v[t]") # "Notice how the time-series data appear to be 'wandering' as well as 'trending' upward. In general, random walk with drift models show definite trends either upward (when the drift a is positive) or downward (when the drift a is negative)." PoE p. 481


# Random walk with trend - the additional term has the effect of strengthening the trend behavior.
y <- ts(rep(0,N))
for (t in 2:N){
  y[t]<- a+l*time(y)[t]+y[t-1]+v[t]  # (12.3c) y_t = \alpha + ∆*t + y_{t-1} + v_t
}
plot(y,type='l', ylab="a+l*time(y)[t]+y[t-1]+v[t]")


#######################################################################################################
# 12.2 SPURIOUS REGRESSIONS 
#######################################################################################################

## "The main reason why it is important to know whether a time series is stationary or non-stationary before one embarks on a regression analysis is that there is a danger of obtaining apparently significant regression results from unrelated data when non-stationary series are used in regression analysis." (PoE p. 482)

# let's generate two series independently
T <- 1000
set.seed(1) # graphs and point estimates won't match those from the book because I don't know what seed they used, but you get the idea.
y <- ts(rep(0,T))
vy <- ts(rnorm(T))
for (t in 2:T){
  y[t] <- y[t-1]+vy[t] # random walk rw1
}

set.seed(3)
x <- ts(rep(0,T))
vx <- ts(rnorm(T))
for (t in 2:T){
  x[t] <- x[t-1]+vx[t]  # random walk rw2
}

y <- ts(y[300:1000])
x <- ts(x[300:1000])
ts.plot(y,x, col=c("red", "blue"), ylab="y and x")

# Moreover, if one estimates a simple regression of series y on series x, one will obtain a highly significant coefficient and may be prone to interpret it as the effect of one variable on the other:
spurious.ols <- lm(y~x)
summary(spurious.ols) 
plot(x, y, type="p", col="grey")
abline(spurious.ols, lty=2) # "The apparent significance of the relationship is false. It results from the fact that we have related one series with a stochastic trend to another series with another stochastic trend." (PoE p. 438)

if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}

bgtest(spurious.ols, order=1, type="Chisq", fill=0) # we reject H0 and conclude that there is serial correlation -> Typically the residuals from such regressions will be HIGHLY correlated. For this example, the LM test value to test for first-order autocorrelation (p-value in parenthesis) is 669.25 (0.000); a sign that there is a problem with the regression.


#######################################################################################################
# 12.3 UNIT ROOT TESTS FOR STATIONARITY
#######################################################################################################

# There are many tests for determining whether a series is stationary or non-stationary. The most popular one, and the one that we discuss, is the Dickey–Fuller test. As noted in our discussion of the autoregressive and random walk models, stochastic processes can include or exclude a constant term and can include or exclude a time trend. There are three variations of the Dickey–Fuller test designed to take account of the role of the constant term and the trend:

## 12.3.1 DICKEY–FULLER TEST 1 (NO CONSTANT AND NO TREND): -     ∆ y_t = \gamma * y_{t-1}  + v_t
## 12.3.2 DICKEY–FULLER TEST 2 (WITH CONSTANT BUT NO TREND) -    ∆ y_t = \alpha + \gamma * y_{t-1}  + v_t
## 12.3.3 DICKEY–FULLER TEST 3 (WITH CONSTANT AND WITH TREND) -  ∆ y_t = \alpha + \gamma * y_{t-1} + \lambda * t + v_t

# H0: \gamma=0 (non-stationarity) vs H1: \gamma<0 (stationarity) - "To test the hypothesis in all three cases, we simply estimate the test equation by least squares and examine the t-statistic for the hypothesis that \gamma = 0. Unfortunately this t-statistic no longer has the t-distribution that we have used previously to test zero null hypotheses for regression coefficients. A problem arises because when the null hypothesis is true, y_t is nonstationary and has a variance that increases as the sample size increases. This increasing variance alters the distribution of the usual t-statistic when H0 is true. To recognize this fact, the statistic is often called a \tau (tau) statistic, and its value must be compared to specially generated critical values. NOTE THAT CRITICAL VALUES ARE GENERATED FOR THE THREE DIFFERENT TESTS because, as we have seen in Section 12.1, the addition of the constant term and the time-trend term changes the behavior of the time series."  (PoE p 485)

if (!require(urca)) {
  install.packages("urca")
  library(urca)
} # Augmented Dickey-Fuller unit root test: function ur.df() in package urca. Failing to reject H0 means non-stationarity.

# let's review the graphs to choose the proper specification for the Dickey-Fuller test
plot(f, main="Federal funds rate") # drift
plot(b, main= "Three-year bond rate") # drift

# Dickey-Fuller -> only one lag, as in PoE section 12.3.6
summary(ur.df(f, type ="drift")) # \tau = -2.5048 > -2.88 -> can't reject H0 -> non-stationarity
summary(ur.df(b, type ="drift")) # \tau = -2.7028 > -2.88 -> can't reject H0 at 5% level -> non-stationarity at 5% level

# AUGMENTED Dickey-Fuller -> as in PoE equation (12.6), but up to 10 lags, to be decided based on BIC
summary(ur.df(f, type ="drift", lags=10, selectlags = "BIC")) # \tau = -1.8868 > -2.88 -> can't reject H0 -> non-stationarity
summary(ur.df(b, type ="drift", lags=10, selectlags = "BIC")) # \tau = -1.0843 > -2.88 -> can't reject H0 -> non-stationarity



#################################################
## 12.3.7 ORDER OF INTEGRATION 
#################################################

# The order of integration of a series is the minimum number of times it must be differenced to make it stationary. If y_t is non-stationary but the transformed series ∆y_t= y_t - y_{t-1} is stationary, we say y_t is "integrated of order 1" or "I(1)" - because we just needed to take 1 difference in order to reach a stationary process. On the other hand, an already stationary series is said "integrated of order 0" or "I(0)". 

# is the first difference ∆F_t = F_t - F_{t-1} stationary? Is the first difference ∆B_t = B_t - B_{t-1} stationary?

plot(D1_f, main="Change in the federal funds rate")
summary(ur.df(D1_f, type ="none")) # \tau = -4.9341 < -1.95 -> reject H0 -> stationarity  
summary(ur.df(D1_f, lags = 10, type ="none", selectlags = "BIC")) # \tau = -3.7749  < -1.95 -> reject H0 -> stationarity  

plot(D1_b, main="Change in the federal funds rate")
summary(ur.df(D1_b, type ="none")) # \tau = -6.4406 < -1.95 -> reject H0 -> stationarity  
summary(ur.df(D1_b, lags = 10, type ="none", selectlags = "BIC")) # \tau = -6.4342 < -1.95 -> reject H0 -> stationarity  

# f is I(1) and b is I(1)


#######################################################################################################
# 12.4 COINTEGRATION
#######################################################################################################

# "If y_t and x_t are nonstationary I(1) variables, then we expect their difference, or any linear combination of them, such as e_t = y_t - \beta_1 - \beta_2 x_t, to be I(1) as well. However, there is an important case when e_t = y_t - \beta_1 - \beta_2 x_t is a stationary I(0) process. In this case y_t and x_t are said to be cointegrated. Cointegration implies that y_t and x_t share similar stochastic trends, and, since the difference e_t is stationary, they never diverge too far from each other. A natural way to test whether y_t and x_t are cointegrated is to test whether the errors e_t = y_t - \beta_1 - \beta_2 x_t are stationary. Since we cannot observe e_t, we test the stationarity of the least squares residuals, \hat{e}_t = y_t - \beta_1 - \beta_2 x_t using a Dickey–Fuller test. The test for cointegration is effectively a test of the stationarity of the residuals. If the residuals are stationary, then y_t and x_t are said to be cointegrated; if the residuals are nonstationary, then y_t and x_t are not cointegrated, and any apparent regression relationship between them is said to be spurious." PoE p. 488-489


#################################################
## 12.4.1 AN EXAMPLE OF A COINTEGRATION TEST
#################################################

# So far we know that f and b are I(1). Let's see whether e_t = b_t - \beta_1 - \beta_2 f_t is I(0)

if (!require(dynlm)) {
  install.packages("dynlm")
  library(dynlm)
}

fb.dyn <- dynlm(b ~ f)
summary(fb.dyn)  # (12.9) PoE
ehat.fb <- residuals(fb.dyn)

summary(dynlm(d(ehat.fb) ~ L(ehat.fb) + L(d(ehat.fb)))) # (Intercept) isn't significant, so we remove it in the regression below

coint.fb <- dynlm(d(ehat.fb) ~ L(ehat.fb) + L(d(ehat.fb)) - 1)
summary(coint.fb) # \tau = -4.196 < -3.37 (from table 12.4 PoE) -> reject H0 -> residuals are stationary -> f and b are cointegrated

# NOTICE THAT THE CRITICAL VALUES TO PERFORM THE DICKEY-FULLER TEST ON RESIDUALS ARE DIFFERENT FROM THOSE IN SECTION 12.3 UNIT ROOT TESTS FOR STATIONARITY. THUS WE WON'T USE THE ur.df() FUNCTION FROM THE urca PACKAGE TO TEST FOR COINTEGRATION. INSTEAD WE WILL PERFORM THE PHILLIPS-OULIARIS TEST, WHICH DOES THE STEPS ABOVE AUTOMATICALLY.

if (!require(tseries)) {
  install.packages("tseries")
  library(tseries)
} # Phillips-Ouliaris test for the null hypothesis that series are not cointegrated: function po.test() in package tseries

# H0 :the series are not cointegrated, residuals are nonstationary
# H1 :the series are cointegrated, residuals are stationary

bfx <- cbind(b,f) # bind both series you are testing for cointegration in a single object
po.test(bfx) # Phillips-Ouliaris test p-value = 0.04986 < 0.05 -> reject H0 -> f and b are cointegrated at 5% signif. level

# In other words, there is a fundamental relationship between these two variables (the estimated regression relationship between them is valid and not spurious).

#################################################
## 12.4.2 THE ERROR CORRECTION MODEL
#################################################
# "A relationship between I(1) variables is also often referred to as a long run relationship while a relationship between I(0) variables is often referred to as a short run relationship...if y and x are cointegrated, it means that there is a long-run relationship between them.The error correction model is a very popular model because it allows for the existence of an underlying or fundamental link between variables (the long-run relationship) as well as for short-run adjustments (i.e. changes) between variables, including adjustments to achieve the cointegrating relationship. It also shows that we can work with I(1) variables (y_{t-1}, x_{t-1}) and I(0) variables (∆y_t, ∆x_t) in the same equation provided that (y, x) are cointegrated, meaning that the term (y_{t-1} - \beta_0 - \beta_1 x_{t-1}) contains stationary residuals." (PoE p. 490-491)

# f and b are I(1)
# e = b - \alpha - \beta*f is I(0)
# b and f are cointegrated -> there is a long run relationship

ecm <- dynlm(d(b) ~ L(ehat.fb) + d(f) + L(d(f)) - 1)
summary(ecm)
 

#######################################################################################################
# If you have questions or suggestions regarding this script, write me: robson.tigre@unibo.it 
# For an alternative replication structure of Chapter 12, check the companion: 
# https://bookdown.org/ccolonescu/RPoE4/time-series-non-stationarity.html (Accessed on Nov 30 2017)
#######################################################################################################