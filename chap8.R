################################################################
## PoE - Chapter 8 - Heteroskedasticity
## Authors: Prof. Sergio Pastorello/T.A Robson Tigre
## Created: Oct 29 2017
## Last modified: Nov 7 2017 -> MINOR COMMENT EDITS ON Nov 8
## R version 3.4.0 (2017-04-21) and RStudio 1.0.143 on Windows 10
################################################################

# Below I provide comments on R code chunks to cover the main topics seen during the lecture. 
# Try to run line by line to get what is happening. Main sections below are:
# 8.1 The nature of heteroskedasticity
# 8.2 Detecting heteroskedasticity
# 8.3 HC robust standard errors

rm(list = ls()) #this frees the "Environment" panel from other objects by removing them.

# Loading the data
if (!require(devtools)) {         # if R can't load and attach the devtools package (i.e. !require)...
  install.packages("devtools")      # ...then R must install this package...
  library(devtools)                   # ...and then try to load and attach it again. Same for other packages below.
}
if (!require(PoEdata)) {
  devtools::install_git("https://github.com/ccolonescu/PoEdata")
  library(PoEdata)
}



data("food") #loads food data from the package PoEdata
attach(food) #attaching the data saves us from having to refer to variables as food$food_exp and food$income

################################################################
# 8.1 The nature of heteroskedasticity
################################################################

## Estimate the regression model so we can later plot the regression line on the income vs food_exp plot 
lm.food = lm(food_exp ~ income); summary(lm.food)
e.hat = lm.food$residuals #alternatively you can obtain the vector of residuals through residuals(lm.food)

## Plotting food_exp vs. income
plot(income, food_exp,
     col = "dimgrey", 
     ylim = c(0, max(food_exp)),
     xlim = c(0, max(income)),
     xlab = "weekly income in $100",
     ylab = "weekly food expenditure in $",
     main = "food data",
     type = "p")
abline(lm.food, col = "blue") #let's just include the regression line obtained from lm.food in our plot
pre <- predict(lm.food)
segments(income,food_exp,income, pre, col="grey") #the greater the income the farther the data points get from the blue regression line: "when the variances for all observations are not the same, we say that heteroskedasticity exists" (PoE, p. 299)


## Now let's check it from the perspective of residuals as the y-axis
# Plotting residuals (i.e., e_i = y_i - E[y_i|x_i] = y_i - beta_1 - beta_2 x_i) vs. income
plot(income, e.hat,
     col = "dimgrey", 
     ylim = c(min(e.hat), max(e.hat)),
     xlim = c(0, max(income)),
     xlab = "weekly income in $100",
     ylab = "LS residuals",
     main = "food data",
     type = "p")
abline(h = 0, lty = 2, lwd = 2, col = "red")
segments(income,e.hat,income, 0, col="grey") 

################################################################
# 8.2 Detecting Heteroskedasticity
################################################################

# 8.2.2 Lagrange multiplier tests -> H0: Homoskedasticity, so rejecting H0 means evidence of heteroskedasticity
## Manual computation - try understand the logic behind the manual computation, later we will see the automatic way

### Auxiliary regression
e.hat2 = e.hat^2 # the auxiliary regression will use squared residuals as the dependent variable (PoE p. 304) 
z2 = income
z3 = log(income) 

lm.bp = lm(e.hat2 ~ z2 + z3) # attention: e.hat2 ~ z2 + z3 is just an example, and is not the "must-be" form of the auxiliary regression. So far we presupposes that we have knowledge of the variables appearing in the variance function if the alternative hypothesis of heteroskedasticity is true. In other words, so far we are assuming we are able to specify z_2,z_3,...,z_s, and in this case we are assuming z2 = income and z3 = log(income). Therefore it could have other forms according to different assumptions. 
summary(lm.bp)

# Breusch and Pagan test statistic
N = length(income)
R2 = summary(lm.bp)$r.squared #Since the R2 measures the proportion of variation in e.hat2 explained by the z's, it is a natural candidate to compose Lagrange multiplier test statistics
BP.test = N * R2 #it turns out the R2 multiplied by the sample size, N, follows the chi-squared distribution with degrees of freedom equals the number of coefficients in the auxiliary regression minus 1 (in this case, 3 - 1 )
p.value = 1-pchisq(BP.test, df = 2)
c(BP.test, p.value) #we reject H0 and conclude that the residuals variance depends on income.

#or just use the function bptest, from the package lmtest to do the robust Breusch-Pagan (results vary a little) 
## Using the 'bptest' command in the 'bstats' package
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}


?bptest
bptest(lm.food, ~ z2 + z3, studentize = FALSE)


# 8.2.2a The White test: "White suggested defining the z's in the auxiliary regression as equal to the x's in the main regression, the squares of the x's, and possibly their cross-products. Frequently, the variables that affect the variance are the same as those in the mean function. Also, by using a quadratic function we can approximate a number of other possible variance functions." (PoE, p 305)

## Manual computation
lm.white = lm(e.hat2 ~ income + I(income^2))
R2 = summary(lm.white)$r.squared
White.test = N * R2
p.value = 1-pchisq(White.test, df = 2)
c(White.test, p.value) #we reject H0 and conclude that the variance depends on income.

## Using the 'bptest' command 
bptest(lm.food, ~ income + I(income^2))


# 8.2.3 The Goldfeld-Quandt test - "It is designed for two groups of data with possibly different variances. To introduce this case, consider a wage equation where earnings per hour (WAGE) depends on years of education (EDUC), years of experience (EXPER) and a dummy variable METRO that is equal to one for workers who live in a metropolitan area and zero for workers who live outside a metropolitan area." (PoE p. 307)

## Loading the data
data("cps2")
attach(cps2)

## Manual computation
### Subsample sizes
N = dim(cps2)[[1]]
cps2.M = subset(cps2, metro == 1) #subset of individuals in the metropolitan area (808 obs)
N.M = dim(cps2.M)[[1]]
cps2.R = subset(cps2, metro == 0) #subset of individuals in the rural area (192 obs)
N.R = dim(cps2.R)[[1]]
c(N, N.M, N.R)

### Metropolitan area model
#notice that since we are running one separate regression for metro == 1 and another for metro == 0, it makes no sense to include the variable metro in the regressions (otherwise it would result in multicollinearity)
lm.cps2.M = lm(wage ~ educ + exper, data = cps2.M) 
sigma2.M = summary(lm.cps2.M)$sigma^2

### Rural area model
lm.cps2.R = lm(wage ~ educ + exper, data = cps2.R)
sigma2.R = summary(lm.cps2.R)$sigma^2

###This is the two tail application of the test, meaning we are not hypothesizing whether the variance in the metropolitan area is greater (right tail) or lesser (left tail) than that in the rural area. We just want to test whether they are equal or not or not -> H0: sigma2_metro = sigma2_rural vs HA: sigma2_metro != sigma2_rural
GQ.test = sigma2.M/sigma2.R #Test statistic under H0: s^2_M = s^2_R (see PoE p. 307 for further explanation)

### Critical values for a two-tail 5% test
alpha = 0.05
lower.Fc = qf(alpha/2, N.M-3, N.R-3); lower.Fc # a two-tail test with alpha=5% -> 2.5% on each tail (thus alpha/2)
upper.Fc = qf(1-alpha/2, N.M-3, N.R-3); upper.Fc
c(GQ.test, lower.Fc, upper.Fc) #since the critical values lower.Fc(0,025; 805; 189)= 0.805 and upper.Fc(0,975; 805; 189)=1.262, and 1.262 < 2.088, we reject H0

## Using the 'gqtest' command in the 'lmtest' package
?gqtest

### Estimating the full model (without METRO)
lm.cps2 = lm(wage ~ educ + exper, data = cps2)

### Computing the test statistic for a two-tail alternative
gqtest(lm.cps2, point = N.R, alternative = "two.sided",
       order.by = ~ metro)

# 8.2.3a The food expenditure example
N = dim(food)[[1]]
gqtest(lm.food, point = N/2, alternative = "two.sided",
       order.by = ~ income)

################################################################
# 8.3 Heteroskedasticity-Consistent Standard Errors
################################################################

## Loading the 'matrix 'sandwich' package 
## and computing the robust variance matrix
if(!require(sandwich)){
  install.packages("sandwich")
  library(sandwich)
}  #the sandwich package implements heteroskedasticity-consistent standard error estimators or simply 'robust' standard error estimators. It is called 'sandwich' because of its algebraic expression - (X'X)^{-1} X' Omega X (X'X)^{-1}, which to some resembles a sandwich (?!). The term 'robust' is used because they are valid in large samples for both heteroskedastic and homoskedastic errors. Given this characteristic it is a good practice in empirical research to always report robust standard errors

?vcovHC

var.HC = vcovHC(lm.food, type = "HC")

## Classical standard errors and t tests
coeftest(lm.food)

## HC robust standard errors and t tests
coeftest(lm.food, vcov = var.HC) #notice that point estimates in column 'Estimate' are the same between coeftest(lm.food) and coeftest(lm.food, vcov = var.HC). This is because coefficient estimates are not affected by heteroskedasticity. Standard errors however are. Since confidence intervals and tests are based on standard errors, not using HC robust standard errors may result in invalid inference about the significance of the estimated effect.

## Classical interval estimates
confint(lm.food)

## HC robust interval estimates
se.HC = sqrt(diag(var.HC))
alpha = 0.05 #95% confidence intervals (1-0.05)
tc = qt(c(alpha/2, 1-alpha/2), df = lm.food$df.residual); tc #those are the critical values
coef(lm.food)                       # Estimates
coef(lm.food)+tc[[1]]*se.HC         # Lower bound
coef(lm.food)+tc[[2]]*se.HC         # Upper bound

## Classical F test for H0: beta2 = 8
if(!require(car)){
  install.packages("car")
  library(car)
}
?linearHypothesis
linearHypothesis(lm.food, c("income = 8"))

## HC robust F test for H0: beta2 = 10
linearHypothesis(lm.food, c("income = 8"), vcov. = var.HC)

##############################################################################################
# If you have questions or suggestions regarding this script, write me: robson.tigre0@gmail.com 
# If you would like to the replicate Chapter 8 in its full extent, check the companion: 
# https://bookdown.org/ccolonescu/RPoE4/heteroskedasticity.html (Accessed on Nov 07 2017)
###############################################################################################
