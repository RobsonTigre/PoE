############################################################################
## PoE - Chapter 10 - Random Regressors and Moment-Based Estimation
## Author: Robson Tigre (T.A)
## Created: Nov 15 2017
## Last modified: Nov 21 2017
## R version 3.4.0 (2017-04-21) and RStudio 1.0.143 on Windows 10
############################################################################

#Before sending an email or posting a question to Piazza regarding the content of this script, please:
  ## Try to run line by line to get what is going on in a context; 
  ## Pay attention to the comments. Answers to your question may be a couple of lines below or above your current line;
  ## Check the help documentation for the function you are trying to understand.
#If the question persists after you've followed those steps, please don't hesitate to discuss it.

# Main sections in this script are:
	#10.2 Cases in Which x and e Are Correlated
	#10.3 Estimators Based on the Method of Moments
  #10.4 Specification Tests

# Loading the data
if (!require(devtools)) {         # if R can't load and attach the devtools package (i.e. !require)...
  install.packages("devtools")      # ...then R must install this package...
  library(devtools)                   # ...and then try to load and attach it again. Same for other packages below.
}
if (!require(PoEdata)) {
  devtools::install_git("https://github.com/ccolonescu/PoEdata")
  library(PoEdata)
}


data("mroz") # data on married women who are in the labor force
data.frame(attr(mroz, "names"), attr(mroz, "var.labels")) # variables and labels -> lfp: dummy variable = 1 if woman worked in 1975, else 0
mroz.1975 <- mroz[mroz$lfp==1,] #sub-setting the data on women who worked in 1975
str(mroz.1975) # 428 obs. of  25 variables -> 428 women in the original data set worked in 1975
attach(mroz.1975)

################################################################
#10.2 Cases in Which x and e Are Correlated
################################################################
# Useful statistical terminology:
  # Population: The population is the set of entities under study. For example, the height of men. This is a hypothetical population because it includes all men in the world!

  # Sample: Even if it is possible to measure the height of the entire male population worldwide it would be very costly and would take a great deal of time. Instead, we take a subset of this population, called a sample, and use this sample to draw inferences about the whole male population.

  # Parameter: A parameter is an unknown value that represents a characteristic of the population. In the example above, the parameter is the mean height of men. It is unknown because our sample doesn't include every man in the population. Since it is unknown, we are left with the task of estimating the true mean height of men (of the whole male population) based on the mean height of men in our sample.

  # Estimator: An estimator is the way or formula we use to estimate the unknown parameter. For instance, based on the observations in our sample, we estimate the mean height of men by summing the height of every man in the sample and dividing it by the sample size N. i.e., x.bar = (x_1 + x_2 + x_3 + ... + x_N)/N

  # An estimate: Is the specific value we obtain once we apply the estimator x.bar to our sample. While the estimator is the formula, the estimate would be a value of the mean height of men, for instance, 1.763 m. 

  # Unbiasedness (finite-sample property): An estimator is said to be unbiased when the expected value of this estimator is equal to the value of the true parameter. 

  # Consistency (asymptotic property): An estimator is said to be consistent if as the sample size increases (goes to infinity) the expected value of the estimator goes to the value of the true parameter (i.e. converges to the true parameter).

  # Efficiency:  Among all unbiased estimators in a given category, an estimator is said to be efficient if its variance is smaller than the variance of any other estimator in that category. 

# Now we can move on to the next review:

  # First let's recall that for the OLS estimator to be unbiased we assumed the correctly specified model y = \beta_1 + \beta_2 x + e with the expected value of the error term e conditional on any value of x being zero, i.e. E(e|x)=0. "This assumption implies that we have:   (i) omitted no important variables, (ii) that we have used the correct functional form, and (iii) that there exist no factors that cause the error term e to be correlated with x." (PoE p. 402)

  # If E(e|x) = 0, then we can show that it is also true that x and e are uncorrelated, i.e. that cov(x, e) = 0. Explanatory variables that are not correlated with the error term are called EXOGENOUS VARIABLES. Conversely, if x and e are correlated, then cov(x,e) != 0 and we can show that E(e|x)!= 0. Explanatory variables that are correlated with the error term are called ENDOGENOUS VARIABLES. (PoE p. 402)

  # "Throughout this chapter we use the relation between wages and years of education as an example. In this case the omitted variable "intelligence" is in the regression error, and it is likely to be positively correlated with the years of education a person receives, with more intelligent individuals choosing to obtain more years of education. When regressing wage on years of education, increases in wages are all attributed to increases in education by the least squares estimator. The effect of education is overstated because some of the increase in wages is also due to higher intelligence. The statistical consequences of correlation between x and e is that the least squares estimator is biased — and this bias will not disappear no matter how large the sample. Consequently the least squares estimator is inconsistent when there is correlation between x and e." (PoE, p. 405)

  # The algebra: log(wage) = \beta_1 + \beta_2 EDUC + \beta_3 EXPER + \beta_4 EXPER^2 + e, but we omitted intelligence, which is correlated with higher wages and is also correlated with higher levels of education. Therefore, e can be seen as e = \alpha_1 + \alpha_2 INTEL + u. And since cov(EDUC, INTEL)!=0, because "more intelligent individuals choose to obtain more years of education", then cov(EDUC, e)!=0, meaning EDUC is an endogenous variable. 

##10.2.4 LEAST SQUARES ESTIMATION OF A WAGE EQUATION
wage.ols <- lm(log(wage) ~ educ + exper + I(exper^2))
summary(wage.ols) #equation (10.6) -> 1 additional year of educ "causes" a 10.74896% increase in wage


################################################################
#10.3 Estimators Based on the Method of Moments
################################################################

if (!require(car)) {
  install.packages("car")
  library(car)
} # car (Companion to Applied Regression) package to perform the F test on the "excluded instruments" using linearHypothesis()

if (!require(AER)) {
  install.packages("AER")
  library(AER)
} # AER (Applied Econometrics with R) package to perform IV 2SLS the using ivreg()

?ivreg # Please take a look. In RStudio you may click the icon on the right of the printer icon to enlarge the documentation

## 10.3.6 INSTRUMENTAL VARIABLES (IV) ESTIMATION OF THE WAGE EQUATION
### Conditions for a good instrumental variable z: (1) z does not have a direct effect on y, and thus it does not belong on the right-hand side of the model as an explanatory variable [this is called redundancy or exclusion restriction]. (2) z is not correlated with the regression error term e. It is exogenous. (3) z is strongly correlated with x, the endogenous explanatory variable. (PoE p. 410)

### In the original model y = \beta_0 + \beta_1 x + e, where cov(x,e)!=0, The estimation using the instrumental variable z follows:
  
  # (1) In what we call "first stage regression", use z to obtain predicted values of the endogenous variable x, i.e. \hat{x}. The first stage regression has the following form: x = \theta_0 + \theta_1 z + u. The estimate \hat{\theta}_1 in this model should be significant, indicating that z is strongly correlated with x, the endogenous explanatory variable.
  
  # (2) the "second stage regression", means estimating the original model but replacing x by \hat{x}. This means that instead of estimating y = \beta_0 + \beta_1 x + e, we will estimate y = \beta_0 + \beta_1 \hat{x} + e, and voilà, we are done! The \hat{\beta}_1 obtained in the second stage regression is your IV estimate of \beta_1.


### Using one instrumental variable - mothereduc: mother's education level
first.stage1 <- lm(educ ~ exper + I(exper^2) + mothereduc) #This line is just informative, ivreg() computes the first stage automatically 
summary(first.stage1) #equation (10.26) - first stage with one IV - note that the coefficient of MOTHEREDUC is very significant. This is important, as it indicates that our instrument is correlated with the variable we suspect to be endogenous, even after accounting for the other exogenous variables (i.e, considering exper and exper^2 are exogenous) in the model. 

#### When instrumental variables are weak, estimates and tests based on the resulting IV estimator are unreliable, so we test...
linearHypothesis(first.stage1, "mothereduc=0") # F=73.946, so our instrument isn't weak. This conclusion comes from the rule-of-thumb that F on excluded instruments should be greater than 10. This rule has been refined by econometric researchers Stock and Yogo, and their analysis is discussed in Appendix 10E. 

second.stage1 <- ivreg(log(wage) ~ educ + exper + I(exper^2) | exper + I(exper^2) + mothereduc)
summary(second.stage1) # were the educ coefficient significant, we would say that 1 additional year of educ causes a 4.92630% increase in wage. This reduction from 10.74896% (OLS) to 4.92630% (IV1) "is consistent with the fact that the least squares estimator tends to overestimate the effect of education if EDUC is positively correlated with the omitted factors in the error term. Also notice that the standard error on the coefficient of education (0.0374) is over 2.5 times larger than the standard error reported with the least squares estimates (0.0141). This reflects the fact that even with a good instrumental variable, THE IV/2SLS ESTIMATOR IS NOT EFFICIENT, as discussed in Section 10.3.3a. How can we improve the efficiency of the instrumental variables estimator? We can obtain a larger sample, if possible, or we can obtain more and stronger instrumental variables." (PoE p. 415-416)

### Using two instrumental variables - mothereduc and fathereduc: both parents' education level
first.stage2 <- lm(educ ~ exper + I(exper^2) + mothereduc + fathereduc) # Informative line, ivreg() computes the first stage automatically 
summary(first.stage2) #In the first stage regression, the estimated coefficients of MOTHEREDUC and FATHEREDUC are highly significant. This is important, as it indicates that INDIVIDUALLY our instruments are correlated with the variable we suspect to be endogenous, even after accounting for the other exogenous variables (i.e, considering exper and exper^2 are exogenous) in the model. 

linearHypothesis(first.stage2, c("mothereduc=0", "fathereduc=0")) # F=55.4 > 10 "The F-statistic value for the null hypothesis that both these coefficients are zero is 55.40. This value is greater than the rule-of-thumb threshold of 10." It indicates that JOINTLY our instruments are correlated with the variable we suspect to be endogenous.

second.stage2 <- ivreg(log(wage) ~ educ + exper + I(exper^2) |  exper + I(exper^2) + mothereduc + fathereduc)
summary(second.stage2) # The educ coefficient is now significant at 10% level: 1 additional year of educ causes a 6.13966% increase in wage. Compared to the previous result using only MOTHEREDUC as an instrument, we see that there is an increase in the estimate of the return to education to 6.14%, and a slight REDUCTION IN THE STANDARD ERROR.

#In sum, if an explanatory variable is correlated with the regression error term, the least squares estimator fails. If a strong instrumental variable is available, the IV estimator is consistent and approximately normally distributed in large samples. But if we use a weak instrument, or an instrument that is invalid in the sense that it is not uncorrelated with the regression error, then IV estimation can be as bad as, or worse than, using the least squares estimator. 

## The mechanics behind the 2 stages 
#------------------------------------------------------------------------------------------------------------------------------------
# ATTENTION: this part is just so you can understand each step of the 2SLS machinery. Although coefficients obtained this way are the
# same as those obtained with ivreg(), STANDARD ERRORS FOLLOWING THE MANUAL ESTIMATION BELOW ARE BIASED, SO THIS IS NOT THE WAY TO GO 
# WHEN PERFORMING 2SLS ON THE EXERCISES, ON THE 2ND EXAM OR IN YOUR PRIVATE LIFE! "While this two-step process yields proper IV/2SLS
# estimates, as we discussed in section 10.3.4, the accompanying standard errors and t-values are not correct." (POE P. 415) #------------------------------------------------------------------------------------------------------------------------------------


### In the original model y = \beta_0 + \beta_1 x + e, where cov(x,e)!=0, The estimation using the instrumental variable z follows:

 # "(1) In what we call 'first stage regression', use z to obtain predicted values of the endogenous variable x, i.e. \hat{x}. The first stage regression has the following form: x = \theta_0 + \theta_1 z + u. The estimate \hat{\theta}_1 in this model should be significant, indicating that z is strongly correlated with x, the endogenous explanatory variable."
  
  #first you obtain fitted (i.e. "predicted") values using the coefficients estimated on the first stage...
  first.stage2 <- lm(educ ~ exper + I(exper^2) + mothereduc + fathereduc); summary(first.stage2)
  educ.hat <- fitted(first.stage2)

 # (2) "the 'second stage regression', means estimating the original model but replacing x by \hat{x}. This means that instead of estimating y = \beta_0 + \beta_1 x + e, we will estimate y = \beta_0 + \beta_1 \hat{x} + e..."
  
  #...then estimate the original model plugging the predicted values educ.hat where you would otherwise have the endogenous variable educ.
  second.stage2.MANUAL <- lm(log(wage) ~ educ.hat + exper + I(exper^2))
  summary(second.stage2.MANUAL) ##... and voilà: \hat{\beta}_1 obtained in the second stage regression is your IV estimate of \beta_1.
  summary(second.stage2) # second.stage2.MANUAL and second.stage2 have the same "Estimate", BUT NOT "Std. Error", "t value" and "Pr(>|t|)"


################################################################
#10.4 Specification Tests
################################################################
# In this section we ask two other important questions that must be answered in each situation in which instrumental variables estimation is considered: (1) Can we test for whether x is correlated with the error term? This might give us a guide for when to use least squares and when to use IV estimators; (2) Can we test if our instrument is valid, and uncorrelated with the regression error, as required?

# Question (1) is important because if cov(x,e)!=0 OLS estimator will be biased and inconsistent, meaning the bias wont get any smaller as the sample size increases, while the IV estimator will be consistent. HOWEVER, if cov(x,e)=0 you don't want to use the IV estimator, because in this scenario the OLS estimator will be unbiased and have lower variance than the IV estimator (remember that both OLS and IV will be consistent under cov(x,e)=0).

## 10.4.1 THE HAUSMAN TEST FOR ENDOGENEITY - H0: cov(x,e)=0 HA: cov(x,e)!=0
###The idea of the test is to compare the performance of the least squares estimator to an instrumental variables estimator

#Manual computation - try understand the logic behind the manual computation, later we will see the automatic way
first.stage2 <- lm(educ ~ exper + I(exper^2) + mothereduc + fathereduc) # estimate the first stage model as usual
v.hat <- residuals(first.stage2) # obtain the residuals of the first stage model above
artif1 <- lm(log(wage) ~ educ + exper + I(exper^2) + v.hat) # include the residuals from the first stage as an explanatory variable in the original regression. let's call this "artificial regression". 
summary(artif1) # Employ the usual t-test for the hypothesis of significance between x and e -> p-value = 0.095441, meaning we reject  H0: cov(x,e)=0 at 10% significance level. This is evidence that educ is endogenous.

#The automatic computation: set the parameter diagnostics to TRUE in summary() and look at Wu-Hausman 
summary(second.stage2, diagnostics=TRUE) # p-value = 0.0954 < 0.1 (just as above) -> reject H0 at 10% -> evidence of endogeneity 


## 10.4.2 TESTING INSTRUMENT VALIDITY - H0: all the surplus moment conditions are valid
## In order to compute the IV estimator for an equation with B possibly endogenous variables, we must have at least B instruments [when we have B endogenous variables and B ivs, we refer to that model as just identified or exactly identified]. The validity of this minimum number of required instruments cannot be tested. Nevertheless, in the case in which we have L > B instruments available [when we have B endogenous variables and L ivs, with L > B, we refer to that model as overidentified], we can test the validity of the L - B extra, or surplus, moment conditions.

#Manual computation - later we will see the automatic way
second.stage2 <- ivreg(log(wage)~educ+exper+I(exper^2) | exper+I(exper^2)+mothereduc+fathereduc) # compute the iv estimates: \hat{\beta}_k
e.hat <- residuals(second.stage2) # obtain the residuals of the second stage
artif2 <- lm(e.hat ~ mothereduc + fathereduc) # Regress e.hat only on the L available instruments used in second.stage2.
R2 <- summary(artif2)$r.squared # get the R2 of the artificial regression above
N <- nrow(mroz.1975) # get the sample size
LM.test  <- N * R2 # looks like the Lagrange multiplier tests from chapter 8, right?!
p.value = 1-pchisq(LM.test, df = 1) # degrees of freedom equals L-B: L=2 (2 ivs), B=1 (1 endog.) -> L-B=1
c(LM.test, p.value) # p-value = 0.5386398 > 0.1 -> we can't reject H0 that all the surplus moment conditions are valid

#The automatic computation: set the parameter diagnostics to TRUE in summary() 
summary(second.stage2, diagnostics=TRUE)  #  p-value=0.5386 -> same conclusion as from the manual way

#CONCLUSION: taken together, result above suggest that our IV estimator for the wage equation is consistent. (PoE p 423)


#######################################################################################################
# If you have questions or suggestions regarding this script, write me: robson.tigre0@gmail.com 
# For an alternative replication structure of Chapter 10, check the companion: 
# https://bookdown.org/ccolonescu/RPoE4/random-regressors.html (Accessed on Nov 15 2017)
#######################################################################################################
