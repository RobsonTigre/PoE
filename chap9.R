############################################################################
## PoE - Chapter 9 - Regression With Time-Series Data: Stationary Variables
## Author: Robson Tigre (T.A)
## Created: Nov 08 2017
## Last modified: Nov 16 2017 -> MINOR COMMENT EDITS ON Nov 18
## R version 3.4.0 (2017-04-21) and RStudio 1.0.143 on Windows 10
############################################################################

#Before sending an email or posting a question to Piazza regarding the content of this script, please:
	## Try to run line by line to get what is going on in a context; 
	## Pay attention to the comments. Answers to your question may be a couple of lines below or above your current line;
	## Check the help documentation for the function you are trying to understand.
#If the question persists after you've followed the steps above, please don't hesitate to discuss it. 

# Main sections in this script are:
	# 9.2 Finite Distributed Lags
	# 9.3 Serial Correlation
	# 9.4 Other Tests For Serially Correlated Errors
	# 9.5 Estimation With Serially Correlated Errors
	# 9.6 Autoregressive Distributed Lag Models
	# 9.7 Forecasting

# Loading the data
if (!require(devtools)) {         # if R can't load and attach the devtools package (i.e. !require)...
  install.packages("devtools")      # ...then R must install this package...
  library(devtools)                   # ...and then try to load and attach it again. Same for other packages below.
}
if (!require(PoEdata)) {
  devtools::install_git("https://github.com/ccolonescu/PoEdata")
  library(PoEdata)
}

#"An assumption that we maintain throughout this chapter is that the variables in our equations are stationary. This assumption will take on more meaning in Chapter 12 when it is relaxed. For the moment we note that a stationary variable is one that is not explosive, nor trending, and nor wandering aimlessly without returning to its mean." (PoE p. 339)

data("okun") #loads okun data from the package PoEdata
str(okun) # g -> "percentage change in U.S. Gross Domestic Product, seasonally adjusted"; u -> "U.S. Civilian Unemployment Rate  (Seasonally adjusted)". Nowhere here it is explicit that the okun data set is understood by R as a time-series data, right?
is.ts(okun) # we use is.ts() to check whether R understands okun as a time-series data, and we see it doesn't -> FALSE

#...this is because we must tell R that this data-set should be understood and processed as a time-series data, along with the time of the first observation (1985Q2 or second quarter of 1985) and the frequency of observations (an observation for every 3-month period, i.e. every quarter of a year).
okun.ts <- ts(okun, start=c(1985,2), end=c(2009,3), frequency=4)
is.ts(okun.ts) #TRUE

################################################################
#9.2 Finite Distributed Lags
################################################################

## 9.2.2 AN EXAMPLE: OKUN'S LAW 
###The economic model known as Okun's Law: "In this model the change in the unemployment rate from one period to the next depends on the rate of growth of output in the economy: U_t - U_{t-1} = -\gamma*(G_t - G_N)" (PoE p.343)

if (!require(dynlm)) {
	install.packages("dynlm")
	library(dynlm)
} # dynlm -> Dynamic Linear Models and Time Series Regression.

#Have in mind the assumption TSMR2 of the distributed lag model: y and x are stationary random variables, and e_t is independent of current, past and future values of x
okunL3.dyn <- dynlm(d(u) ~ g + L(g, 1:3), data=okun.ts); summary(okunL3.dyn)  # as in PoE Table 9.2 - this is the same as entering the formula as d(u) ~ L(g, 0:3) instead of d(u) ~ g + L(g, 1:3). L(g, 0:3) will include the "lag 0" of g_t, which is g_t, the lag 1 of g_t, which is g_{t-1} and so on. You may want to try and compare using the two formulas. 
okunL2.dyn <- dynlm(d(u) ~ L(g, 0:2), data=okun.ts); summary(okunL2.dyn) # as in PoE Table 9.2
c(AIC(okunL3.dyn), BIC(okunL3.dyn), AIC(okunL2.dyn), BIC(okunL2.dyn)) #-55.43179 > -58.95107 ; -40.10853 > -46.12932, meaning okunL2.dyn is preferable to okunL3.dyn acording to Akaike and Bayesian criteria

plot(okun.ts[,"g"],
	col = "dimgrey",
	xlab = "Year",
	ylab = "Growth (g)",
	main = "U.S. GDP growth: 1985Q2 to 2009Q3",
	type="o")

plot(okun.ts[,"u"], 
	col = "dimgrey",
	xlab = "Year",
	ylab = "Unemployment (u)",
	main = "U.S. unemployment rate: 1985Q3 to 2009Q3",
	type="o")  #ATTENTION: the model on the book involves u_t - u_{t-1}, while the okun data provides us with u_t. you may compare this graph with FIGURE (9.4) A and notice that they are not the same...

du <- diff(okun.ts[,"u"]) #use the function diff() to create u_t - u_{t-1}
plot(du,
	col = "dimgrey",
	ylim = c(-0.50, 1.25),
	xlab = "Year",
	ylab = "Change in unemployment (du)",
	main = "Change in the U.S. unemployment rate: 1985Q3 to 2009Q3",
	type="o")


################################################################
#9.3 Serial Correlation
################################################################
#Now recall we've saw that TSMR2, of the distributed lag model, assumes e_t is independent of current, past and future values of x. on the other hand, the existence of serial correlation implies violation of this assumption since current, past, and future values of some (or all) of the variables on the right hand side of the regression function (i.e., right-hand side variables) may be correlated with e_t.

##9.3.1 SERIAL CORRELATION IN OUTPUT GROWTH
ggL1 <- data.frame(cbind(okun.ts[,"g"], lag(okun.ts[,"g"],-1)))
names(ggL1) <- c("g","gL1")
plot(ggL1,
	col = "dimgrey",
	xlab = "Growth (g)",
	ylab = "1st Lag of Growth (g_{t-1})",
	main = expression("Scatter diagram for g and g"[t-1]))
meang <- mean(ggL1$g, na.rm=TRUE)
abline(v=meang, lty=2)
abline(h=mean(ggL1$gL1, na.rm=TRUE), lty=2)
reg1 <- lm(ggL1$gL1 ~ ggL1$g)
abline(reg1, col = "blue") 
summary(reg1) #there is a positive and significant correlation between growth in time t and growth of time t-1 

ggL2 <- data.frame(cbind(okun.ts[,"g"], lag(okun.ts[,"g"],-2)))
names(ggL2) <- c("g","gL2")
plot(ggL2,
	col = "dimgrey",
	xlab = "Growth (g)",
	ylab = "2nd Lag of Growth (g_{t-2})",
	main = expression("Scatter diagram for g and g"[t-2]))
meang <- mean(ggL2$g, na.rm=TRUE)
abline(v=meang, lty=2)
abline(h=mean(ggL2$gL2, na.rm=TRUE), lty=2)
reg2 <- lm(ggL2$gL2 ~ ggL2$g)
abline(reg2, col = "blue")
summary(reg2) #there is a positive and significant correlation between growth in time t and growth of time t-2 

#a more convenient way of checking the existence of this correlation is by using a correlogram: 
growth_rate <- okun.ts[,"g"]
acf.growth <- acf(growth_rate, lag.max=12, ylab = "Correlation", # Auto- and Cross- Covariance and -Correlation Function Estimation
						 main = "Correlogram for Growth (g)") #Notice that since the frequency our data is quarters, the x-axis is marked into quarters (0; 0.25; 0.5; 0.75), but have in mind that each bar in the correlogram corresponds to one lag, starting with lag 0. The blue dashed lines indicate confidence intervals for the existence of serial correlation. When bars cross these lines it is evidence of existence of serial correlation. 
acf.growth[c(0.25,0.50,0.75,1)] #to see the numerical values on the y-axis of the correlogram

###9.3.2a A Phillips Curve -> describes the relationship between inflation and unemployment. INF_t = INF^E_t - \gamma*(U_t - U_{t-1}), where INF^E_t  denotes inflationary expectations for period t
data("phillips_aus")
str(phillips_aus) # inf -> "Australian Inflation Rate" u -> "Australian Unemployment Rate (Seasonally adjusted)"
phill.ts <- ts(phillips_aus, start=c(1987,1), end=c(2009,3), frequency=4)
inflation <- phill.ts[,"inf"]

plot(inflation,
	col = "dimgrey",
	ylim = c(-0.5, 3),
	xlab = "Year",
	ylab = "Australian Inflation Rate (inf)",
	main = "Time series for Australian price inflation: 1987Q1 to 2009Q3",
	type="o")

Du <- diff(phill.ts[,"u"]) #ATTENTION the model on the book involves u_t - u_{t-1}
plot(Du,
	col = "dimgrey",
	ylim = c(-0.8, 1),
	xlab = "Year",
	ylab = "Australian Unemployment Rate (Du)",
	main = "Time series for the quarterly change in the Australian unemployment rate: 1987Q1 to 2009Q3",
	type="o")

phill.dyn <- dynlm(inf ~ diff(u), data=phill.ts); summary(phill.dyn) #notice that you may use diff() instead of d() to denote (u_t - u_{t-1}), but the converse is not true. d() and L() are only valid inside dynlm, while diff() and lag() are always valid.  
ehat <- residuals(phill.dyn) #now we want to test time/serial correlation in the residuals of the regression phill.dyn which is not the same thing as we did acf.growth, since g in okun's model is an explanatory variable.
acf.ehat <- acf(ehat, lag.max=12, ylim = c(-0.3, 1), ylab = "Correlation", main = "Correlogram for residuals from least-squares estimated Phillips curve") # each bar corresponding to one lag, starting with lag 0. 
acf.ehat[c(0.25,0.50,0.75,1)] #what are your conclusions from this graph?


################################################################
#9.4 Other Tests For Serially Correlated Errors
################################################################

##9.4.1 A LAGRANGE MULTIPLIER (LM) TEST
##as convenient as correlograms are, there are more formal ways to investigate the existence of serial correlation
## Manual computation - try understand the logic behind the manual computation to 
LM1 <- dynlm(inf ~ diff(u) + L(ehat), data=phill.ts); summary(LM1) # in y_t = \beta_1 + \beta_2 x_t + \rho e_{t-1} + v_t we find that \rho is significant, meaning there is a correlation between y_t and the residual from t-1, e_{t-1}. 
LM2 <- dynlm(ehat ~ diff(u) + L(ehat), data=phill.ts); summary(LM2) # e_t = \gamma_1 + \gamma_2 x_t + \rho e_{t-1} + v_t. look at the coefficient for L(ehat) in LM1 and LM2. Yes.
R2 <- summary(LM2)$r.squared
T <- length(LM2$residuals)
LM.test  <- T * R2 #looks like the lagrange multiplier tests from chapter 8, right?!
p.value = 1-pchisq(LM.test, df = 1) #degrees of freedom equals the number of lagged residuals included in LM2 (only e_{t-1} included, so df=1)
c(LM.test, p.value) # 2.760881e+01 = 27.60881. we reject H0 and conclude that there is serial correlation

#The automatic computation: use the function bgtest from the package lmtest to perform the Breusch-Godfrey Test 
if (!require(lmtest)) {
  install.packages("lmtest")
  library(lmtest)
}

bgtest(phill.dyn, order=1, type="Chisq", fill=NA)
bgtest(phill.dyn, order=1, type="Chisq", fill=0) 
bgtest(phill.dyn, order=4, type="Chisq", fill=NA) #9.4.1a Testing Correlation at Longer Lags - testing for four lags of the residuals
bgtest(phill.dyn, order=4, type="Chisq", fill=0) #9.4.1a Testing Correlation at Longer Lags - testing for four lags of the residuals


##9.4.2 THE DURBIN-WATSON TEST - "It is used less frequently today because its critical values are not available in all software packages, and one has to examine upper and lower critical bounds instead. Also, unlike the LM and correlogram tests, its distribution no longer holds when the equation contains a lagged dependent variable." (PoE p. 355)
dwtest(phill.dyn) #also rejects H0

################################################################
#9.5 Estimation With Serially Correlated Errors
################################################################
## In cases such y_t = \delta + \theta_1 y_{t-1} + \delta_0 x_t + \delta_1 x_{t-1} + v_t "the time-series assumption TSMR2 introduced in Section 9.2.1 is no longer valid. In the context of the above equation, this assumption says that v_t is not correlated with CURRENT, PAST AND FUTURE values of y_{t-1}; x_t and x_{t-1}. Since y_t is a future value of y_{t-1} and y_t depends directly on v_t, the assumption will be violated. Under this assumption, the least squares estimator is no longer unbiased, but it does have the desirable large sample property of consistency, and, if the errors are normally distributed, it is best in a large sample sense. Thus, we replace TSMR2 with the following assumption: Assumption for models with a lagged dependent variable TSMR2A - In the multiple regression model y_t = \beta_1 + \beta_2 x_{t2} + ... + \beta_K x_{tK} + v_t where some of the x_{tk} may be lagged values of y, v_t is uncorrelated with ALL x_{tk} AND THEIR PAST VALUES" (PoE p. 356)

if(!require(sandwich)){
  install.packages("sandwich")
  library(sandwich)
}  ##you surely remember this from chap8.R (or at least you should). the sandwich package also brings Heteroskedasticity and Autocorrelation Consistent (HAC) Covariance Matrix estimators that we can use to correct estimated standard errors, since under violation of TMSR2 the usual formula for standard errors is no longer correct. "The HAC variance estimate is equal to the HC variance estimate multiplied by an extra term that depends on the serial correlation in the errors." (PoE p 357)

coeftest(phill.dyn)
coeftest(phill.dyn, vcov.=vcovHAC(phill.dyn)) #Heteroskedasticity and Autocorrelation Consistent (HAC) Covariance Matrix Estimation - The HAC standard errors are larger than those from least squares, implying that if we ignore the autocorrelation, we will overstate the reliability of the least squares estimates.

##"Using least squares with HAC standard errors overcomes the negative consequences that autocorrelated errors have FOR LEAST SQUARES STANDARD ERRORS. However, it does not address the issue of finding an estimator that is better, in the sense that it has a lower variance. One way to proceed is to make an assumption about the model that generates the autocorrelated errors, and to derive an estimator compatible with this assumption" (PoE, p. 359)

## 9.5.2 ESTIMATING AN AR(1) ERROR MODEL
### 9.5.2b Nonlinear Least Squares Estimation
#### this model has the structure y_t = \beta_1 + \beta_2 x_t + e_t with e_t = \rho e_{t-1} + v_t. Since errors are unobservable, this means that the regression model for y_t is omitting one relevant term in explaining y_t. Our role now is to transform the model y_t = \beta_1 + \beta_2 x_t + e_t with the autocorrelated error e_t = \rho e_{t-1} + v_t into a new model that has an error term v_t uncorrelated over time. Through the transformations explained in PoE p. 361, we depart from equation (9.38) y_t = \beta_1 + \beta_2 x_t + e_t with e_t = \rho e_{t-1} + v_t to the new model that has an error term v_t uncorrelated over time in equation (9.43) y_t = \beta_1*(1-\rho) + \beta_2 x_t + \rho y_{t-1} + \rho*\betha_2 x_{t-1} + v_t.

#### Before we proceed, make sure you understand that in our example y_t = \beta_1 + \beta_2 x_t + e_t is represented by phill.dyn <- dynlm(inf ~ diff(u), data=phill.ts), which we showed to be serially correlated by using the Lagrange multiplier tests. Now we are going to depart from phill.dyn <- dynlm(inf ~ diff(u), data=phill.ts) to build a model that includes the AR(1) process and can represent the Phillips model with a new equation that has an error term v_t uncorrelated over time. This new model is estimated by phill.nls <- nls(inf ~ b1*(1-rho) + b2*Du + rho*Linf - rho*b2*LDu, data=phill.dfr, start=list(rho=0.5, b1=0.5, b2=-0.5)) below. "The coefficient of x_{t-1} is equal to \rho*\beta2 which is the negative product of \rho (the coefficient of y_{t-1}) and \beta_2 (the coefficient of x_t). This fact means that although (9.43) is a linear function of the variables x_t, y_{t-1} and x_{t-1}, it is not a linear function of the parameters (\beta_1, \beta_2, \rho)." (PoE p. 361)

summary(phill.dyn) #\beta_1 =  0.77762; \beta_2 = -0.52786
acf.ehat[0.25] #\rho=0.549

# Non-linear AR(1) from equation (9.45)
phill.ts.tab <- cbind(phill.ts[,"inf"],
                      phill.ts[,"u"],
                      lag(phill.ts[,"inf"], -1), 
                      diff(phill.ts[,"u"], lag=1),
                      lag(diff(phill.ts[,2],lag=1), -1))
phill.dfr <- data.frame(phill.ts.tab)
names(phill.dfr) <- c("inf", "u", "Linf", "Du", "LDu")
phill.nls <- nls(inf ~ b1*(1-rho) + b2*Du + rho*Linf - rho*b2*LDu,
						data=phill.dfr, start=list(rho=0.549, b1=0.77762, b2=-0.52786)) #nls() estimates nonlinear least squares
summary(phill.nls)

##9.5.3 ESTIMATING A MORE GENERAL MODEL -> equation (9.49)
if(!require(nlWaldTest)){
  install.packages("nlWaldTest")
  library(nlWaldTest)
}  

s.nls <- summary(phill.nls)
phill.gen <- dynlm(inf ~ L(inf) + d(u) + L(d(u)), data=phill.ts)
s.gen <- summary(phill.gen)
nlW <- nlWaldtest(phill.gen, texts="b[4]=-b[2]*b[3]") #testing the null hypothesis that H0 : \delta_1 = -\theta*\delta_0. 
nlW #p-value = 0.7376, so we can't reject H0 : \delta_1 = -\theta*\delta_0, meaning that in this specific case the more general model in phill.gen is as good as the somewhat more restrictive model in phill.nls regarding coefficient estimation. HOWEVER, "Specification and estimation of the more general model does have some advantages. It makes the dependence of y_t on its lag and that of x more explicit, and it can often provide a useful economic interpretation" (PoE p.364)
 
summary(dynlm(inf ~ L(inf) + d(u), data=phill.ts)) #equation (9.51)

# "While including one lag of y and one lag of x will correct for serially correlated errors if they follow an AR(1) model, it might not solve the problem if the form of serial correlation is more complex. How do we check whether some serial correlation still remains? If we include a lagged y and a lagged x and the errors are still serially correlated, how do we proceed? Checking for serial correlation proceeds along the same lines as we have described in Section 9.4. We apply the same tests to the errors from the new model with lags. Also, if we have doubts about whether the errors in the new model are correlated, we can use HAC standard errors with this model. Alternatively, including more lags of y on the right side of the equation can have the effect of eliminating any remaining serial correlation in the errors. Models with a general number of lags of y and x are called autoregressive distributed lag models; we consider them in the next section."(PoE p. 365)


################################################################
#9.6 Autoregressive Distributed Lag Models
################################################################
## Autoregressive distributed lag models are referred as ARDL(p,q) where p refers to the number of lags of y and q refers to the number of lags of x in the model. "(...) the main concern for estimation is choice of the lag lengths p and q. There are a number of different criteria for choosing p and q. Because they all do not necessarily lead to the same choice, there is a degree of subjective judgment that must be used. Four possible criteria are [listed on PoE p. 366]".

##9.6.1 THE PHILLIPS CURVE
phill1.gen <- dynlm(inf ~ L(inf) + d(u) + L(d(u)), data=phill.ts); summary(phill1.gen) #ARDL(1,1) - but the coefficient for L(d(u)) is no significant
phill2.gen <- dynlm(inf~L(inf)+d(u), data=phill.ts); summary(phill2.gen) #ARDL(1,0)
ehat.phill2.gen <- residuals(phill2.gen)
acf.ehat.phill2.gen <- acf(ehat.phill2.gen, lag.max=12, ylim = c(-0.3, 1), ylab = "Correlation",
						main = "Correlogram for residuals from least-squares estimated Phillips curve") 
acf.ehat.phill2.gen[c(0.25,0.50,0.75,1)]

##from the correlogram acf.ehat.phill2.gen autocorrelations are not significantly different from zero, they provide no evidence of serial correlation HOWEVER check the results of the Lagrange multiplier tests below for more objective information on autocorrelation
bgtest(phill2.gen, order=1, type="Chisq", fill=0) #reject H0 with p-value = 0.04213
bgtest(phill2.gen, order=2, type="Chisq", fill=0) #reject H0 with p-value = 0.07718
bgtest(phill2.gen, order=3, type="Chisq", fill=0) #can't reject H0 because p-value = 0.1563 (i.e. > 0.1)
bgtest(phill2.gen, order=4, type="Chisq", fill=0) #reject H0 with p-value = 0.04865
bgtest(phill2.gen, order=5, type="Chisq", fill=0) #reject H0 with p-value = 0.02872
###"Using a 5% significance level, tests for orders one, four, and five reject a null hypothesis of no autocorrelation. Taken together, these test results provide some evidence, but not overwhelming evidence, that serial correlation in the errors still exists; one lag of the dependent variable INF has not been sufficient to eliminate the autocorrelation." Let's try including more lags

#See PoE page 368 for the elimination criteria used to conclude that ARDL(4,0) seems to be the most adequate model 
phill3.gen <- dynlm(inf ~ L(inf, 1:4) + d(u), data=phill.ts); summary(phill3.gen) #equation (9.57)

aics <- rep(0,6)
bics <- rep(0,6)
for (i in 1:6){
  ari <- dynlm(inf ~ L(inf, 1:i) + d(u), start=i, data=phill.ts)
  aics[i] <- AIC(ari)
  bics[i] <- BIC(ari)
}
tbl1 <- data.frame(rbind(aics, bics))
names(tbl1) <- c("1","2","3","4","5")
row.names(tbl1) <- c("AIC","BIC")
tbl1  #ARDL(4,0) has minimum AIC BIC

aics <- rep(0,6)
bics <- rep(0,6)
for (i in 1:6){
  ari <- dynlm(inf~ L(inf, 1:i) + d(u) + L(d(u)), start=i, data=phill.ts)
  aics[i] <- AIC(ari)
  bics[i] <- BIC(ari)
}
tbl2 <- data.frame(rbind(aics, bics))
names(tbl2) <- c("1","2","3","4","5")
row.names(tbl2) <- c("AIC","BIC")
tbl2 #ARDL(4,1) has minimum AIC BIC, but AIC BIC for ARDL(4,0) are even lower


###9.6.2 OKUN'S LAW - PoE pages 369-370
summary(okunL2.dyn) # recall that okunL2.dyn <- dynlm(d(u) ~ L(g, 0:2), data=okun.ts) ; summary(okunL2.dyn)
ehat.okun <- residuals(okunL2.dyn)
acf(ehat.okun, lag.max=12, ylim = c(-0.4, 1), ylab = "Correlation", 
		main = "Correlogram for residuals from Okun's law ARDL(0,2) model.")
bgtest(okunL2.dyn, order=1, type="Chisq", fill=0) 

okunL1D1 <- dynlm(d(u) ~ L(g, 0:1) + L(d(u)), data=okun.ts); summary(okunL1D1) #(1,1) is the best according to table Table 9.5 PoE
AIC(okunL1D1); BIC(okunL1D1)
acf(residuals(okunL1D1), lag.max=12, ylim = c(-0.4, 1), ylab = "Correlation", 
		main = "Correlogram for residuals from Okun's law ARDL(1,1) model.")
bgtest(okunL1D1, order=1, type="Chisq", fill=0) ##can't reject H0 because p-value = 0.6804 (i.e. > 0.1) -> no serial correlation

################################################################
#9.7 Forecasting
################################################################
###9.7.2 FORECASTING WITH AN ARDL MODEL
####CHOOSING THE MODEL WE ARE GOING TO USE FOR PREDICTION
aics <- rep(0,5)
bics <- rep(0,5)
for (i in 1:5){
  ari <- dynlm(g ~ L(g,1:i), start=i, data=okun.ts)
  aics[i] <- AIC(ari)
  bics[i] <- BIC(ari)
}
tbl <- data.frame(rbind(aics, bics))
names(tbl) <- c("1","2","3","4","5")
row.names(tbl) <- c("AIC","BIC")
tbl #based on the BIC, which penalizes additional lags more heavily than does the AIC, we choose AR(2)

###PERFORMING THE PREDICTION BASED ON THE MODEL WE CHOSE ABOVE, I.E. AR(2)
if(!require(forecast)){
  install.packages("forecast")
  library(forecast)
} 

gdp <- okun[,"g"]
gdp.ts<-ts(gdp, start=c(1985,2), frequency=4)

ar2g <- ar(gdp.ts, aic=FALSE, order.max=2, method="ols") #aic = TRUE instead of FALSE selects the complexity of the model automatically from AIC
fcst <- data.frame(forecast(ar2g, 3)); fcst
plot(forecast(ar2g,3))

###9.7.3 EXPONENTIAL SMOOTHING - the name comes from the fact that the weights attributed to each point in the series decline exponentially as the observations get older. "The smaller the value of \alpha, the greater the contribution of past observations to a forecast, and the smoother the series of within-sample forecasts is. With large values of \alpha, the most recent observation is the major contributor to a forecast, and the series of forecasts more closely mimics the actual series." PoE p. 377
okun.HW1 <- HoltWinters(gdp.ts, beta=FALSE, gamma=FALSE); okun.HW1
plot.okun.HW1 <- plot(okun.HW1, main = "Holt-Winters smoothing for alpha=0.3805126")
okunf.HW1 <- forecast(okun.HW1, 1, main = "T+1 Holt-Winters forecast for alpha=0.3805126"); okunf.HW1
plot.okunf.HW1 <- plot(okunf.HW1)

okun.HW2 <- HoltWinters(gdp.ts, alpha=0.8, beta=FALSE, gamma=FALSE); okun.HW2
plot.okun.HW2 <- plot(okun.HW2, main = "Holt-Winters smoothing for alpha=0.8")
okunf.HW2 <- forecast(okun.HW2,1, main = "T+1 Holt-Winters forecast for alpha=0.8"); okunf.HW2
plot.okunf.HW2 <- plot(okunf.HW2)

# "In Figure 9.12(a), where \alpha 0.38, the smoothed series is much less volatile than the actual series. It retains a jagged appearance, but the peaks and troughs are much less extreme. In Figure 9.12(b), where \alpha = 0.8, the peaks and troughs of the smoothed series are only slightly less pronounced, and the forecasts closely follow the actual series by one period reflecting the high weight placed on the most recent value." (PoE p. 377)


#######################################################################################################
# If you have questions or suggestions regarding this script, write me: robson.tigre@unibo.it 
# For an alternative replication structure of Chapter 9, check the companion: 
# https://bookdown.org/ccolonescu/RPoE4/time-series-stationary-variables.html (Accessed on Nov 13 2017)
#######################################################################################################