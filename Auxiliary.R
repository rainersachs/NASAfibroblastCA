#An R version of parts of NASAfibroCA.Rmd
#All the libraries needed for this script
library(deSolve) # library for solving differential equations
library(minpack.lm) #for non-linear regression package
rm(list=ls())
# Create dataframe that stores the data
Oxygen = data.frame(d = c(0, .0125, .02, .025, .05, .075, .1, .2, .4), 
                    CA = c(.24, 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91))

Si = data.frame(d = c(0, .02, .04, .06, .08, .1, .12, .2, .4, .8, 1.2), 
                    CA = c(.11, 1.26, 1.14, 1.58, 1.22, 1.89, 3.47, 4.6, 9.79, 27.01, 38.84))

Fe600 = data.frame(d = c(0, .01, .02, .04, .06, .08, .1, .12, .2, .4, .8), 
                     CA = c(.13, .76, .99, 1.2, 1.74, 1.28, 1.2, 1.7, 3.02, 5.52, 12.42))

Fe450 = data.frame(d = c(0, .02, .04, .06, .08, .1, .2, .4), 
                   CA = c(0, .86, .6, .8, 1.22, 2.02, 2.3, 4.77))

Fe300 = data.frame(d = c(0, .005, .01,  0.02, .04, .07, .1, .2, .4, .8), 
                   CA = c(0.41, 1.23, 1.47, 1.22, .97, 1.46, 1.21, 4.38, 6.22, 13.6))

Ti = data.frame(d = c(0,  0.02, .04, .06, .08, .1, .15, .3, .6), 
                   CA = c(0, 1.99, 1.88, 1.44, 2.67, 2.57, 2.50, 5.64, 11.19))

param = data.frame(ion = c("O", "Si", "Ti", "Fe600", "Fe450", "Fe300"),
                   Z = c(8, 14, 22, 26, 26, 26), L = c(75, 100, 125, 175, 195, 240), 
                   Z.beta = c(595, 690, 770, 1075, 1245, 1585))

#putting it in one big data frame
big_df = rbind(Oxygen, Si, Ti, Fe600, Fe450, Fe300)
big_df$Z = rep(param$Z, times = c(9, 11, 9, 11, 8, 10))
big_df$Z.beta = rep(param$Z.beta, times = c(9, 11, 9, 11, 8, 10))
big_df$L = rep(param$L, times = c(9, 11, 9, 11, 8, 10))
big_df$error = c(0.24, 0.63, 0.77, 0.75, 0.52, 0.82, 0.78, 1.31, 1.59, 0.12, 0.05, 0.07, 0.56, 0.18, 0.60, 1.23, 1.60, 1.55, 4.27, 7.21, 0, 0.70, 0.66, 0.59, 0.80, 0.78, 0.48, 1.15, 2.39, 0.16, 0.38, 0.24, 0.21, 0.02, 0.37, 0.54, 0.17, 0.55, 1.75, 2.59, 0, 0.43, 0.34, 0.40, 0.50, 0.64, 0.73, 1.09, 0.29, 0.55, 0.60, 0.55, 0.49, 0.60, 0.54, 1.03, 1.22, 3.62)
big_df$ion = rep(param$ion, times = c(9, 11, 9, 11, 8, 10))

#will modify the data frame to get rid of the zero dose points irrelevant to the parameter estimation
modified_df = big_df[big_df$d != 0, ]
modified_df$CA = modified_df$CA*0.01
modified_df$error = modified_df$error*0.01
big_df$CA = big_df$CA * 0.01
big_df$error = big_df$error * 0.01
big_df$errorbar_lower = big_df$CA - big_df$error
big_df$errorbar_upper = big_df$CA + big_df$error

#NTE1 and NTE2 models in Cacao's report using their parameters
#NTE1 function
NTE1_function = function(d, L, Z.beta, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796) {
0.0017 + eta0*L*exp(-eta1*L)*(d != 0) + 
    (6.242*(d/L))*(sig0*(1-exp(-Z.beta/kap))^2 + 0.041/6.24*L*(1 - (1-exp(-Z.beta/kap))^2))
} 

#NTE2 function
NTE2_function = function(d, L, Z.beta, eta0 = 0.00047, eta1 = 0.011, sig0 = 6.75, kap = 590) {
0.0017 + eta0*L*exp(-eta1*L)*exp(-(1012*(d/L)))*(d != 0) + #what is .0017? mistake?
    (6.242*(d/L))*(1-exp(-(1012*(d/L))))*
    (sig0*(1-exp(-Z.beta/kap))^2 + 0.041/6.24*L*(1 - (1-exp(-Z.beta/kap))^2))
} 

#Our IDER
#our proposed function#RKS double check this and NTE1, NTE2
IDER = function(d, L, Z.beta, eta0, eta1, sig0, kap) {
  P = (1-exp(-Z.beta/kap))^2
  sig = sig0*P + 0.041/6.24*L*P
  eta = eta0*L*exp(-eta1*L)
  0.00071+ sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))
} 

#nls method to get the parameters needed (4 parameter estimation)
IDER_model = nlsLM(CA ~ IDER(d, L, Z.beta, eta0, eta1, sig0, kap), data = modified_df, start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5, kap = 500), 
weights = (1/(modified_df$error)^2))
coef(IDER_model)
#eta0 = 1.479796e-04, eta1 = 3.480127e-03, sig0 = 2.467512e+00, kap = 2.523256e+02 #RKS checks with MS
vcov(IDER_model)# checked MS  RKS
summary(IDER_model, correlation=T) #not same as cor(vcov(IDER_model))RKS
#When creating any sort of figures we use this IDER now without the background effect then add the background effect at the end.
IDER = function(d, L, Z.beta, eta0 = 1.479796e-04, eta1 = 3.480127e-03, sig0 = 2.467512e+00, kap = 2.523256e+02) {
  P = (1-exp(-Z.beta/kap))^2
  sig = sig0*P + 0.041/6.24*L*P
  eta = eta0*L*exp(-eta1*L)
  sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))#+0.00071 RKS
} 
#Information crieteria (AIC and BIC)#RKS spelling
L_function = function(func, eta0, eta1, sig0, kap) {
  a = vector(length = 0)
  for (i in 1:length(modified_df[, 1])) {
    a = c(a, log(func(d = modified_df$d[i], L = modified_df$L[i], Z.beta = modified_df$Z.beta[i], eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap)))
  }
  return(sum(a))
}
L_NTE1 = L_function(NTE1_function, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796)
L_NTE2 = L_function(NTE2_function, eta0 = 0.00047, eta1 = 0.011, sig0 = 6.75, kap = 590)
L_IDER = L_function(IDER, eta0 = 1.479796e-04, eta1 = 3.480127e-03, sig0 = 2.467512e+00, kap = 2.523256e+02 )

AIC_function = function(L, k = 2) {
  2*k - 2*L
}

BIC_function = function(n = length(modified_df[, 1]), k = 2, L) {
  log(n)*k - 2*L
}

NTE1_AIC = AIC_function(L = L_NTE1)
NTE2_AIC = AIC_function(L = L_NTE2)
IDER_AIC = AIC_function(L = L_IDER)
NTE1_BIC = BIC_function(L = L_NTE1)
NTE2_BIC = BIC_function(L = L_NTE2)
IDER_BIC = BIC_function(L = L_IDER)
information_critera_df = data.frame(AIC = c(NTE1_AIC, NTE2_AIC, IDER_AIC), BIC = c(NTE1_BIC, NTE2_BIC, IDER_BIC), row.names = c("NTE1 model", "NTE2 model", "IDER model"))
information_critera_df

#Creating a general function that will give mix_IDER based on any input
MIXDER_function = function(r, L, Z.beta, d = seq(0, 0.2, by = 0.001), eta0 = 1.300771e-04, eta1 = 3.164156e-03, sig0 = 2.481817e+00, kap = 2.565276e+02) {
  dE=function(yini,State,Pars){
  eta0 = eta0; eta1 = eta1; sig0 = sig0; kap = kap
  with(as.list(c(State, Pars)), {
    P = vector(length = length(L))
    sig = vector(length = length(L))
    etaa = vector(length = length(L))
    u = vector(length = length(L))
    for (i in 1:length(L)) {
      P[i] = (1-exp(-Z.beta[i]/kap))^2
      sig[i] = sig0*P[i] + 0.041/6.24*L[i]*P[i]
      etaa[i] = eta0*L[i]*exp(-eta1*L[i])
      u[i] = uniroot(function(d) sig[i]*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + etaa[i]*(1-exp(-10^5*d)) - I, lower = 0, upper = 1, tol = 10^-10)$root
    }
    dI = vector(length = length(L))
    for (i in 1:length(L)) {
      dI[i] = r[i]*(sig[i]*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + etaa[i]*10^5*exp(-10^5*u[i]))
    }
    dI = sum(dI)
    return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini,times = d, dE, pars, method = "radau")
  return(out)
}
#Note I had to make ode solver and uniroot to work together in this algorithmn. The alogirthmn is a recursive loop 
# that finds the root at every turn of the inputted dose, quite computationally taxing. #Note: Had to use method = "radau" or else
#the small dose points with high concavity will not have accurate effect and causes massive accumulated errors later.#RKS false?
r=rep(1/6,6);L = c(75, 100, 125, 175, 195, 240); Z.beta = c(595, 690, 770, 1075, 1245, 1585)
ttrial=MIXDER_function(r,L,Z.beta,d=seq(0,.4,by=0.001))#for this mixture can't go much above 0.6#RKS  MIXIDER -->MIXDER
plot(ttrial[,1],ttrial[,2],type='l')


#Graphs for NTE1 and IDER model 12_ion figure


#Comparing their NTE1 model and our IDER Model

#Use these dose points toe generate the curves necessary. Had to include a lot more dose points at lower doses when it is very sensitive to change
# d_O = c(1e-6*(1:9), 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:9), 0.01*(10:40))
# d_1 = c(1e-6*(1:9), 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:9), 0.01*(10:15))
# 
# #setEPS()
# postscript("NTE1_model test.eps", width = 5, height = 7)
# par(mfrow=c(2,1))
# 
# plot(x = big_df$d[1:9]*100, y = big_df$CA[1:9]*100, main = "Oxygen", ylim = c(0, 10),pch=19, ylab = "CA", xlab = "Dose (cGy)")
# arrows(big_df$d[1:9]*100, big_df$errorbar_lower[1:9]*100, big_df$d[1:9]*100, big_df$errorbar_upper[1:9]*100, length=0.02, angle=90, code=3)
# lines(x = d_O*100, y = 100*(0.00071 + IDER(d = d_O, L = 75, Z.beta = 595)), col = "red")
# lines(x = d_O*100, y = 100 * NTE1_function(d = d_O, L= 75, Z.beta = 595), col = "blue", lty = 2)
# 
# plot(x = big_df$d[1:7]*100, y = big_df$CA[1:7]*100, main = "Oxygen",pch=19, ylab = "CA", xlab = "Dose (cGy)", ylim = c(0, 4), xlim = c(0, 15))
# arrows(big_df$d[1:7]*100, big_df$errorbar_lower[1:7]*100, big_df$d[1:7]*100, big_df$errorbar_upper[1:7]*100, length=0.02, angle=90, code=3)
# lines(x = d_1*100, y = 100*(0.00071 + IDER(d = d_1, L = 75, Z.beta = 595)), col = "red")
# lines(x = d_1*100, y = 100 * NTE1_function(d = d_1, L= 75, Z.beta = 595), col = "blue", lty = 2)
# 
# dev.off()