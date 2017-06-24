#Calibration
rm(list=ls())
library(deSolve) # library for solving differential equations
library(minpack.lm) #for non-linear regression package
d = seq(0, 0.1, 0.01)

#First IDER
IDER_1 = function(d, a, b) {
  1 - exp(-(a*d + b*d^2))
} 

#Second IDER
IDER_2 = function(d, k) {
  1 - exp(-(exp(k*d) - 1))
}

#third IDER
dE = function(t, E, parms) {
    c = parms$c;
    a = parms$a;
    dH = c*(a-log(1-E))*(0.5+log(1-E))*(1-E)
    return(list(dH))
}

IDER_3 = function(d, c, a) {
  yini = c(H= 0)
  out = ode(y = yini, times = d, func = dE, parms = list(c = c, a = a))
  return(out[, 2])
}


#Mixture

dE_1 = function(d, a, b) {
  (2*b*d + a)*(exp(-b*d^2 - a*d))
}

dE_2 = function(d, k) {
  k*exp(-exp(k*d) + k*d + 1)
}

dE_3 = function(E, a, c) {
  c*(a-log(1-E))*(0.5+log(1-E))*(1-E)
}

MIXDER_function = function(r , a1 = 1, b = 2, k = 4, a2 = 0.2, c = 10, d)  {
  dE=function(yini,State,Pars){
  a1 = a1; b = b; a2 = a2; k = k; c = c;
  with(as.list(c(State, Pars)), {
    u = vector(length = 3)
    u[1] = uniroot(function(d) (1 - exp(-(a1*d + b*d^2))) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root 
    u[2] = uniroot(function(d) (1 - exp(-(exp(k*d) - 1))) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
    dI = vector(length = 3)
    dI[1] = r[1]*dE_1(d = u[1], a = a1, b = b)
    dI[2] = r[2]*dE_2(d = u[2], k = k)
    dI[3] = r[3]*dE_3(E = I, a = a2, c = c)
    dI = sum(dI)
    return(list(c(dI)))
      })
    }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini,times = d, dE, pars, method = "radau")
  return(out)
}
d = seq(0, 1.5, 0.01)

#setEPS()
#postscript("mouse_MIXDER_one_third.eps", width = 10, height = 7)
r3=.9; r2=(1-r3)/2; r1=r2
plot(x = d, y = IDER_2(d = d, k = 4), col = "green", type = "l",lwd=2, ann="F")
lines(x = d, y = MIXDER_function(r = c(r1,r2,r3), d = d)[, 2], col = "red",lwd=2)
lines(x = d, y = IDER_1(d = d, a = 1, b = 2), col = "blue")
lines(x = d, y = IDER_3(d = d, c = 10, a = 0.1))
#dev.off()


