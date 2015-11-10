
library("Rcpp")
library("RcppArmadillo")
library("MASS")
library("R.matlab")
require(R.utils)

setwd("C:\\Users\\david\\Documents\\GitHub\\Didris")

rm(list=ls())
set.seed(90)

sourceCpp("general_eq.cpp")

A = 1;
beta = 0.971;
psi = 2;
eta = 0.5;
Abar = 0;
gamma = 0.66;
alpha_h = 0.5;
z = 1.5;
phi = 2;

param   = array(0,dim=c(4));
param[1] = beta;
param[2] = eta;
param[3] = psi;
param[4] = Abar;
param[5] = A;
param[6] = gamma;
param[7] = alpha_h;
param[8] = z;
param[9] = phi;

minb = 0;
maxb = 5;
numb = 51;

minthe = 0.01;
maxthe = 0.99;
numthe = 51;

wl = 0.2955;
wh = 1.8346;
r = 0.1754;
Ph = 1;

# 1. First estimation of utility for each type (b, theta)
start = proc.time()[3];

Result = utility(param, minb, maxb, numb, minthe, maxthe, numthe, wl, wh, r, Ph);

finish = proc.time()[3] - start;
print(finish)

H = array(0,dim=c(numb, numthe));
A = array(0,dim=c(numb, numthe));
U = array(0,dim=c(numb, numthe));
C = array(0,dim=c(numb, numthe));

for (ib in 1:numb){
  for (it in 1:numthe){
    H[ib, it] = Result[[1]][[ib]][[it]];
    A[ib, it] = Result[[2]][[ib]][[it]];
    U[ib, it] = Result[[3]][[ib]][[it]];
    C[ib, it] = Result[[4]][[ib]][[it]];
  }
}


# 2. Estimation of general equilibrium prices: w_l, w_h, r

size = 10;

HH = array(0,dim=c(size, size));
LL = array(0,dim=c(size, size));
KK = array(0,dim=c(size, size));
WL = array(0,dim=c(size, size));
WH = array(0,dim=c(size, size));
RR = array(0,dim=c(size, size));
codes = array(0,dim=c(size, size));
minim = array(0,dim=c(size, size));

zz = seq(1.1, 2.8, length=size);
PP = seq(0.2, 1, length=size);

# Starting point
w_l = 0.3; 
w_h = 1.94; 
r = 0.08;

for(iz in 1:size){
  for(ip in 1:size){
    z = zz[iz]; 
    Ph = PP[ip];
    
    param[8] = z;
    
    f = function(x) residual(param, minb, maxb, numb, minthe, maxthe, numthe, exp(x[1]), exp(x[2]), exp(x[3]), Ph);
    
    tryCatch(
      
      expr = {
        evalWithTimeout( {
          
          x = nlm(f, c(w_l, w_h, r), iterlim = 1);
          
          sol = x$estimate;
          w_l = sol[1];
          w_h = sol[2];
          r = sol[3];
          
          x = nlm(f, c(w_l, w_h, r), iterlim = 100);
          
          sol = x$estimate;
          w_l = exp(sol[1]);
          w_h = min(exp(sol[2]), 10^10);
          r   = exp(sol[3]);
          
          caps = general_eq(param, minb, maxb, numb, minthe, maxthe, numthe, w_l, w_h, r, Ph);
          
          codes[iz, ip] = x$code;
          minim[iz, ip] = x$minimum;
          HH[iz, ip] = caps[2];
          LL[iz, ip] = caps[3];
          KK[iz, ip] = caps[4];
          WL[iz, ip] = sol[1];
          WH[iz, ip] = sol[2];
          RR[iz, ip] = sol[3]; }, timeout = 20)
        },
      TimeoutException = function(ex) cat("Timeout. Skipping.\n")
    )
    
    print(paste("Zeta: ", toString(iz), "Price: ", toString(ip)));
  }
}
