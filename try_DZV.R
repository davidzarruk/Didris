
library("Rcpp")
library("RcppArmadillo")
library("MASS")
library("R.matlab")

setwd("C:\\Users\\david\\Documents\\GitHub\\Didris")

rm(list=ls())
set.seed(90)

sourceCpp("compall_sinGSL.cpp")
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
numb = 301;

minthe = 0.01;
maxthe = 0.99;
numthe = 301;

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

HH = array(0,dim=c(25));
LL = array(0,dim=c(25));
KK = array(0,dim=c(25));
WL = array(0,dim=c(25));
WH = array(0,dim=c(25));
RR = array(0,dim=c(25));

zz = seq(1.1, 1.7, length=5);

for(iz in 1:5){

  z = zz[iz]; 
  param[8] = z;
  
  f = function(x) residual(param, minb, maxb, numb, minthe, maxthe, numthe, x[1], x[2], x[3], Ph);
  
  x = nlm(f, c(0.3, 1.94, 0.08));
  
  sol = x$estimate;
  w_l = sol[1];
  w_h = sol[2];
  r = sol[3];
  
  caps = general_eq(param, minb, maxb, numb, minthe, maxthe, numthe, w_l, w_h, r, Ph);
  
  HH[iz] = caps[2];
  LL[iz] = caps[3];
  KK[iz] = caps[4];
  WL[iz] = sol[1];
  WH[iz] = sol[2];
  RR[iz] = sol[3];
  
  print(iz);
  
}
