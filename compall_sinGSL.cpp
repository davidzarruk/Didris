// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
using std::vector;

using namespace arma;
using namespace Rcpp;


/*** Absolute value of elements of a vector ***/
double abs(double a) {
  
  double x = a;

  if(x<0){
    x = -a;
  }
  return x;
}

//I. Solving the optimal level of a_{t+1}. First I define a c-structure of
//parameters for the solver

double resid(double b, double a, double wl, double wh, double i, double Ph, double ppsi, double ttheta, double eeta, int h, double bbeta){
  double residual=pow(b-a+wl*(1-h)-Ph*h,-ppsi)-bbeta*(1+i)*(pow(h*ttheta,eeta)*pow((wh+(1+i)*a),-ppsi)+(1-pow(h*ttheta,eeta))*pow(wl+(1+i)*a,-ppsi));
  return residual;
}

double bisection(double b, double amin, double amax, double wl, double wh, double i, double Ph, double ppsi, double ttheta, double eeta, int h, double bbeta){
  
  double p = 0;
  double fa = resid(b, amin, wl, wh, i, Ph, ppsi, ttheta, eeta, h, bbeta);
  double fb = resid(b, amax, wl, wh, i, Ph, ppsi, ttheta, eeta, h, bbeta);
  double mult = fa*fb;
  
  if(mult <= 0){
    p = (amax + amin)/2;
    double fp = resid(b, p, wl, wh, i, Ph, ppsi, ttheta, eeta, h, bbeta);
    double err = abs(fp);
    while(err>pow(10,-7)){
      fp = resid(b, p, wl, wh, i, Ph, ppsi, ttheta, eeta, h, bbeta);
      mult = fa*fp;
      if(mult < 0){
        amax = p;
      }
      else{
        amin = p;
      }
      p = (amax + amin)/2;
      err = abs(fp);
    }
  }
  else{
    Rcout << "La cagaste weeeeyyyy" << std::endl;
  }
  return p;
}

vector<double> grid(double min, double max, int num){
  vector<double> ans;
  
  ans.resize(num);
  
  for(int i=0; i<num; i++){
    ans[i] = min + ((max - min)/(num-1))*i;
  }
  return ans;
}


// [[Rcpp::export]]
vector<vector<vector<double> > > utility(vec par, double minb, double maxb, int numb, double minthe, double maxthe, int numthe, double wl, double wh, double r, double Ph){

    /*** Initializing matRices  ***/
  vector<vector<vector<double> > > ans;
  ans.resize(4);
  for(int k=0; k<4; k++){
    ans[k].resize(numb);
    for(int j=0; j<numb; j++){
      ans[k][j].resize(numthe);
      for(int l=0; l<numthe; l++){
        ans[k][j][l] = 0;
      }
    }
  }
  //0. Initializing utility if education or not
  double ht = 0;
  double at = 0;
  double ut = 0;
  double ct = 0;
  
  //1. Loading parameters
  double bbeta=par[0];
  double eeta=par[1];
  double ppsi=par[2];
  double Abar = par[3];
  vector<double> bgr = grid(minb, maxb, numb);
  vector<double> tthetagr = grid(minthe, maxthe, numthe);

  double b = 0.0;
  double ttheta = 0.0;
  
  for(int k=0; k<numb; k++){
    for(int j=0; j<numthe; j++){
      
      //0. Initializing utility if education or not
      ht = 0;
      at = 0;
      ut = 0;
      ct = 0;
      
      b = bgr[k];
      ttheta = tthetagr[j];
      
      double amin = -wl/(1+r) + pow(10, -10);
      double amax = b + wl - pow(10, -10);
      
      double a0 = bisection(b, amin, amax, wl, wh, r, Ph, ppsi, ttheta, eeta, 0, bbeta);
      
      if(a0 < -Abar){
        a0 = -Abar;
      } else if(a0 > b + wl){
        a0 = b + wl - pow(10,-10);
      }
      
      if(-Abar < b - Ph){
        if(ttheta != 1){
          amin = -wl/(1+r) + pow(10, -10);
        } else{
          amin = -wh/(1+r) + pow(10, -10);
        }
        
        amax = b - Ph - pow(10, -20);
        
        double a1 = bisection(b, amin, amax, wl, wh, r, Ph, ppsi, ttheta, eeta, 1, bbeta);
        
        if(a1 < -Abar){
          a1 = -Abar;
        } else if(a1 > b - Ph){
          a1 = b - Ph - pow(10,-10);
        }
        
        double u0 = pow(b - a0 + wl, 1 - ppsi)/(1 - ppsi) + bbeta*pow(wl + (1+r)*a0, 1 - ppsi)/(1 - ppsi);
        double u1 = pow(b - a1 - Ph, 1 - ppsi)/(1 - ppsi) + bbeta*(pow(ttheta, eeta)*pow(wh + (1+r)*a1, 1 - ppsi)/(1-ppsi) + (1-pow(ttheta, eeta))*pow(wl + (1+r)*a1, 1 - ppsi)/(1 - ppsi));
        
        if(u0 >= u1){
          ht = 0;
          at = a0;
          ut = u0;
          ct = b - a0 + wl;
        } else{
          ht = 1;
          at = a1;
          ut = u1;
          ct = b - a1 - Ph;
        }
        
      } else{
        ht = 0;
        at = a0;
        ut = pow(b - a0 + wl,1 - ppsi)/(1 - ppsi) + bbeta*pow(wl + (1+r)*a0, 1-ppsi)/(1-ppsi);
        ct = b - a0 + wl;
      }
      ans[0][k][j] = ht;
      ans[1][k][j] = at;
      ans[2][k][j] = ut;
      ans[3][k][j] = ct;
    }
  }
  return ans; 
}
