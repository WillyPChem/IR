#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double pi = 4.*atan(1.0);
double Eigenfunction(double x, double alpha, int quantumnumber);



int main(){

  double psi;
  double x = 1.0;
  double alpha = 8.7e7;
  int n = 0;

  psi = Eigenfunction( 0.001, 8.7e8, 40);

  

}

double Eigenfunction(double x, double alpha, int n) {

  double term1;
  double term2;
  double term12;

  term1 = pow(alpha/pi,1./4); 
  term2 = exp(-alpha*x*x/2);
  term12 = term1*term2; 
  
  printf(" %e",term12);
  return term12;
}
