#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
// I = imaginary unit
#include<malloc.h>

double pi = 3.14159265359;
void dfdt(int dim,double complex *psivec,double complex *dpsi,double dx);
// Function for taking 2nd derivative
// Arguments: arg1 - dim, number of pts in wavefxn
//            arg2 - psivec is array of wavefxn values
//            arg3 - dpsi is array of d/dt wavefxn values
//            arg4 - dx is increment along x axis
//            void dftdt(int dim,double *psivec,double *dpsi,double dx);

int main(){

  int dim = 10;
  double *x;
  double complex *wfn, *dpsi;
  double L = 1.;

  double dx;
  int i;
 
  dx = L/dim;
  x = (double *)malloc(dim*sizeof(double));
  wfn = (double complex *)malloc(dim*sizeof(double complex));
  dpsi  = (double complex *)malloc(dim*sizeof(double complex));

  for (i=0; i<=dim; i++) {

    x[i] = i*dx;
    wfn[i] = sqrt(2./L)*sin(pi*x[i]/L) + 0.*I;
    
    printf(" %f %f %f\n",x[i],creal(wfn[i]),cimag(wfn[i]));
    
  }

  dfdt (dim,wfn,dpsi,dx);
  
  for (i=0; i<=dim; i++) {

    printf(" %f (%f, %f) (%f, %f)\n",x[i],creal(dpsi[i]),cimag(wfn[i])*pi*pi/2.,cimag(dpsi[i]),-1*creal(wfn[i])*pi*pi/2.);

  } 

}


void dfdt(int dim,double complex *psivec,double complex *dpsi,double dx) {

  int j;
  dpsi[0] = 0. + 0.*I;
  dpsi[dim] = 0. + 0.*I;

  for (j=1; j<dim; j++) {

    dpsi[j] = (I/2.)*(psivec[j+1] - 2*psivec[j] + psivec[j-1])/(dx*dx);

  }

}

