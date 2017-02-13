#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){

  int i, iMAX;
  double x, dx, fx, fpx, fxf, xf;

  iMAX = 8;
  dx = 10./iMAX;

  // rise = fxf-fx
  // run = xf - x
  // slope = (fxf-fx)/(xf-x)

  double pi = 3.14;
	for(i=0; i<=iMAX; i++){

	  x = i*dx;
	  //fx = -(x-4)*(x-4) + 4;
	  //fx = sqrt(2/10.)*sin(pi*x/10.);


	  // Derivative below thus line
	  fx = sqrt(2/10.)*sin(pi*x/10.);
	  xf = (i+1)*dx;
	  fxf = sqrt(2/10.)*sin(pi*fxf/10.);
	
	  double ror = (fxf-fx)/(xf-x); 

	  fpx = sqrt(2/10.)*(pi/10.)*cos(pi*x/10.);	

	  double px = fx*fx;
	  printf(" %i %f %f %f \n",i,x,ror,fpx);	

	}

}
