#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){

  int i, iMAX;
  double x, dx, fx;

  //double pi = 4*atan(1.0);
  double  pi = 3.141592653589793;
  
  dx = 10./iMAX;
  iMAX = 8;


	for(i=0; i<=iMAX; i++){

	  x = i*dx;
	  fx = sqrt(2/10.)*sin(pi*x/10.);

	  printf(" %i %f %f \n",i,x,fx);


	}
}
