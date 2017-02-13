#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){

  int i, iMAX;
  double x, dx, fx;

  iMAX = 8;
  dx = 10./iMAX;

  double pi = 3.14;
	for(i=0; i<=iMAX; i++){

	  x = i*dx;
	  //fx = -(x-4)*(x-4) + 4;
	  fx = sqrt(2/10.)*sin(pi*x/10.);

	  double px = fx*fx;
	  printf(" %i %f %f \n",i,x,px);	

	}

}
