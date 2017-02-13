#include<stdio.h>
#include<stdlib.h>

int main(){

  int i, iMAX;
  double x, dx, fx;

  iMAX = 8;
  dx = 4./iMAX;

	for(i=0; i<=iMAX; i++){

	  x = i*dx;
	  fx = -(x-4)*(x-4) + 4;

	  printf(" %i %f %f \n",i,x,fx);	

	}

}
