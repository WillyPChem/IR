#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>
// I = imaginary unit
#include<malloc.h>

/*
Step 1: Take 2nd derivative of Psi(x) = ground-state harmonic oscillator
Step 2: Multiply by -hbar^2/2mu
Step 3: Add 1/2kx^2 term
Step 4: Subtract -qE(t)x term
*/

// Global Constants
double pi,alpha,mu,hbar,q,k;
  pi = 4*arctan(1);
  alpha = 8.386*pow(10,7);
  mu = 1.626*pow(10,-27);
  hbar = 1.0546*pow(10,-34);
  q = 1.60*pow(10,-19);
  k = 481.;

// Function Prototypes

// 2nd derivative of psi_j/dt
void dfdt(int dim, double complex *psivec, double complex *dpsi, double dx);
	// dim = number of points in wavefxn
	// psivec = array of wavefxn values
	// d2psi = array of dpsi^2/dx^2 wavefxn values
	// dx = increment along x axis

// H hat applied to psi
void Hpsi(int dim, double complex *dpsi, double complex *psivec, double complex *dpsij, double dx, double *x);
	// dim = number of points in wavefxn
	// psivec = array of wavefxn values
	// d2psi = array of dpsi^2/dx^2 wavefxn values
	// dx = increment along x axis
	// x = location along curve

// Electronic Field,E(x)
double E(double x);
	// x = location along curve

// 3rd Order Runge–Kutta :: Advances wavefxn forward in time
void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);
	// dim = number of points in the wavefunction
	// xvec = vector of x-values that the wfn is evaluated at
	// wfn = vector that stores the wfn values at time t_i
	// 		when the function is called, and when the function has been executed,
	//		this vector will hold the wfn values at time t_f = t_i + dt
	// dx = differential step along x between subsequent points in *xvec
	// dt = differential step along time between subsequent times (i.e. t_f - t_i)

/* **************************************************************************************************** */

int main()
{
	// Initialize & Define Variables
	int dim = 10;
	double *x;
	double complex *wfn, *dpsi, *dpsij;
	double L = 1.;
	double dx = L/dim;
	int i;

	// Arrays via malloc()
	x = (double *)malloc(dim*sizeof(double));
	wfn = (double complex *)malloc(dim*sizeof(double complex));
	dpsi = (double complex *)malloc(dim*sizeof(double complex));
	dpsij = (double complex *)malloc((dim+1)*sizeof(double complex));

	// Psi(x)
	for (i=0; i<dim; i++)
	{
		x[i] = i*dx;
		wfn[i] = pow(alpha/pi,1/4)*exp((-alpha/2)*pow(x[i],2));

		// Particle in a box
		//wfn[i] = sqrt(2./L)*sin(pi*x[i]/L) + 0.*I;
		printf("wfn[%d] = %f %f*i\n",i,creal(wfn[i]),cimag(wfn[i]));
	}
	printf("\n\n");
	// Euler Step
	dfdt(dim,wfn,dpsij,dx);

    // Apply -hbar^2/2mu to Euler Step
    Hpsi(dim,dpsi,wfn,dpsij,dx,x);

    // print dpsi/dt
    for (i=0; i<dim; i++)
	{
		printf("dpsi[%d] = %f %f*i\n",i,creal(dpsi[i]),cimag(dpsi[i]));
	}

    // Free memory
    free(x);
    free(wfn);
    free(dpsi);
    free(dpsij);

	return 0;
}

/* **************************************************************************************************** */

void dfdt(int dim,double complex *psivec,double complex *dpsij,double dx) {
	int j;
	dpsij[0] = 0. + 0.*I;
	dpsij[dim] = 0. + 0.*I;

	for (j=1; j<dim; j++)
	{
		dpsij[j] = (1.*I/2)*((psivec[j+1] - (2*psivec[j]) + psivec[j-1])/(dx*dx));
	}
}

/* **************************************************************************************************** */

void Hpsi(int dim, double complex *dpsi, double complex *psivec, double complex *dpsij, double dx, double *x)
{
	int j;
	for (j=0; j<dim; j++)
	{
		// (1) is placeholder for d^2/dx^2 aka dpsi[i]
		dpsi[j] = - ((pow(hbar,2)/(2*mu))*(dpsi[i])) + ((0.5*k*pow(x[j],2))*(psivec[j])) - ((q*E(x[j])*x[j])*(psivec[j]));
	}
}
/* **************************************************************************************************** */

double E(double x)
{
	return 0;
}

/* **************************************************************************************************** */

void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt) {
	double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;
	int i;

	// Temporary arrays for computing derivatives of wfns and approximate updates to wfns
	wfn_dot = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn2 = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn3 = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn_np1 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k1 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k2 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k3 = (double complex *)malloc((dim+1)*sizeof(double complex));

  //Initialize all (real and imaginary parts) to zero
	for (i=0; i<=dim; i++)
	{
		wfn_dot[i] = 0. + 0.*I;
		wfn2[i] = 0. + 0.*I;
		wfn3[i] = 0. + 0.*I;
		wfn_np1[i] = 0. + 0.*I;
		k1[i] = 0. + 0.*I;
		k2[i] = 0. + 0.*I;
		k3[i] = 0. + 0.*I;
	}

	// Get dPsi(n)/dt at initial time!
	dfdt(dim, wfn, wfn_dot, dx);

	// Compute approximate wfn update with Euler step
	for (i=0; i<=dim; i++)
	{
		k1[i] = dt*wfn_dot[i];
		wfn2[i] = wfn[i] + k1[i]/2.;
	}

	// Get dPsi(n+k1/2)/dt
	dfdt(dim, wfn2, wfn_dot, dx);

	// Compute approximate wfn update with Euler step
	for (i=0; i<=dim; i++)
	{
		k2[i] = dt*wfn_dot[i];
		wfn3[i] = wfn[i] + k2[i]/2.;
	}

	// Get dPsi(n+k2/2)/dt
	dfdt(dim, wfn3, wfn_dot, dx);

	// Compute approximate update with Euler step

	// Then update actual wfn

	for (i=0; i<=dim; i++)
	{
		k3[i] = dt*wfn_dot[i];
		wfn_np1[i] = wfn[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
		wfn[i] = wfn_np1[i];
	}

	// wfn vector has now been updated!

	// Now free memory associated with temporary vectors
	free(wfn_dot);
	free(wfn2);
	free(wfn3);
	free(wfn_np1);
	free(k1);
	free(k2);
	free(k3);
}
