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
That generates psi_dot
*/

// Global Constants
double pi = 4*atan(1.);
// double alpha = 0.5502;
// Reduced mass in SI units
// double mu = 1.626*pow(10,-27);
// Reduced mass in atomic units - verified for HCl
double mu = 1784 ;
// double mu = 200.;
// double hbar = 1.0546*pow(10,-34);
double hbar = 1.;
double q = 1.;
// Force constant in atomic units - verified for HCl
double k = 0.309;
// double k = 1.;
double alpha = sqrt(mu*k);
// Timestep
double dt = 0.001;
int MAX_TIME;

// Function Prototypes

// 2nd derivative of psi_j/dt
void dfdt(int dim, double complex *psivec, double complex *dpsi, double dx);
	// dim = number of points in wavefxn
	// psivec = array of wavefxn values
	// d2psi = array of dpsi^2/dx^2 wavefxn values
	// dx = increment along x axis

// H hat applied to psi
double Hpsi(double t,int dim, double complex *psivec, double complex *dpsij, double dx, double *x);
	// dim = number of points in wavefxn
	// psivec = array of wavefxn values
	// d2psi = array of dpsi^2/dx^2 wavefxn values
	// dx = increment along x axis
	// x = location along curve
	// t = time
        // Returns energy from <Psi | H_0 | Psi>
// Electronic Field,E(t)
double E(double t);
	// t = current time

// 3rd Order Runge–Kutta :: Advances wavefxn forward in time
double RK3(double t,int dim, double *xvec, double complex *wfn, double dx);
	// dim = number of points in the wavefunction
	// xvec = vector of x-values that the wfn is evaluated at
	// wfn = vector that stores the wfn values at time t_i
	// 		when the function is called, and when the function has been executed,
	//		this vector will hold the wfn values at time t_f = t_i + dt
	// dx = differential step along x between subsequent points in *xvec
	// dt = differential step along time between subsequent times (i.e. t_f - t_i)
	// t = time
        // Returns Energy at current time

// Function to compute dipole moment expectation value
double complex DipoleMoment(int dim, double *xvec, double complex *wfn, double dx); 
void Fourier (double complex *dm, int n);
//void Fourier (const double inreal [], double outreal [], double outimag [], int n);
/* **************************************************************************************************** */

int main()
{
	// Initialize & Define Variables
	int dim = 400;
	double *x;
	double complex dpm, *dm, *wfn, *dpsi, *dpsij;
	double L = 16.;
	double dx = L/dim;
	int i, pidx;
        double time, Et, Ei, SlowPeriod;
        FILE *dpfp;        
 
        // Open Dipole Moment File
        dpfp = fopen("DipoleMoment.txt","w");

	// Arrays via malloc()
	x = (double *)malloc(dim*sizeof(double));
	wfn = (double complex *)malloc(dim*sizeof(double complex));
	dpsi = (double complex *)malloc(dim*sizeof(double complex));
	dpsij = (double complex *)malloc((dim+1)*sizeof(double complex));

	// Psi(x) - initialize as ground state wf
	for (i=0; i<dim; i++)
	{
		x[i] = (i-dim/2)*dx;
		wfn[i] = pow( (alpha/pi) , (1./4))*exp((-alpha/2.)*pow(x[i],2.));

	}
        // Ground state energy - aka slowest frequency in the system
        Ei = 0.5*sqrt(k/mu);
        // Period of slowest oscillation
        SlowPeriod = 1./Ei;
        // Must run for at least 1 period of slowest oscillation, but will run a bit longer
        MAX_TIME = (int)(10*SlowPeriod/dt); 
        printf("  MAX_TIME is %i\n",MAX_TIME);

        // Dipole moment array must be at least MAX_TIME long - will also zero pad for good measure
        dm = (double complex *)malloc((MAX_TIME+10000)*sizeof(double complex));

        fprintf(dpfp, " Time (a.u.) \t Real(mu)  \t Im(mu) \t  E(t=0) (a.u.)  \t E(t) (a.u.)\n");
        pidx=1;
        for (int j=0; j<MAX_TIME; j++){
          time = j*dt;
          Et = RK3(time, dim, x, wfn, dx); 
          dm[j] = DipoleMoment(dim, x, wfn, dx);

          if (j%500==0) {

            fprintf(dpfp, " %12.10f  %12.10e  %12.10e  %12.10e  %12.10e  \n",time,creal(dm[j]),cimag(dm[j]),Ei, Et);
            printf("\n\n#%i\n",pidx);

            for (i=0; i<dim; i++) {
                  printf("%f %e %e %e  %e\n",x[i],creal(wfn[i]),creal( conj(wfn[i])*wfn[i]), (k/2.)*x[i]*x[i], -q*E(time)*x[i]);

	    }
            pidx++;
          }
        }
        // Dipole moment information collected - now zero pad
        for (int j=0; j<10000; j++) {
          dm[MAX_TIME+j] = 0. + 0.*I;
        }


        Fourier(dm, MAX_TIME+10000);

        free(x);
        free(wfn);
        free(dpsi);
        free(dpsij);
        fclose(dpfp);
        free(dm);
        return 0;
}

/* **************************************************************************************************** */
// Just finite difference second derivative of psi -> dpsi
void dfdt(int dim,double complex *psivec,double complex *dpsij,double dx) {
	int j;
	dpsij[0] = 0. + 0.*I;
	dpsij[dim] = 0. + 0.*I;

	for (j=1; j<dim; j++)
	{
		dpsij[j] = (psivec[j+1] - 2.*psivec[j] + psivec[j-1])/(dx*dx);
	}
}

/* **************************************************************************************************** */


double Hpsi(double t,int dim, double complex *psivec, double complex *dpsij, double dx, double *x)
{

   double complex Et;
   // A temporary vector for the second derivative of psivec
   double complex *temp;
   temp = (double complex *)malloc((dim+1)*sizeof(double complex));

   // Call dfdt, which calculates the second derivative of psivec and multiplies it by i/2
   dfdt(dim, psivec, temp, dx);
   // dpsiij is the time derivative!
   // It results from the action of the Hamiltonian on the current wavefunction!


	int j;
        Et = 0. + 0.*I;
	for (j=0; j<dim; j++)
	{
		// (1) is placeholder for d^2/dx^2 aka dpsi[i]
		dpsij[j] = I*temp[j]/(2*mu) - 0.5*I*k*pow(x[j],2)*psivec[j] + I*q*E(t)*x[j]*psivec[j];
                // sum psi^* H psi dx
                Et += conj(psivec[j])*(0.5*k*pow(x[j],2)*psivec[j] - temp[j]/(2*mu))*dx;
                dpsij[j] += I*q*E(t)*x[j]*psivec[j];
	}


   free(temp);
   return creal(Et);
}
/* **************************************************************************************************** */

double E(double x)
{
   double Emax= 0.005;
   double sigma = 200;
   double field;
   double f1 = 0.5*sqrt(k/mu);
   double f2 = sqrt(k/mu);	
   double amplitude;
   if (x<sigma) {

     amplitude = Emax*sin(pi*x/(2*sigma))*sin(pi*x/(2*sigma));
   }
   else {
 
     amplitude = 0.;

   }

   field = amplitude*(sin(f1*x) + sin(f2*x));
   return field;
}

/* **************************************************************************************************** */

double RK3(double t,int dim, double *xvec, double complex *wfn, double dx) {
	double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;
	int i;
        double E0, E1, E2;
	// Temporary arrays for computing derivatives of wfns and approximate updates to wfns
	wfn_dot = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn2 = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn3 = (double complex *)malloc((dim+1)*sizeof(double complex));
	wfn_np1 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k1 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k2 = (double complex *)malloc((dim+1)*sizeof(double complex));
	k3 = (double complex *)malloc((dim+1)*sizeof(double complex));

  //Initialize all (real and imaginary parts) to zero
        //printf("  dx is %f\n",dx);
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
	//dfdt(dim, wfn, wfn_dot, dx);
        E0 = Hpsi(t, dim, wfn, wfn_dot, dx, xvec);
	// Compute approximate wfn update with Euler step
	for (i=0; i<=dim; i++)
	{
		k1[i] = dt*wfn_dot[i];
		wfn2[i] = wfn[i] + k1[i]/2.;
	}

	// Get dPsi(n+k1/2)/dt
	//dfdt(dim, wfn2, wfn_dot, dx);
        E1 = Hpsi((t+dt/2.), dim, wfn2, wfn_dot, dx, xvec);
	
	// Compute approximate wfn update with Euler step
	for (i=0; i<=dim; i++)
	{
		k2[i] = dt*wfn_dot[i];
		wfn3[i] = wfn[i] + k2[i]/2.;
	}

	// Get dPsi(n+k2/2)/dt
	//dfdt(dim, wfn3, wfn_dot, dx);
        E2 = Hpsi((t+dt/2.), dim, wfn3, wfn_dot, dx, xvec);
	
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

        return E0;

}
// Compute dipole moment given current wavefunction and dipole operator:
// mu = <Psi | mu | Psi> = <Psi | -q x | Psi>
double complex DipoleMoment(int dim, double *xvec, double complex *wfn, double dx) {
  double complex dm;
  dm = 0. + 0.*I;

  for (int i=0; i<dim; i++) {

    dm -= conj(wfn[i])*q*xvec[i]*wfn[i]*dx;

  }

  return dm;

}
/***************************************************************************/
void Fourier (double complex *dm, int n){
  FILE *fp;
  fp = fopen("Absorption_Spectrum.txt","w");
  double wmin=0.1*sqrt(k/mu);
  double wmax=4*sqrt(k/mu);
  int maxk = 500;
  double dw = (wmax-wmin)/maxk;
  double time;
 
  for (int k = 0; k <=maxk; k++) {
    double sumreal = 0;
    double sumimag = 0;
    double  w = wmin+k*dw;

    for (int t = 0; t < n; t++){
      time = dt*t;
      double angle = time*w;
      sumreal += creal(dm[t]) * cos(angle) + 0. * sin(angle);
      sumimag  += -creal(dm[t]) * sin(angle) + 0. * cos(angle);
    }
    fprintf(fp," %12.10e  %12.10e\n",w*(21947482.58263),sumreal*sumreal+sumimag*sumimag);
  }
  fclose(fp);
}
