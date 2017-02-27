# include <cstdlib>
# include <cmath>
# include <iostream>
# include <sstream>
# include <iomanip>
# include <ctime>
# include <cstring>
# includte_polynomial.hppm <stdio.h>
using namespace std;

# include "hermite_polynomial.hpp"

int main ( );
void hermite_polynomial_test01 ( );
string i4_to_string ( int i4 );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HERMITE_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    HERMITE_POLYNOMIAL_PRB tests the HERMITE_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int e;
  int p;

  timestamp ( );
  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the HERMITE_POLYNOMIAL library.\n";

  hermite_polynomial_test01 ( );

  double *fx2_vec;
  double h1, h2;
  double x_vec[10];
  for (int i=0; i<10; i++) {

    x_vec[i] = 0.1*i;
  }

  // going to get a vector of length 2*10, which 
  // contains the the first 2 Hermite polynomials
  // evaluated at the 10 x values contained in x_vec
  fx2_vec = h_polynomial_value ( 10, 3, x_vec );
    
  for (int i=0; i<10; i++) {

    // Get the first Hermite polynomial
    h1 = fx2_vec[i];
    h2 = fx2_vec[3*10+i];
    printf("  %f  %f  %f\n",x_vec[i], h1, h2);
  }


//
//  Terminate.
//
  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void hermite_polynomial_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLYNOMIAL_TEST01 tests H_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int n;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "HERMITE_POLYNOMIAL_TEST01:\n";
  cout << "  H_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the physicist's Hermite polynomials.\n";
  cout << "  H_POLYNOMIAL_VALUE evaluates the polynomial.\n";
  cout << "\n";
  cout << "                        Tabulated                 Computed\n";
  cout << "     N        X           H(N,X)                    H(N,X)                     Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    // Look-up table
    h_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    // function
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double H_POLYNOMIAL_VALUE[M*(N+1)], the values of the first
//    N+1 Hermite polynomials at the evaluation points.

    fx2_vec = h_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    delete [] fx2_vec;

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setw(8) << e << "\n";
  }
  return;
}
//****************************************************************************80

//****************************************************************************80

string i4_to_string ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}

