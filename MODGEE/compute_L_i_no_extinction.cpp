#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <time.h>
#include <Rcpp.h>
#include <math.h>
#include <float.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double compute_L_i_no_extinction_Cpp(const NumericVector &M_i, const NumericVector &x_i)
  {
	int T = M_i.size();  
	const double *M_i_pointer = &M_i[0];
	const double *x_i_pointer = &x_i[0];

	double LL_i = 0;

	for ( int t=0; t<T; t++, M_i_pointer++, x_i_pointer++ )
		{
		if ( *M_i_pointer == 1 )
			LL_i += log((1-exp(-(*x_i_pointer))));
		else
			LL_i += -(*x_i_pointer);
		}
	
  return(exp(LL_i));
  }
