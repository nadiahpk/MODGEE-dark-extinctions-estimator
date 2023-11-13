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

const int T_MAX = 215;
double g_log_p_M_and_extant[T_MAX];

// [[Rcpp::export]]
double compute_L_i_Cpp(const NumericVector &M_i,
	const NumericVector &x_i, const NumericVector &m)
  {
	int t;
	int T = M_i.size();
	if (T>T_MAX) return 0;
	
	int last_t_1 = -1;
	double *log_p_M_and_extant = &g_log_p_M_and_extant[0];

	const double *M_i_pointer = &M_i[0];
	const double *x_i_pointer = &x_i[0];
	const double *m_pointer = &m[0];
	double *log_p_M_and_extant_t = log_p_M_and_extant;

	for ( t=0; t<T; t++, M_i_pointer++, x_i_pointer++, m_pointer++, log_p_M_and_extant_t++ )
		{
		if ( *M_i_pointer == 1 )
			{
			last_t_1 = t;
			(*log_p_M_and_extant_t) = log((1-exp(-(*x_i_pointer)))*exp(-(*m_pointer)));
			}
		else
			{
			(*log_p_M_and_extant_t) = -(*x_i_pointer)-(*m_pointer);
			}
		}
	
	double temp = 0, L_i = 0;
	m_pointer = &m[0];
	log_p_M_and_extant_t = log_p_M_and_extant;

	for ( t=0; t<T; t++, m_pointer++, log_p_M_and_extant_t++ )
		{
		if ( t > last_t_1 ) L_i += (1-exp(-*m_pointer))*exp(temp);
		temp += (*log_p_M_and_extant_t);
		}

	L_i += exp(temp);

  return(L_i);
  }
