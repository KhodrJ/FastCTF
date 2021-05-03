#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <tuple>
#include <iomanip>

const double TIME_STEP = 3600.00;

using cdouble = std::complex<double>;
using VCD = std::vector<cdouble>;

extern "C"
{
	extern int dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
	extern int zhseqr_(char *JOB, char *COMPZ, int *N, int *ILO, int *IHI, cdouble *H, int *LDH, cdouble *W, cdouble *Z, int *LDZ, cdouble *WORK, int *LWORK, int *INFO);
}

inline double factorial(int n)
{
	double f = 1.0;

	for (int i = 1; i <= n; i++)
		f *= i;

	return f;
}

inline cdouble c(double a, double b)
{
	return cdouble(a,b);
}

inline double max(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;
}

inline double min(double a, double b)
{
	if (a <= b)
		return a;
	else
		return b;
}

