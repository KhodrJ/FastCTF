#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <tuple>

// const double TIME_STEP = 60.00;
const double TIME_STEP = 3600.00;

using cdouble = std::complex<double>;
using VCD = std::vector<cdouble>;

extern "C"
{
	extern int dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
	extern int zhseqr_(char *JOB, char *COMPZ, int *N, int *ILO, int *IHI, cdouble *H, int *LDH, cdouble *W, cdouble *Z, int *LDZ, cdouble *WORK, int *LWORK, int *INFO);
	extern int zgetrf_(int *M, int *N, cdouble *A, int *LDA, int *IPIV, int *INFO);
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

inline cdouble Determinant(int N, cdouble *A)
{
	cdouble det = 1.0;
	int *IPIV = new int[N];
	int INFO;

	zgetrf_(&N, &N, A, &N, IPIV, &INFO);

	/*
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			std::cout << A[i+j*N] << " ";
		std::cout << std::endl;
	}
	*/

	for (int i = 0; i < N; i++)
	{
		det *= A[i*(N+1)];
		if (i+1 != IPIV[i])
			det = -det;
	}

	delete[] IPIV;

	return det;
}
