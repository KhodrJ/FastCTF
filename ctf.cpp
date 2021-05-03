/*
 * ==================================================================================
 * 
 *                   ______        _    _____ _______ ______ 
 *                  |  ____|      | |  / ____|__   __|  ____|
 *                  | |__ __ _ ___| |_| |       | |  | |__   
 *                  |  __/ _` / __| __| |       | |  |  __|  
 *                  | | | (_| \__ \ |_| |____   | |  | |     
 *                  |_|  \__,_|___/\__|\_____|  |_|  |_|     
 *                                         
 * 
 *    FastCTF: A Robust Solver for Conduction Transfer Function Coefficients
 * 
 *    Date: 28 March, 2021   
 *    Author: Khodr Jaber
 * 
 *    
 * 
 *    A simple code that computes both conduction transfer functions and response factors based on Taylor series expansions of the Laplace domain solution to the one-dimensional heat equation. The inversion of this approximated solution is performed by generating a polynomial s-transfer function and applying the residue theorem. The z-transfer function is then applied to recover the CTFs and response factors.
 * 
 * ==================================================================================
*/
 
#include "c_poly.h"

Poly GetSeriesT_00(int ord, double L, double alpha)
{
	VCD coeffs = {};
	for (int i = 0; i <= ord; i++)
		coeffs.push_back( pow(L / sqrt(alpha), 2*i) / factorial(2*i) );

	return Poly(coeffs);
}

Poly GetSeriesT_01(int ord, double L, double R, double alpha)
{
	VCD coeffs = {};
	for (int i = 0; i <= ord; i++)
		coeffs.push_back( (R/L) * pow(L, 2*i+1) * pow(1.0 / sqrt(alpha), 2*i) / factorial(2*i+1) );

	return Poly(coeffs);
}

Poly GetSeriesT_10(int ord, double L, double k, double alpha)
{
	VCD coeffs = {0.0};
	for (int i = 0; i < ord; i++)
		coeffs.push_back( k * pow(L, 2*i+1) * pow(1.0 / sqrt(alpha), 2*(i+1)) / factorial(2*i+1) );

	return Poly(coeffs);
}

Poly TFtoTaylor(Poly Num, Poly Den, int length)
{
	VCD coeffs = {};

	coeffs.push_back(Num.coeffs[0] / Den.coeffs[0]);
	for (int i = 1; i < Den.coeffs.size(); i++)
	{
		cdouble sum = 0;
		for (int j = 0; j < i; j++)
			sum += coeffs[j] * Den.coeffs[i-j];

		coeffs.push_back( (Num.coeffs[i] - sum) / Den.coeffs[0] );
	}
	for (int i = 0; i < length; i++)
	{
		cdouble sum = 0;
		for (int j = 0; j < Den.coeffs.size()-1; j++)
			sum += coeffs[j+i+1] * Den.coeffs[Den.coeffs.size()-1-j];

		coeffs.push_back( -sum / Den.coeffs[0] );
	}

	return Poly(coeffs);
}

std::pair<Poly,Poly> ZTransform(cdouble res, cdouble root, int mult)
{
	double A = -TIME_STEP * root.real();
	double B = -TIME_STEP * root.imag();
	double C = 2.0 * res.real();
	double D = 2.0 * res.imag();
	Poly Num(VCD{0});
	Poly Den(VCD{0});
	Poly z(VCD{0.0, -1.0});

	if (root.imag() != 0.0)
	{
		Num = Poly(VCD {0.0, exp(A)*(D*sin(B) - C*cos(B)), C*exp(2*A)});
		Den = Poly(VCD {1.0, -2*exp(A)*cos(B), exp(2*A)});
	}
	if (root.imag() == 0.0)
	{
		if (root.real() != 0.0)
		{
			Num = Poly(VCD {0.0, res.real()*exp(A)});
			Den = Poly(VCD {-1.0, exp(A)});
		}
		else
			std::cout << "Oops..." << std::endl;
	}

	return std::make_pair(Num, Den);
}

std::tuple<Poly, Poly, Poly> GetZTF(Poly Num, Poly Den, double U, int N_k)
{
	// First:	Numerator of z-transfer function.
	// Second:	Denominator of z-transfer function.
	// Third:	Response factors.
    
	std::pair<VCD, std::vector<int>> RaR = RootsPoly(Den);
	VCD roots = RaR.first;
	std::vector<int> reps = RaR.second;
    
	int N = roots.size();
	Poly ZNum(VCD{0.0});
	Poly ZDen(VCD{1.0});
	std::vector<Poly> M_bar;
	for (int i = 0; i < N; i++)
		M_bar.push_back(Poly(VCD {1.0}));

	// Get residues of linear and constant term.
	cdouble alpha_1 = 0.0;
	cdouble alpha_2 = TIME_STEP * U;

	// Compute transfer function.
	for (int i = 0; i < N; i++)
	{
		cdouble res = Residue(roots[i], 1, 1, Num, Den);
		alpha_1 -= res;
		if (fabs(res.imag()) > 0)
			alpha_1 -= c(res.real(), -res.imag());

		std::pair<Poly, Poly> ztf = ZTransform(res, roots[i], 1);

		M_bar[i] = M_bar[i]*ztf.first;
		for (int j = 0; j < N; j++)
			if (j != i)
				M_bar[j] = M_bar[j]*ztf.second;
		ZDen = ZDen * ztf.second;
	}

	// Assemble.
	Poly ZM_bar = Poly(VCD {0.0});
	for (int i = 0; i < N; i++)
		ZM_bar = ZM_bar + M_bar[i];
	alpha_1 = -ZM_bar.coeffs[ZM_bar.Degree()] / ZDen.coeffs[ZDen.Degree()];

	Poly M_bar_poly = ZM_bar;
	Poly M_bar_z_poly = ZM_bar * Poly(VCD {0.0, 1.0});
	Poly M_bar_zm_poly = ZM_bar; M_bar_zm_poly.coeffs.erase(M_bar_zm_poly.coeffs.begin());

	ZNum = ZDen * Poly(VCD {alpha_2-alpha_1, alpha_1}) + (M_bar_z_poly + Poly(VCD {-2.0})*ZM_bar + M_bar_zm_poly);
	if (ZNum.Degree() > ZDen.Degree())
		ZNum.coeffs.erase(ZNum.coeffs.end()-1);

	cdouble normal = 1.0 / ZDen.coeffs[ZDen.Degree()];
	ZNum = ZNum * Poly(VCD {normal / TIME_STEP});
	ZDen = ZDen * Poly(VCD {normal});

	// Compute response factors.
	Poly Yk = Poly(VCD {0});
	Poly Ykin = Poly(VCD {0});
	Poly ZNum_F = ZNum; ZNum_F.Flip();
	Poly ZDen_F = ZDen; ZDen_F.Flip();
	Yk = TFtoTaylor(ZNum_F, ZDen_F, N_k);

	return std::make_tuple(ZNum, ZDen, Yk);
}

std::pair<Poly, Poly> Pade(Poly taylor, int m, int n)
{
	VCD vec_num = {};
	VCD vec_den = {1.0};

	int		N = m+n+1;
	int		NRHS = 1;
	double		*A = new double[N*N]{0.0};
	int		LDA = N;
	int		*IPIV = new int[N];
	double		*B = new double[N];
	int		LDB = N;
	int		INFO;

	for (int i = 0; i < n+1; i++)
	{
		A[i+N*i] = -1.0;
		B[i] = -taylor.coeffs[i].real();
		for (int j = n+1; j < (n+1)+i; j++)
			A[i+N*j] = taylor.coeffs[i - (j-(n+1)+1)].real();
	}
	for (int i = n+1; i < N; i++)
	{
		B[i] = -taylor.coeffs[m+1 + (i-(n+1))].real();
		for (int j = n+1; j < N; j++)
			A[i+N*j] = taylor.coeffs[(i-(n+1)) + (m+1)-1 - (j-(n+1))].real();
	}

	dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);

	for (int i = 0; i < n+1; i++)
		vec_num.push_back(B[i]);
	for (int i = n+1; i < N; i++)
		vec_den.push_back(B[i]);

	delete[] A, B, IPIV;

	return std::make_pair(Poly(vec_num), Poly(vec_den));
}

Poly TaylorRecover(Poly Num, Poly Den, Poly Q)
{
    Poly P(VCD {0});
    
    Poly NQ = Num*Q;
    for (int i = 0; i < Q.coeffs.size(); i++)
    {
        cdouble sum = 0.0;
        for (int j = 0; j < i; j++)
            sum += P.coeffs[j] * Den.coeffs[i-j];
        
        P.coeffs.push_back( (NQ.coeffs[i] - sum) / Den.coeffs[0] );
    }
    
    return P;
}

int main(int argc, char *argv[])
{
	int			i, j, N, N_layers, ord, p_ord, N_k;
	double			R_o, R_i, **layer_dat;
	std::ifstream		input;
	std::ofstream		output_x, output_y, output_z;

    
    
	// Read input.
	input.open("input.txt");
	input >> N_layers >> R_o >> R_i;
	layer_dat = new double*[N_layers];
	for (i = 0; i < N_layers; i++)
	{
		layer_dat[i] = new double[5];
		for (j = 0; j < 5; j++)
			input >> layer_dat[i][j];
	}
	input.close();
	for (i = 0; i < N_layers; i++)
	{
		if (layer_dat[i][0] != 0)
			layer_dat[i][4] = 0.001*layer_dat[i][0] / layer_dat[i][1];
	}
    
    
    
   	// Compute U value.
	double U = 0.0;
	for (i = 0; i < N_layers; i++)
		U += layer_dat[i][4];
	U = 1.0/(U+R_i+R_o);

    
    
	// Parameters.
	ord = 20;
	p_ord = 5;
    	N_k = 144;
	N = 50;


	// STEP ONE: Compute transmission matrix polynomials.

	Poly T_00 = Poly(VCD {1.0});
	Poly T_01 = Poly(VCD {R_o});
	Poly T_10 = Poly(VCD {0.0});
	Poly T_11 = Poly(VCD {1.0});
	Poly T_00_t = Poly(VCD {0});
	Poly T_01_t = Poly(VCD {0});
	Poly T_10_t = Poly(VCD {0});
	Poly T_11_t = Poly(VCD {0});

	for (i = 0; i < N_layers+1; i++)
	{
		if (i < N_layers)
		{
			double L_i = layer_dat[i][0]*0.001;
			double k_i = layer_dat[i][1];
			double rho_i = layer_dat[i][2];
			double C_pi = layer_dat[i][3];
			double R_i_ = layer_dat[i][4];
			double alpha_i = k_i / (rho_i*C_pi);

			if (L_i != 0)
			{
				T_00_t = GetSeriesT_00(ord, L_i, alpha_i);
				T_01_t = GetSeriesT_01(ord, L_i, R_i_, alpha_i);
				T_10_t = GetSeriesT_10(ord, L_i, k_i, alpha_i);
				T_11_t = GetSeriesT_00(ord, L_i, alpha_i);
			}
			else
			{
				T_00_t = Poly(VCD {1.0});
				T_01_t = Poly(VCD {R_i_});
				T_10_t = Poly(VCD {0.0});
				T_11_t = Poly(VCD {1.0});
			}
		}
		else
		{	
			T_00_t = Poly(VCD {1.0});
			T_01_t = Poly(VCD {R_i});
			T_10_t = Poly(VCD {0.0});
			T_11_t = Poly(VCD {1.0});
		}
		Poly T_00_tmp = T_00_t*T_00 + T_01_t*T_10;
		Poly T_01_tmp = T_00_t*T_01 + T_01_t*T_11;
		Poly T_10_tmp = T_10_t*T_00 + T_11_t*T_10;
		Poly T_11_tmp = T_10_t*T_01 + T_11_t*T_11;

		T_00 = T_00_tmp;
		T_01 = T_01_tmp;
		T_10 = T_10_tmp;
		T_11 = T_11_tmp;
	}
	
	
	
	// STEP TWO: Compute inversion and assemble transfer functions.
    
		// Cross.
	Poly NumS_Y(VCD {1.0});
	Poly DenS_Y = T_01;
	std::pair<Poly, Poly> PadeS_Y = Pade(DenS_Y, p_ord, p_ord);
	std::tuple<Poly, Poly, Poly> PadeZ_Y = GetZTF(PadeS_Y.second, Poly(VCD {0, 0, 1}) * PadeS_Y.first, U, N_k);
	Poly Yk_Y = std::get<2>(PadeZ_Y);
    std::cout << "Verification of Y: " << std::get<0>(PadeZ_Y).SumCoeffs() / std::get<1>(PadeZ_Y).SumCoeffs() << std::endl;
    
    		// External.
	Poly NumS_X(VCD {0}), DenS_X(VCD {0});
	for (int i = 0; i < p_ord+1; i++)
	{
		//NumS_X.coeffs.push_back(T_00.coeffs[i]);
		//DenS_X.coeffs.push_back(T_01.coeffs[i]);
	}
	DenS_X = PadeS_Y.first;
	NumS_X = TaylorRecover(T_00, T_01, DenS_X);
	std::tuple<Poly, Poly, Poly> PadeZ_X = GetZTF(NumS_X, Poly(VCD {0, 0, 1}) * DenS_X, U, N_k);
	Poly Yk_X = std::get<2>(PadeZ_X);
	std::cout << "Verification of X: " << std::get<0>(PadeZ_X).SumCoeffs() / std::get<1>(PadeZ_X).SumCoeffs() << std::endl;
    
		//Internal.
	Poly NumS_Z(VCD {0}), DenS_Z(VCD {0});
	for (int i = 0; i < p_ord+1; i++)
	{
		//NumS_Z.coeffs.push_back(T_11.coeffs[i]);
		//DenS_Z.coeffs.push_back(T_01.coeffs[i]);
	}
    DenS_Z = PadeS_Y.first;
	NumS_Z = TaylorRecover(T_11, T_01, DenS_Z);
	std::tuple<Poly, Poly, Poly> PadeZ_Z = GetZTF(NumS_Z, Poly(VCD {0, 0, 1}) * DenS_Z, U, N_k);
	Poly Yk_Z = std::get<2>(PadeZ_Z);
	std::cout << "Verification of Z: " << std::get<0>(PadeZ_Z).SumCoeffs() / std::get<1>(PadeZ_Z).SumCoeffs() << std::endl;
	std::cout << std::endl;


    
	// PRINT OUTPUT.
	std::cout << "In descending order...\n" << std::endl;
	std::cout << std::setprecision(10) << "a_k: " << std::endl;
	std::get<0>(PadeZ_X).Print(p_ord, 1); std::cout << "Sum of a_k: " << std::get<0>(PadeZ_X).SumCoeffs() << 	std::endl;
	std::cout << "\nb_k: " << std::endl;
	std::get<0>(PadeZ_Y).Print(p_ord, 1);  std::cout << "Sum of b_k: " << std::get<0>(PadeZ_Y).SumCoeffs() << 	std::endl;
	std::cout << "\nc_k: " << std::endl;
	std::get<0>(PadeZ_Z).Print(p_ord, 1); std::cout << "Sum of c_k: " << std::get<0>(PadeZ_Z).SumCoeffs() << 	std::endl;
	std::cout << "\nd_k: " << std::endl;
	std::get<1>(PadeZ_Z).Print(p_ord, 1);  std::cout << "Sum of d_k: " << std::get<1>(PadeZ_Y).SumCoeffs() << std::endl << std::endl;

		// Uncomment to print response factors.
	//std::get<2>(PadeZ_X).Print(std::get<2>(PadeZ_X).Degree(), 1); std::cout << std::endl;
	//std::get<2>(PadeZ_Y).Print(std::get<2>(PadeZ_Y).Degree(), 1); std::cout << std::endl;
	//std::get<2>(PadeZ_Z).Print(std::get<2>(PadeZ_Z).Degree(), 1); dstd::cout << std::endl;

	std::cout << "U: " << U << std::endl;

	return 0;
}
