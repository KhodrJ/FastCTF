#include "cppspec.h"

class Poly
{
	public:
		// coeffs	Vector of coefficients in ascending order of x^N.

		Poly(double z)
		{
			coeffs.clear();
			coeffs.resize(1);
			coeffs[0] = z;

			coeffs.resize(Degree()+1);

		}
		Poly(int m, double z)
		{
			coeffs.clear();
			coeffs.resize(m);

			for (int i = 0; i < m; i++)
				coeffs[i] = z;

			coeffs.resize(Degree()+1);
		}
		Poly(VCD coeffs_)
		{
			coeffs.clear();
			coeffs.resize(coeffs_.size());

			for (int i = 0; i < coeffs_.size(); i++)
				coeffs[i] = coeffs_[i];

			coeffs.resize(Degree()+1);
		}


		VCD coeffs;
		static int DividePolyT(Poly *q, Poly *r, Poly n, Poly d);
		friend std::ostream& operator<< (std::ostream &os, const Poly &poly);

		int Degree()
		{
			int d = coeffs.size()-1;
			bool sent = true;
			while (sent)
			{
				if (coeffs[d] == cdouble(0,0))
					d -= 1;
				else
					sent = false;
			}

			return d;
		}
		int Normalize(int piv)
		{
			cdouble normal = coeffs[piv];
			for (int i = 0; i < coeffs.size(); i++)
				coeffs[i] /= normal;
			
			return 0;
		}

		int Copy(Poly poly)
		{
			int m = poly.coeffs.size();
			coeffs.clear();
			coeffs.resize(m);

			for (int i = 0; i < m; i++)
				coeffs[i] = poly.coeffs[i];

			return 0;
		}

		int Diff(int k)
		{
			if (k > 0)
			{
				if (coeffs.size() != 0)
				{
					for (int i = 1; i < coeffs.size(); i++)
						coeffs[i] *= (double)i;

					coeffs.erase(coeffs.begin());
					Diff(k-1);
				}
			}

			return 0;
		}
		int Print(int ord, int var)
		{
			for (int i = 0; i <= ord; i++)
			{
				if (var == 0)
				{
					if (i < ord)
						std::cout << coeffs[i].real() << "x^" << i << " + ";
					else
						std::cout << coeffs[i].real() << "x^" << i << std::endl;
				}
				else
				{
					if (i < ord)
						std::cout << coeffs[i].real() << ",";
					else
						std::cout << coeffs[i].real() << std::endl;

				}
			}

			return 0;
		}
		double SumCoeffs()
        {
            double sum = 0.0;
            for (int i = 0; i < coeffs.size(); i++)
                sum += coeffs[i].real();
            
            return sum;
        }

		cdouble Eval(cdouble z)
		{
			cdouble sum = 0.0;
			if (z == c(0,0))
				return coeffs[0];
			else
			{
				for (int i = 0; i < coeffs.size(); i++)
					sum += coeffs[i]*pow(z, i);
			}

			return sum;
		}
};

inline std::ostream& operator<< (std::ostream &os, const Poly &poly)
{
	os << "(";
	for (int i = 0; i < poly.coeffs.size(); i++)
		os << poly.coeffs[i] << " ";
	os << ")";

	return os;
}

inline Poly operator+ (const Poly &poly1, const Poly &poly2)
{
	if (poly1.coeffs.size() == 0 && poly2.coeffs.size() == 0)
		return Poly(VCD {0});

	VCD new_coeffs;
	int m = poly1.coeffs.size();
	int n = poly2.coeffs.size();
	new_coeffs.resize(max(m,n));

	for (int i = 0; i < new_coeffs.size(); i++)
	{
		if (i < m)
			new_coeffs[i] += poly1.coeffs[i];
		if (i < n)
			new_coeffs[i] += poly2.coeffs[i];
	}
	
	return Poly(new_coeffs);
}

inline Poly operator- (const Poly &poly1, const Poly &poly2)
{
	if (poly1.coeffs.size() == 0 && poly2.coeffs.size() == 0)
		return Poly(VCD {0});

	VCD new_coeffs;
	int m = poly1.coeffs.size();
	int n = poly2.coeffs.size();
	new_coeffs.resize(max(m,n));

	for (int i = 0; i < m; i++)
		new_coeffs[i] = poly1.coeffs[i];

	for (int i = 0; i < new_coeffs.size(); i++)
	{
		if (i < n)
			new_coeffs[i] -= poly2.coeffs[i];
	}
	
	return Poly(new_coeffs);
}


inline Poly operator* (const Poly &poly1, const Poly &poly2)
{
	if (poly1.coeffs.size() == 0 || poly2.coeffs.size() == 0)
		return Poly(VCD {0});
	else
	{
		VCD new_coeffs;
		int m = poly1.coeffs.size();
		int n = poly2.coeffs.size();
		new_coeffs.resize(m+n-1);

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				new_coeffs[i+j] += poly1.coeffs[i]*poly2.coeffs[j];

		return Poly(new_coeffs);
	}
}

inline std::pair<Poly,Poly> DividePoly(Poly n, Poly d) 
{
    // Currently unused.
    
	Poly q_(1,0);
	Poly r_(1,0);

        if (n.Degree() < d.Degree()) 
        { 
                q_.coeffs.clear(); 
                r_.Copy(n); 
 
                return std::make_pair(Poly(1,0), n); 
        } 
        else 
        { 
                q_.coeffs.clear(); 
                r_.Copy(n); 

		int S = 10000;

                while (S >= 0 && r_.coeffs.size() != 0 && r_.Degree() >= d.Degree()) 
                { 
                        int d_r = r_.Degree(); 
                        int d_d = d.Degree(); 
                        VCD coeffs_t; 
                        coeffs_t.resize(d_r - d_d+1);
		       	coeffs_t[d_r - d_d] = r_.coeffs[d_r] / d.coeffs[d_d];
                        Poly t = Poly(coeffs_t); 
                       
                        q_ = q_ + t; 
                        r_ = r_ - (t*d);

			S--;
                }

		if (S < 0)
			std::cout << "Something went wrong during polynomial division..." << std::endl;
        } 
 
        return std::make_pair(q_, r_); 
}

inline cdouble Residue(cdouble root, int j, int mult, Poly P, Poly Q)
{
	// Computes simple residues using rational function polynomials.
    // Note: j starts at 1.

	VCD roots;
	roots.resize(mult);
	Poly Q_prime = Q; Q_prime.Diff(1);

	return P.Eval(root) / Q_prime.Eval(root);
}

inline std::pair<VCD,std::vector<int>> RootsPoly(Poly poly)
{
	int zero_counter = 0;
	if (poly.coeffs[0] == c(0,0))
	{
		int i = 0;
		while (poly.coeffs[i] == c(0,0))
		{
			zero_counter++;
			i++;
		}
	}
	for (int i = 0; i < zero_counter; i++)
		poly.coeffs.erase(poly.coeffs.begin());
	
	char 			JOB = 'E';
	char 			COMPZ = 'N';
	int 			N = poly.coeffs.size()-1;
	int			ILO = 1;
	int			IHI = N;
	cdouble			*H = new cdouble[N*N]{c(0,0)};
	cdouble 		*W = new cdouble[N]{c(0,0)};
	int			LDH = N;
	int 			LDZ = 1;
	cdouble  		*WORK = new cdouble[N]{c(0,0)};
	int			LWORK = N;
	int			INFO;
	VCD			roots;
	std::vector<int>	reps;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			int rc = i*N+j;
			if (j == 0)
				H[rc] = -(poly.coeffs[i+1]/poly.coeffs[0]);
			if (j == i+1)
				H[rc] = 1.0;
		}
	}

	zhseqr_(&JOB, &COMPZ, &N, &ILO, &IHI, H, &LDH, W, 0, &LDZ, WORK, &LWORK, &INFO);

	// Clean roots, ignore small relativality small imaginary components.
	for (int i = 0; i < N; i++)
	{
		W[i] = 1.0 / W[i];
		if ( std::abs(W[i].real() / W[i].imag()) > 1e6 )
			W[i] = c(W[i].real(), 0.0);
		if ( std::abs(W[i].imag() / W[i].real()) > 1e6 )
			W[i] = c(0.0, W[i].imag());
	}

	int counter = 1;
	for (int i = 0; i < N; i++)
	{
		if (W[i].imag() >= 0)
		{
			roots.push_back(W[i]);
			reps.push_back(1);
		}
	} 

	delete[] H, W, WORK;

	return std::make_pair(roots, reps);
}
