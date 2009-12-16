#pragma once
#include <veclib/cblas.h>
#include <veclib/clapack.h>


template<int n> class fixedVector;

///////////////////////////////////////////////////////////////////////////////////////
//!
//! fixedMatrix
//! fixed n x n - Matrix. 
//! supports 
//! +, -, *, +=, -=, = with other fixedMatrix<n>
//! * with fixedVector<n>, * with double
//! =, == double (== id*d)
//! getAt(int row, int column), operator() (row, column)
//! ostream operator <<
//! norm(), norm2(), print(), p(){print();}

template<int n>
class fixedMatrix
{
public:
	fixedMatrix() {}
	fixedMatrix(double d)
	{
		operator = (d);
	}
	double &getAt(int r, int c)
	{
		return values[c + r*n];
	}
	double getAt(int r, int c) const
	{
		return values[c + r*n];
	}
	double operator ()(int r, int c) const
	{
		return getAt(r, c);
	}
	double &operator ()(int r, int c)
	{
		return getAt(r, c);
	}
	
	friend ostream &operator << (ostream &out, const fixedMatrix<n> &s)
	{
		//out << "fixed mat " << rows << "x" << cols << " : "<< endl;
		out << "[ ";
		for(int r=0; r < n; r++)
		{
			for(int c=0; c< n; c++)
				out << s(r, c) << " ";			
			if(r != n-1) out	<< "| ";
		}
		out << "] ";
		return out;
	}
	

	double operator = (double d)
	{
		memset(values, 0, sizeof(double)*n*n);
		if(d == 0.0) return d;
		for(int i=0; i<n; i++)
			getAt(i, i) = d;
		return d;
	}
	
	void operator = (const fixedMatrix<n> &other)
	{
		memcpy(values, other.values, sizeof(double)*n*n);	
	}
	
	fixedMatrix<n> operator + (const fixedMatrix<n> &other ) const
	{
		fixedMatrix<n> erg;
		for(int i=0; i<n*n; i++)
			erg.values[i] = values[i] + other.values[i];
		return erg;
	}
	
	void operator += (const fixedMatrix<n> &other )
	{
		for(int i=0; i<n*n; i++)
			values[i] += other.values[i];
	}	
	
	fixedMatrix<n> operator - (const fixedMatrix<n> &other ) const
	{
		fixedMatrix<n> erg;
		for(int i=0; i<n*n; i++)
			erg.values[i] = values[i] - other.values[i];
		return erg;
	}

	void operator -= (const fixedMatrix<n> &other )
	{
		for(int i=0; i<n*n; i++)
			values[i] -= other.values[i];
	}	
	
	fixedMatrix<n> operator * (const fixedMatrix<n> &other ) const
	{
		fixedMatrix<n> erg;
		for(int r=0; r<n; r++)
			for(int c=0; c<n; c++)
			{
				double s = 0;
				for(int i=0; i<n; i++)
					s += getAt(r, i) * getAt(i, c);
				erg(r, c) = s;
			}
		return erg;
	}
	
	fixedMatrix<n> operator * (double d) const
	{
		fixedMatrix<n> erg;
		for(int r=0; r<n; r++)
			for(int c=0; c<n; c++)
				erg(r, c) = getAt(r, c)*d;
		return erg;
	}
	
	fixedVector<n> operator * (const fixedVector<n> &vec ) const
	{
		fixedVector<n> erg;
		for(int r=0; r<n; r++)
		{
			double s = 0;
			for(int c=0; c<n; c++)
				s += getAt(r, c) * vec(c);
			erg(r) = s;
		}			
		return erg;
	}
	
	bool operator == (double d) const
	{
		for(int i=0; i<n*n; i++)
		{
			if(values[i] != d) return false;
		}
		return true;
	}
	bool operator != (double d) const
	{
		return ! operator == (d);
	}
	double norm() const
	{
		double s = 0;
		for(int i=0; i<n*n; i++)
			s += values[i]*values[i];
		return sqrt(s);			
	}
	double norm2() const
	{
		double s = 0;
		for(int i=0; i<n*n; i++)
			s += values[i]*values[i];
		return s;
	}
	
	//inline void setAsInverseOf(const fixedMatrix<n> &mat );
	
	void p();
	void print() { p(); }
	
private:
	double values[n*n];
};

///////////////////////////////////////////////////////////////////////////////////////
//!
//! smallInverse<int n>
//! A class to hold a inverse of a smallMatrix<n>
//! implemented with LAPACKs LU-Decomposition dgetrf
//! (uses double[n*n] for LU and interchange[n] for pivoting
//! functions:
//! setAsInverseOf(const fixedMatrix<n> &mat) : init as inverse of mat
//! fixedVector<n> * smallInverse<n> = smallInverse<n> * fixedVector<n>
//! = A^{-1} b
template<int n>
class smallInverse
{
public:
	double densemat[n*n];
	__CLPK_integer interchange[n];
	
	void setAsInverseOf(const fixedMatrix<n> &mat)
	{
		for(int r=0; r<n; r++)
			for(int c=0; c<n; c++)
				densemat[c + r*n] = mat(r, c);
		__CLPK_integer info = 0;
		__CLPK_integer dim = n;
		dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
		ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}
	
	fixedVector<n> operator * (const fixedVector<n> &vec) const
	{
		fixedVector<n> erg =vec;
		char trans ='N';
		__CLPK_integer nrhs = 1;
		__CLPK_integer dim = n;
		__CLPK_integer info = 0;
		dgetrs_(&trans, &dim, &nrhs, const_cast<double*> (densemat), &dim, 
					const_cast<__CLPK_integer*> (interchange), erg.values, &dim, &info);	
		return erg;
	}
};

template<int n>
fixedVector<n> operator * (const fixedVector<n> &vec, const smallInverse<n> &mat)
{
	return mat * vec;
}


template<int n>
class fixedVector
{
public:
	double &getAt(int i)
	{
		return values[i];
	}
	double getAt(int i) const
	{
		return values[i];
	}
	double &operator ()(int i)
	{
		return values[i];
	}
	double operator () (int i) const
	{
		return values[i];
	}
	

	friend ostream &operator << (ostream &out, const fixedVector<n> &v)
	{
		out << "( ";
		for(int i=0; i<n; i++)
			out << v(i) << " ";			
		out << ") ";
		return out;
	}
	
	double operator = (double d)
	{
		for(int i=0; i<n; i++)
			getAt(i) = d;
		return d;
	}
	
	void operator = (const fixedVector<n> &other)
	{
		memcpy(values, other.values, sizeof(double)*n);	
	}
	
	fixedVector<n> operator + (const fixedVector<n> &other ) const
	{
		fixedVector<n> erg;
		for(int i=0; i<n; i++)
			erg.values[i] = values[i] + other.values[i];
		return erg;
	}
	
	void operator += (const fixedVector<n> &other )
	{
		for(int i=0; i<n; i++)
			values[i] += other.values[i];		
	}
	
	fixedVector<n> operator - (const fixedVector<n> &other ) const
	{
		fixedVector<n> erg;
		for(int i=0; i<n; i++)
			erg.values[i] = values[i] - other.values[i];
		return erg;
	}

	void operator -= (const fixedVector<n> &other )
	{
		for(int i=0; i<n; i++)
			values[i] -= other.values[i];		
	}
	
	double operator * (const fixedVector<n> &other ) const
	{
		double s=0;
		for(int i=0; i<n; i++)
			s += values[i] * other.values[i];
		return s;
	}
	
	fixedVector<n> operator * (double alpha) const
	{
		fixedVector<n> erg;
		for(int i=0; i<n; i++)
			erg(i) = getAt(i) * alpha;
		return erg;
	}

	fixedVector<n> operator * (const fixedMatrix<n> &mat )
	{
		return (mat * (*this));
	}
	
	
	
	inline void operator /= (const fixedMatrix<n> &mat);
	
	fixedVector<n> operator / (const fixedMatrix<n> &mat )
	{
		fixedVector<n> erg = *this;
		erg /= mat;
		return erg;
	}
	
	double norm2() const
	{
		double s =0;
		for(int i=0; i<n; i++) s += getAt(i)*getAt(i);
		return s;
	}
	
	void p();
	void print() { p(); }
	
	double values[n];
};

template<int n>
inline fixedVector<n> operator * (double alpha, const fixedVector<n> &v) 
{
	return v * alpha;
}


template<>
inline void fixedVector<1>::operator /= (const fixedMatrix<1> &mat )
{
	getAt(0) = getAt(0) / mat(0, 0);
}
template<>
inline void fixedVector<2>::operator /= (const fixedMatrix<2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	ASSERT(invD != 0.0);

	double a = invD*(getAt(0) * mat(1,1) - getAt(1) * mat(0,1));
	double b = invD*(getAt(1) * mat(0,0) - getAt(0) * mat(1,0));
	getAt(0) = a;
	getAt(1) = b;
}

template<int n>
inline void fixedVector<n>::operator /= (const fixedMatrix<n> &mat )
{
	double densemat[n*n];	
	
	for(int r=0; r<n; r++)
		for(int c=0; c<n; c++)
			densemat[c + r*n] = mat(r, c);
	
	__CLPK_integer interchange[n];
	
	__CLPK_integer info = 0;
	__CLPK_integer dim = n;
	dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
	ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
		
	char trans ='N';
	__CLPK_integer nrhs = 1;
	dgetrs_(&trans, &dim, &nrhs, densemat, &dim, interchange, values, &dim, &info);	
}

/*template<>
inline void fixedMatrix<2>::setAsInverseOf(const fixedMatrix<2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	ASSERT(invD != 0.0);
	getAt(0,0) = invD * mat(1,1);
	getAt(0,1) = -invD * mat(0,1);
	getAt(1,0) = -invD * mat(1,0);
	getAt(1,1) = invD * mat(0,0);
}*/

/*
a b   1 0
c d   0 1

c    bc/a  ,   c/a      0
c    d         0        1

c    bc/a  ,     c/a      0
0    d-bc/a      -c/a     1


1    b/a  ,     1/a         0
0    1          -c/(ad-bc)  a/(ad-bc)

1    0    ,     d/(ad-bc)      -b/(ad-bc)
0    1          -c/(ad-bc)      a/(ad-bc)
*/




template<int n>
fixedVector<n> operator * (double alpha, fixedVector<n> &vec)
{
	return vec*alpha;
}
/*
template<int n>
inline double mnorm2(const  fixedVector<n> &v)
{
	return v.norm2();
}*/

template<typename TYPE>
inline double mnorm2(const TYPE &v)
{
	return v.norm2();
}

template<typename TYPE>
inline double mnorm(const TYPE &v)
{
	return sqrt(v.norm2());
}

template <>
inline double mnorm(const double &a)
{
	return dabs(a);
}

template <>
inline double mnorm2(const double &a)
{
	return a*a;
}

template<typename M> inline double getAt(const M &m, int i)
{
	return m(i);
}
template<> inline double getAt(const double &m, int i) { return m; }

template<typename M> inline double getAt(const M &m, int i, int j)
{
	return m(i, j);
}
template<> inline double getAt(const double &m, int i, int j) { return m; }


template<typename M> inline double setAt(M &m, int i, double a)
{
	return m(i) = a;
}
template<> inline double setAt(double &m, int i, double a) { m = a; return a; }

template<typename M> inline double setAt(M &m, int i, int j, double a)
{
	return m(i, j) = a;
}
template<> inline double setAt(double &m, int i, int j, double a) { m = a; return a; }


///////////////////////////////////////////////////////////////////

template <typename t> class matrix_trait;

template<>
struct matrix_trait<double>
{
	typedef double vec_type;
	typedef double inverse_type;
	enum { nrOfUnknowns = 1 } ;
};
template <>
struct matrix_trait<fixedMatrix<2> >
{
	typedef fixedVector<2> vec_type;
	typedef fixedVector<2> inverse_type;
	enum { nrOfUnknowns = 2 } ;
};

template <int n>
struct matrix_trait<fixedMatrix<n> >
{
	typedef fixedVector<n> vec_type;
	typedef smallInverse<n> inverse_type;
	enum { nrOfUnknowns = n } ;
};

////////////////////////////////////////////////////////////////////
template <typename t> class vec_traits;

template<>
struct vec_traits<double>
{
	enum { nrOfUnknowns = 1 };
};
template <int n>
struct vec_traits<fixedVector<n> >
{
	enum { nrOfUnknowns = n };
};


////////////////////////////////////////////////////////////////////

template<typename mat_type, typename vec_type> struct Mult_Traits;
template<> struct Mult_Traits<double, double>
{
	typedef double ReturnType;
};
template<int n> struct Mult_Traits<fixedMatrix<n>, fixedVector<n> >
{
	typedef fixedVector<n> ReturnType;
};
template<int n> struct Mult_Traits<double, fixedVector<n> >
{
	typedef fixedVector<n> ReturnType;
};



