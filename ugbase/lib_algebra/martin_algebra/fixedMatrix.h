/*
 *  fixedMatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#include <veclib/cblas.h>
#include <veclib/clapack.h>

template<typename storage_type, int n> class fixedVector;

#include "arrayStorage.h"
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

template<typename storage_type, int rows_=0, int cols_=0>
class fixedMatrix
{
// storage specific
//--------------------
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array_type;
	typedef fixedMatrix<storage_type, rows_, cols_> matrix_type;
	typedef fixedVector<storage_type, rows_> vector_type;
	enum { fixed_rows=rows_, fixed_cols = cols_ };
public:
	inline void setSize(int rows, int cols, bool bZero=true)
	{
		values.setSize(rows, cols, bZero);
	}

	inline int getCols() const
	{
		return values.getCols();
	}	

	inline int getRows() const
	{
		return values.getRows();
	}
private:
	array_type values;


public:
	fixedMatrix() : values() 
	{		
	}
	fixedMatrix(int rows, int cols) : values(rows, cols)
	{			
	}
	

	
	inline double &getAt(int r, int c)
	{
		return values(r, c);
	}
	inline double getAt(int r, int c) const
	{
		return values(r, c);
	}
	inline double operator ()(int r, int c) const
	{
		return getAt(r, c);
	}
	inline double &operator ()(int r, int c)
	{
		return getAt(r, c);
	}
	
	friend ostream &operator << (ostream &out, const matrix_type &s)
	{
		out <<  storage_type::getType() << " mat " << s.getRows() << "x" << s.getCols() << " : ";
		out << "[ ";
		for(int r=0; r < s.getRows(); r++)
		{
			for(int c=0; c< s.getCols(); c++)
				out << s(r, c) << " ";			
			if(r != s.getRows() -1) out << "| ";
		}
		out << "] ";
		return out;
	}	
	
	double operator = (double d)
	{
		for(int i=0; i < values.size(); i++)
			values[i] = 0;
		if(d == 0.0) return d;
		int n = min(getRows(), getCols());
		for(int i=0; i < n; i++)
			getAt(i, i) = d;
		return d;
	}
	
	void operator = (const matrix_type &other)
	{
		values.setSize(other.getRows(), other.getCols(), false);
		for(int i=0; i < other.values.size(); i++) values[i] = other.values[i];
		//memcpy(values, other.values, sizeof(double)*n*n);	
	}
	
	matrix_type operator + (const matrix_type &other ) const
	{
		ASSERT(getRows() == other.getRows() && getCols() == other.getCols());
		matrix_type erg;
		erg.setSize(getRows(), getCols(), false);
		for(int i=0; i < values.size(); i++)
			erg.values[i] = values[i] + other.values[i];
		return erg;
	}
	
	void operator += (const matrix_type &other )
	{
		if(getRows() == 0 && getCols() == 0)
			setSize(other.getRows(), other.getCols());
		else 
		{ ASSERT(getRows() == other.getRows() && getCols() == other.getCols()); }
		for(int i=0; i < values.size(); i++)
			values[i] += other.values[i];
	}	
	
	matrix_type operator - (const matrix_type &other ) const
	{
		ASSERT(getRows() == other.getRows() && getCols() == other.getCols());
		matrix_type erg;
		erg.setSize(getRows(), getCols(), false);

		for(int i=0; i<values.size(); i++)
			erg.values[i] = values[i] - other.values[i];
		return erg;
	}
	
	void operator -= (const matrix_type &other )
	{
		if(getRows() == 0 && getCols() == 0)
			setSize(other.getRows(), other.getCols());
		else 
		{ ASSERT(getRows() == other.getRows() && getCols() == other.getCols()); }
		
		for(int i=0; i<values.size(); i++)
			values[i] -= other.values[i];
	}	
	
	template<int other_rows, int other_cols>
	fixedMatrix<storage_type, rows_, other_cols> operator * (const fixedMatrix<storage_type, other_rows, other_cols> &other ) const
	{
		ASSERT(getCols() == other.getRows());
		
		fixedMatrix<storage_type, rows_, other_cols> erg;
		erg.setSize(getRows(), other.getCols(), false);
		
		for(int r=0; r < getRows(); r++)
			for(int c=0; c < other.getCols(); c++)
			{
				double s = 0;
				for(int i=0; i < getCols(); i++)
					s += getAt(r, i) * other.getAt(i, c);
				erg(r, c) = s;
			}
		return erg;
	}
	
	matrix_type operator * (double d) const
	{
		matrix_type erg(getRows(), getCols());
		for(int i=0; i<values.size(); i++)
				erg.values[i] = values[i]*d;
		return erg;
	}
	
	vector_type operator * (const vector_type &vec ) const
	{
		ASSERT(getRows() == vec.getSize());
		vector_type erg(getRows());
		for(int r=0; r < getRows(); r++)
		{
			double s = 0;
			for(int c=0; c < getCols(); c++)
				s += getAt(r, c) * vec(c);
			erg(r) = s;
		}			
		return erg;
	}
	
	bool operator == (double d) const
	{
		for(int i=0; i<values.size(); i++)
			if(values[i] != d) 
				return false;
		return true;
	}
	bool operator != (double d) const
	{
		return ! operator == (d);
	}
	double norm() const
	{
		double s = 0;
		for(int i=0; i<values.size(); i++)
			s += values[i]*values[i];
		return sqrt(s);			
	}
	double norm2() const
	{
		double s = 0;
		for(int i=0; i<values.size(); i++)
			s += values[i]*values[i];
		return s;
	}
	
	inline void setAsInverseOf(const matrix_type &mat );
	
	void p();
	void print() { p(); }
	

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
template<typename storage_type, int rows_, int cols_>
class smallInverse
{
	//storage
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array2_type;
	typedef typename storage_traits<storage_type, __CLPK_integer, rows_, 0>::array_type interchange_array_type;
	
	typedef fixedVector<storage_type, rows_> vector_type;
public:
	array2_type densemat;
	interchange_array_type interchange;
	
	inline int getCols() const
	{
		return densemat.getCols();
	}	
	
	inline int getRows() const
	{
		return densemat.getRows();
	}
	
///
public:
	
	void setAsInverseOf(const fixedMatrix<storage_type, rows_, cols_> &mat)
	{
		ASSERT(mat.getRows() == mat.getCols());
		__CLPK_integer rows = mat.getRows();
		__CLPK_integer cols = mat.getCols();
		
		densemat.setSize(rows, cols);
		for(int r=0; r < rows; r++)
			for(int c=0; c < cols; c++)
				densemat[c + r*cols] = mat(r, c);
		
		interchange.setSize(rows);
		
		__CLPK_integer info = 0;
		
		dgetrf_(&rows, &cols, &densemat[0], &rows, &interchange[0], &info);
		ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}
	
	vector_type operator * (const vector_type &vec) 
	{
		ASSERT(getRows() == getCols() && getCols() == vec.getSize());
		vector_type erg =vec;
		char trans ='N';
		__CLPK_integer nrhs = 1;
		__CLPK_integer dim = getRows();
		__CLPK_integer info = 0;
		dgetrs_(&trans, &dim, &nrhs,  const_cast<double*> (&densemat(0,0)), &dim, const_cast<__CLPK_integer*> (&interchange[0]), &erg.values[0], &dim, &info);	
		return erg;
	}
};

template<typename storage_type, int rows, int cols>
fixedVector<storage_type, rows> operator * (const fixedVector<storage_type, rows> &vec, const smallInverse<storage_type, rows, cols> &mat)
{
	return mat * vec;
}


template<typename storage_type, int n_>
class fixedVector
{
	//storage
	typedef typename storage_traits<storage_type, double, n_, 0>::array_type array_type;
	typedef fixedVector<storage_type, n_> vector_type;
	enum { fixed_n=n_};
public:
	inline void setSize(int n, bool bZero=true)
	{
		values.setSize(n, bZero);
	}
	
	inline int getSize() const
	{
		return values.size();
	}	
//private:
	array_type values;

	///
	
	
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
	
	
	fixedVector() : values()
	{		
	}
	
	fixedVector(int n) : values(n)
	{
	}
	
	fixedVector(const vector_type &other)
	{
		values = other.values;
	}
	
	friend ostream &operator << (ostream &out, const vector_type &v)
	{
		out << "( ";
		for(int i=0; i < v.getSize(); i++)
			out << v(i) << " ";			
		out << ") ";
		return out;
	}
	
	double operator = (double d)
	{
		for(int i=0; i<getSize(); i++)
			getAt(i) = d;
		return d;
	}
	
	void operator = (const vector_type &other)
	{
		values = other.values;
	//	memcpy(values.values, other.values, sizeof(double)*getSize());	
	}
	
	vector_type operator + (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			ASSERT(getSize() == other.getSize());
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] + other.values[i];
			return erg;
		}
	}
	
	void operator += (const vector_type &other )
	{
		if(other.getSize() == 0) return;
		ASSERT(getSize() == other.getSize());
		for(int i=0; i<getSize(); i++)
			values[i] += other.values[i];		
	}
	
	vector_type operator - (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			ASSERT(getSize() == other.getSize());
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] - other.values[i];
			return erg;
		}
	}
	
	void operator -= (const vector_type &other )
	{
		if(other.getSize() == 0) return;
			
		ASSERT(getSize() == other.getSize());
		for(int i=0; i<getSize(); i++)
			values[i] -= other.values[i];		
	}
	
	double operator * (const vector_type &other ) const
	{
		if(other.getSize() == 0) return 0.0;			
		ASSERT(getSize() == other.getSize());
		double s=0;
		for(int i=0; i<getSize(); i++)
			s += values[i] * other.values[i];
		return s;
	}
	
	vector_type operator * (double alpha) const
	{
		vector_type erg(getSize());
		for(int i=0; i<getSize(); i++)
			erg(i) = getAt(i) * alpha;
		return erg;
	}
	
	template<typename array_type>
	vector_type operator * (const fixedMatrix<array_type> &mat )
	{
		return (mat * (*this));
	}
	
	
	template<typename array_type>
	inline void operator /= (const fixedMatrix<array_type> &mat);
	
	template<typename array_type>
	vector_type operator / (const fixedMatrix<array_type> &mat )
	{
		vector_type erg = *this;
		erg /= mat;
		return erg;
	}
	
	double norm2() const
	{
		double s =0;
		for(int i=0; i<getSize(); i++) s += getAt(i)*getAt(i);
		return s;
	}
	
	void p();
	void print() { p(); }
};

template<typename storage_type, int n>
inline fixedVector<storage_type, n> operator * (double alpha, const fixedVector<storage_type, n> &v) 
{
	return v * alpha;
}

/*
template<typename storage_type>
inline void fixedVector<storage_type, 1>::operator /= (const fixedMatrix<storage_type, 1, 1> &mat )
{
	getAt(0) = getAt(0) / mat(0, 0);
}

template<typename storage_type>
inline void fixedVector<storage_type, 2>::operator /= (const fixedMatrix<storage_type, 2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	ASSERT(invD != 0.0);
	
	double a = invD*(getAt(0) * mat(1,1) - getAt(1) * mat(0,1));
	double b = invD*(getAt(1) * mat(0,0) - getAt(0) * mat(1,0));
	getAt(0) = a;
	getAt(1) = b;
}


template<typename storage_type>
inline void fixedVector<storage_type, rows_>::operator /= (const fixedMatrix<storage_type, rows_, rows_> &mat )
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
*/
/*template<>
template<typename array_type>
inline void fixedMatrix<array_type, 2>::setAsInverseOf(const fixedMatrix<array_type, 2> &mat )
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



template<typename storage_type, int n>
fixedVector<storage_type, n> operator * (double alpha, fixedVector<storage_type, n> &vec)
{
	return vec*alpha;
}

/////////////////////////////////////////////////////////////////////

template <typename t> class matrix_trait;
template <typename t> class vec_traits;
template <typename A, typename B> class Mult_Traits;

template <typename storage_type, int rows_, int cols_>
struct matrix_trait< fixedMatrix<storage_type, rows_, cols_> >
{
	typedef fixedVector< storage_type, rows_ > vec_type;
	typedef smallInverse<storage_type, rows_, cols_> inverse_type;
	enum { nrOfUnknowns = rows_ } ;
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, int n>
struct vec_traits<fixedVector<storage_type, n> >
{
	enum { nrOfUnknowns = n };
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, int rows_, int cols_>
struct Mult_Traits<fixedMatrix<storage_type, rows_, cols_>, fixedVector<storage_type, cols_> >
{
	typedef fixedVector<storage_type, rows_> ReturnType;
};

template <typename storage_type, int rows_, int cr_, int cols2_>
struct Mult_Traits<fixedMatrix<storage_type, rows_, cr_>, fixedMatrix<storage_type, cr_, cols2_> >
{
	typedef fixedMatrix<storage_type, rows_, cols2_> ReturnType;
};

template <typename storage_type, int n_>
struct Mult_Traits<double, fixedVector<storage_type, n_> >
{
	typedef fixedVector<storage_type, n_> ReturnType;
};
//#endif

//////////////////////////////////////////////////////
// wrapper for using doubles
//setSize(t, a, b)
template<typename T>
inline void setSize(T &t, int a, int b)
{
	t.setSize(a, b);
}

template<>
inline void setSize(double &d, int a, int b)
{
	return;
}
//setSize(t, a)
template<typename T>
inline void setSize(T &t, int a)
{
	t.setSize(a);
}

template<>
inline void setSize(double &d, int a)
{
	return;
}
// getSize
template<typename T>
inline int getSize(T &t)
{
	return t.getSize();
}

template<>
inline int getSize(double &t)
{
	return 1;
}

//getRows
template<typename T>
inline int getRows(const T &t)
{
	return t.getRows();
}

template<typename T>
inline int getCols(const T &t)
{
	return t.getCols();
}

template<>
inline int getRows(const double &t)
{
	return 1;
}

template<>
inline int getCols(const double &t)
{
	return 1;
}
