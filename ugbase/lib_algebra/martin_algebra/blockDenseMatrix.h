/*
 *  blockDenseMatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#include <veclib/cblas.h>
#include <veclib/clapack.h>

template<typename storage_type, int n> class blockVector;

#include "arrayStorage.h"
///////////////////////////////////////////////////////////////////////////////////////
//!
//! blockDenseMatrix
//! fixed n x n - Matrix. 
//! supports 
//! +, -, *, +=, -=, = with other blockDenseMatrix<n>
//! * with blockVector<n>, * with double
//! =, == double (== id*d)
//! getAt(int row, int column), operator() (row, column)
//! ostream operator <<
//! norm(), norm2(), print(), p(){print();}

template<typename storage_type, int rows_=0, int cols_=0>
class blockDenseMatrix
{
// storage specific
//--------------------
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array_type;
	typedef blockDenseMatrix<storage_type, rows_, cols_> matrix_type;
	typedef blockVector<storage_type, rows_> vector_type;
	enum { fixed_rows=rows_, fixed_cols = cols_ };
public:
	inline void setSize(int rows, int cols, bool bZero=true)
	{
		values.setSize(rows, cols, bZero);
	}

	void swap(matrix_type &other)
	{
		values.swap(other.values);
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
	blockDenseMatrix() : values() 
	{		
	}
	blockDenseMatrix(int rows, int cols) : values(rows, cols)
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
		{ ASSERT1(getRows() == other.getRows() && getCols() == other.getCols()); }
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
	blockDenseMatrix<storage_type, rows_, other_cols> operator * (const blockDenseMatrix<storage_type, other_rows, other_cols> &other ) const
	{
		ASSERT(getCols() == other.getRows());
		
		blockDenseMatrix<storage_type, rows_, other_cols> erg;
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
		ASSERT1(getRows() == vec.getSize());
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
//! setAsInverseOf(const blockDenseMatrix<n> &mat) : init as inverse of mat
//! blockVector<n> * smallInverse<n> = smallInverse<n> * blockVector<n>
//! = A^{-1} b
template<typename storage_type, int rows_, int cols_>
class smallInverse
{
	//storage
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array2_type;
	typedef typename storage_traits<storage_type, __CLPK_integer, rows_, 0>::array_type interchange_array_type;
	
	typedef blockVector<storage_type, rows_> vector_type;
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
	
	void setAsInverseOf(const blockDenseMatrix<storage_type, rows_, cols_> &mat)
	{
		ASSERT1(mat.getRows() == mat.getCols());
		__CLPK_integer rows = mat.getRows();
		__CLPK_integer cols = mat.getCols();
		
		densemat.setSize(rows, cols);
		for(int r=0; r < rows; r++)
			for(int c=0; c < cols; c++)
				densemat[c + r*cols] = mat(r, c);
		
		interchange.setSize(rows);
		
		__CLPK_integer info = 0;
		
		dgetrf_(&rows, &cols, &densemat[0], &rows, &interchange[0], &info);
		ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}
	
	vector_type operator * (const vector_type &vec) 
	{
		ASSERT1(getRows() == getCols() && getCols() == vec.getSize());
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
blockVector<storage_type, rows> operator * (const blockVector<storage_type, rows> &vec, const smallInverse<storage_type, rows, cols> &mat)
{
	return mat * vec;
}


/////////////////////////////////////////////////////////////////////

template <typename t> class matrix_trait;
template <typename t> class vec_traits;
template <typename A, typename B> class Mult_Traits;

template <typename storage_type, int rows_, int cols_>
struct matrix_trait< blockDenseMatrix<storage_type, rows_, cols_> >
{
	typedef blockVector< storage_type, rows_ > vec_type;
	typedef smallInverse<storage_type, rows_, cols_> inverse_type;
	enum { nrOfUnknowns = rows_ } ;
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, int n>
struct vec_traits<blockVector<storage_type, n> >
{
	enum { nrOfUnknowns = n };
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, int rows_, int cols_>
struct Mult_Traits<blockDenseMatrix<storage_type, rows_, cols_>, blockVector<storage_type, cols_> >
{
	typedef blockVector<storage_type, rows_> ReturnType;
};

template <typename storage_type, int rows_, int cr_, int cols2_>
struct Mult_Traits<blockDenseMatrix<storage_type, rows_, cr_>, blockDenseMatrix<storage_type, cr_, cols2_> >
{
	typedef blockDenseMatrix<storage_type, rows_, cols2_> ReturnType;
};

template <typename storage_type, int n_>
struct Mult_Traits<double, blockVector<storage_type, n_> >
{
	typedef blockVector<storage_type, n_> ReturnType;
};
//#endif
