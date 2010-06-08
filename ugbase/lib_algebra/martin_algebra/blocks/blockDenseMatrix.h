/*
 *  blockDenseMatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  Dense matrices mainly for use inside of SparseMatrix
 *  optimized for small memory consumption
 *  SmallInverse: Inverse for blockDenseMatrices using LAPACK
 */
#pragma once

#include "math.h"

#include "blocks.h"
#include "double.h"

#include <vector>
#include <veclib/cblas.h>
#include <veclib/clapack.h>

#include "arrayStorage.h"

#include "../algebra_misc.h"

namespace ug {
template<typename value_type, typename storage_type, int n> class blockVector;

///////////////////////////////////////////////////////////////////////////////////////
//!
//! blockDenseMatrix
//! template parameters:
//! 1. value_type: i.e. float, double or blockDenseMatrix (recursive)
//! 2. storage_type: fixedStorage, arrayStorage
//! 3. rows, 4. cols: with storage_type=fixedStorage, size of the fixed matrix
//! if storage_type=variableStorage, rows and cols are ignored
//! supports 
//! +, -, *, +=, -=, = with other blockDenseMatrix<n>
//! * with blockVector<n>, * with double
//! =, == double (== id*d)
//! getAt(int row, int column), operator() (row, column)
//! ostream operator <<
//! norm(), norm2(), print(), p(){print();}
template<typename value_type, typename storage_type, int rows_=0, int cols_=0>
class blockDenseMatrix
{	
private: 
// storage specific
	typedef typename storage_traits<storage_type, value_type, rows_, cols_>::array2_type array_type;
	typedef blockDenseMatrix<value_type, storage_type, rows_, cols_> matrix_type;
	typedef blockVector<value_type, storage_type, rows_> vector_type;	
	enum { fixed_rows=rows_, fixed_cols = cols_ };
	
	array_type values;
	
public:
	inline void setSize(int rows, int cols, bool bZero=true)
	{
		if(rows != getRows() || cols != getCols())
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


public: 
// create
	blockDenseMatrix() : values() 
	{
		FORCE_CREATION { p(); print(); }
	}
	blockDenseMatrix(int rows, int cols) : values(rows, cols)
	{	
		FORCE_CREATION { p(); print(); }
	}

	blockDenseMatrix(const blockDenseMatrix &other) : values(other.values)
	{
		//ASSERT2(0, "thou shall not use the copy constructor for it is slow");
		// copy constructor inevitable in vector<connection> cons;
	}

// access functions
	inline value_type &getAt(int r, int c)
	{
		return values(r, c);
	}
	inline const value_type &getAt(int r, int c) const
	{
		return values(r, c);
	}
	inline value_type &operator ()(int r, int c)
	{
		return getAt(r, c);
	}
	inline const value_type &operator ()(int r, int c) const
	{
		return getAt(r, c);
	}


// algebra functions
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
	
// add / substract
	matrix_type operator + (const matrix_type &other ) const
	{
		UG_ASSERT(getRows() == other.getRows() && getCols() == other.getCols(), "");
		matrix_type erg;
		erg.setSize(getRows(), getCols(), false);
		for(int i=0; i< values.size(); i++)
			erg.values[i] = values[i] + other.values[i];
		return erg;
	}
	
	void operator += (const matrix_type &other )
	{
		if(getRows() == 0 && getCols() == 0)
			setSize(other.getRows(), other.getCols());
		else 
		{ UG_ASSERT(getRows() == other.getRows() && getCols() == other.getCols(), ""); }
		
		for(int i=0; i< values.size(); i++)
			values[i] += other.values[i];
	}	
	
	matrix_type operator - (const matrix_type &other ) const
	{
		UG_ASSERT(getRows() == other.getRows() && getCols() == other.getCols(), "");
		matrix_type erg;
		erg.setSize(getRows(), getCols(), false);

		for(int i=0; i< values.size(); i++)
			erg.values[i] = values[i] - other.values[i];
		return erg;
	}
	
	void operator -= (const matrix_type &other )
	{
		if(getRows() == 0 && getCols() == 0)
			setSize(other.getRows(), other.getCols());
		else 
		{ UG_ASSERT(getRows() == other.getRows() && getCols() == other.getCols(), ""); }
		
		for(int i=0; i< values.size(); i++)
			values[i] -= other.values[i];
	}	
	
// multiply
	template<int other_rows, int other_cols>
	blockDenseMatrix<value_type, storage_type, rows_, other_cols> operator * (const blockDenseMatrix<value_type, storage_type, other_rows, other_cols> &other ) const
	{
		UG_ASSERT(getCols() == other.getRows(), "");
		
		blockDenseMatrix<value_type, storage_type, rows_, other_cols> erg;
		erg.setSize(getRows(), other.getCols(), false);
		
		for(int r=0; r < getRows(); r++)
			for(int c=0; c < other.getCols(); c++)
			{
				for(int i=0; i < getCols(); i++)
					add_mult(erg(r,c), getAt(r, i), other.getAt(i, c));
			}
		return erg;
	}


	matrix_type operator * (double d) const
	{
		matrix_type erg(getRows(), getCols());
		
		for(int i=0; i< values.size(); i++)
				erg.values[i] = values[i]*d;
		return erg;
	}
	
	vector_type operator * (const vector_type &vec ) const
	{
		vector_type erg(getRows());
		assign_mult(erg, vec);
		return erg;
	}

// compare
	bool operator == (double d) const
	{		
		for(int i=0; i< values.size(); i++)
			if(values[i] != d) 
				return false;
		return true;
	}
	bool operator != (double d) const
	{
		return ! operator == (d);
	}

// temporary prevention
	//! dest -= this*vec . use this to prevent temporary variables	
	void sub_mult(vector_type &dest, const vector_type &vec) const
	{
		UG_ASSERT1(getRows() == vec.getSize());
		for(int r=0; r < getRows(); r++)
		{
			for(int c=0; c < getCols(); c++)
				SubMult(dest(r), getAt(r, c), vec(c));
		}			
	}

	//! dest += this*vec . use this to prevent temporary variables	
	void add_mult(vector_type &dest, const vector_type &vec) const
	{
		UG_ASSERT(getRows() == vec.getSize(), "");
		for(int r=0; r < getRows(); r++)
		{
			for(int c=0; c < getCols(); c++)
				AddMult(dest(r), getAt(r, c), vec(c));
		}			
	}

	//! dest = this*vec . use this to prevent temporary variables	
	void assign_mult(vector_type &dest, const vector_type &vec) const
	{
		UG_ASSERT(getRows() == vec.getSize(), "");
		for(int r=0; r < getRows(); r++)
		{
			dest(r) = 0.0;
			for(int c=0; c < getCols(); c++)
				AddMult(dest(r), getAt(r, c), vec(c));
		}			
	}


	//! this -= alpha*mat . use this to prevent temporary variables	
	void sub_mult(double alpha, const matrix_type &mat) 
	{
		setSize(mat.getRows(), mat.getCols());
		for(int r=0; r < getRows(); r++)
		{
			for(int c=0; c < getCols(); c++)
				SubMult(getAt(r, c), alpha, mat.getAt(r, c));
		}
	}

	//! this += alpha*mat . use this to prevent temporary variables	
	void add_mult(const double alpha, const matrix_type &mat) 
	{
		setSize(mat.getRows(), mat.getCols());
		for(int r=0; r < getRows(); r++)
		{
			for(int c=0; c < getCols(); c++)
				AddMult(getAt(r, c), alpha, mat.getAt(r, c));
		}
	}

	//! this = alpha*mat . use this to prevent temporary variables	
	void assign_mult(const double alpha, const matrix_type &mat) 
	{
		setSize(mat.getRows(), mat.getCols());
		for(int r=0; r < getRows(); r++)
		{
			for(int c=0; c < getCols(); c++)
				AssignMult(getAt(r, c), alpha, mat.getAt(r, c));
		}			
	}
	void assign_mult(const matrix_type &mat, double alpha)
	{
		assign_mult(alpha, mat);
	}
	
// other
	double norm() const
	{
		double s = 0;
		for(int i=0; i< values.size(); i++)
			s += mnorm2(values[i]);
		return sqrt(s);			
	}
	double norm2() const
	{
		double s = 0;
		for(int i=0; i< values.size(); i++)
			s += mnorm2(values[i]);
		return s;
	}


	void invert()
	{
		static vector<__CLPK_integer> interchange;
		interchange.resize(max(getCols(), getRows()));
		__CLPK_integer info = 0;
		__CLPK_integer rows = getRows();
		__CLPK_integer cols = getCols();
		
		double *ptr = &getAt(0,0); 
		
		dgetrf_(&rows, &cols, ptr, &rows, &interchange[0], &info);
		UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
		
		// calc work size
		double worksize; __CLPK_integer iWorksize = -1;
		dgetri_(&rows, ptr, &rows, &interchange[0], &worksize, &iWorksize, &info);
		UG_ASSERT(info == 0, "");
		iWorksize = worksize;
		
		static vector<double> work;
		work.resize(iWorksize);

		dgetri_(&rows, ptr, &rows, &interchange[0], &work[0], &iWorksize, &info);
		UG_ASSERT(info == 0, "");
	}
	
	//inline void setAsInverseOf(const matrix_type &mat ); // deprecated
	
// print
	void p() { cout << *this; cout.flush();}
	void print() { p(); }

	friend ostream &operator << (ostream &out, const matrix_type &s)
	{
		//out <<  storage_type::getType() << " mat " << s.getRows() << "x" << s.getCols() << " : ";
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
};

template<typename value_type, typename storage_type, int rows, int cols>
blockDenseMatrix<value_type, storage_type, rows, cols> operator * (double d, const blockDenseMatrix<value_type, storage_type, rows, cols> &mat)
{
	return mat * d;
}

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
private: // storage
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array2_type;
	typedef typename storage_traits<storage_type, __CLPK_integer, rows_, 0>::array_type interchange_array_type;
	
	typedef blockVector<double, storage_type, rows_> vector_type;
	
	array2_type densemat;
	interchange_array_type interchange;
	
public:	
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
	//! initializes this object as inverse of mat
	void setAsInverseOf(const blockDenseMatrix<double, storage_type, rows_, cols_> &mat)
	{
		UG_ASSERT(mat.getRows() == mat.getCols(), "");
		__CLPK_integer rows = mat.getRows();
		__CLPK_integer cols = mat.getCols();
		
		densemat.setSize(rows, cols);
		for(int r=0; r < rows; r++)
			for(int c=0; c < cols; c++)
				densemat[c + r*cols] = mat(r, c);
		
		interchange.setSize(rows);
		
		__CLPK_integer info = 0;
		
		dgetrf_(&rows, &cols, &densemat[0], &rows, &interchange[0], &info);
		UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}
	
	//! calc dest = mat^{-1} * vec
	void apply(double *dest, const vector_type &vec) const
	{
		UG_ASSERT(getRows() == getCols() && getCols() == vec.getSize(), "");
		for(int i=0; i<vec.getSize(); i++)
			dest[i] = vec(i);
		char trans ='N';
		__CLPK_integer nrhs = 1;
		__CLPK_integer dim = getRows();
		__CLPK_integer info = 0;
		
		dgetrs_(&trans, 
				&dim, 
				&nrhs,  
				&(*const_cast<array2_type*> (&densemat))(0,0), 
				&dim, 
				&(*const_cast<interchange_array_type*> (&interchange))[0], 
				dest, 
				&dim, 
				&info);	
	}
	
	vector_type operator * (const vector_type &vec) const
	{
		vector_type erg(vec.getSize());
		apply(&erg(0), vec);
		return erg;
	}

// temporary prevention
	//! dest = this*vec . use this to prevent temporary variables
	void assign_mult(vector_type &dest, const vector_type &vec) const
	{
		apply(&dest(0), vec);
	}
	//! dest += this*vec . use this to prevent temporary variables
	void add_mult(vector_type &dest, const vector_type &vec) const
	{
		// we need one temporary variable
		// keep static so it gets reallocated only once or twice
		static vector<double> erg;
		erg.resize(vec.getSize());
		
		apply(&erg[0], vec);
		for(int i=0; i<vec.getSize(); i++)
			dest(i) += erg[i];
	}
	//! dest -= this*vec . use this to prevent temporary variables
	void sub_mult(vector_type &dest, const vector_type &vec) const
	{
		// we need one temporary variable
		// keep static so it gets reallocated only once or twice
		static vector<double> erg;
		erg.resize(vec.getSize());
		
		apply(&erg[0], vec);
		for(int i=0; i<vec.getSize(); i++)
			dest(i) -= erg[i];
	}
};

template<typename storage_type, int rows, int cols>
blockVector<double, storage_type, rows> operator * (const blockVector<double, storage_type, rows> &vec, const smallInverse<storage_type, rows, cols> &mat)
{
	return mat * vec;
}

/////////////////////////////////////////////////////////////////////
//
//						MATRIX TRAITS
//
/////////////////////////////////////////////////////////////////////

template <typename t> class block_matrix_traits;
template <typename t> class block_vector_traits;
template <typename A, typename B> class block_multiply_Traits;

template <typename storage_type, int rows_, int cols_>
struct block_matrix_traits< blockDenseMatrix<double, storage_type, rows_, cols_> >
{
	typedef blockVector<double, storage_type, rows_ > vec_type;
	typedef smallInverse<storage_type, rows_, cols_> inverse_type;
	enum { nrOfUnknowns = rows_ } ;
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, int n>
struct block_vector_traits<blockVector<double, storage_type, n> >
{
	enum { nrOfUnknowns = n };
};

/////////////////////////////////////////////////////////////////

//! mat * vec
template <typename storage_type, int rows_, int cols_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cols_>, blockVector<double, storage_type, cols_> >
{
	typedef blockVector<double, storage_type, rows_> ReturnType;
};

//! mat * mat
template <typename storage_type, int rows_, int cr_, int cols2_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cr_>, blockDenseMatrix<double, storage_type, cr_, cols2_> >
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols2_> ReturnType;
};

//! double * mat
template <typename storage_type, int rows_, int cols_>
struct block_multiply_traits<double, blockDenseMatrix<double, storage_type, rows_, cols_> >
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols_> ReturnType;
};

//! mat * double
template <typename storage_type, int rows_, int cols_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cols_>, double>
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols_> ReturnType;
};

//! double * vec
template <typename storage_type, int n_>
struct block_multiply_traits<double, blockVector<double, storage_type, n_> >
{
	typedef blockVector<double, storage_type, n_> ReturnType;
};

//! vec * double
template <typename storage_type, int n_>
struct block_multiply_traits<blockVector<double, storage_type, n_>, double >
{
	typedef blockVector<double, storage_type, n_> ReturnType;
};
//#endif

} // namespace ug
