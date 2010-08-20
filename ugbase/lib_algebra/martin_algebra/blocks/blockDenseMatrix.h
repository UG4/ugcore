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
#ifndef __H__UG__MARTIN_ALGEBRA__BLOCK_DENSE_MATRIX__
#define __H__UG__MARTIN_ALGEBRA__BLOCK_DENSE_MATRIX__

#include "math.h"

#include "blocks.h"
#include "double.h"

#include <vector>

#include "../lapack/lapack.h"

#include "arrayStorage.h"

#include "../algebra_misc.h"

namespace ug {
template<typename value_type, typename storage_type, size_t n> class blockVector;

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
//! at(size_t row, size_t column), operator() (row, column)
//! ostream operator <<
//! norm(), norm2(), print(), p(){print();}
template<typename value_type, typename storage_type, size_t rows_=0, size_t cols_=0>
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
	inline void resize(size_t rows, size_t cols, bool bZero=true)
	{
		if(rows != num_rows() || cols != num_cols())
			values.resize(rows, cols, bZero);
	}

	void swap(matrix_type &other)
	{
		values.swap(other.values);
	}


	inline size_t num_cols() const
	{
		return values.num_cols();
	}

	inline size_t num_rows() const
	{
		return values.num_rows();
	}


public:
// create
	blockDenseMatrix() : values()
	{
		FORCE_CREATION { p(); print(); }
	}
	blockDenseMatrix(size_t rows, size_t cols) : values(rows, cols)
	{
		FORCE_CREATION { p(); print(); }
	}

	blockDenseMatrix(const blockDenseMatrix &other) : values(other.values)
	{
		//ASSERT2(0, "thou shall not use the copy constructor for it is slow");
		// copy constructor inevitable in vector<connection> cons;
	}

// access functions
	inline value_type &at(size_t r, size_t c)
	{
		return values(r, c);
	}
	inline const value_type &at(size_t r, size_t c) const
	{
		return values(r, c);
	}
	inline value_type &operator ()(size_t r, size_t c)
	{
		return at(r, c);
	}
	inline const value_type &operator ()(size_t r, size_t c) const
	{
		return at(r, c);
	}


// algebra functions
	double operator = (double d)
	{
		for(size_t i=0; i < values.size(); i++)
			values[i] = 0;
		if(d == 0.0) return d;
		size_t n = min(num_rows(), num_cols());
		for(size_t i=0; i < n; i++)
			at(i, i) = d;
		return d;
	}

	void operator = (const matrix_type &other)
	{
		values.resize(other.num_rows(), other.num_cols(), false);
		for(size_t i=0; i < other.values.size(); i++) values[i] = other.values[i];
		//memcpy(values, other.values, sizeof(double)*n*n);
	}

// add / substract
	matrix_type operator + (const matrix_type &other ) const
	{
		UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
		matrix_type erg;
		erg.resize(num_rows(), num_cols(), false);
		for(size_t i=0; i< values.size(); i++)
			erg.values[i] = values[i] + other.values[i];
		return erg;
	}

	void operator += (const matrix_type &other )
	{
		if(num_rows() == 0 && num_cols() == 0)
			resize(other.num_rows(), other.num_cols());
		else
		{ UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), ""); }

		for(size_t i=0; i< values.size(); i++)
			values[i] += other.values[i];
	}

	matrix_type operator - (const matrix_type &other ) const
	{
		UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
		matrix_type erg;
		erg.resize(num_rows(), num_cols(), false);

		for(size_t i=0; i< values.size(); i++)
			erg.values[i] = values[i] - other.values[i];
		return erg;
	}

	void operator -= (const matrix_type &other )
	{
		if(num_rows() == 0 && num_cols() == 0)
			resize(other.num_rows(), other.num_cols());
		else
		{ UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), ""); }

		for(size_t i=0; i< values.size(); i++)
			values[i] -= other.values[i];
	}

// multiply
	template<size_t other_rows, size_t other_cols>
	blockDenseMatrix<value_type, storage_type, rows_, other_cols> operator * (const blockDenseMatrix<value_type, storage_type, other_rows, other_cols> &other ) const
	{
		UG_ASSERT(num_cols() == other.num_rows(), "");

		blockDenseMatrix<value_type, storage_type, rows_, other_cols> erg;
		erg.resize(num_rows(), other.num_cols(), false);

		for(size_t r=0; r < num_rows(); r++)
			for(size_t c=0; c < other.num_cols(); c++)
			{
				for(size_t i=0; i < num_cols(); i++)
					AddMult(erg(r,c), at(r, i), other.at(i, c));
			}
		return erg;
	}


	matrix_type operator * (double d) const
	{
		matrix_type erg(num_rows(), num_cols());

		for(size_t i=0; i< values.size(); i++)
				erg.values[i] = values[i]*d;
		return erg;
	}

	void operator *= (double d)
	{
		for(size_t i=0; i< values.size(); i++)
			values[i] *= d;
	}

	void operator /= (double d)
	{
		for(size_t i=0; i< values.size(); i++)
			values[i] /= d;
	}

	matrix_type &operator /= (matrix_type &other)
	{
		matrix_type tmp = other;
		tmp.invert();

		(*this) = (*this) * tmp;

		return *this;
	}

	vector_type operator * (const vector_type &vec ) const
	{
		vector_type erg(num_rows());
		assign_mult(erg, vec);
		return erg;
	}

// compare
	bool operator == (double d) const
	{
		for(size_t i=0; i< values.size(); i++)
			if(values[i] != d)
				return false;
		return true;
	}
	bool operator != (double d) const
	{
		return ! operator == (d);
	}

// TODO: I inserted this, since it is still used in operator*. Otherwise we get compile errors. Andreas Vogel
	//! dest = this*vec . use this to prevent temporary variables
	void assign_mult(vector_type &dest, const vector_type &vec) const
	{
		UG_ASSERT(num_rows() == vec.size(), "");
		for(size_t r=0; r < num_rows(); r++)
		{
			dest(r) = 0.0;
			for(size_t c=0; c < num_cols(); c++)
				AddMult(dest(r), at(r, c), vec(c));
		}
	}


// other
	double norm() const
	{
		double s = 0;
		for(size_t i=0; i< values.size(); i++)
			s += BlockNorm2(values[i]);
		return sqrt(s);
	}
	double norm2() const
	{
		double s = 0;
		for(size_t i=0; i< values.size(); i++)
			s += BlockNorm2(values[i]);
		return s;
	}


	void invert()
	{
		std::vector<int> interchange;
		interchange.resize(max(num_cols(), num_rows()));

		int info = getrf(num_rows(), num_cols(), &at(0,0), num_rows(), &interchange[0]);
		UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));

		// calc work size
		double worksize; int iWorksize = -1;
		info = getri(num_rows(), &at(0,0), num_rows(), &interchange[0], &worksize, iWorksize);
		UG_ASSERT(info == 0, "");
		iWorksize = worksize;

		std::vector<double> work;
		work.resize(iWorksize);

		info = getri(num_rows(), &at(0,0), num_rows(), &interchange[0], &work[0], iWorksize);
		UG_ASSERT(info == 0, "");
	}

	//inline void setAsInverseOf(const matrix_type &mat ); // deprecated

// print
	void p() { std::cout << *this << std::endl; }
	void print() { p(); }

	friend std::ostream &operator << (std::ostream &out, const matrix_type &s)
	{
		//out <<  storage_type::getType() << " mat " << s.num_rows() << "x" << s.num_cols() << " : ";
		out << "[ ";
		for(size_t r=0; r < s.num_rows(); r++)
		{
			for(size_t c=0; c< s.num_cols(); c++)
				out << s(r, c) << " ";
			if(r != s.num_rows() -1) out << "| ";
		}
		out << "] ";
		return out;
	}
};

template<typename value_type, typename storage_type, size_t rows, size_t cols>
blockDenseMatrix<value_type, storage_type, rows, cols> operator * (double d, const blockDenseMatrix<value_type, storage_type, rows, cols> &mat)
{
	return mat * d;
}

///////////////////////////////////////////////////////////////////////////////////////
//!
//! smallInverse<size_t n>
//! A class to hold a inverse of a smallMatrix<n>
//! implemented with LAPACKs LU-Decomposition dgetrf
//! (uses double[n*n] for LU and interchange[n] for pivoting
//! functions:
//! setAsInverseOf(const blockDenseMatrix<n> &mat) : init as inverse of mat
//! blockVector<n> * smallInverse<n> = smallInverse<n> * blockVector<n>
//! = A^{-1} b
template<typename storage_type, size_t rows_, size_t cols_>
class smallInverse
{
private: // storage
	typedef typename storage_traits<storage_type, double, rows_, cols_>::array2_type array2_type;
	typedef typename storage_traits<storage_type, lapack_int, rows_, 0>::array_type interchange_array_type;

	typedef blockVector<double, storage_type, rows_> vector_type;

	array2_type densemat;
	interchange_array_type interchange;

public:
	inline size_t num_cols() const
	{
		return densemat.num_cols();
	}

	inline size_t num_rows() const
	{
		return densemat.num_rows();
	}

///
public:
	//! initializes this object as inverse of mat
	void setAsInverseOf(const blockDenseMatrix<double, storage_type, rows_, cols_> &mat)
	{
		int rows = mat.num_rows(), cols = mat.num_cols();
		UG_ASSERT(rows == cols, "only for square matrices");

		densemat.resize(rows, cols);
		for(size_t r=0; r < rows; r++)
			for(size_t c=0; c < cols; c++)
				densemat[r + c*rows] = mat(r, c);

		interchange.resize(rows);

		int info = getrf(rows, cols, &densemat(0,0), rows, &interchange[0]);
		UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
	}

	//! calc dest = mat^{-1} * vec
	void apply(double *dest, const vector_type &vec) const
	{
		UG_ASSERT(num_rows() == num_cols(), "only for square matrices");
		UG_ASSERT(num_cols() == vec.size(), "size mismatch");
		for(size_t i=0; i<vec.size(); i++)
			dest[i] = vec(i);

		int info = getrs(ModeNoTrans, num_rows(), 1, &densemat(0,0), num_rows(), &interchange[0], dest,	num_rows());
	}

	vector_type operator * (const vector_type &vec) const
	{
		vector_type erg(vec.size());
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
		static std::vector<double> erg;
		erg.resize(vec.size());

		apply(&erg[0], vec);
		for(size_t i=0; i<vec.size(); i++)
			dest(i) += erg[i];
	}
	//! dest -= this*vec . use this to prevent temporary variables
	void sub_mult(vector_type &dest, const vector_type &vec) const
	{
		// we need one temporary variable
		// keep static so it gets reallocated only once or twice
		static std::vector<double> erg;
		erg.resize(vec.size());

		apply(&erg[0], vec);
		for(size_t i=0; i<vec.size(); i++)
			dest(i) -= erg[i];
	}
};

template<typename storage_type, size_t rows, size_t cols>
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

template <typename storage_type, size_t rows_, size_t cols_>
struct block_matrix_traits< blockDenseMatrix<double, storage_type, rows_, cols_> >
{
	typedef blockVector<double, storage_type, rows_ > vec_type;
	typedef smallInverse<storage_type, rows_, cols_> inverse_type;
	enum { nrOfRows = rows_ } ;
	enum { nrOfCols = cols_ } ;
};

/////////////////////////////////////////////////////////////////
template <typename storage_type, size_t n>
struct block_vector_traits<blockVector<double, storage_type, n> >
{
	enum { nrOfUnknowns = n };
};

/////////////////////////////////////////////////////////////////

//! mat * vec
template <typename storage_type, size_t rows_, size_t cols_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cols_>, blockVector<double, storage_type, cols_> >
{
	typedef blockVector<double, storage_type, rows_> ReturnType;
};

//! mat * mat
template <typename storage_type, size_t rows_, size_t cr_, size_t cols2_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cr_>, blockDenseMatrix<double, storage_type, cr_, cols2_> >
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols2_> ReturnType;
};

//! double * mat
template <typename storage_type, size_t rows_, size_t cols_>
struct block_multiply_traits<double, blockDenseMatrix<double, storage_type, rows_, cols_> >
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols_> ReturnType;
};

//! mat * double
template <typename storage_type, size_t rows_, size_t cols_>
struct block_multiply_traits<blockDenseMatrix<double, storage_type, rows_, cols_>, double>
{
	typedef blockDenseMatrix<double, storage_type, rows_, cols_> ReturnType;
};

//! double * vec
template <typename storage_type, size_t n_>
struct block_multiply_traits<double, blockVector<double, storage_type, n_> >
{
	typedef blockVector<double, storage_type, n_> ReturnType;
};

//! vec * double
template <typename storage_type, size_t n_>
struct block_multiply_traits<blockVector<double, storage_type, n_>, double >
{
	typedef blockVector<double, storage_type, n_> ReturnType;
};
//#endif


template<typename M>
inline void GetInverse(typename block_matrix_traits<M>::inverse_type &inv, const M &m)
{
	inv.setAsInverseOf(m);
}





} // namespace ug

#endif

