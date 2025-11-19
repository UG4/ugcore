/*
 * Copyright (c) 2010-2013:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


#ifndef __H__UG__SMALL_ALGEBRA__BLOCK_DENSE__
#define __H__UG__SMALL_ALGEBRA__BLOCK_DENSE__

#include "densematrix.h"
#include "densevector.h"
#include <algorithm>

namespace ug{

template<typename A>
inline double BlockNorm2(const DenseMatrix<A> &mat)
{
	double sum=0;
	for(size_t r=0; r < mat.num_rows(); r++)
		for(size_t c=0; c < mat.num_cols(); c++)
			sum += BlockNorm2(mat(r, c));

	return sum;
}


template<typename A>
inline double BlockNorm2(const DenseVector<A> &v)
{
	double sum=0;
	for(size_t i=0; i < v.size(); i++)
		sum += BlockNorm2(v[i]);
	return sum;
}



template<typename A>
inline double BlockMaxNorm(const DenseVector<A> &v)
{
	double max=0;
	for(size_t i=0; i < v.size(); i++)
		max = std::max(max, BlockMaxNorm(v[i]));
	return max;
}




//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables

// vec = mat * vec
// todo: replace add_mult etc. with template expressions
// dest = b*vec
template<typename A, typename B, typename C>
inline void AssignMult(DenseVector<A> &dest, const DenseMatrix<B> &mat, const DenseVector<C> &vec)
{
	UG_ASSERT(mat.num_cols() == vec.size(), "");
	dest.resize(mat.num_rows());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		AssignMult(dest(r), mat(r, 0), vec(0));
		for(size_t c=1; c < mat.num_cols(); c++)
			AddMult(dest(r), mat(r, c), vec(c));
	}
}

template<typename A, typename B, typename C>
inline void AssignMult(DenseMatrix<A> &dest, const DenseMatrix<B> &mA, const DenseMatrix<C> &mB)
{
	UG_ASSERT(mA.num_cols() == mB.num_rows(), "");
	dest.resize(mA.num_rows(), mB.num_cols());
	for(size_t r=0; r < mA.num_rows(); r++)
		for(size_t c=0; c < mB.num_cols(); c++)
		{
			AssignMult(dest(r, c), mA(r, 0), mB(0, c));
			for(size_t k=1; k < mB.num_rows(); k++)
				AddMult(dest(r, c), mA(r, k), mB(k, c));
		}
}

// dest += b*vec
template<typename A, typename B, typename C>
inline void AddMult(DenseVector<A> &dest, const DenseMatrix<B> &mat, const DenseVector<C> &vec)
{
	UG_ASSERT(mat.num_cols() == vec.size(), "");
	dest.resize(mat.num_rows());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		for(size_t c=0; c < mat.num_cols(); c++)
			AddMult(dest(r), mat(r, c), vec(c));
	}
}


// dest += b*vec
template<typename A, typename B, typename C>
inline void AddMult(DenseMatrix<A> &dest, const DenseMatrix<B> &mA, const DenseMatrix<C> &mB)
{
	UG_ASSERT(mA.num_cols() == mB.num_rows(), "");
	UG_ASSERT(dest.num_rows()==mA.num_rows() && dest.num_cols()==mB.num_cols(), "");
	for(size_t r=0; r < mA.num_rows(); r++)
		for(size_t c=0; c < mB.num_cols(); c++)
		{
			for(size_t k=0; k < mB.num_rows(); k++)
				AddMult(dest(r, c), mA(r, k), mB(k, c));
		}
}



// dest -= b*vec
template<typename A, typename B, typename C>
inline void SubMult(DenseVector<A> &dest, const DenseMatrix<B> &mat, const DenseVector<C> &vec)
{
	UG_ASSERT(mat.num_cols() == vec.size(), "");
	dest.resize(mat.num_rows());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		for(size_t c=0; c < mat.num_cols(); c++)
			SubMult(dest(r), mat(r, c), vec(c));
	}
}


// mat = alpha * mat

// dest = b*vec
template<typename A, typename B>
inline void AssignMult(DenseMatrix<A> &dest, const double &alpha, const DenseMatrix<B> &mat)
{
	dest.resize(mat.num_rows(), mat.num_cols());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		for(size_t c=0; c < mat.num_cols(); c++)
			AssignMult(dest(r, c), alpha, mat(r, c));
	}
}

template<typename A, typename B>
inline void AddMult(DenseMatrix<A> &dest, const double &alpha, const DenseMatrix<B> &mat)
{
	dest.resize(mat.num_rows(), mat.num_cols());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		for(size_t c=0; c < mat.num_cols(); c++)
			AddMult(dest(r, c), alpha, mat(r, c));
	}
}

template<typename A, typename B>
inline void SubMult(DenseMatrix<A> &dest, const double &alpha, const DenseMatrix<B> &mat)
{
	dest.resize(mat.num_rows(), mat.num_cols());
	for(size_t r=0; r < mat.num_rows(); r++)
	{
		for(size_t c=0; c < mat.num_cols(); c++)
			SubMult(dest(r, c), alpha, mat(r, c));
	}
}


// VECTORs

// dest = b*vec
template<typename A, typename B>
inline void AssignMult(DenseVector<A> &dest, const double &b, const DenseVector<B> &vec)
{
	dest.resize(vec.size());
	for(size_t i=0; i<vec.size(); i++)
		AssignMult(dest[i], b, vec[i]);
}

// dest += b*vec
template<typename A, typename B>
inline void AddMult(DenseVector<A> &dest, const double &b, const A &vec)
{
	dest.resize(vec.size());
	for(size_t i=0; i<vec.size(); i++)
		AddMult(dest[i], b, vec[i]);
}

// dest -= b*vec
template<typename A, typename B>
inline void SubMult(DenseVector<A> &dest, const double &b, const DenseVector<B> &vec)
{
	dest.resize(vec.size());
	for(size_t i=0; i<vec.size(); i++)
		SubMult(dest[i], b, vec[i]);
}

// dest = vec*b
template<typename A> inline void AssignMult(A &dest, const A &vec, const double &b)
{
	AssignMult(dest, b, vec);
}
// dest += vec*b
template<typename A> inline void AddMult(A &dest, const A &vec, const double &b)
{
	AddMult(dest, b, vec);
}
// dest -= vec*b
template<typename A> inline void SubMult(A &dest, const A &vec, const double &b)
{
	SubMult(dest, b, vec);
}


template<typename T>
inline void SetSize(DenseMatrix<T> &t, size_t a, size_t b)
{
	t.resize(a, b);
}

//setSize(t, a) for vectors
template<typename T>
inline void SetSize(DenseVector<T> &t, size_t a)
{
	t.resize(a);
}

// getSize
template<typename T>
inline size_t GetSize(const DenseVector<T> &t)
{
	return t.size();
}

//getRows
template<typename T>
inline size_t GetRows(const DenseMatrix<T> &t)
{
	return t.num_rows();
}

template<typename T>
inline size_t GetCols(const DenseMatrix<T> &t)
{
	return t.num_cols();
}

template<typename T>
struct block_traits;

template<typename T>
struct block_traits< DenseMatrix<T> >
{
	enum { ordering = DenseMatrix<T>::ordering };
	enum { is_static = DenseMatrix<T>::is_static};
	enum { static_num_rows = DenseMatrix<T>::static_num_rows};
	enum { static_num_cols = DenseMatrix<T>::static_num_cols};

	// todo: to be implemented
	// using inverse_type = DenseMatrixInverse;
};

template<typename T>
struct block_traits< DenseVector<T> >
{
	enum { is_static = DenseVector<T>::is_static};
	enum { static_size = DenseVector<T>::static_size};
};

template<typename T1, typename T2>
struct block_multiply_traits;

template<typename T1, typename T2>
struct block_multiply_traits<DenseMatrix<T1>, DenseVector<T2> >
{
	using ReturnType = DenseVector<T2>;
};

//////////////////////////////////////////////////////////////////////////////////////////////
// block_traits

template<typename T> class DenseMatrixInverse;

//////////////////////////////////////////////////////////////////////////////////////////////
// variable matrix
template<eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< VariableArray2<number, TOrdering> > >
{
	enum { ordering = TOrdering };
	enum { is_static = false};
	enum { static_num_rows = 0};
	enum { static_num_cols = 0};
	enum { depth = 1 };

	using inverse_type = DenseMatrixInverse< VariableArray2<number, TOrdering> >;
};
//////////////////////////////////////////////////////////////////////////////////////////////
// fixed matrix
template<size_t TBlockSize, eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< FixedArray2<number, TBlockSize, TBlockSize, TOrdering> > >
{
	enum { ordering = TOrdering };
	enum { is_static = true};
	enum { static_num_rows = TBlockSize};
	enum { static_num_cols = TBlockSize};
	enum { depth = 1 };

	using inverse_type = DenseMatrixInverse< FixedArray2<number, TBlockSize, TBlockSize, TOrdering>  >;
};

//////////////////////////////////////////////////////////////////////////////////////////////
// variable block matrix
template<typename TValue, eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< VariableArray2<TValue, TOrdering> > >
{
	enum { ordering = TOrdering };
	enum { is_static = false};
	enum { static_num_rows = 0};
	enum { static_num_cols = 0};
	enum { depth = block_traits<TValue>::depth+1 };

	using inverse_type = DenseMatrixInverse< VariableArray2<number, TOrdering> >;
};

template<typename TValue, size_t TBlockSize, eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< FixedArray2<TValue, TBlockSize, TBlockSize, TOrdering> > >
{
	enum { ordering = TOrdering };
	enum { is_static = false};
	enum { static_num_rows = 0};
	enum { static_num_cols = 0};
	enum { depth = block_traits<TValue>::depth+1 };

	using inverse_type = DenseMatrixInverse< VariableArray2<number, TOrdering> >;
};

//////////////////////////////////////////////////////////////////////////////////////////////
// fixed 1x1 to 3x3 : inverse is matrix
template<eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< FixedArray2<number, 1, 1, TOrdering> > >
{
	enum { ordering = DenseMatrix< FixedArray2<number, 1, 1> >::ordering };
	enum { is_static = true};
	enum { static_num_rows = 1};
	enum { static_num_cols = 1};

	using inverse_type = DenseMatrix< FixedArray2<number, 1, 1, TOrdering> >;
};

template<eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< FixedArray2<number, 2, 2, TOrdering> > >
{
	enum { ordering = DenseMatrix< FixedArray2<number, 2, 2> >::ordering };
	enum { is_static = true};
	enum { static_num_rows = 2};
	enum { static_num_cols = 2};

	using inverse_type = DenseMatrix< FixedArray2<number, 2, 2, TOrdering> >;
};

template<eMatrixOrdering TOrdering>
struct block_traits< DenseMatrix< FixedArray2<number, 3, 3, TOrdering> > >
{
	enum { ordering = DenseMatrix< FixedArray2<number, 3, 3> >::ordering };
	enum { is_static = true};
	enum { static_num_rows = 3};
	enum { static_num_cols = 3};

	using inverse_type = DenseMatrix< FixedArray2<number, 3, 3, TOrdering> >;
};


template<typename T> struct block_multiply_traits<DenseMatrix<T>, DenseMatrix<T> >
{
	using ReturnType = DenseMatrix<T>;
};


//////////////////////////////////////////////////////////////////////////////////////////////




} // namespace ug

#endif