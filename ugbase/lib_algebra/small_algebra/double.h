/*
 *  blockMatrix.h
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  Header File for general block matrix / double accessing 
 *  i.e. setAt(mat, i, j, d) -> mat(i,j) = d
 *	and setAt(f, 0, 0, d) -> f = d.
 *  This means doubles and block matrices can be accessed
 *	by the same methods.
 */

// todo: also with float / complex<float> / complex<double>

#ifndef __H__UG__SMALL_ALGEBRA__DOUBLE__
#define __H__UG__SMALL_ALGEBRA__DOUBLE__

#include "blocks.h"
#include "common/common.h"

namespace ug{


//////////////////////////////////////////////////////
template<typename T> number BlockNorm(const T &t);
template <>
inline number BlockNorm(const number &a)
{
	return a>0 ? a : -a;
}

template<typename T> number BlockNorm2(const T &t);
template <>
inline number BlockNorm2(const number &a)
{
	return a*a;
}

//////////////////////////////////////////////////////
// get/set specialization for numbers

template<> inline number &BlockRef(number &m, size_t i)
{
	UG_ASSERT(i == 0, "block is number, doesnt have component (" << i << ").");
	return m;
}
template<> inline const number &BlockRef(const number &m, size_t i)
{
	UG_ASSERT(i == 0, "block is number, doesnt have component (" << i << ").");
	return m;
}

template<> inline number &BlockRef(number &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is number, doesnt have component (" << i << ", " << j << ").");
	return m;
}
template<> inline const number &BlockRef(const number &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is number, doesnt have component (" << i << ", " << j << ").");
	return m;
}

//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables 

inline void AssignMult(number &dest, const number &b, const number &vec)
{
	dest = b*vec;
}
// dest += vec*b
inline void AddMult(number &dest, const number &b, const number &vec)
{
	dest += b*vec;
}


// dest -= vec*b
inline void SubMult(number &dest, const number &b, const number &vec)
{
	dest -= b*vec;
}


//////////////////////////////////////////////////////
//setSize(t, a, b) for numbers
template<>
inline void SetSize(number &d, size_t a)
{
	UG_ASSERT(a == 1, "block is number, cannot change size to " << a << ".");
	return;
}

template<>
inline void SetSize(number &d, size_t a, size_t b)
{
	UG_ASSERT(a == 1 && b == 1, "block is number, cannot change size to (" << a << ", " << b << ").");
	return;
}

template<>
inline size_t GetSize(const number &t)
{
	return 1;
}

template<>
inline size_t GetRows(const number &t)
{
	return 1;
}

template<>
inline size_t GetCols(const number &t)
{
	return 1;
}
///////////////////////////////////////////////////////////////////

inline bool InverseMatMult(number &dest, const double &beta, const number &mat, const number &vec)
{
	dest = beta*vec/mat;
	return true;
}

///////////////////////////////////////////////////////////////////
// traits: information for numbers


template<>
struct block_traits<number>
{
	typedef number vec_type;
	typedef number inverse_type;

	enum { is_static = true};
	enum { static_num_rows = 1};
	enum { static_num_cols = 1};
	enum { static_size = 1 };
	enum { depth = 0 };
};

template<> struct block_multiply_traits<number, number>
{
	typedef number ReturnType;
};

template<typename T>
struct block_multiply_traits<number, T>
{
	typedef T ReturnType;
};

template<typename T>
struct block_multiply_traits<T, number>
{
	typedef T ReturnType;
};

inline bool GetInverse(number &inv, const number &m)
{
	inv = 1.0/m;
	return (m != 0.0);
}

inline bool Invert(number &m)
{
	bool b = (m != 0.0);
	m = 1/m;
	return b;
}

} // namespace ug

#endif



