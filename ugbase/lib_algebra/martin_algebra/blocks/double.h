/*
 *  blockMatrix.h
 *  flexamg
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

#ifndef __H__UG__MARTIN_ALGEBRA__DOUBLE__
#define __H__UG__MARTIN_ALGEBRA__DOUBLE__
#include "blocks.h"
namespace ug{

template <typename t> struct block_matrix_traits;
template <typename t> struct block_vector_traits;
template<typename entry_type, typename vec_type> struct block_multiply_traits;


//////////////////////////////////////////////////////
template<typename T> double BlockNorm(const T &t);
template <>
inline double BlockNorm(const double &a)
{
	return a>0 ? a : -a;
}

template<typename T> double BlockNorm2(const T &t);
template <>
inline double BlockNorm2(const double &a)
{
	return a*a;
}

//////////////////////////////////////////////////////
// get/set specialization for doubles

template<> inline double &BlockRef(double &m, size_t i)
{
	UG_ASSERT(i == 0, "block is double, doesnt have component (" << i << ").");
	return m;
}
template<> inline const double &BlockRef(const double &m, size_t i)
{
	UG_ASSERT(i == 0, "block is double, doesnt have component (" << i << ").");
	return m;
}

template<> inline double &BlockRef(double &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is double, doesnt have component (" << i << ", " << j << ").");
	return m;
}
template<> inline const double &BlockRef(const double &m, size_t i, size_t j)
{
	UG_ASSERT(i == 0 && j == 0, "block is double, doesnt have component (" << i << ", " << j << ").");
	return m;
}

//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables

// with double, it remains the question if this is really fast,
// because we are writing into memory with double &dest (do we???)
// temporary variables would be faster, does gcc know???
// dest = vec*b
inline void AssignMult(double &dest, const double &b, const double &vec)
{
	dest = b*vec;
}
// dest += vec*b
inline void AddMult(double &dest, const double &b, const double &vec)
{
	dest += b*vec;
}


// dest -= vec*b
inline void SubMult(double &dest, const double &b, const double &vec)
{
	dest -= b*vec;
}


//////////////////////////////////////////////////////
//setSize(t, a, b) for doubles
template<>
inline void SetSize(double &d, size_t a)
{
	UG_ASSERT(a == 1, "block is double, cannot change size to " << a << ".");
	return;
}

template<>
inline void SetSize(double &d, size_t a, size_t b)
{
	UG_ASSERT(a == 1 && b == 1, "block is double, cannot change size to (" << a << ", " << b << ").");
	return;
}

template<>
inline size_t GetSize(const double &t)
{
	return 1;
}

template<>
inline size_t GetRows(const double &t)
{
	return 1;
}

template<>
inline size_t GetCols(const double &t)
{
	return 1;
}

///////////////////////////////////////////////////////////////////
// traits: information for doubles



template<>
struct block_matrix_traits<double>
{
	typedef double vec_type;
	typedef double inverse_type;
	enum { nrOfUnknowns = 1 } ;

	//TODO: Decide weather to use nrOfUnknows or nrOfRows/nrOfCols. Currently mixed. Therefore I added the two below here. Andreas Vogel
	enum { nrOfRows = 1 } ;
	enum { nrOfCols = 1 } ;
};

template<>
struct block_vector_traits<double>
{
	enum { nrOfUnknowns = 1 };
};

template<> struct block_multiply_traits<double, double>
{
	typedef double ReturnType;
};

template<typename M>
inline void GetInverse(typename block_matrix_traits<M>::inverse_type &inv, const M &m);

template<>
inline void GetInverse(double &inv, const double &m)
{
	inv = 1.0/m;
}

template<>
inline void Invert(double &m)
{
	m = 1/m;
}

} // namespace ug

#endif



