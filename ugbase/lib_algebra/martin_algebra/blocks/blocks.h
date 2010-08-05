/*
 *  blocks.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  Header File for general block matrix ccessing 
 *  i.e. setAt(mat, i, j, d) -> mat(i,j) = d
 *	and setAt(f, 0, 0, d) -> f = d.
 *  This means doubles and block matrices can be accessed
 *	by the same methods.
 */

#ifndef __H__UG__MARTIN_ALGEBRA__BLOCKS__
#define __H__UG__MARTIN_ALGEBRA__BLOCKS__

namespace ug{
	
inline double dabs(double a) { return a > 0 ? a : -a; }

//////////////////////////////////////////////////////
/*
template<size_t n>
inline double mnorm2(const  blockVector<n> &v)
{
	return v.norm2();
}*/

template<typename TYPE>
inline double BlockNorm2(const TYPE &v)
{
	return v.norm2();
}

template<typename TYPE>
inline double BlockNorm(const TYPE &v)
{
	return sqrt(v.norm2());
}


//////////////////////////////////////////////////////

// get/set vector
template<typename M> inline double &BlockRef(M &m, size_t i)
{
	return m(i);
}

template<typename M> inline const double &BlockRef(const M &m, size_t i)
{
	return m(i);
}

// get/set matrix
template<typename M> inline double &BlockRef(M &m, size_t i, size_t j)
{
	return m(i, j);
}

template<typename M> inline const double &BlockRef(const M &m, size_t i, size_t j)
{
	return m(i, j);
}

//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables 

	
// MATRICES

// todo: replace add_mult etc. with template expressions
// dest = b*vec
template<typename A, typename B, typename C> inline void AssignMult(A &dest, const B &b, const C &vec)
{
	b.assign_mult(dest, vec);
}

// dest += b*vec
template<typename A, typename B, typename C> inline void AddMult(A &dest, const B &b, const C &vec)
{
	b.add_mult(dest, vec);
}

// dest -= b*vec
template<typename A, typename B, typename C> inline void SubMult(A &dest, const B &b, const C &vec)
{
	b.sub_mult(dest, vec);
}


// VECTORs

// dest = b*vec
template<typename A> inline void AssignMult(A &dest, const double &b, const A &vec)
{
	dest.assign_mult(b, vec);
}
// dest += b*vec
template<typename A> inline void AddMult(A &dest, const double &b, const A &vec)
{
	dest.add_mult(b, vec);
}
// dest -= b*vec
template<typename A> inline void SubMult(A &dest, const double &b, const A &vec)
{
	dest.sub_mult(b, vec);
}

// dest = vec*b
template<typename A> inline void AssignMult(A &dest, const A &vec, const double &b)
{
	dest.assign_mult(b, vec);
}
// dest += vec*b
template<typename A> inline void AddMult(A &dest, const A &vec, const double &b)
{
	dest.add_mult(b, vec);
}
// dest -= vec*b
template<typename A> inline void SubMult(A &dest, const A &vec, const double &b)
{
	dest.sub_mult(b, vec);
}

/*inline void AssignMult(double &dest, const double &b, const double &vec)
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
}*/

	

//////////////////////////////////////////////////////
//setSize(t, a, b) for matrices
template<typename T>
inline void SetSize(T &t, size_t a, size_t b)
{
	t.resize(a, b);
}

//setSize(t, a) for vectors
template<typename T>
inline void SetSize(T &t, size_t a)
{
	t.resize(a);
}

// getSize
template<typename T>
inline size_t GetSize(T &t)
{
	return t.size();
}

//getRows
template<typename T>
inline size_t GetRows(const T &t)
{
	return t.num_rows();
}

template<typename T>
inline size_t GetCols(const T &t)
{
	return t.num_cols();
}

} // namespace ug

#include "double.h"
#endif
