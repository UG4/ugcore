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

#pragma once

namespace ug{
	
template <typename t> class matrix_trait;
template <typename t> class vec_traits;
template<typename entry_type, typename vec_type> struct Mult_Traits;

inline double dabs(double a) { return a > 0 ? a : -a; }

//////////////////////////////////////////////////////
/*
template<int n>
inline double mnorm2(const  blockVector<n> &v)
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

//////////////////////////////////////////////////////

// get/set vector
template<typename M> inline double getAt(const M &m, int i)
{
	return m(i);
}
template<typename M> inline double setAt(M &m, int i, double a)
{
	return m(i) = a;
}

// get/set matrix
template<typename M> inline double getAt(const M &m, int i, int j)
{
	return m(i, j);
}

template<typename M> inline double setAt(M &m, int i, int j, double a)
{
	return m(i, j) = a;
}

//////////////////////////////////////////////////////
// get/set specialization for doubles

template<> inline double getAt(const double &m, int i) { return m; }
template<> inline double setAt(double &m, int i, double a) { m = a; return a; }

template<> inline double getAt(const double &m, int i, int j) { return m; }
template<> inline double setAt(double &m, int i, int j, double a) { m = a; return a; }



//////////////////////////////////////////////////////
// algebra stuff to avoid temporary variables 

// todo: replace add_mult etc. with template expressions
// dest = b*vec
template<typename A, typename B, typename C> inline void assign_mult(A &dest, const B &b, const C &vec)
{
	b.assign_mult(dest, vec);
}

// dest += b*vec
template<typename A, typename B, typename C> inline void add_mult(A &dest, const B &b, const C &vec)
{
	b.add_mult(dest, vec);
}

// dest -= b*vec
template<typename A, typename B, typename C> inline void sub_mult(A &dest, const B &b, const C &vec)
{
	b.sub_mult(dest, vec);
}

// dest = b*vec
template<typename A> inline void assign_mult(A &dest, const double &b, const A &vec)
{
	dest.assign_mult(b, vec);
}
// dest += b*vec
template<typename A> inline void add_mult(A &dest, const double &b, const A &vec)
{
	dest.add_mult(b, vec);
}
// dest -= b*vec
template<typename A> inline void sub_mult(A &dest, const double &b, const A &vec)
{
	dest.sub_mult(b, vec);
}

// dest = vec*b
template<typename A> inline void assign_mult(A &dest, const A &vec, const double &b)
{
	dest.assign_mult(b, vec);
}
// dest += vec*b
template<typename A> inline void add_mult(A &dest, const A &vec, const double &b)
{
	dest.add_mult(b, vec);
}
// dest -= vec*b
template<typename A> inline void sub_mult(A &dest, const A &vec, const double &b)
{
	dest.sub_mult(b, vec);
}

// with double, it remains the question if this is really fast,
// because we are writing into memory with double &dest (do we???) 
// temporary variables would be faster, does gcc know???
// dest = vec*b
inline void assign_mult(double &dest, const double &b, const double &vec)
{
	dest = b*vec;
}
// dest += vec*b
inline void add_mult(double &dest, const double &b, const double &vec)
{
	dest += b*vec;
}
// dest -= vec*b
inline void sub_mult(double &dest, const double &b, const double &vec)
{
	dest -= b*vec;
}

//////////////////////////////////////////////////////
//setSize(t, a, b) for matrices
template<typename T>
inline void setSize(T &t, int a, int b)
{
	t.setSize(a, b);
}

//setSize(t, a) for vectors
template<typename T>
inline void setSize(T &t, int a)
{
	t.setSize(a);
}

// getSize
template<typename T>
inline int getSize(T &t)
{
	return t.getSize();
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
//////////////////////////////////////////////////////
// wrapper for using doubles
template<>
inline void setSize(double &d, int a)
{
	return;
}
template<>
inline void setSize(double &d, int a, int b)
{
	return;
}
template<>
inline int getSize(double &t)
{
	return 1;
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

///////////////////////////////////////////////////////////////////
// traits: information for doubles
template<>
struct matrix_trait<double>
{
	typedef double vec_type;
	typedef double inverse_type;
	enum { nrOfUnknowns = 1 } ;
};

template<>
struct vec_traits<double>
{
	enum { nrOfUnknowns = 1 };
};

template<> struct Mult_Traits<double, double>
{
	typedef double ReturnType;
};


} // namespace ug

#include "blockDenseMatrix.h"
#include "blockVector.h"
