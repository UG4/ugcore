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
template<typename T> double mnorm(const T &t);

template <>
inline double mnorm(const double &a)
{
	return a>0 ? a : -a;
}

template<typename T> double mnorm2(const T &t);
template <>
inline double mnorm2(const double &a)
{
	return a*a;
}

//////////////////////////////////////////////////////
// get/set specialization for doubles

template<> inline double getAt(const double &m, int i) { return m; }
template<> inline double setAt(double &m, int i, double a) { m = a; return a; }

template<> inline double getAt(const double &m, int i, int j) { return m; }
template<> inline double setAt(double &m, int i, int j, double a) { m = a; return a; }


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
inline int getSize(const double &t)
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
struct block_matrix_traits<double>
{
	typedef double vec_type;
	typedef double inverse_type;
	enum { nrOfUnknowns = 1 } ;
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


} // namespace ug

#endif



