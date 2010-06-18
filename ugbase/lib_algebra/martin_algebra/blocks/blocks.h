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

} // namespace ug

#endif
