/*
 *  matrixRow.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__CPU_ALGEBRA__MATRIXROW_IMPL__
#define __H__UG__CPU_ALGEBRA__MATRIXROW_IMPL__


#include "algebra_misc.h"

namespace ug{

//#pragma mark -
////////////////////////////////////////////////////////////////////////////////

// operator *
//-------------
//! constructs temporary variable
//! @param	x
//! @return A[row] * x 
template<typename value_type>
template<typename vec_type>
inline vec_type matrixrow<value_type>::operator *(const Vector<vec_type> &x) const
{
	vec_type d=0.0;
	for(cRowIterator it = beginRow(); !it.isEnd(); ++it)
		AddMult(d, (*it).dValue, x[(*it).iIndex]);
	return d;
}

// the following functions are necessary when one will prevent
// temporary variables in x[i] = A[i] * y or x[i] += A[i] *y.

// copyToMult
//-------------
//! @param	d	out: d = A[row] * x
//! @param	x	in
//! @return d = A[row] * x 
template<typename value_type>
template<typename vec_type>
inline void matrixrow<value_type>::assign_mult(vec_type &d, const Vector<vec_type> &x) const
{
	d = 0.0;
	add_mult(d, x);
}

// sub_mult
//-------------
//! @param	d	out: d -= A[row] * x
//! @param	x	in
template<typename value_type>
template<typename vec_type>
inline void matrixrow<value_type>::sub_mult(vec_type &d, const Vector<vec_type> &x) const
{
	for(cRowIterator it = beginRow(); !it.isEnd(); ++it)
		SubMult(d, (*it).dValue, x[(*it).iIndex]);
}


// add_mult
//-------------
//! @param	d	out: d += A[row] * x
//! @param	x	in
template<typename value_type>
template<typename vec_type>
inline void matrixrow<value_type>::add_mult(vec_type &d, const Vector<vec_type> &x) const
{
	for(cRowIterator it = beginRow(); !it.isEnd(); ++it)
		AddMult(d, (*it).dValue, x[(*it).iIndex]);
}


//#pragma mark -
////////////////////////////////////////////////////////////////////////////////

template<typename value_type, typename vec_type>
inline void multiplyCopyTo(vec_type &d, const matrixrow<value_type> &r, const Vector<vec_type> &x)
{
	r.assign_mult(d, x);
}

template<typename value_type, typename vec_type>
inline void multiplyAddTo(vec_type &d, const matrixrow<value_type> &r, const Vector<vec_type> &x)
{
	r.add_mult(d, x);
}

template<typename value_type, typename vec_type>
inline void multiplySubstractFrom(vec_type &d, const matrixrow<value_type> &r, const Vector<vec_type> &x)
{
	r.sub_mult(d, x);
}

} // namespace ug

#endif
