/*
 *  matrixRow.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__MATRIXROW_IMPL__
#define __H__UG__MARTIN_ALGEBRA__MATRIXROW_IMPL__

#ifndef FLEXAMG
#include "algebra_misc.h"
#endif

namespace ug{
// operator []
//-------------
//! access connection nr i
//! @param	i	connection nr i 
//! @return A[row] * x 
template<typename entry_type>
inline const typename matrixrow<entry_type>::connection &matrixrow<entry_type>::operator [] (size_t i) const
{
	UG_ASSERT(i < A.getNrOfConnections(row) && i >= 0, *this << " has no connection nr. " << i);
	return A.pRowStart[row][i];
}

/*
template<typename entry_type>
inline size_t matrixrow<entry_type>::getConNr(size_t index) const
{
	for(size_t i=0; i< getNrOfConnections(); i++)
	{
		if(A.cons[row][i].iIndex == index)
			return i;
	}
	return -1;
}*/

 //#pragma mark -
////////////////////////////////////////////////////////////////////////////////

// operator *
//-------------
//! constructs temporary variable
//! @param	x
//! @return A[row] * x 
template<typename entry_type>
template<typename vec_type>
inline vec_type matrixrow<entry_type>::operator *(const Vector<vec_type> &x) const
{
	vec_type d=0;
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
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::assign_mult(vec_type &d, const Vector<vec_type> &x) const
{
	d = 0.0;
	add_mult(d, x);
}

// sub_mult
//-------------
//! @param	d	out: d -= A[row] * x
//! @param	x	in
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::sub_mult(vec_type &d, const Vector<vec_type> &x) const
{
	for(cRowIterator it = beginRow(); !it.isEnd(); ++it)
		SubMult(d, (*it).dValue, x[(*it).iIndex]);
}


// add_mult
//-------------
//! @param	d	out: d += A[row] * x
//! @param	x	in
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::add_mult(vec_type &d, const Vector<vec_type> &x) const
{
	for(cRowIterator it = beginRow(); !it.isEnd(); ++it)
		AddMult(d, (*it).dValue, x[(*it).iIndex]);
}


//#pragma mark -
////////////////////////////////////////////////////////////////////////////////

template<typename entry_type, typename vec_type>
inline void multiplyCopyTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.assign_mult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplyAddTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.add_mult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplySubstractFrom(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.sub_mult(d, x);
}

} // namespace ug

#endif
