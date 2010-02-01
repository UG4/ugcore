
/*
 *  matrixRow.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */


// operator []
//-------------
//! access connection nr i
//! @param	i	connection nr i 
//! @return A[row] * x 
template<typename entry_type>
inline const typename matrixrow<entry_type>::connection &matrixrow<entry_type>::operator [] (int i) const
{
	ASSERT2(i < A.iNrOfConnections[row] && i >= 0, *this << " has no connection nr. " << i);
	/*ASSERT(A.cons[row]+i < A.consmem+A.iMaxTotalNrOfConnections
	 && A.cons[row]+i >= A.consmem);*/
	//ASSERT(A.consmem.isMemInChunk(A.cons[row][i]));		
	return A.cons[row][i];
}



template<typename entry_type>
inline int matrixrow<entry_type>::getConNr(int index) const
{
	for(int i=0; i< getNrOfConnections(); i++)
	{
		if(A.cons[row][i].iIndex == index)
			return i;
	}
	return -1;
}

#pragma mark -
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
	vec_type d;
	cRowIterator it = beginRow();
	
	if(row >= x.getLength()) ++it; // skip diag		
	if(it.isEnd()) return 0.0;
	
	d = (*it).dValue * x[(*it).iIndex];
	++it;	
	for(; !it.isEnd(); ++it)
		d += (*it).dValue * x[(*it).iIndex];
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
inline void matrixrow<entry_type>::copyToMult(vec_type &d, const Vector<vec_type> &x) const
{
	d = 0.0;
	addToMult(d, x);
}

// substractFromMult
//-------------
//! @param	d	out: d -= A[row] * x
//! @param	x	in
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::substractFromMult(vec_type &d, const Vector<vec_type> &x) const
{
	cRowIterator it = beginRow();
	if(row >= x.getLength()) ++it; // skip diag

	for(; !it.isEnd(); ++it)
		d -= (*it).dValue * x[(*it).iIndex];
}

// addToMult
//-------------
//! @param	d	out: d += A[row] * x
//! @param	x	in
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::addToMult(vec_type &d, const Vector<vec_type> &x) const
{
	cRowIterator it = beginRow();
	if(row >= x.getLength()) ++it; // skip diag
	
	for(; !it.isEnd(); ++it)
		d += (*it).dValue * x[(*it).iIndex];
}


#pragma mark -
////////////////////////////////////////////////////////////////////////////////

template<typename entry_type, typename vec_type>
inline void multiplyCopyTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.copyToMult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplyAddTo(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.addToMult(d, x);
}

template<typename entry_type, typename vec_type>
inline void multiplySubstractFrom(vec_type &d, const matrixrow<entry_type> &r, const Vector<vec_type> &x)
{
	r.substractFromMult(d, x);
}

