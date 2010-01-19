
/*
 *  matrixRow.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 18.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

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

// the != 0.0 is bad but we need this for restriction, since A.cons[i][0].iIndex = i.
template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::copyToMult(vec_type &d, const Vector<vec_type> &x) const
{
	d = 0.0;
	addToMult(d, x);
}

template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::substractFromMult(vec_type &d, const Vector<vec_type> &x) const
{
	cRowIterator it = beginRow();
	if(!it.isEnd() && (*it).dValue == 0.0) 
		++it;
	for(; !it.isEnd(); ++it)
		d -= (*it).dValue * x[(*it).iIndex];
}

template<typename entry_type>
template<typename vec_type>
inline void matrixrow<entry_type>::addToMult(vec_type &d, const Vector<vec_type> &x) const
{
	cRowIterator it = beginRow();
	if(!it.isEnd() && (*it).dValue == 0.0) 
		++it;
	
	for(; !it.isEnd(); ++it)
		d += (*it).dValue * x[(*it).iIndex];
}

template<typename entry_type>
template<typename vec_type>
inline vec_type matrixrow<entry_type>::operator *(const Vector<vec_type> &x) const
{
	vec_type d;
	cRowIterator it = beginRow();
	if(!it.isEnd() && (*it).dValue == 0.0) 
		++it;
	
	for(; !it.isEnd(); ++it)
	{
		// otherwise we dont know how big d is.
		if((*it).dValue != 0.0)
		{
			d = (*it).dValue * x[(*it).iIndex];
			++it;
			break;
		}
	}
	
	for(; !it.isEnd(); ++it)
		d += (*it).dValue * x[(*it).iIndex];
	return d;
}


template<typename entry_type>
inline const typename matrixrow<entry_type>::connection &matrixrow<entry_type>::operator [] (int i) const
{
	ASSERT2(i < A.iNrOfConnections[row] && i >= 0, *this << " has no connection nr. " << i);
	/*ASSERT(A.cons[row]+i < A.consmem+A.iMaxTotalNrOfConnections
	 && A.cons[row]+i >= A.consmem);*/
	//ASSERT(A.consmem.isMemInChunk(A.cons[row][i]));		
	return A.cons[row][i];
}

