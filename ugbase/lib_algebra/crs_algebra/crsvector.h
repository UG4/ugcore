/*
 *  CRSVector.h
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__CRS_ALGEBRA__VECTOR__
#define __H__UG__CRS_ALGEBRA__VECTOR__

#include "../cpu_algebra/vector.h"

namespace ug{
///////////////////////////////////////////////////////////////////
//							CRSVector
///////////////////////////////////////////////////////////////////

/// \addtogroup lib_algebra
///	@{

//!
template <typename TValueType>
class CRSVector : public Vector<TValueType>
{
public:
	typedef TValueType value_type;
	typedef CRSVector<TValueType> vector_type;

	//! constructor
	CRSVector() : Vector<TValueType>() {}

	//! constructor with length
	CRSVector(size_t _length) : Vector<TValueType>(_length) {}

};

// @}

} // namespace ug

#endif
