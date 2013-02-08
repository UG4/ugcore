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

	//! clones the vector (deep-copy) including values
	SmartPtr<vector_type> clone() const;

	//! clones the vector (deep-copy) excluding values
	SmartPtr<vector_type> clone_without_values() const;

protected:
	//! virtual clone using covariant return type
	virtual vector_type* virtual_clone() const;

	//! virtual clone using covariant return type excluding values
	virtual vector_type* virtual_clone_without_values() const;
};

template<typename value_type>
CRSVector<value_type>* CRSVector<value_type>::virtual_clone() const
{
	return new CRSVector<value_type>(*this);
}

template<typename value_type>
SmartPtr<CRSVector<value_type> > CRSVector<value_type>::clone() const
{
	return SmartPtr<CRSVector<value_type> >(this->virtual_clone());
}

template<typename value_type>
CRSVector<value_type>* CRSVector<value_type>::virtual_clone_without_values() const
{
	return new CRSVector<value_type>(this->size());
}

template<typename value_type>
SmartPtr<CRSVector<value_type> > CRSVector<value_type>::clone_without_values() const
{
	return SmartPtr<CRSVector<value_type> >(this->virtual_clone_without_values());
}

// @}

} // namespace ug

#endif
