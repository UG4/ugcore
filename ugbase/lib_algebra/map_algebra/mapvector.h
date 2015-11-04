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
class MapVector : public Vector<TValueType>
{
public:
	typedef TValueType value_type;
	typedef MapVector<TValueType> vector_type;

	//! constructor
	MapVector() : Vector<TValueType>() {}

	//! constructor with length
	MapVector(size_t _length) : Vector<TValueType>(_length) {}

};

// @}

} // namespace ug

#endif
