/*
 * ScalarVectorAdapter.hh
 *
 *  Created on: Jun 10, 2013
 *      Author: anaegel
 */

#ifndef SCALAR_VECTOR_ADAPTER_HH_
#define SCALAR_VECTOR_ADAPTER_HH_

#include <cstdlib>

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/lib_algebra_impl.h"
#include "lib_algebra/cpu_algebra_types.h"


using namespace ug;

// provides an interface for matrix of algebra type B for matrices originally of algebra type A
// allows to access a CPUBlockAlgebra (AT) as a scalar CPUAlgebra (ST)

template<class AT, class ST>
class ScalarVectorAdapter{

public:
	typedef typename AT::vector_type encapsulated_vector_type;
	typedef typename ST::vector_type::value_type value_type;
	static const int blockSize = AT::blockSize;

	//typedef typename ST::vector_type::const_row_iterator const_row_iterator;

	ScalarVectorAdapter(encapsulated_vector_type& vec) : m_src(vec) {};

	inline value_type &operator [] (size_t i)
	{ return BlockRef(m_src[i/blockSize], i%blockSize);}

	inline const value_type &operator [] (size_t i) const { return (i);}

	void resize(size_t newSize, bool bCopyValues=true)
	{
		m_src.resize_exactly(newSize/blockSize, bCopyValues);
	}
	void reserve(size_t newCapacity, bool bCopyValues=true)
	{
		m_src.reserve_exactly(newCapacity/blockSize, bCopyValues);
	}
	void print(const char * const text = NULL) const
	{
		m_src.print(text);
	}
private:
	encapsulated_vector_type &m_src;
};


// partielle Spezialisierung fuer Block Algebra 1,2,3
template<>
ScalarVectorAdapter<CPUAlgebra, CPUAlgebra>::value_type&
ScalarVectorAdapter<CPUAlgebra, CPUAlgebra>::operator [] (size_t i)
{ return m_src[i]; }


#endif /* SPARSEMATRIXPROXY_HH_ */
