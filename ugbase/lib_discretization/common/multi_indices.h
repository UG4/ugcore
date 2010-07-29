/*
 * multi_indices.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__MULTI_INDICES__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__MULTI_INDICES__

#include <vector>
#include <iostream>
#include <assert.h>

namespace ug{

template<int N, typename TSingleIndexType = size_t>
class MultiIndex
{
	public:
		typedef TSingleIndexType single_index_type;

	public:
		/// number of indices in multi index
		inline size_t size() const {return N;}

		/// access to index component
		inline single_index_type& operator[] (size_t i)
		{
			UG_ASSERT(i < N, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator[] (size_t i) const
		{
			UG_ASSERT(i < N, "Index invalid");
			return m_indices[i];
		}

	private:
		single_index_type m_indices[N];
};

template <int N>
std::ostream& operator<< (std::ostream& outStream, const ug::MultiIndex<N>& v)
{
	outStream << "[" ;
	for(size_t i = 0; i < N; ++i)
		outStream << v[i];
	outStream << "]";
	return outStream;
}

}


#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__MULTI_INDICES__ */
