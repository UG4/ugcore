/*
 * multi_index.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__COMMON__MULTI_INDEX__
#define __H__UG__LIB_DISC__COMMON__MULTI_INDEX__

#include <vector>
#include <iostream>
#include <assert.h>

#include "common/common.h"

namespace ug{

/**
 * A MultiIndex is just a vector of integers.
 */
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

		///	comparison operator
		bool operator==(const MultiIndex& o) const
		{
			for(size_t i=0; i < N; ++i)
				if(m_indices[i] != o[i]) return false;
			return true;
		}

		bool operator!=(const MultiIndex& o) const
		{
			return !(*this==o);
		}

	private:
		single_index_type m_indices[N];
};

// specialization of 1
template <>
class MultiIndex<1, size_t>
{
	public:
		typedef size_t single_index_type;

	public:
	///	Default constructor
		MultiIndex(){};

	///	Constructor with values
		MultiIndex(single_index_type a)
			: m_indices(a)
		{};

		/// number of indices in multi index
		inline size_t size() const {return 1;}

		/// access to index component
		inline single_index_type& operator[] (size_t i)
		{
			UG_ASSERT(i == 0, "Index invalid");
			return m_indices;
		}

		/// const access to index component
		inline const single_index_type& operator[] (size_t i) const
		{
			UG_ASSERT(i == 0, "Index invalid");
			return m_indices;
		}

		///	comparison operator
		bool operator==(const MultiIndex& o) const
		{
			return m_indices == o[0];
		}

		bool operator!=(const MultiIndex& o) const
		{
			return !(*this==o);
		}

	private:
		single_index_type m_indices;
};

// specialization of 2
template <>
class MultiIndex<2, size_t>
{
	public:
		typedef size_t single_index_type;

	public:
	///	Default constructor
		MultiIndex(){};

	///	Constructor with values
		MultiIndex(single_index_type a, single_index_type b){
			m_indices[0] = a;
			m_indices[1] = b;
		}

		/// number of indices in multi index
		inline size_t size() const {return 2;}

		/// access to index component
		inline single_index_type& operator[] (size_t i)
		{
			UG_ASSERT(i < 2, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator[] (size_t i) const
		{
			UG_ASSERT(i < 2, "Index invalid");
			return m_indices[i];
		}

		///	comparison operator
		bool operator==(const MultiIndex& o) const
		{
			return (m_indices[0] == o[0]) && (m_indices[1] == o[1]);
		}

		bool operator!=(const MultiIndex& o) const
		{
			return !(*this==o);
		}


	private:
		single_index_type m_indices[2];
};

// specialization of 3
template <>
class MultiIndex<3, size_t>
{
	public:
		typedef size_t single_index_type;

	public:
	///	Default constructor
		MultiIndex(){};

	///	Constructor with values
		MultiIndex(single_index_type a, single_index_type b, single_index_type c){
			m_indices[0] = a;
			m_indices[1] = b;
			m_indices[2] = c;
		}

		/// number of indices in multi index
		inline size_t size() const {return 3;}

		/// access to index component
		inline single_index_type& operator[] (size_t i)
		{
			UG_ASSERT(i < 3, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator[] (size_t i) const
		{
			UG_ASSERT(i < 3, "Index invalid");
			return m_indices[i];
		}

		///	comparison operator
		bool operator==(const MultiIndex& o) const
		{
			return 	(m_indices[0] == o[0]) &&
					(m_indices[1] == o[1]) &&
					(m_indices[2] == o[2]);
		}

		bool operator!=(const MultiIndex& o) const
		{
			return !(*this==o);
		}

	private:
		single_index_type m_indices[3];
};

template <int N>
std::ostream& operator<< (std::ostream& outStream, const ug::MultiIndex<N>& v)
{
	outStream << "[" ;
	for(size_t i = 0; i < N; ++i)
	{
		outStream << v[i];
		if(i != N-1) outStream << ",";
	}
	outStream << "]";
	return outStream;
}

////////////////////////////////////////////////////////////////////////////////
//	degree of freedom access using multi indices
////////////////////////////////////////////////////////////////////////////////

/// type of DoF-Index used to identify an DoF in the Algebra
typedef MultiIndex<2> DoFIndex;

template <typename TMatrix>
inline number&
DoFRef(TMatrix& mat, const DoFIndex& iInd, const DoFIndex& jInd)
{
	return BlockRef(mat(iInd[0], jInd[0]), iInd[1], jInd[1]);
}

template <typename TMatrix>
inline const number&
DoFRef(const TMatrix& mat, const DoFIndex& iInd, const DoFIndex& jInd)
{
	return BlockRef(mat(iInd[0], jInd[0]), iInd[1], jInd[1]);
}

template <typename TVector>
inline number&
DoFRef(TVector& vec, const DoFIndex& ind)
{
	return BlockRef(vec[ind[0]], ind[1]);
}

template <typename TVector>
inline const number&
DoFRef(const TVector& vec, const DoFIndex& ind)
{
	return BlockRef(vec[ind[0]], ind[1]);
}

template <typename TMatrix>
void SetDirichletRow(TMatrix& mat, const DoFIndex& ind)
{
	SetDirichletRow(mat, ind[0], ind[1]);
}

}


#endif /* __H__UG__LIB_DISC__COMMON__MULTI_INDEX__ */
