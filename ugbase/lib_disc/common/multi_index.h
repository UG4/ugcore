/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__COMMON__MULTI_INDEX__
#define __H__UG__LIB_DISC__COMMON__MULTI_INDEX__


#include <iostream>
#include <cassert>

//#include "common/common.h"
#include "lib_algebra/small_algebra/blocks.h"	// needed for BlockRef

namespace ug {

/**
 * A MultiIndex is just a vector of integers.
 */
template<int N, typename TSingleIndexType = size_t>
class MultiIndex
{
	public:
		using single_index_type = TSingleIndexType;

	public:
		/// number of indices in multi index
		static constexpr inline size_t size() {return N;}

		/// access to index component
		inline single_index_type& operator [] (size_t i)
		{
			UG_ASSERT(i < N, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator [] (size_t i) const
		{
			UG_ASSERT(i < N, "Index invalid");
			return m_indices[i];
		}

		///	comparison operator
		bool operator == (const MultiIndex& o) const
		{
			for(size_t i=0; i < N; ++i)
				if(m_indices[i] != o[i]) return false;
			return true;
		}

		bool operator != (const MultiIndex& o) const
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
		using single_index_type = size_t;

	public:
	///	Default constructor
		MultiIndex() = default;

	///	Constructor with values
		explicit MultiIndex(single_index_type a)
			: m_indices(a)
		{};

		/// number of indices in multi index
		static constexpr inline size_t size() {return 1;}

		/// access to index component
		inline single_index_type& operator [] (size_t i)
		{
			UG_ASSERT(i == 0, "Index invalid");
			return m_indices;
		}

		/// const access to index component
		inline const single_index_type& operator [] (size_t i) const
		{
			UG_ASSERT(i == 0, "Index invalid");
			return m_indices;
		}

		///	comparison operator
		bool operator == (const MultiIndex& o) const
		{
			return m_indices == o[0];
		}

		bool operator != (const MultiIndex& o) const
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
	using single_index_type = size_t;

	public:
	///	Default constructor
		MultiIndex() = default;

	///	Constructor with values
		MultiIndex(single_index_type a, single_index_type b){
			m_indices[0] = a;
			m_indices[1] = b;
		}

		/// number of indices in multi index
		static constexpr inline size_t size() {return 2;}

		/// access to index component
		inline single_index_type& operator [] (size_t i)
		{
			UG_ASSERT(i < 2, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator [] (size_t i) const
		{
			UG_ASSERT(i < 2, "Index invalid");
			return m_indices[i];
		}

		///	comparison operator
		bool operator == (const MultiIndex& o) const
		{
			return (m_indices[0] == o[0]) && (m_indices[1] == o[1]);
		}

		bool operator != (const MultiIndex& o) const
		{
			return (m_indices[0] != o[0]) || (m_indices[1] != o[1]);
		}

		bool operator < (const MultiIndex& o) const
		{
			if(m_indices[0] < o[0]) return true;
			if(m_indices[0] == o[0])
				if(m_indices[1] < o[1]) return true;
			return false;
		}

		bool operator > (const MultiIndex& o) const
		{
			if(m_indices[0] > o[0]) return true;
			if(m_indices[0] == o[0])
				if(m_indices[1] > o[1]) return true;
			return false;
		}

	private:
		single_index_type m_indices[2];
};

// specialization of 3
template <>
class MultiIndex<3, size_t>
{
	public:
	using single_index_type = size_t;

	public:
	///	Default constructor
		MultiIndex() = default;

	///	Constructor with values
		MultiIndex(single_index_type a, single_index_type b, single_index_type c){
			m_indices[0] = a;
			m_indices[1] = b;
			m_indices[2] = c;
		}

		/// number of indices in multi index
		static constexpr inline size_t size() {return 3;}

		/// access to index component
		inline single_index_type& operator [] (size_t i)
		{
			UG_ASSERT(i < 3, "Index invalid");
			return m_indices[i];
		}

		/// const access to index component
		inline const single_index_type& operator [] (size_t i) const
		{
			UG_ASSERT(i < 3, "Index invalid");
			return m_indices[i];
		}

		///	comparison operator
		bool operator == (const MultiIndex& o) const
		{
			return 	(m_indices[0] == o[0]) &&
					(m_indices[1] == o[1]) &&
					(m_indices[2] == o[2]);
		}

		bool operator != (const MultiIndex& o) const
		{
			return !(*this==o);
		}

	private:
		single_index_type m_indices[3];
};

template <int N>
std::ostream& operator << (std::ostream& outStream, const ug::MultiIndex<N>& v)
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
using DoFIndex = MultiIndex<2>;

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

template <typename TMatrix>
void SetRow(TMatrix& mat, const DoFIndex& ind, number val = 0.0)
{
	SetRow(mat, ind[0], ind[1], val);
}

}


#endif