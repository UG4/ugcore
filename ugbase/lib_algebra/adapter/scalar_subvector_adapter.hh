/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef SCALAR_SUBVECTOR_ADAPTER_HH_
#define SCALAR_SUBVECTOR_ADAPTER_HH_

#include <cstdlib>

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/lib_algebra_impl.h"
#include "lib_algebra/cpu_algebra_types.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_vector.h"
#endif

#include "lib_algebra/small_algebra/blocks.h"
namespace ug{

// provides an interface for matrix of algebra type B for matrices originally of algebra type A
// allows to access a CPUBlockAlgebra (AT) as a scalar CPUAlgebra (ST)

template<class InVT, class ST>
class ConstScalarSubVectorAdapter{
public:
	using encapsulated_vector_type = InVT;
	using value_type = typename ST::vector_type::value_type;
	static constexpr int blockSize = block_traits<typename InVT::value_type>::static_size;

	ConstScalarSubVectorAdapter(const encapsulated_vector_type& vec, size_t alpha) : m_src(vec), m_alpha(alpha) {};

	inline const value_type &operator [] (size_t i) const
	{ return BlockRef(m_src[i], m_alpha); }

	void print(const char * const text = nullptr) const
	{ m_src.print(text);}

	size_t size() const
	{ return (m_src.size());}

private:
	const encapsulated_vector_type &m_src;
	const size_t m_alpha;
};

template<class InVT, class ST=CPUAlgebra::vector_type>
class ScalarSubVectorAdapter{

public:
	using encapsulated_vector_type = InVT ;
	using value_type = typename ST::vector_type::value_type;
	static constexpr int blockSize = block_traits<typename InVT::value_type>::static_size;
	// using const_row_iterator= typename ST::vector_type::const_row_iterator

	ScalarSubVectorAdapter(encapsulated_vector_type& vec, size_t alpha) : m_src(vec), m_alpha(alpha) {};

	inline value_type &operator [] (size_t i)
	{ return BlockRef(m_src[i], m_alpha); }

	inline const value_type &operator [] (size_t i) const
	{ return BlockRef(m_src[i], m_alpha); }

	void resize(size_t newSize, bool bCopyValues=true)
	{
		m_src.resize_exactly(newSize, bCopyValues);
	}
	void reserve(size_t newCapacity, bool bCopyValues=true)
	{
		m_src.reserve_exactly(newCapacity, bCopyValues);
	}
	void print(const char * const text = nullptr) const
	{
		m_src.print(text);
	}

	void set_random(double from, double to)
	{
		for(size_t i=0; i<size(); i++)
			BlockRef(m_src[i], m_alpha) = urand(from, to);
	}

	size_t size() const {return (m_src.size());}
private:
	encapsulated_vector_type &m_src;
	const size_t m_alpha;
};


// partielle Spezialisierung fuer CPUAlgebra
template<>
inline ScalarSubVectorAdapter<CPUAlgebra::vector_type, CPUAlgebra>::value_type&
ScalarSubVectorAdapter<CPUAlgebra::vector_type, CPUAlgebra>::operator [] (size_t i)
{ return m_src[i]; }

#ifdef UG_PARALLEL
template<>
inline ScalarSubVectorAdapter<ParallelVector<CPUAlgebra::vector_type>, CPUAlgebra>::value_type&
ScalarSubVectorAdapter<ParallelVector<CPUAlgebra::vector_type> , CPUAlgebra>::operator [] (size_t i)
{ return m_src[i]; }
#endif

} // namespace ug

#endif /* SPARSEMATRIXPROXY_HH_ */
