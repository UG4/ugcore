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
	typedef InVT encapsulated_vector_type;
	typedef typename ST::vector_type::value_type value_type;
	static const int blockSize = block_traits<typename InVT::value_type>::static_size;

	ConstScalarSubVectorAdapter(const encapsulated_vector_type& vec, size_t alpha) : m_src(vec), m_alpha(alpha) {};

	inline const value_type &operator [] (size_t i) const
	{ return BlockRef(m_src[i], m_alpha); }

	void print(const char * const text = NULL) const
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
	typedef InVT encapsulated_vector_type;
	typedef typename ST::vector_type::value_type value_type;
	static const int blockSize = block_traits<typename InVT::value_type>::static_size;
	//typedef typename ST::vector_type::const_row_iterator const_row_iterator

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
	void print(const char * const text = NULL) const
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
