
#ifndef __H__UG__LIB_ALGEBRA__OPERATIONS_VEC_ON_INDEX_SET__
#define __H__UG__LIB_ALGEBRA__OPERATIONS_VEC_ON_INDEX_SET__

#include "operations_vec.h"
#include <vector>

namespace ug
{

/// sets dest = alpha on a given index set
template<typename vector_t>
inline void VecSet(vector_t &dest, number alpha,
                   const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < dest.size(), "Index to large in index set.");
		dest[i] = alpha;
	}
}

/// calculates dest = alpha1*v1 on a given index set
template<typename vector_t>
inline void VecScaleAssign(vector_t &dest,
                           number alpha1, const vector_t &v1,
                           const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < dest.size(), "Index to large in index set.");
		UG_ASSERT(i < v1.size(), "Index to large in index set.");
		VecScaleAssign(dest[i], alpha1, v1[i]);
	}
}

/// calculates dest = alpha1*v1 + alpha2*v2 on a given index set
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest,
                         number alpha1, const vector_t &v1,
                         number alpha2, const vector_t &v2,
                         const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < dest.size(), "Index to large in index set.");
		UG_ASSERT(i < v1.size(), "Index to large in index set.");
		UG_ASSERT(i < v2.size(), "Index to large in index set.");
		VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i]);
	}
}


/// calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3 on a given index set
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest,
                        number alpha1, const vector_t &v1,
                        number alpha2, const vector_t &v2,
                        number alpha3, const vector_t &v3,
                        const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < dest.size(), "Index to large in index set.");
		UG_ASSERT(i < v1.size(), "Index to large in index set.");
		UG_ASSERT(i < v2.size(), "Index to large in index set.");
		UG_ASSERT(i < v3.size(), "Index to large in index set.");
		VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i]);
	}
}


// VecProd

/// calculates s += scal<a, b> on a given index set
template<typename vector_t>
inline void VecProd(const vector_t &a, const vector_t &b,
                    number &sum, const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < a.size(), "Index to large in index set.");
		UG_ASSERT(i < b.size(), "Index to large in index set.");
		VecProdAdd(a[i], b[i], sum);
	}
}

/// returns scal<a, b> on a given index set
template<typename vector_t>
inline number VecProd(const vector_t &a, const vector_t &b,
                      const std::vector<size_t> vIndex)
{
	number sum=0;
	VecProd(a, b, sum, vIndex);
	return sum;
}


/// calculates s += norm_2^2(a) on a given index set
template<typename vector_t>
inline void VecNormSquaredAdd(const vector_t &a, const vector_t &b,
                              number &sum, const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < a.size(), "Index to large in index set.");
		VecNormSquaredAdd(a[i], sum);
	}
}

/// returns norm_2^2(a) on a given index set
template<typename vector_t>
inline number VecNormSquared(const vector_t &a, const vector_t &b,
                             const std::vector<size_t> vIndex)
{
	number sum=0;
	VecNormSquaredAdd(a, sum, vIndex);
	return sum;
}



}

#endif /* __H__UG__LIB_ALGEBRA__OPERATIONS_VEC_ON_INDEX_SET__ */
