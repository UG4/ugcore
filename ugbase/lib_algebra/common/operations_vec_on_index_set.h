/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_ALGEBRA__OPERATIONS_VEC_ON_INDEX_SET__
#define __H__UG__LIB_ALGEBRA__OPERATIONS_VEC_ON_INDEX_SET__

#include <vector>
#include "operations_vec.h"
#include "common/types.h"

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

/* unused variable "const vector_t &b" maybe unused functions

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

*/

}

#endif