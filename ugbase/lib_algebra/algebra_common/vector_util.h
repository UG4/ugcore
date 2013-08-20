/*
 * vector_util.h
 *
 *  Created on: 18.06.2013
 *      Author: mrupp
 */

#ifndef VECTOR_UTIL_H_
#define VECTOR_UTIL_H_

#include "lib_algebra/small_algebra/blocks.h"
#include <boost/math/special_functions/fpclassify.hpp>
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug{

template<typename TBlock>
bool BlockVectorFiniteAndNotTooBig(TBlock &v, double tooBigValue=1e24)
{
	for(size_t j=0; j< GetSize(v); j++)
	{
		double d = BlockRef(v, j);
		if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
				return false;
	}
	return true;
}

template<typename TBlock>
bool BlockMatrixFiniteAndNotTooBig(TBlock &m, double tooBigValue=1e24)
{
	for(size_t r=0; r<GetRows(m); r++)
		for(size_t c=0; c< GetCols(m); c++)
		{
			double d = BlockRef(m, r, c);
			if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
					return false;
		}

	return true;
}


template<typename TVector>
bool IsFiniteAndNotTooBig(const TVector &v, double tooBigValue=1e24)
{
	for(size_t i=0; i<v.size(); i++)
		if(BlockVectorFiniteAndNotTooBig(v[i], tooBigValue) == false)
			return false;
	return true;
}

/*template<typename TVector>
bool IsFinite(const TVector &v)
{
	for(size_t i=0; i<v.size(); i++)
	{
		for(size_t j=0; j< GetSize(v[i]); j++)
		{
			double d = BlockRef(v[i], j);
			if(isfinite(d) == false)
				return false;
		}
	}
	return true;
}*/

#ifdef UG_PARALLEL
template<typename TVector>
bool IsFiniteAndNotTooBig(const ParallelVector<TVector> &v)
{
	return AllProcsTrue(IsFiniteAndNotTooBig((TVector&)v), v.layouts()->proc_comm());
}
template<typename TVector>
bool IsFinite(const ParallelVector<TVector> &v)
{
	return AllProcsTrue(IsFiniteAndNotTooBig((TVector&)v), v.layouts()->proc_comm());
}
#endif

}
#endif /* VECTOR_UTIL_H_ */
