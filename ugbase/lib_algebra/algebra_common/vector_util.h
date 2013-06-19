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

namespace ug{
template<typename TVector>
bool IsFiniteAndNotTooBig(TVector &v, double tooBigValue=1e24)
{
	for(size_t i=0; i<v.size(); i++)
	{
		for(size_t j=0; j< GetSize(v[i]); j++)
		{
			double d = BlockRef(v[i], j);
			if(d > tooBigValue || d < -tooBigValue || isfinite(d) == false)
				return false;
		}
	}
	return true;
}

template<typename TVector>
bool IsFinite(TVector &v)
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
}

template<typename TVector>
bool IsNormalAndNotTooBig(TVector &v, double tooBigValue=1e24)
{
	return IsFiniteAndNotTooBig(v, tooBigValue);
}

}
#endif /* VECTOR_UTIL_H_ */
