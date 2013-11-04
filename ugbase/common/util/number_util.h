/*
 * number_util.h
 *
 *  Created on: 01.11.2013
 *      Author: mrupp
 */

#ifndef NUMBER_UTIL_H_
#define NUMBER_UTIL_H_
#include <boost/math/special_functions/fpclassify.hpp>

#include "common/math/ugmath.h"
namespace ug{
inline bool IsFiniteAndNotTooBig(double d)
{
	const double tooBigValue = 1e30;
	if(d > tooBigValue || d < -tooBigValue || std::isfinite(d) == false)
//	if(std::isfinite(d) == false)
		return false;
	return true;
}

inline bool CloseToZero(double d)
{
	if(d > -1e-100 && d < +1e-100)
		return true;
	else
		return false;
}

//inline bool IsFiniteAndNotTooBig(float d)
//{
//	return IsFiniteAndNotTooBig((double) d);
//}


template <std::size_t N, std::size_t M, typename T>
inline bool IsFiniteAndNotTooBig(const MathMatrix<N, M, T> &m)
{
	for(size_t r=0; r<m.num_rows(); r++)
		for(size_t c=0; c<m.num_cols(); c++)
			if(IsFiniteAndNotTooBig(m(r,c)) == false) return false;
	return true;
}

template <std::size_t N, typename T>
inline bool IsFiniteAndNotTooBig(const MathVector<N, T> &m)
{
	for(size_t r=0; r<m.size(); r++)
			if(IsFiniteAndNotTooBig(m[r]) == false) return false;
	return true;
}



template <size_t TRank, size_t N, typename T>
inline bool IsFiniteAndNotTooBig(const MathTensor<TRank, N, T> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}


template<typename TData, size_t N>
inline bool IsFiniteAndNotTooBig(const MathTensorX<TData, N> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}

template<typename TData>
inline bool IsFiniteAndNotTooBig(const std::vector<TData> &t)
{
	for(size_t i=0; i<t.size(); i++)
		if(IsFiniteAndNotTooBig(t[i]) == false) return false;

	return true;
}



}


#endif /* NUMBER_UTIL_H_ */
