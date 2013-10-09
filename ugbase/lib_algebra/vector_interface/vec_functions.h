/*
 * vec_vec_functions.h
 *
 *  Created on: 26.09.2013
 *      Author: mrupp
 */

#ifndef VEC_VEC_FUNCTIONS_H_
#define VEC_VEC_FUNCTIONS_H_
#include "common/error.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// T, T functions
template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2)
{
	ASSERT_EQUAL(v1.size(), v2.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i];
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2, double a3, const T &v3)
{
	ASSERT_EQUAL_3(v1.size(), v2.size(), v3.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i] + a3*v3[i];
}

template<typename T>
inline void VecAdd(double a1, T &v1, double a2, const T &v2, double a3, const T &v3, double a4, const T &v4)
{
	ASSERT_EQUAL_4(v1.size(), v2.size(), v3.size(), v4.size());
	for(size_t i=0; i<v1.size(); i++)
		v1[i] = a1*v1[i] + a2*v2[i] + a3*v3[i] + a4*v4[i];
}

template<typename T>
inline double VecProd(T &v1, T &v2)
{
	ASSERT_EQUAL(v1.size(), v2.size());
	double s=0;
	for(size_t i=0; i<v1.size(); i++)
		s += v1[i]*v2[i];
	return s;
}

template<typename T>
inline double VecNorm2(T &v)
{
	double s=0;
	for(size_t i=0; i<v.size(); i++)
		s += v[i]*v[i];
	return s;
}


#endif /* VEC_VEC_FUNCTIONS_H_ */
