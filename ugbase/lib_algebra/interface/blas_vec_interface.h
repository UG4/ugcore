/*
 * blas_vec_interface.h
 *
 *  Created on: 22.03.2011
 *      Author: mrupp
 */

#ifndef BLAS_VEC_INTERFACE_H_
#define BLAS_VEC_INTERFACE_H_


// specializations for double
//! calculates dest = alpha1*v1. for doubles
inline void VecScaleAssign(double &dest, double alpha1, const double &v1)
{
	dest = alpha1*v1;
}
// and so on.


//! calculates dest = alpha1*v1
//
template<typename vector_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAssign(dest[i], alpha1, v1[i]);
}

//! calculates dest = v1
//
template<typename vector_t>
inline void VecAssign(vector_t &dest, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
	    dest[i]=v1[i];
}


//! calculates dest = alpha1*v1 + alpha2*v2
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2);

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t>
inline void VecScaleAdd(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const vector_t &v2, double alpha3, const vector_t &v3);


//! calculates dest = dest + v
template<typename vector_t>
inline void VecAdd(vector_t &dest, vector_t &v);

//! calculates dest = dest - v
template<typename vector_t>
inline void VecSubstract(vector_t &dest, vector_t &v);


// VecProd

//! calculates s += scal<a, b>
template<typename vector_t>
inline void VecProdAdd(const vector_t &a, const vector_t &b, double &sum)
{
	for(size_t i=0; i<a.size(); i++) VecProdAdd(a[i], b[i], sum);
}

//! returns scal<a, b>
template<typename vector_t>
inline double VecProd(const vector_t &a, const vector_t &b)
{
	double sum=0;
	VecProdAdd(a, b, sum);
	return sum;
}

// VecNorm

//! calculates s += norm_2^2(a)
template<typename vector_t>
inline void VecNormSquaredAdd(const vector_t &a, const vector_t &b, double &sum)
{
	for(int i=0; i<a.size(); i++) VecNormSquaredAdd(a[i], sum);
}

//! returns norm_2^2(a)
template<typename vector_t>
inline double VecNormSquared(const vector_t &a, const vector_t &b)
{
	double sum=0;
	VecNormSquaredAdd(a, sum);
	return sum;
}

#endif /* BLAS_VEC_INTERFACE_H_ */
