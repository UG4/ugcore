/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__LIB_ALGEBRA__OPERATIONS_VEC__
#define __H__UG__LIB_ALGEBRA__OPERATIONS_VEC__

#include <stddef.h> // size_t
#include <cmath> // log, exp

#include "../small_algebra/storage/fixed_array.h"
#include "../small_algebra/small_matrix/densevector.h"


namespace ug
{

// operations for doubles
//-----------------------------------------------------------------------------
// todo: check if we might replace double with template<T>

// VecScale: These function calculate dest = sum_i alpha_i v_i

//! calculates dest = alpha1*v1. for doubles
inline void VecScaleAssign(double &dest, double alpha1, const double &v1)
{
	dest = alpha1*v1;
}


//! calculates dest = alpha1*v1 + dest for doubles (FMA?)
inline void VecScaleAdd(double &dest, double alpha1, const double &v1)
{

	//#ifdef FP_FAST_FMA
	//			y[i] =std::fma(alpha, x[i], y[i]);
	//#else
	dest = alpha1*v1 + dest;
}

//! calculates dest = alpha1*v1 + alpha2*v2. for doubles
inline void VecScaleAdd(double &dest, double alpha1, const double &v1, double alpha2, const double &v2)
{
	//#ifdef FP_FAST_FMA
	//			dest[i] = std::fma(alpha1, v1[i], alpha2 *v2[i]);
	//#else
	dest = alpha1*v1 + alpha2*v2;
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3. for doubles
inline void VecScaleAdd(double &dest, double alpha1, const double &v1, double alpha2, const double &v2, double alpha3, const double &v3)
{
	dest = alpha1*v1 + alpha2*v2 + alpha3*v3;
}


// VecProd

//! calculates s += scal<a, b>
inline void VecProdAdd(const double &a, const double &b, double &s)
{
	s += a*b;
}

template<typename vector_t>
inline void VecProdAdd(const vector_t &a, const vector_t &b, double &s)
{
	// was: #pragma omp simd reduction(+:s)
	const size_t N= a.size();
	// #pragma omp parallel for simd shared(a,b,N) schedule(static) reduction(+:s)
	for(size_t i=0; i<N; i++)
	{
		double dot=0.0;
		VecProdAdd(a[i], b[i], dot);
		s += dot;
	}
}


//! returns scal<a, b>
inline double VecProd(const double &a, const double &b)
{
	return a*b;
}


//! computes scal<a, b>
inline void VecProd(const double &a, const double &b, double &s)
{
	s = a*b;
}


// VecNorm

//! returns norm_2^2(a)
inline double VecNormSquared(const double &a)
{
	return a*a;
}

//! calculates s += norm_2^2(a)
inline void VecNormSquaredAdd(const double &a, double &s)
{
	s += a*a;
}

// Elementwise (Hadamard) product of two vectors

//! calculates s = a * b (the Hadamard product)
inline void VecHadamardProd(double &dest, const double &v1, const double &v2)
{
	dest = v1 * v2;
}

// Some unary mathematical elementwise operations on vectors

//! calculates elementwise exp
inline void VecExp(double &dest, const double &v)
{
	dest = exp (v);
}

//! calculates elementwise log (natural logarithm)
inline void VecLog(double &dest, const double &v)
{
	dest = log (v);
}

// templated

// operations for vectors
//-----------------------------------------------------------------------------
// these functions execute vector operations by using the operations on the elements of the vector

// todo: change vector_t to TE_VEC<vector_t>


// VecScale: These function calculate dest = sum_i alpha_i v_i

//! calculates dest = alpha1*v1
template<typename vector_t>
inline void VecScaleAssign(vector_t &dest, double alpha1, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAssign(dest[i], alpha1, v1[i]);
}

//! sets dest = v1 entrywise
template<typename vector_t>
inline void VecAssign(vector_t &dest, const vector_t &v1)
{
	for(size_t i=0; i<dest.size(); i++)
		dest[i] = v1[i];
}


//! calculates dest = alpha1*v1 + alpha2*v2
template <int N>
inline void VecScaleAdd(DenseVector<FixedArray1<double, N>> &dest,
			double alpha1, const DenseVector<FixedArray1<double, 2>> &v1,
			double alpha2, const DenseVector<FixedArray1<double, 2>> &v2)
{
	// #pragma unroll N
	for (size_t i = 0; i<N; ++i)
	{ dest[i] = alpha1*v1[i] + alpha2*v2[i]; }
}

//! calculates dest = alpha1*v1 + alpha2*v2
/*inline void VecScaleAdd(DenseVector<FixedArray1<double, 3>> &dest,
			double alpha1, const DenseVector<FixedArray1<double, 3>> &v1,
			double alpha2, const DenseVector<FixedArray1<double, 3>> &v2)
{
	dest[0] = alpha1*v1[0] + alpha2*v2[0];
	dest[1] = alpha1*v1[1] + alpha2*v2[1];
	dest[2] = alpha1*v1[2] + alpha2*v2[2];
}

//! calculates dest = alpha1*v1 + alpha2*v2
inline void VecScaleAdd(DenseVector<FixedArray1<double, 4>> &dest,
			double alpha1, const DenseVector<FixedArray1<double, 4>> &v1,
			double alpha2, const DenseVector<FixedArray1<double, 4>> &v2)
{

	dest[0] = alpha1*v1[0] + alpha2*v2[0];
	dest[1] = alpha1*v1[1] + alpha2*v2[1];
	dest[2] = alpha1*v1[2] + alpha2*v2[2];
	dest[3] = alpha1*v1[3] + alpha2*v2[3];

}

*/

//! calculates dest = alpha1*v1 + dest
template<typename vector_t, template <class T> class TE_VEC>
inline void VecScaleAdd(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1)
{
	// #pragma omp parallel for simd
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAdd(dest[i], alpha1, v1[i]);
}

//! calculates dest = alpha1*v1 + alpha2*v2
template<typename vector_t, template <class T> class TE_VEC>
inline void VecScaleAdd(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2)
{
	// #pragma omp parallel for simd
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i]);
}

//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename vector_t, template <class T> class TE_VEC>
inline void VecScaleAdd(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2, double alpha3, const TE_VEC<vector_t> &v3)
{
	for(size_t i=0; i<dest.size(); i++)
		VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i]);
}


// VecProd

//! calculates s += scal<a, b>
template<typename vector_t>
inline void VecProd(const vector_t &a, const vector_t &b, double &sum)
{
	// was: #pragma omp simd reduction(+:sum)
	// #pragma omp parallel for simd shared(a,b) schedule(static) reduction(+:sum)
	for(size_t i=0; i<a.size(); i++)
	{
		// double dot = 0.0;
		VecProdAdd(a[i], b[i], sum);
		// sum +=dot;
	}
}

//! returns scal<a, b>
template<typename vector_t>
constexpr double VecProd(const vector_t &a, const vector_t &b)
{
	double sum=0;
	VecProdAdd(a, b, sum);
	return sum;
}


//! calculates s += norm_2^2(a)
template<typename vector_t>
inline void VecNormSquaredAdd(const vector_t &a, double &sum)
{
	const size_t N =a.size();
	// #pragma omp parallel for simd shared(a) reduction(+:sum)
	for(size_t i=0; i<N; i++)
	{
		// double norm2=0.0;
		VecNormSquaredAdd(a[i], sum);
		// sum +=norm2;
	}
}

//! returns norm_2^2(a)
template<typename vector_t>
inline double VecNormSquared(const vector_t &a)
{
	double sum=0;
	VecNormSquaredAdd(a, sum);
	return sum;
}

// Elementwise (Hadamard) product of two vectors
template<typename vector_t>
inline void VecHadamardProd(vector_t &dest, const vector_t &v1, const vector_t &v2)
{
	for(size_t i=0; i<dest.size(); i++)
		VecHadamardProd(dest[i], v1[i], v2[i]);
}

// Elementwise exp on a vector
template<typename vector_t>
inline void VecExp(vector_t &dest, const vector_t &v)
{
	for(size_t i=0; i<dest.size(); i++)
		VecExp(dest[i], v[i]);
}

// Elementwise log (natural logarithm) on a vector
template<typename vector_t>
inline void VecLog(vector_t &dest, const vector_t &v)
{
	for(size_t i=0; i<dest.size(); i++)
		VecLog(dest[i], v[i]);
}



//! We also provide this as a struct. This uses thread-parallel execution.
template<typename TVector>
struct vector_operations
{

	static double dot(const TVector &a, const TVector &b, double sum=0.0)
	{
		const size_t N =a.size();
		#pragma omp parallel for simd shared(a,b) schedule(static) reduction(+:sum)
		for(size_t i=0; i<N; i++)
		{ VecProdAdd(a[i], b[i], sum); }
		return sum;
	}

	static double norm2(const TVector &a, double sum=0.0)
	{
		const size_t N =a.size();
		#pragma omp parallel for simd shared(a) reduction(+:sum)
		for(size_t i=0; i<N; i++)
		{ VecNormSquaredAdd(a[i], sum); }
		return sum;
	}

	// y = a*x + y
	static void axpy(double alpha, const TVector &x, TVector &y)
	{
		const size_t N =y.size();
		#pragma omp parallel for simd
		for(size_t i=0; i<N; i++)
		{ VecScaleAdd(y[i], alpha, x[i]); }
	}

	// z = alpha1 * v1 + alpha2 *v2
	static void vec_scale_add(TVector &dest, double alpha1, const TVector &v1, double alpha2, const TVector &v2)
	{
		const size_t N =dest.size();
		#pragma omp parallel for simd
		for(size_t i=0; i<N; i++)
		{ VecScaleAdd(dest[i], alpha1, v1[i], alpha2, v2[i]); }
	}

};


} // namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATIONS_VEC__ */
