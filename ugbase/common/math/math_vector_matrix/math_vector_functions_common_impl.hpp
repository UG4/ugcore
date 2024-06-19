/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 ��7):
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
 * "Vogel, A., Reiter, S., Rupp, M., N��gel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This file supplies a general implementation of vector_functions.
//	The methods implemented here work on every vector-type, that supports
//	the vector-description in lgmath_vector_descriptor.h.
//	Since this implementation is quite general it is not the fastest.
//	One should consider implementing specializations for commonly used types
//	like vector3 etc.
////////////////////////////////////////////////////////////////////////

#ifndef __H__COMMON__VECTOR_FUNCTIONS_COMMON_IMPL__
#define __H__COMMON__VECTOR_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include <algorithm>
#include "common/error.h"
#include "common/assert.h"
#include "common/static_assert.h"

namespace ug
{

template <typename vector_target_t, typename vector_source_t>
void VecCopy(vector_target_t& target, const vector_source_t& source,
			 typename vector_target_t::value_type fill)
{
	using std::min;
	size_t minSize = min(target.size(), source.size());
	for(size_t i = 0; i < minSize; ++i)
		target[i] = source[i];

	for(size_t i = minSize; i < target.size(); ++i)
		target[i] = fill;
}

///	adds a MathVector<N>s to a second one
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += v1[i];
	}
}

///	adds two MathVector<N>s and adds the result to a third one
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += v1[i] + v2[i];
	}
}

///	adds three Vectors and adds the result to a fourth one
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += v1[i] + v2[i] + v3[i];
	}
}

///	adds four Vectors and adds the result to a fifth one
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2,
			const vector_t& v3, const vector_t& v4)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += v1[i] + v2[i] + v3[i] + v4[i];
	}
}

/// Scales a Vector and adds it to a second vector
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += s1 * v1[i];
	}
}

/// Scales two Vectors, adds the sum to a third vector
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += s1 * v1[i] + s2 * v2[i];
	}
}

/// Scales three Vectors, adds the sum to a fourth vector
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += s1 * v1[i] + s2 * v2[i] + s3 * v3[i];
	}
}

/// Scales four Vectors, adds the sum to a fifth vector
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3,
								 typename vector_t::value_type s4, const vector_t& v4)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] += s1 * v1[i] + s2 * v2[i] + s3 * v3[i] + s4 * v4[i];
	}
}


///	adds two MathVector<N>s and stores the result in a third one
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v1[i] + v2[i];
	}
}

///	adds three Vectors and stores the result in a fourth one
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v1[i] + v2[i] + v3[i];
	}
}

///	adds four Vectors and stores the result in a fifth one
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2,
			const vector_t& v3, const vector_t& v4)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v1[i] + v2[i] + v3[i] + v4[i];
	}
}

///	subtracts v2 from v1 and stores the result in a vOut
template <typename vector_t>
inline
void
VecSubtract(vector_t& vOut, const vector_t& v1, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v1[i] - v2[i];
	}
}

/// component-wise exponentiation of a vector
template <typename vector_t>
inline
void
VecPow(vector_t& vOut, const vector_t& v1, typename vector_t::value_type s)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = std::pow(v1[i], s);
	}
}

///	scales a MathVector<N>
template <typename vector_t>
inline
void
VecScale(vector_t& vOut, const vector_t& v, typename vector_t::value_type s)
{
	typedef typename vector_t::size_type size_type;
	const size_type N = vOut.size();
	#pragma omp parallel for simd
	for(size_type i = 0; i < N; ++i)
	{
		vOut[i] = s * v[i];
	}
}

/// Scales two Vectors, adds them and returns the sum in a third vector
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	const size_type N = vOut.size();
	#pragma omp parallel for simd
	for(size_type i = 0; i < N; ++i)
	{
		vOut[i] = s1 * v1[i] + s2 * v2[i];
	}
}

/// Scales three Vectors, adds them and returns the sum in a fourth vector
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3)
{
	typedef typename vector_t::size_type size_type;
	const size_type N = vOut.size();
	#pragma omp parallel for simd
	for(size_type i = 0; i < N; ++i)
	{
		vOut[i] = s1 * v1[i] + s2 * v2[i] + s3 * v3[i];
	}
}

/// Scales four Vectors, adds them and returns the sum in a fifth vector
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3,
								 typename vector_t::value_type s4, const vector_t& v4)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = s1 * v1[i] + s2 * v2[i] + s3 * v3[i] + s4 * v4[i];
	}
}

//	performs linear interpolation between two MathVector<N>s.
template <typename vector_t>
inline
void
VecInterpolateLinear(vector_t& vOut, const vector_t& v1, const vector_t& v2,
						typename vector_t::value_type interpAmount)
{
	typedef typename vector_t::size_type size_type;
	const size_type N =vOut.size();

	#pragma omp parallel for simd
	for(size_type i = 0; i < N; ++i)
	{
		vOut[i] = (1. - interpAmount) * v1[i] + interpAmount * v2[i];
	}
}

///	returns the squared length of v. Faster than VecLength.
template <typename vector_t>
inline
typename vector_t::value_type
VecLengthSq(const vector_t& v)
{
	typename vector_t::value_type len = 0;

	typedef typename vector_t::size_type size_type;
	const size_type N =v.size();

	#pragma omp parallel for simd shared(v, N) reduction(+:len)
	for(size_type i = 0; i < N; ++i)
	{
		len += v[i] * v[i];
	}

	return len;
}

///	returns the length of v. Slower than VecLengthSq.
template <typename vector_t>
inline
typename vector_t::value_type
VecLength(const vector_t& v)
{
	return static_cast<typename vector_t::value_type>(
			std::sqrt(VecLengthSq(v)));
}

///	returns the squared distance of two MathVector<N>s.
template <typename vector_t>
inline
typename vector_t::value_type
VecDistanceSq(const vector_t& v1, const vector_t& v2)
{
	vector_t v;
	VecSubtract(v, v1, v2);
	return VecLengthSq(v);
}

///	returns the squared distance of two MathVector<N>s.
template <typename TVector, typename TMatrix>
inline
typename TVector::value_type
VecDistanceSq(const TVector& v1, const TVector& v2, const TMatrix& M)
{
	TVector delta;
	TVector deltaM;
	VecSubtract(delta, v1, v2);
	MatVecMult(deltaM, M, delta);
	return VecDot(deltaM, delta);
}

///	returns the distance of two MathVector<N>s.
template <typename vector_t>
inline
typename vector_t::value_type
VecDistance(const vector_t& v1, const vector_t& v2)
{
	return static_cast<typename vector_t::value_type>
			(std::sqrt(VecDistanceSq(v1, v2)));
}

///	calculates the dot-product of two MathVector<N>s
template <typename vector_t>
inline
typename vector_t::value_type
VecDot(const vector_t& v1, const vector_t& v2)
{
	typename vector_t::value_type dp = 0;
	typedef typename vector_t::size_type size_type;

	for(size_type i = 0; i < v1.size(); ++i)
	{
		dp += v1[i] * v2[i];
	}

	return dp;
}

template <typename vector_t>
inline
typename vector_t::value_type
VecAngle(const vector_t& v1, const vector_t& v2)
{
	typedef typename vector_t::value_type value_t;

	value_t l = sqrt(VecLengthSq(v1) * VecLengthSq(v2));
	if(l > 0){
		number a = VecDot(v1, v2) / l;
		if(a >= 1)
			return 0;
		else if(a <= -1)
			return PI;
		return acos(a);
	}

	return 0;
}

template <typename vector_t>
inline
typename vector_t::value_type
VecAngleNorm(const vector_t& v1, const vector_t& v2)
{
	number a = VecDot(v1, v2);
	if(a >= 1)
		return 0;
	else if(a <= -1)
		return PI;
	return acos(a);
}

///	calculates the cross product of two MathVector<N>s of dimension 3. It makes no sense to use VecCross for
/// MathVector<N>s that have not dimension 3.
/// todo: Create compile-time check for dimensionality in case of MathVector type.
template <typename vector_t>
inline
void
VecCross(vector_t& vOut, const vector_t& v1, const vector_t& v2)
{
	if(&vOut != &v1 && &vOut != &v2)
	{
		vOut[0] = v1[1] * v2[2] - v2[1] * v1[2];
		vOut[1] = v1[2] * v2[0] - v2[2] * v1[0];
		vOut[2] = v1[0] * v2[1] - v2[0] * v1[1];
	}
	else
	{
		vector_t _temp;
		_temp[0] = v1[1] * v2[2] - v2[1] * v1[2];
		_temp[1] = v1[2] * v2[0] - v2[2] * v1[0];
		_temp[2] = v1[0] * v2[1] - v2[0] * v1[1];
		vOut = _temp;
	}
}

/** computes the "generalized vector product" of two vectors
 *
 * In 3d, this is a usual cross-product. In 2d, the first component
 * of the result is the determinant of the 2x2-matrix build of the operands,
 * and the second component is always 0.
 *
 * \tparam dim	dimensionality of the vectors
 */
template <size_t dim>
inline void GenVecCross
(
	MathVector<dim> & result, ///< = v_1 X v_2
	const MathVector<dim> & v_1, const MathVector<dim> & v_2
)
{
	/* Cf. the specializations below! */
	UG_THROW ("The generalized vector product is defined only in 2 and 3 dimensions");
};

/// specialization of the "generalized vector product" in 2d. \see GenVecCross
template <>
inline void GenVecCross<2>
(
	MathVector<2> & result,
	const MathVector<2> & v_1, const MathVector<2> & v_2
)
{
	result[0] = v_1[0] * v_2[1] - v_1[1] * v_2[0];
	result[1] = 0;
};

/// specialization of the "generalized vector product" in 3d. \see GenVecCross
template <>
inline void GenVecCross<3>
(
	MathVector<3> & result,
	const MathVector<3> & v_1, const MathVector<3> & v_2
)
{
	VecCross (result, v_1, v_2);
};

///	scales a MathVector<N> to unit length
template <typename vector_t>
inline
void
VecNormalize(vector_t& vOut, const vector_t& v)
{
	typename vector_t::value_type len = VecLength(v);
	if(len > 0)
		VecScale(vOut, v, (typename vector_t::value_type) 1 / len);
	else
		vOut = v;
}

///	Calculates a triangle-normal in 3d.
template <typename vector_t>
inline
void
CalculateTriangleNormalNoNormalize(vector_t& vOut, const vector_t& v1,
						const vector_t& v2, const vector_t& v3)
{
	vector_t e1, e2;
	VecSubtract(e1, v2, v1);
	VecSubtract(e2, v3, v1);

	VecCross(vOut, e1, e2);
}

///	Calculates a triangle-normal in 3d.
template <typename vector_t>
inline
void
CalculateTriangleNormal(vector_t& vOut, const vector_t& v1,
						const vector_t& v2, const vector_t& v3)
{
	CalculateTriangleNormalNoNormalize(vOut, v1, v2, v3);
	VecNormalize(vOut, vOut);
}

/// Set each vector component to scalar
template <typename vector_t>
inline
void
VecSet(vector_t& vInOut, typename vector_t::value_type s)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vInOut.size(); ++i)
	{
		vInOut[i] = s;
	}
}

/// Add a scalar to a vector (componentwise)
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v, typename vector_t::value_type s)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v[i] + s;
	}
}

/// Subtract a scalar from a vector (componentwise)
template <typename vector_t>
inline
void
VecSubtract(vector_t& vOut, const vector_t& v, typename vector_t::value_type s)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v[i] - s;
	}
}

template <typename vector_t>
inline
typename vector_t::value_type
VecTwoNorm(const vector_t& v)
{
	return VecLength(v);
}

template <typename vector_t>
inline
typename vector_t::value_type
VecTwoNormSq(const vector_t& v)
{
	return VecLengthSq(v);
}

template <typename vector_t>
inline
typename vector_t::value_type
VecOneNorm(const vector_t& v)
{
	typename vector_t::value_type len = 0;
	typedef typename vector_t::size_type size_type;

	for(size_type i = 0; i < v.size(); ++i)
	{
		len += std::abs(v[i]);
	}

	return len;
}

template <typename vector_t>
inline
typename vector_t::value_type
VecPNorm(const vector_t& v, unsigned int p)
{
	typename vector_t::value_type len = 0;
	typedef typename vector_t::size_type size_type;

	for(size_type i = 0; i < v.size(); ++i)
	{
		len += std::pow(v[i], p);
	}

	return std::pow(len, (typename vector_t::value_type) 1/ p);
}

template <typename vector_t>
inline
typename vector_t::value_type
VecMaxNorm(const vector_t& v)
{
	typename vector_t::value_type m = 0;
	typedef typename vector_t::size_type size_type;

	for(size_type i = 0; i < v.size(); ++i)
	{
		m = std::max(m, v[i]);
	}

	return m;
}

template <typename vector_t>
inline
typename vector_t::value_type
VecInftyNorm(const vector_t& v)
{
	return VecMaxNorm(v);
}


///	component-wise product
template <typename vector_t>
inline
void
VecElemProd(vector_t& vOut, const vector_t& v1, const vector_t& v2)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = v1[i] * v2[i];
	}
}

///	component-wise square root
template <typename vector_t>
inline
void
VecElemSqrt(vector_t& vOut, const vector_t& v1)
{
	typedef typename vector_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = sqrt(v1[i]);
	}
}

///	component-wise comparison of two vectors (in the absolute values)
template <typename vector_t>
inline
bool
VecAbsIsLess(const vector_t& v1, const vector_t& v2)
{
	for(typename vector_t::size_type i = 0; i < v1.size(); ++i)
		if (std::abs (v1[i]) >= std::abs (v2[i]))
			return false;
	return true;
}

///	component-wise comparison of a vector (in the absolute values) with a given number
template <typename vector_t>
inline
bool
VecAbsIsLess(const vector_t& v1, const typename vector_t::value_type s)
{
	for(typename vector_t::size_type i = 0; i < v1.size(); ++i)
		if (std::abs (v1[i]) >= s)
			return false;
	return true;
}

}//	end of namespace

#endif /* __H__COMMON__MathVector_FUNCTIONS_COMMON_IMPL__ */

