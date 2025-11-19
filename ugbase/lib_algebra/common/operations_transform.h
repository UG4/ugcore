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


#ifndef __H__UG__LIB_ALGEBRA__OPERATIONS_TRANSFORM__
#define __H__UG__LIB_ALGEBRA__OPERATIONS_TRANSFORM__

// here we transform a Template Expression x = X1 [+/- X2 [+/- X3]]  into
// a function like VecScaleAdd (in operations_vec.h) or MatMultAdd (in operations_mat.h),
// where Xi can be Matrix*Vector, double*Matrix*vector, double*vector or vector.
// x += X1 [+/- X2]  and x -= X1 [+/- X2] also possible.
// todo: be careful with x left and right. check.
///////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace ug
{
// 3 Operants
//--------------

//! v = Mv + Mv + v
template<typename vector_t, typename matrix_t, typename vector_t2>
inline void VectorAssign(vector_t &dest, double alpha1, const MatVec_Expression<matrix_t, vector_t> &m1,
					   double alpha2, const MatVec_Expression<matrix_t, vector_t> &m2,
					   double alpha3, const vector_t2 &v3)
{
	MatMultAdd(dest, alpha1*m1.alpha, m1.l, m1.r, alpha2*m2.alpha, m2.l, m2.r, alpha3*getScaling(v3), getVector(v3));
}

//! v = Mv + v + v
template<typename vector_t, typename matrix_t, typename vector_t2, typename vector_t3>
inline void VectorAssign(vector_t &dest, double alpha1, const MatVec_Expression<matrix_t, vector_t> &m1,
					   double alpha2, const vector_t2 &v2,
					   double alpha3, const vector_t3 &v3)
{
	MatMultAdd(dest, alpha1*m1.alpha, m1.l, m1.r, getScaling(v2)*alpha2, getVector(v2), getScaling(v3)*alpha3, getVector(v3));
}

//!
//! v = x + v + Mv  ->  v = x + Mv + v
//! v = v + v + Mv
//! v = Mv + v + Mv
template<typename vector_t, typename T0, typename matrix_t, typename vector_t2>
inline void VectorAssign(vector_t &dest,
					   double alpha0, const T0 &t0,
					   double alpha1, const vector_t2 &v1,
					   double alpha2, const MatVec_Expression<matrix_t, vector_t> &m2)
{
	VectorAssign(dest, alpha0, t0, alpha2, m2, alpha1, v1);
}


//!
//! v = v + Mv + X  ->  v = Mv + v + X
//! v = v + Mv + v
//! v = v + Mv + Mv
template<typename vector_t, typename matrix_t, typename vector_t2, typename T3>
inline void VectorAssign(vector_t &dest,
					   double alpha1, const vector_t2 &v1,
					   double alpha2, const MatVec_Expression<matrix_t, vector_t> &m2,
					   double alpha3, const T3 &t3)
{
	VectorAssign(dest, alpha2, m2, alpha1, v1, alpha3, t3);
}

//! v = v + v + v
template<typename vector_t, typename vector_t1, typename vector_t2, typename vector_t3>
inline void VectorAssign(vector_t &dest,
						double alpha1, const vector_t1 &v1,
						double alpha2, const vector_t2 &v2,
						double alpha3, const vector_t3 &v3)

{
	VecScaleAdd(dest, getScaling(v1)*alpha1, getVector(v1), getScaling(v2)*alpha2, getVector(v2), getScaling(v3)*alpha3, getVector(v3));
}



// 2 Operants
//--------------

//! v = Mv + Mv
template<typename vector_t, typename matrix_t>
inline void VectorAssign(vector_t &dest,	double alpha1, const MatVec_Expression<matrix_t, vector_t> &m1,
										double alpha2, const MatVec_Expression<matrix_t, vector_t> &m2)
{
	MatMultAdd(dest, alpha1*m1.alpha, m1.l, m1.r, alpha2*m2.alpha, m2.l, m2.r);
}


//! v = Mv + v
template<typename vector_t, typename matrix_t, typename vector_t2>
inline void VectorAssign(vector_t &dest, double alpha1, const MatVec_Expression<matrix_t, vector_t> &m1, double alpha2, const vector_t2 &v2)
{
	MatMultAdd(dest, alpha1*m1.alpha, m1.l, m1.r, getScaling(v2)*alpha2, getVector(v2));
}

//! v = v + Mv
template<typename vector_t, typename matrix_t>
inline void VectorAssign(vector_t &dest, double alpha1, const vector_t &v1, double alpha2, const MatVec_Expression<matrix_t, vector_t> &m2)
{
	VectorAssign(dest, alpha2, m2, alpha1, v1);
}

//! v = v + v
template<typename vector_t, typename vector_t2, typename vector_t3>
inline void VectorAssign(vector_t &dest, double alpha1, const vector_t2 &v1, double alpha2, const vector_t3 &v2)
{
	VecScaleAdd(dest, getScaling(v1)*alpha1, getVector(v1), getScaling(v2)*alpha2, getVector(v2));
}




///
// 1 Operant

//! v = Mv
template<typename vector_t, typename matrix_t>
inline void VectorAssign(vector_t &dest, const MatVec_Expression<matrix_t, vector_t> &m1)
{
	MatMult(dest, m1.alpha, m1.l, m1.r);
}


//! v = v
template<typename vector_t, typename vector_t2, typename vector_t3>
inline void VectorAssign(vector_t &dest, const vector_t2 &v1)
{
	VecScaleAssign(dest, getScaling(v1), getVector(v1));
}


/////
// transform all AlphaMatVec_X_Expression of the form a+b, (a+b)+c or a+(b+c)
//
template<typename vector_t, typename T1, typename operation, typename T2, typename T3>
inline void VectorAssign(vector_t &dest, double alpha1, const AlphaMatVec_X_Expression<T1, operation, T2> &t1, double alpha2, const T3 &t2)
{
	VectorAssign(dest, alpha1, t1.cast().l, operation::is_add() ? alpha1 : -alpha1, t1.cast().r, alpha2, t2);
}

template<typename vector_t, typename T1, typename T2, typename operation, typename T3>
inline void VectorAssign(vector_t &dest, double alpha1, const T1 &t1, double alpha2, const AlphaMatVec_X_Expression<T2, operation, T3> &t2)
{
	VectorAssign(dest, alpha1, t1, alpha2, t2.l, operation::is_add() ? alpha2 : -alpha2, t2.r);
}

template<typename vector_t, typename T1, typename operation, typename T2>
inline void VectorAssign(vector_t &dest, const AlphaMatVec_X_Expression<T1, operation, T2 > &t)
{
	VectorAssign(dest, 1.0, t.l, operation::is_add() ? 1.0 : -1.0, t.r);
}


/////////////////////////////////////////////////////////////////////////////////
//! transforms x += X1 into x = (1.0)*X1 + (1.0)*x
template<typename vector_t, typename T1>
inline void VectorAdd(vector_t &dest, const T1 &t1)
{
	VectorAssign<vector_t, T1, vector_t>(dest, 1.0, t1, 1.0, dest);
}

//! transforms x -= X1 into x = (-1.0)*X1 + (1.0)*x
template<typename vector_t, typename T1>
inline void VectorSub(vector_t &dest, const T1 &t1)
{
	VectorAssign(dest, -1.0, t1, 1.0, dest);
}


} // namespace ug

#endif