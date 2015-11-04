/*
 * \file operations_transform.h
 *
 * \author Martin Rupp
 *
 * \date 29.09.2010
 *
 * Goethe-Center for Scientific Computing 2010.
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

#endif /* __H__UG__LIB_ALGEBRA__OPERATIONS_TRANSFORM__ */
