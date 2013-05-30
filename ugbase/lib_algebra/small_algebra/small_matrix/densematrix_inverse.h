#ifndef __H__UG__CPU_ALGEBRA__DENSEMATRIX_SMALL_INVERSE_H__
#define __H__UG__CPU_ALGEBRA__DENSEMATRIX_SMALL_INVERSE_H__

#include "densematrix.h"
#include "densevector.h"
#include "block_dense.h"
#include "common/common.h"
#include <algorithm>

//
namespace ug {

/// \addtogroup small_algebra
/// \{

//////////////////////////////////////////////////////
// 1x1

template<typename T>
inline bool GetInverse1(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	//UG_ASSERT(mat(0,0)!=0.0, "Determinant zero, cannot invert matrix.");
	if(mat(0,0) == 0.0) return false;
	inv(0,0) = 1/mat(0,0);
	return true;
}

template<typename T>
bool Invert1(DenseMatrix<T> &mat)
{
	//UG_ASSERT(mat(0,0)!=0.0, "Determinant zero, cannot invert matrix.");
	if(mat(0,0) == 0.0) return false;
	mat(0,0) = 1/mat(0,0);
	return true;
};

inline bool GetInverse(DenseMatrix<FixedArray2<double, 1, 1> > &inv, const DenseMatrix<FixedArray2<double, 1, 1> > &mat)
{
	return GetInverse1(inv, mat);
}

inline bool Invert(DenseMatrix< FixedArray2<double, 1, 1> > &mat)
{
	return Invert1(mat);
}

template<typename vector_t, typename matrix_t>
inline bool InverseMatMult1(DenseVector<vector_t> &dest, double beta1,
		const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	if(A1(0,0) == 0.0) return false;
	UG_ASSERT(&dest != &w1, "");
	dest[0] = beta1*w1[0]/A1(0,0);
	return true;
}

inline bool InverseMatMult(DenseVector< FixedArray1<double, 1> > &dest, double beta1,
		const DenseMatrix< FixedArray2<double, 1, 1> > &A1, const DenseVector< FixedArray1<double, 1> > &w1)
{
	return InverseMatMult1(dest, beta1, A1, w1);
}

//////////////////////
// 2x2

template<typename T>
inline double GetDet2(const DenseMatrix<T> &mat)
{
	UG_ASSERT(mat.num_rows() == 2 && mat.num_cols() == 2, "only for 2x2-matrices");
	return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
}

template<typename T>
inline bool GetInverse2(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	UG_ASSERT(&inv != &mat, "inv and mat have to be different. Otherwise use Invert/Invert2");
	double invdet = GetDet2(mat);
	//UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;
	inv(0,0) = mat(1,1) * invdet;
	inv(1,1) = mat(0,0) * invdet;
	inv(0,1) = mat(0,1) * -invdet;
	inv(1,0) = mat(1,0) * -invdet;
	return true;
}

template<typename T>
bool Invert2(DenseMatrix<T> &mat)
{
	double invdet = GetDet2(mat);
	//UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;

	std::swap(mat(0,0), mat(1,1));

	mat(0,0) *= invdet;
	mat(0,1) *= -invdet;
	mat(1,0) *= -invdet;
	mat(1,1) *= invdet;
	return true;
};


inline bool GetInverse(DenseMatrix<FixedArray2<double, 2, 2> > &inv, const DenseMatrix<FixedArray2<double, 2, 2> > &mat)
{
	return GetInverse2(inv, mat);
}

inline bool Invert(DenseMatrix< FixedArray2<double, 2, 2> > &mat)
{
	return Invert2(mat);
}

template<typename vector_t, typename matrix_t>
inline bool InverseMatMult2(DenseVector<vector_t> &dest, double beta,
		const DenseMatrix<matrix_t> &mat, const DenseVector<vector_t> &vec)
{
	number det = GetDet2(mat);
	//UG_ASSERT(det != 0, "Determinant zero, cannot invert matrix.");
	UG_ASSERT(&dest != &vec, "");
	if(det == 0.0) return false;
	dest[0] = beta * (mat(1,1)*vec[0] - mat(0,1)*vec[1]) / det;
	dest[1] = beta * (-mat(1,0)*vec[0] + mat(0,0)*vec[1]) / det;
	return true;
}

template<typename T>
inline bool InverseMatMult(DenseVector< FixedArray1<double, 2> > &dest, double beta,
		const DenseMatrix< FixedArray2<double, 2, 2> > &mat, const DenseVector< FixedArray1<double, 2> > &vec)
{
	return InverseMatMult2(dest, beta, mat, vec);
}

//////////////////////
// 3x3


template<typename T>
inline double GetDet3(const DenseMatrix<T> &mat)
{
	UG_ASSERT(mat.num_rows() == 3 && mat.num_cols() == 3, "only for 3x3-matrices");
	return 	mat(0,0)*mat(1,1)*mat(2,2) + mat(0,1)*mat(1,2)*mat(2,0) + mat(0,2)*mat(1,0)*mat(2,1)
			- mat(0,0)*mat(1,2)*mat(2,1) - mat(0,1)*mat(1,0)*mat(2,2) - mat(0,2)*mat(1,1)*mat(2,0);
}

template<typename T>
inline bool GetInverse3(DenseMatrix<T> &inv, const DenseMatrix<T> &mat)
{
	UG_ASSERT(&inv != &mat, "inv and mat have to be different. Otherwise use Invert/Invert3");
	double invdet = GetDet3(mat);
	//UG_ASSERT(invdet != 0, "Determinant zero, cannot invert matrix.");
	if(invdet == 0.0) return false;
	invdet = 1.0/invdet;

	inv(0,0) = ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) * invdet;
	inv(0,1) = (-mat(0,1)*mat(2,2) + mat(0,2)*mat(2,1)) * invdet;
	inv(0,2) = ( mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1)) * invdet;
	inv(1,0) = (-mat(1,0)*mat(2,2) + mat(1,2)*mat(2,0)) * invdet;
	inv(1,1) = ( mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0)) * invdet;
	inv(1,2) = (-mat(0,0)*mat(1,2) + mat(0,2)*mat(1,0)) * invdet;
	inv(2,0) = ( mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0)) * invdet;
	inv(2,1) = (-mat(0,0)*mat(2,1) + mat(0,1)*mat(2,0)) * invdet;
	inv(2,2) = ( mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0)) * invdet;
	return true;
}

inline bool Invert3(DenseMatrix<FixedArray2<double, 3, 3> > & mat)
{
	DenseMatrix<FixedArray2<double, 3, 3> > inv;
	if(GetInverse3(inv, mat) == false) return false;
	mat = inv;
	return true;
}

inline bool Invert3(DenseMatrix<VariableArray2<double> > & mat)
{
	DenseMatrix<VariableArray2<double> > inv;
	inv.resize(3,3);
	if(GetInverse3(inv, mat) == false) return false;
	mat = inv;
	return true;
}

inline bool GetInverse(DenseMatrix<FixedArray2<double, 3, 3> > &inv, const DenseMatrix<FixedArray2<double, 3, 3> > &mat)
{
	return GetInverse3(inv, mat);
}

inline bool Invert(DenseMatrix< FixedArray2<double, 3, 3> > &mat)
{
	return Invert3(mat);
}

template<typename vector_t, typename matrix_t>
inline bool InverseMatMult3(DenseVector<vector_t> &dest, double beta,
		const DenseMatrix<matrix_t> &mat, const DenseVector<vector_t> &vec)
{
	number det = GetDet3(mat);
	//UG_ASSERT(det != 0, "Determinant zero, cannot invert matrix.");
	UG_ASSERT(&dest != &vec, "");
	if(det == 0.0) return false;
	dest[0] = ( ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) *vec[0] +
				(-mat(0,1)*mat(2,2) + mat(0,2)*mat(2,1)) *vec[1] +
				( mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1)) *vec[2] ) * beta / det;
	dest[1] = ( (-mat(1,0)*mat(2,2) + mat(1,2)*mat(2,0)) * vec[0] +
				( mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0)) * vec[1] +
				(-mat(0,0)*mat(1,2) + mat(0,2)*mat(1,0)) * vec[2] ) * beta / det;
	dest[2] = ( ( mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0)) * vec[0] +
				(-mat(0,0)*mat(2,1) + mat(0,1)*mat(2,0)) * vec[1] +
				( mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0)) * vec[2] ) * beta / det;
	return true;
}

template<typename T>
inline bool InverseMatMult(DenseVector< FixedArray1<double, 3> > &dest, double beta,
		const DenseMatrix< FixedArray2<double, 3, 3> > &mat, const DenseVector< FixedArray1<double, 3> > &vec)
{
	return InverseMatMult3(dest, beta, mat, vec);
}
//////////////////////



//! calculates dest = beta1 * A1;
template<typename vector_t, typename matrix_t>
inline void MatMult(DenseVector<vector_t> &dest,
		const number &beta1, const DenseMatrixInverse<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	if(beta1 == 1.0)
	{
		dest = w1;
		A1.apply(dest);
	}
	else
	{
		DenseVector<vector_t> tmp;
		tmp = w1;
		A1.apply(tmp);
		VecScaleAssign(dest, beta1, dest);
	}
}


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(DenseVector<vector_t> &dest,
		const number &alpha1, const DenseVector<vector_t> &v1,
		const number &beta1, const DenseMatrixInverse<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	// todo: use dynamic stack here
	DenseVector<vector_t> tmp;
	tmp = w1;
	A1.apply(tmp);
	VecScaleAdd(dest, alpha1, v1, beta1, tmp);
}



template<typename T, eMatrixOrdering TOrdering>
inline bool GetInverse(DenseMatrixInverse<VariableArray2<T, TOrdering> > &inv, const DenseMatrix<VariableArray2<T, TOrdering> > &mat)
{
	return inv.set_as_inverse_of(mat);
}

template<typename T, size_t TBlockSize, eMatrixOrdering TOrdering>
inline bool GetInverse(DenseMatrixInverse<FixedArray2<T, TBlockSize, TBlockSize, TOrdering> > &inv, const DenseMatrix<FixedArray2<T, TBlockSize, TBlockSize, TOrdering> > &mat)
{
	return inv.set_as_inverse_of(mat);
}

//////////////////////
template<typename T>
inline bool Invert(DenseMatrix<T> &mat)
{
	switch(mat.num_rows())
	{
		case 1: return Invert1(mat);
		case 2: return Invert2(mat);
		case 3: return Invert3(mat);
		default: return InvertNdyn(mat);
	}
}

template<typename vector_t, typename matrix_t>
inline bool InverseMatMultN(DenseVector<vector_t> &dest, double beta,
		const DenseMatrix<matrix_t> &mat, const DenseVector<vector_t> &vec)
{
	typename block_traits<DenseMatrix<matrix_t> >::inverse_type inv;
	if(!GetInverse(inv, mat)) return false;
	MatMult(dest, beta, inv, vec);
	return true;
}


template<typename vector_t, typename matrix_t>
inline bool InverseMatMult(DenseVector<vector_t> &dest, double beta,
		const DenseMatrix<matrix_t> &mat, const DenseVector<vector_t> &vec)
{
	switch(mat.num_rows())
	{
		case 1: return InverseMatMult1(dest, beta, mat, vec);
		case 2: return InverseMatMult2(dest, beta, mat, vec);
		case 3: return InverseMatMult3(dest, beta, mat, vec);
		default: return InverseMatMultN(dest, beta, mat, vec);
	}
}

///////////////////////////////////////////

// end group small_algebra
/// \}

}

#endif // __H__UG__CPU_ALGEBRA__DENSEMATRIX_SMALL_INVERSE_H__
