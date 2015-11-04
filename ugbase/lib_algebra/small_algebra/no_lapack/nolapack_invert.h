
#ifndef __H__UG__CPU_ALGEBRA__NOLAPACK_INVERT_H_
#define __H__UG__CPU_ALGEBRA__NOLAPACK_INVERT_H_

#include "../small_matrix/densematrix.h"
#include "../small_matrix/densevector.h"
#include "../small_matrix/block_dense.h"

namespace ug{

template<typename TMatrix, typename TInverseMatrixType, typename TVector>
void InverseFromInverseType(TMatrix &mat, TInverseMatrixType &inv, TVector &x, TVector &e)
{
	x.resize(mat.num_rows());
	e.resize(mat.num_rows());
	for(size_t c=0; c<mat.num_rows(); c++)
	{
		e = 0.0;
		e[c] = 1.0;
		MatMult(x, 1.0, inv, e);
		for(size_t r=0; r<mat.num_cols(); r++)
			mat(r, c) = x[r];
	}

}

template<typename TMatrix, typename TVector>
void InverseFromInverseType(TMatrix &mat, TMatrix &inv, TVector &x, TVector &e)
{
	mat = inv;
}

template<typename T>
bool InvertNdyn(DenseMatrix<T> &mat)
{
	typename block_traits<DenseMatrix<T> >::inverse_type inv;
	if(!GetInverse(inv, mat)) return false;

	DenseVector<VariableArray1<typename DenseMatrix<T>::value_type > > e, x;

	InverseFromInverseType(mat, inv, x, e);

	return true;
}

template<typename T, size_t TUnknowns>
bool Invert(DenseMatrix<FixedArray2<T, TUnknowns, TUnknowns> > &mat)
{
	typename block_traits<DenseMatrix<T> >::inverse_type inv;
	if(!GetInverse(inv, mat)) return false;

	DenseVector<FixedArray1<typename DenseMatrix<T>::value_type , TUnknowns > > e, x;
	InverseFromInverseType(mat, inv, x, e);

	return true;
}


}


#endif /* NOLAPACK_INVERT_H_ */
