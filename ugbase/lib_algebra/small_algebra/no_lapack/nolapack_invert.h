/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
