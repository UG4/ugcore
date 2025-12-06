/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__CPU_ALGEBRA__LAPACK_INVERT_H__
#define __H__UG__CPU_ALGEBRA__LAPACK_INVERT_H__

#include "../small_matrix/densematrix.h"
#include "../small_matrix/densevector.h"
#include "../small_matrix/block_dense.h"
#include "lapack.h"
#include <vector>

namespace ug{



template<typename T>
bool InvertNdyn(DenseMatrix<T> &mat)
{
	std::vector<lapack_int> interchange(mat.num_rows());

	int info = getrf(mat.num_rows(), mat.num_cols(), &mat(0,0), mat.num_rows(), &interchange[0]);
	//UG_ASSERT(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had an illegal value"));
	if(info != 0) return false;

	// calc work size
	double worksize; int iWorksize = -1;
	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), &interchange[0], &worksize, iWorksize);
	//UG_ASSERT(info == 0, "");
	iWorksize = static_cast<int>(worksize);

	std::vector<double> work(iWorksize);

	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), &interchange[0], &work[0], iWorksize);
	//UG_ASSERT(info == 0, "");
	if(info != 0) return false;

	return true;
}

template<typename T, size_t TUnknowns>
bool Invert(DenseMatrix<FixedArray2<T, TUnknowns, TUnknowns> > &mat)
{
	lapack_int interchange[TUnknowns];

	int info = getrf(mat.num_rows(), mat.num_cols(), &mat(0,0), mat.num_rows(), interchange);
	UG_COND_THROW(info != 0, "info is " << info << ( info > 0 ? ": Matrix singular in mat(i,i)" : ": i-th argument had an illegal value"));
	if(info != 0) return false;

	// calc work size
	// todo: make this static
	double worksize; int iWorksize = -1;
	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), interchange, &worksize, iWorksize);
	UG_ASSERT(info == 0, "");
	iWorksize = static_cast<int>(worksize);

	std::vector<double> work;
	work.resize(iWorksize);

	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), interchange, &work[0], iWorksize);
	UG_ASSERT(info == 0, "");
	if(info != 0) return false;

	return true;
}


}

#endif
