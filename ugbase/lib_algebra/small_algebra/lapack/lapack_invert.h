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
	iWorksize = (int)worksize;

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
	iWorksize = (int)worksize;

	std::vector<double> work;
	work.resize(iWorksize);

	info = getri(mat.num_rows(), &mat(0,0), mat.num_rows(), interchange, &work[0], iWorksize);
	UG_ASSERT(info == 0, "");
	if(info != 0) return false;

	return true;
}


}

#endif
