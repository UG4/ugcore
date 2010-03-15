/*
 * arnematrix.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__ARNEMATRIX__
#define __H__LIB_ALGEBRA__ARNEMATRIX__

#include <iostream>
#include "common/common.h"
#include "arnevector.h"
#include "../solver/BoostBlock.hh"
#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"

namespace ug{

class ArneMatrix{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		typedef FlexLocalMatrix local_matrix_type;

		typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;

	public:

	ArneMatrix() {};

	bool create_matrix(int nrow, int ncol);

	bool delete_matrix();

	bool set_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool add_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool set_dirichletrows(int nrows, int* rows);

	bool set(number w);

	bool finalize();

	~ArneMatrix();

	// not generic part

	// b := A*x (A = this Object)
	bool apply(ArneVector&b, ArneVector& x)
	{
		ublas::axpy_prod(*_Matrix, *x.getStorage(), *b.getStorage(), true);
		return true;
	}

	// b := A^T * x (A^T = transposed of this object)
	bool applyTransposed(ArneVector&b, ArneVector& x)
	{
		ublas::axpy_prod(*x.getStorage(), *_Matrix, *b.getStorage(), true);
		return true;
	}

	// b := b - A * x (A = this object)
	bool matmul_minus(ArneVector&b, ArneVector& x)
	{
		*x.getStorage() *= -1.0;
		ublas::axpy_prod(*_Matrix, *x.getStorage(), *b.getStorage(), false);
		*x.getStorage() *= -1.0;
		return true;
	}

	bool printToFile(const char* filename);

	ScalarMatrix* getStorage();

	private:
		ScalarMatrix* _Matrix;

};

}

#endif /* __H__LIB_ALGEBRA__ARNEMATRIX__ */
