/*
 * matrix.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__MATRIX__
#define __H__LIB_ALGEBRA__MATRIX__

#include <iostream>
#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include "../../common/types.h"

namespace ug{

class Matrix{

	public:
	bool create_matrix(int nrow, int ncol);

	bool delete_matrix();

	bool set_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool add_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool set_dirichletrows(int nrows, int* rows);

	bool finalize();

	bool printToFile(const char* filename);

	bool set(number w);

	HYPRE_IJMatrix getStorage();

	~Matrix();

	private:
		HYPRE_IJMatrix m_hypreA;

};

}

#endif /* __H__LIB_ALGEBRA__MATRIX__ */
