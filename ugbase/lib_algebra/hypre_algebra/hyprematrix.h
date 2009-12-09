/*
 * hyprematrix.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__HYPREMATRIX__
#define __H__LIB_ALGEBRA__HYPREMATRIX__

#include <iostream>
#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include "../../common/types.h"

namespace ug{

class HypreMatrix{

	public:
	
	HypreMatrix(){};
	
	bool create_matrix(int nrow, int ncol);

	bool delete_matrix();

	bool set_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool add_values(int nrows, int* ncols, int* rows, int* cols, double* values);

	bool set_dirichletrows(int nrows, int* rows);

	bool set(number w);

	bool finalize();

	~HypreMatrix();

	// not generic part

	bool printToFile(const char* filename);

	HYPRE_IJMatrix getStorage();

	private:
		HYPRE_IJMatrix m_hypreA;

};

}

#endif /* __H__LIB_ALGEBRA__HYPREMATRIX__ */
