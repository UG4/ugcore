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
#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"

namespace ug{

class HypreMatrix{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		typedef FlexLocalMatrix local_matrix_type;

		typedef std::vector<index_type> local_index_type;

	public:

		HypreMatrix(){};

		bool create(int nrow, int ncol);
		bool destroy();

		bool set(const local_matrix_type& mat, local_index_type& I, local_index_type& J);
		bool add(const local_matrix_type& mat, local_index_type& I, local_index_type& J);

		bool set_dirichletrows(int nrows, int* rows);

		bool set(number w);

		bool finalize();

		~HypreMatrix();

	// not generic part
	private:
		bool set_values(int nrows, int* ncols, int* rows, int* cols, double* values);
		bool add_values(int nrows, int* ncols, int* rows, int* cols, double* values);

		bool printToFile(const char* filename);

		friend class HYPREboomerAMG;
		HYPRE_IJMatrix getStorage();

	private:
		HYPRE_IJMatrix m_hypreA;

};

}

#endif /* __H__LIB_ALGEBRA__HYPREMATRIX__ */
