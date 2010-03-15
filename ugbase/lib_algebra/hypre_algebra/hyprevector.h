/*
 * vector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__HYPREVECTOR__
#define __H__LIB_ALGEBRA__HYPREVECTOR__

#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <cmath>
#include "assert.h"

// other ug4 modules
#include "../../common/types.h"

// library intern header
#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"


namespace ug{

class HypreVector{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		// local vector type
		typedef  FlexLocalVector local_vector_type;

		typedef std::vector<index_type> local_index_type;

	public:
	bool create_vector(int nentries);

	bool delete_vector();

	bool set_values(int nvalues, int* indices, double* values);
	bool add_values(int nvalues, int* indices, double* values);
	bool get_values(int nvalues, int* indices, double* values) const;

	bool set(const local_vector_type& u, local_index_type& ind);
	bool add(const local_vector_type& u, local_index_type& ind);
	bool get(local_vector_type& u, local_index_type& ind) const;



	bool finalize();

	bool printToFile(const char* filename);

	HYPRE_IJVector getStorage();

	~HypreVector();

	HypreVector& operator+= (const HypreVector& v);

	HypreVector& operator-= (const HypreVector& v);

	number norm2();

	bool set(number w);

	int length();


	private:
		HYPRE_IJVector m_hyprex;


};

}

#endif /* __H__LIB_ALGEBRA__HYPREVECTOR__ */
