/*
 * vector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__HYPRE_ALGEBRA__HYPREVECTOR__
#define __H__LIB_ALGEBRA__HYPRE_ALGEBRA__HYPREVECTOR__

#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <cmath>
#include "assert.h"

// other ug4 modules
#include "common/common.h"

// library intern header
#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"


namespace ug{

class HypreMatrix;
class HYPREboomerAMG;

class HypreVector{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		// local vector type
		typedef  FlexLocalVector local_vector_type;

		typedef std::vector<index_type> local_index_type;

	public:
		bool create(size_t nentries);
		bool destroy();

		bool set(const local_vector_type& u, local_index_type& ind);
		bool add(const local_vector_type& u, local_index_type& ind);
		bool get(local_vector_type& u, local_index_type& ind) const;

		HypreVector& operator+= (const HypreVector& v);
		HypreVector& operator-= (const HypreVector& v);

		bool finalize();

		number one_norm() const;
		number two_norm() const;

		bool set(number w);
		bool operator= (number w);
		bool operator*= (number w);

		size_t size() const;

		~HypreVector();

	private:
		bool set_values(int nvalues, int* indices, double* values);
		bool add_values(int nvalues, int* indices, double* values);
		bool get_values(int nvalues, int* indices, double* values) const;

		bool printToFile(const char* filename);

		friend class HypreMatrix;
		friend class HYPREboomerAMG;
		HYPRE_IJVector getStorage();

	private:
		HYPRE_IJVector m_hyprex;


};

}

#endif /* __H__LIB_ALGEBRA__HYPRE_ALGEBRA__HYPREVECTOR__ */
