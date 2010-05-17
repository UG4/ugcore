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

		typedef std::vector<index_type> local_index_type;

	public:

		ArneMatrix() : _Matrix(NULL) {};

		bool create(uint nrow, uint ncol);
		bool create(const ArneMatrix& v);
		bool destroy();

		// add, set, get
		bool set(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J);
		bool add(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J);
		bool get(local_matrix_type& mat, const local_index_type& I, const local_index_type& J) const;

		// normalize rows
		bool set_dirichlet_rows(const local_index_type& I);

		// set all entries (only currently allocated memory pattern)
		bool operator= (number w);
		bool set(number w);

		// fix memory pattern
		bool finalize();

		// destructor
		~ArneMatrix();

		// b := A*x (A = this Object)
		bool apply(ArneVector&b, const ArneVector& x);

		// b := A^T * x (A^T = transposed of this object)
		bool applyTransposed(ArneVector&b, const ArneVector& x);

		// b := b - A * x (A = this object)
		bool matmul_minus(ArneVector&b, const ArneVector& x);

		// sizes
		uint row_size() const;
		uint col_size() const;

	// not generic part
	private:
		bool printToFile(const char* filename);

		friend class ArneJacobi;
		friend bool diag_step(const ArneMatrix& A, ArneVector& c, ArneVector& d, number damp);

		typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
		ScalarMatrix* getStorage();
		const ScalarMatrix* getStorage() const;

	private:
		ScalarMatrix* _Matrix;

		friend std::ostream& operator<< (std::ostream& outStream, const ug::ArneMatrix& m);

};

std::ostream& operator<< (std::ostream& outStream, const ug::ArneMatrix& m);

}

#endif /* __H__LIB_ALGEBRA__ARNEMATRIX__ */
