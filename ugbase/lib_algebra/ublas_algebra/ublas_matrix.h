/*
 * ublas_matrix.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_MATRIX__
#define __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_MATRIX__

#include <iostream>
#include "common/common.h"

#include "ublas_vector.h"
#include "solver/BoostBlock.hh"

namespace ug{

class UblasMatrix{
	public:

		UblasMatrix() : m_pMatrix(NULL) {};

		bool create(size_t nrow, size_t ncol);
		bool create(const UblasMatrix& v);
		bool destroy();

		// TODO: rework this. Interface has changed
/*		// add, set, get
		bool set(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J);
		bool add(const local_matrix_type& mat, const local_index_type& I, const local_index_type& J);
		bool get(local_matrix_type& mat, const local_index_type& I, const local_index_type& J) const;

		// normalize rows
		bool set_dirichlet_rows(const local_index_type& I);
*/
		// set all entries (only currently allocated memory pattern)
		bool operator= (number w);
		bool set(number w);

		// fix memory pattern
		bool finalize();

		// destructor
		~UblasMatrix();

		// b := A*x (A = this Object)
		bool apply(UblasVector&b, const UblasVector& x);

		// b := A^T * x (A^T = transposed of this object)
		bool apply_transposed(UblasVector&b, const UblasVector& x);

		// b := b - A * x (A = this object)
		bool matmul_minus(UblasVector&b, const UblasVector& x);

		// sizes
		size_t row_size() const;
		size_t col_size() const;

	// not generic part
	private:
		bool printToFile(const char* filename);

		friend class UblasJacobi;
		friend bool diag_step(const UblasMatrix& A, UblasVector& c, UblasVector& d, number damp);

		typedef ublas::compressed_matrix<double, ublas::row_major> ScalarMatrix;
		ScalarMatrix* getStorage();
		const ScalarMatrix* getStorage() const;

	private:
		ScalarMatrix* m_pMatrix;

		friend std::ostream& operator<< (std::ostream& outStream, const ug::UblasMatrix& m);

	public:
		void p() const { std::cout << *this; }

};

std::ostream& operator<< (std::ostream& outStream, const ug::UblasMatrix& m);

}

#endif /* __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_MATRIX__ */
