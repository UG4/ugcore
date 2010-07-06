/*
 * ublas_vector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_VECTOR__
#define __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_VECTOR__

#include "solver/BoostBlock.hh"

#include <cmath>
#include <vector>
#include "assert.h"

#include "common/common.h"

#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"

namespace ug{

class UblasMatrix;
class UblasJacobi;

class UblasVector{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		typedef FlexLocalVector local_vector_type;

		typedef std::vector<index_type> local_index_type;

	public:
		// constructor
		UblasVector() : m_pVector(NULL){};

		// create and destroy
		bool create(size_t nentries);

		// create as copy of other vector
		bool create(const UblasVector& v);

		// destroy
		bool destroy();

		// add, set, get a local function
		bool set(const local_vector_type& u, const local_index_type& ind);
		bool add(const local_vector_type& u, const local_index_type& ind);
		bool get(local_vector_type& u, const local_index_type& ind) const;

		// fix memory pattern
		bool finalize();

		// operations with other vectors
		UblasVector& operator+= (const UblasVector& v);
		UblasVector& operator-= (const UblasVector& v);
		UblasVector& operator= (const UblasVector& v);

		// scalar product
		number operator *(const UblasVector& v);
		number dotprod(const UblasVector &v) { return this->operator*(v);}

		// norms
		number one_norm() const;
		number two_norm() const;

		// set vector to constant value
		bool set(number w);
		bool operator= (number w);
		bool operator*= (number w);

		// number of Blocks
		size_t size() const;

		// write to file
		bool printToFile(const char* filename) const;

		// destructor
		~UblasVector();

	private:
		// disallow copy constructor
		UblasVector(const UblasVector& v);

	private:
		friend class UblasMatrix;
		friend class UblasJacobi;
		friend bool diag_step(const UblasMatrix& A, UblasVector& c, UblasVector& d, number damp);

		typedef ublas::vector<double> ScalarVector;

		UblasVector::ScalarVector* getStorage();
		const UblasVector::ScalarVector* getStorage() const;

	private:
		ScalarVector* m_pVector;

	friend std::ostream& operator<< (std::ostream& outStream, const ug::UblasVector& v);
	public:
		void p() const { std::cout << *this; }

};

std::ostream& operator<< (std::ostream& outStream, const ug::UblasVector& v);

std::ostream& operator<< (std::ostream& outStream, const ug::UblasVector::local_index_type& ind);


}

#endif /* __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_VECTOR__ */
