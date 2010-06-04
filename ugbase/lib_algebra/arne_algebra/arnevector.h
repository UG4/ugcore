/*
 * arnevector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__ARNEVECTOR__
#define __H__LIB_ALGEBRA__ARNEVECTOR__

#include "../solver/BoostBlock.hh"

#include <cmath>
#include <vector>
#include "assert.h"

#include "common/common.h"

#include "lib_algebra/multi_index/multi_indices.h"
#include "lib_algebra/local_matrix_vector/flex_local_matrix_vector.h"

namespace ug{

class ArneMatrix;
class ArneJacobi;

class ArneVector{
	public:
		// index_type
		typedef MultiIndex<1> index_type;

		typedef FlexLocalVector local_vector_type;

		typedef std::vector<index_type> local_index_type;

	public:
		// constructor
		ArneVector() {};

		// create and destroy
		bool create(uint nentries);

		// create as copy of other vector
		bool create(const ArneVector& v);

		// destroy
		bool destroy();

		// add, set, get a local function
		bool set(const local_vector_type& u, const local_index_type& ind);
		bool add(const local_vector_type& u, const local_index_type& ind);
		bool get(local_vector_type& u, const local_index_type& ind) const;

		// fix memory pattern
		bool finalize();

		// operations with other vectors
		ArneVector& operator+= (const ArneVector& v);
		ArneVector& operator-= (const ArneVector& v);
		ArneVector& operator= (const ArneVector& v);

		// scalar product
		number operator *(const ArneVector& v);

		// norms
		number one_norm() const;
		number two_norm() const;

		// set vector to constant value
		bool set(number w);
		bool operator= (number w);
		bool operator*= (number w);

		// number of Blocks
		uint size() const;

		// write to file
		bool printToFile(const char* filename) const;

		// destructor
		~ArneVector();

	private:
		// disallow copy constructor
		ArneVector(const ArneVector& v);

	private:
		friend class ArneMatrix;
		friend class ArneJacobi;
		friend bool diag_step(const ArneMatrix& A, ArneVector& c, ArneVector& d, number damp);

		typedef ublas::vector<double> ScalarVector;

		ArneVector::ScalarVector* getStorage();
		const ArneVector::ScalarVector* getStorage() const;

	private:
		ScalarVector* _Vector;

	friend std::ostream& operator<< (std::ostream& outStream, const ug::ArneVector& v);
	public:
		void p() const { std::cout << *this; }

};

std::ostream& operator<< (std::ostream& outStream, const ug::ArneVector& v);

std::ostream& operator<< (std::ostream& outStream, const ug::ArneVector::local_index_type& ind);


}

#endif /* __H__LIB_ALGEBRA__ARNEVECTOR__ */
