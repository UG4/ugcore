/*
 * arnevector.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__ARNEVECTOR__
#define __H__LIB_ALGEBRA__ARNEVECTOR__

#include "../solver/BoostBlock.hh"

#include "../../common/types.h"
#include <cmath>
#include "assert.h"

namespace ug{

class ArneVector{
	typedef ublas::vector<double> ScalarVector;

	public:
	bool create_vector(int nentries);

	bool delete_vector();

	bool set_values(int nvalues, int* indices, double* values);

	bool add_values(int nvalues, int* indices, double* values);

	bool get_values(int nvalues, int* indices, double* values) const;

	bool finalize();

	bool printToFile(const char* filename);

	~ArneVector();

	ArneVector& operator+= (const ArneVector& v);

	ArneVector& operator-= (const ArneVector& v);

	number norm2();

	bool set(number w);

	int length();

	ArneVector::ScalarVector* getStorage();

	private:
	ScalarVector* _Vector;

};

}

#endif /* __H__LIB_ALGEBRA__ARNEVECTOR__ */
