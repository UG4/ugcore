/*
 * ublas_linearsolver.h
 *
 *  Created on: 04.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_LINEARSOLVER__
#define __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_LINEARSOLVER__

#include "ublas_vector.h"
#include "ublas_matrix.h"

#include "solver/ExtalgSolver.hh"
#include "solver/BoostBlock.hh"


namespace ug{

class UblasJacobi{

	public:
	UblasJacobi(int maxIter = 1e6, double tol = 1e-15, bool verbose = true);

	bool solve(UblasMatrix& A, UblasVector& x, UblasVector& b);

	bool step(UblasMatrix& A, UblasVector& c, UblasVector& d, number damp);

	~UblasJacobi();

	private:
		double m_tol;
		int m_maxIter;
		bool m_verbose;
};

bool diag_step(const UblasMatrix& A, UblasVector& c, UblasVector& d, number damp);

}

#endif /* __H__LIB_ALGEBRA__UBLAS_ALGEBRA__UBLAS_LINEARSOLVER__ */
