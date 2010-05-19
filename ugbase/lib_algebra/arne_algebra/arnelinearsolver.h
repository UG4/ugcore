/*
 * linearsolver.h
 *
 *  Created on: 04.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__ARNELINEARSOLVER__
#define __H__LIB_DISCRETIZATION__ARNELINEARSOLVER__

#include "arnevector.h"
#include "arnematrix.h"
#include "../solver/ExtalgSolver.hh"
#include "../solver/BoostBlock.hh"


namespace ug{

class ArneJacobi{

	public:
	ArneJacobi(int maxIter = 1e6, double tol = 1e-15, bool verbose = true);

	bool solve(ArneMatrix& A, ArneVector& x, ArneVector& b);

	bool step(ArneMatrix& A, ArneVector& c, ArneVector& d, number damp);

	~ArneJacobi();

	private:
		double _tol;
		int _maxIter;
		bool _verbose;
};

bool diag_step(const ArneMatrix& A, ArneVector& c, ArneVector& d, number damp);

}

#endif /* __H__LIB_DISCRETIZATION__ARNELINEARSOLVER__ */
