
#ifndef SOLVE_DEFICIT_H_
#define SOLVE_DEFICIT_H_

#include "lib_algebra/small_algebra/small_matrix/densematrix.h"
#include <vector>

namespace ug{


bool SolveDeficit(DenseMatrix< VariableArray2<double> > &A,
		DenseVector<VariableArray1<double> > &x, DenseVector<VariableArray1<double> > &rhs, double deficitTolerance);

bool Decomp(DenseMatrix< VariableArray2<double> > &A, DenseVector<VariableArray1<double> > &rhs,
		size_t &iNonNullRows, double deficitTolerance);
}

#endif /* SOLVE_DEFICIT_H_ */
