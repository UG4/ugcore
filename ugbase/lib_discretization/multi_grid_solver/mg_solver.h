/*
 * mg_solver.h
 *
 *  Created on: 07.12.2009
 *      Author: andreasvogel
 */

#ifndef MG_SOLVER_H_
#define MG_SOLVER_H_

// extern includes
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/assemble.h"
#include "lib_discretization/function_spaces/grid_function_space.h"

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
class MultiGridSolver : public ILinearSolver {
	public:
		typedef TDomain domain_type;

		typedef TDiscreteFunction discrete_function_type;

		typedef TAlgebra algebra_type;

		typedef typename TAlgebra::matrix_type matrix_type;

		typedef typename TAlgebra::vector_type vector_type;

	public:

		MultiGridSolver(IAssemble<TAlgebra, discrete_function_type>& ass, domain_type& domain, uint surfaceLevel, uint baseLevel, discrete_function_type& u,
						int maxIter, number tol, int cycle, LinearSolver& Smoother, number damp, int nu1, int nu2, LinearSolver& BaseSolver);

		bool solve(matrix_type& A, vector_type& x, vector_type& b);

		bool print_matrices();

	protected:
		bool lmgc(uint l);


	protected:
		matrix_type** _A;
		matrix_type** _I;
		vector_type* _d;
		vector_type* _t;
		vector_type* _c;

		int _num_cycle;

		IAssemble<TAlgebra, discrete_function_type>* _ass;
		domain_type* _domain;
		discrete_function_type* _u;

		LinearSolver* _Smoother;
		number _damp;
		int _nu1;
		int _nu2;

		LinearSolver* _BaseSolver;
		uint _baseLevel;
		uint _surfaceLevel;

		int _maxIter;
		number _tol;
};

}

#include "mg_solver_impl.hpp"


#endif /* MG_SOLVER_H_ */
