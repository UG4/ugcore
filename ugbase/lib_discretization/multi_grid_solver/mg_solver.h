/*
 * mg_solver.h
 *
 *  Created on: 07.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__

// extern includes
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/assemble.h"
#include "lib_discretization/function_spaces/grid_function_space.h"
#include "lib_discretization/linear_operator/transfer_operator.h"

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
class MultiGridCycle : public IIterativeStep<TAlgebra> {
	public:
		typedef TDomain domain_type;

		typedef TDiscreteFunction discrete_function_type;

		typedef TAlgebra algebra_type;

		typedef typename TAlgebra::matrix_type matrix_type;

		typedef typename TAlgebra::vector_type vector_type;

	public:
		// constructore
		MultiGridCycle(	IAssemble<algebra_type, discrete_function_type>& ass, domain_type& domain, discrete_function_type& u,
						uint surfaceLevel, uint baseLevel, int cycle_type,
						TransferOperator<TDomain, TAlgebra, typename TDiscreteFunction::dof_manager_type>& transferOperator,
						IIterativeStep<TAlgebra>& smoother, int nu1, int nu2, ILinearSolver<TAlgebra>& baseSolver);

		// This functions allocates the Memory for the solver
		// and assembles coarse grid matrices using 'ass'
		bool prepare();

		// This function performes one multi-grid cycle step
		// A correction c is returned as well as the updated defect d := d - A*c
		// The Matrix A remains unchanged
		bool step(matrix_type& A, vector_type& c, vector_type &d);

		// This functions deallocates the Memory for the solver
		bool finish();

 	protected:
		bool lmgc(matrix_type* A[], discrete_function_type& c, discrete_function_type& d, uint l);


	protected:
		IAssemble<TAlgebra, discrete_function_type>& m_ass;
		domain_type& m_domain;
		discrete_function_type& m_u;

		uint m_surfaceLevel;
		uint m_baseLevel;
		int m_cycle_type;

		IIterativeStep<TAlgebra>& m_smoother;
		int m_nu1;
		int m_nu2;

		ILinearSolver<TAlgebra>& m_baseSolver;

		matrix_type** m_A;
		matrix_type** m_I;
		discrete_function_type& m_c;
		discrete_function_type& m_d;
		discrete_function_type& m_t;

		TransferOperator<TDomain, TAlgebra, typename  TDiscreteFunction::dof_manager_type>& m_trans;
};

}

#include "mg_solver_impl.hpp"


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__ */
