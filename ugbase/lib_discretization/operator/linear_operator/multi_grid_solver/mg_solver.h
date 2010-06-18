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
#include "lib_discretization/lib_discretization.h"

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
class AssembledMultiGridCycle :
	public ILinearizedIteratorOperator<typename TApproximationSpace::surface_function_type, typename TApproximationSpace::surface_function_type> {
	public:
		typedef typename TApproximationSpace::domain_type phys_domain_type;

		typedef typename TApproximationSpace::surface_function_type surface_function_type;

		typedef typename TApproximationSpace::level_function_type level_function_type;

		typedef TAlgebra algebra_type;

		typedef typename algebra_type::matrix_type matrix_type;

		typedef typename algebra_type::vector_type vector_type;

		typedef AssembledLinearizedOperator<level_function_type> level_operator_type;

		typedef ProlongationOperator<level_function_type> prolongation_operator_type;

		typedef ProjectionOperator<level_function_type> projection_operator_type;

	private:
		typedef TApproximationSpace approximation_space_type;

		typedef ILinearizedOperatorInverse<level_function_type, level_function_type> base_solver_type;

		typedef ILinearizedIteratorOperator<level_function_type, level_function_type> smoother_type;

	public:
		// constructore
		AssembledMultiGridCycle(	IAssemble<level_function_type, algebra_type>& ass, approximation_space_type& approxSpace,
									uint surfaceLevel, uint baseLevel, int cycle_type,
									smoother_type& smoother, int nu1, int nu2, base_solver_type& baseSolver, bool grid_changes = true);

		bool init(ILinearizedOperator<surface_function_type,surface_function_type>& A);

		// This functions allocates the Memory for the solver
		// and assembles coarse grid matrices using 'ass'
		bool prepare(surface_function_type &u, surface_function_type& d, surface_function_type &c);

		// This function performes one multi-grid cycle step
		// A correction c is returned as well as the updated defect d := d - A*c
		bool apply(surface_function_type& d, surface_function_type &c);

		// This functions deallocates the Memory for the solver
		~AssembledMultiGridCycle();

 	protected:
 		// smooth on level l, restrict defect, call lmgc (..., l-1) and interpolate correction
		bool lmgc(uint l);

		// smmoth help function: perform smoothing on level l, nu times
		bool smooth(level_function_type& d, level_function_type& c, uint l, int nu);

		bool allocate_memory();
		bool free_memory();

	protected:
		// operator to invert (surface level)
		AssembledLinearizedOperator<surface_function_type>* m_Op;

		IAssemble<level_function_type, algebra_type>& m_ass;
		approximation_space_type& m_approxSpace;
		phys_domain_type& m_domain;

		uint m_surfaceLevel;
		uint m_baseLevel;
		int m_cycle_type;

		smoother_type& m_smoother;
		int m_nu1;
		int m_nu2;

		base_solver_type& m_baseSolver;

		level_operator_type** m_A;
		projection_operator_type** m_P;
		prolongation_operator_type** m_I;

		level_function_type** m_u;
		level_function_type** m_c;
		level_function_type** m_d;
		level_function_type** m_t;

		// true -> allocate new matrices on every prepare
		bool m_grid_changes;
		bool m_allocated;
};

}

#include "mg_solver_impl.hpp"


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__ */
