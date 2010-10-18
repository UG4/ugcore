/*
 * mg_solver.h
 *
 *  Created on: 07.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__

// extern includes
#include <vector>
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/lib_discretization.h"

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
class AssembledMultiGridCycle :
	virtual public ILinearIterator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	//	Approximation Space
		typedef TApproximationSpace approximation_space_type;

	//	Function type this preconditioner works on
		typedef typename TApproximationSpace::function_type function_type;

	//	DoFDistribution Type
		typedef typename TApproximationSpace::dof_distribution_type dof_distribution_type;

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename algebra_type::vector_type vector_type;

	//	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	//	Level Operator Type
		typedef AssembledLinearOperator<dof_distribution_type, algebra_type> operator_type;

	//	Prolongation Operator
		typedef IProlongationOperator<vector_type, vector_type> prolongation_operator_type;

	//	Projection Operator
		typedef IProjectionOperator<vector_type, vector_type> projection_operator_type;

	private:
	//	Base Solver type
		typedef ILinearOperatorInverse<vector_type, vector_type> base_solver_type;

	//	Smoother type
		typedef ILinearIterator<vector_type, vector_type> smoother_type;

	public:
	// 	Constructor
		AssembledMultiGridCycle(	IAssemble<dof_distribution_type, algebra_type>& ass, approximation_space_type& approxSpace,
									size_t surfaceLevel, size_t baseLevel, int cycle_type,
									smoother_type& smoother, int nu1, int nu2, base_solver_type& baseSolver, bool grid_changes = true) :
			m_pAss(&ass), m_pApproxSpace(&approxSpace),
			m_surfaceLevel(surfaceLevel), m_baseLevel(baseLevel), m_cycle_type(cycle_type),
			m_nu1(nu1), m_nu2(nu2), m_pBaseSolver(&baseSolver),
			m_grid_changes(grid_changes), m_allocated(false)
		{
			m_vSmoother.resize(1);
			m_vSmoother[0] = &smoother;
		};

	// 	Constructor
		AssembledMultiGridCycle() :
			m_pAss(NULL), m_pApproxSpace(NULL),
			m_surfaceLevel(0), m_baseLevel(0), m_cycle_type(1),
			m_nu1(1), m_nu2(1), m_pBaseSolver(NULL),
			m_grid_changes(false), m_allocated(false)
		{
			m_vSmoother.resize(1);
			m_vSmoother[0] = NULL;

			m_vProlongation.resize(1);
			m_vProlongation[0] = NULL;

			m_vProjection.resize(1);
			m_vProjection[0] = NULL;
		};

	// 	Setup
		void set_discretization(IAssemble<dof_distribution_type, algebra_type>& ass) {m_pAss = &ass;}
		void set_approximation_space(approximation_space_type& approxSpace) {m_pApproxSpace = &approxSpace;}
		void set_surface_level(int surfLevel) {m_surfaceLevel = surfLevel;}
		void set_base_level(int baseLevel) {m_baseLevel = baseLevel;}
		void set_base_solver(base_solver_type& baseSolver) {m_pBaseSolver = &baseSolver;}
		void set_smoother(smoother_type& smoother) {m_vSmoother[0] = & smoother;}
		void set_cycle_type(int type) {m_cycle_type = type;}
		void set_num_presmooth(int num) {m_nu1 = num;}
		void set_num_postsmooth(int num) {m_nu2 = num;}
		void set_prolongation_operator(IProlongationOperator<vector_type, vector_type>& P) {m_vProlongation[0] = &P;}
		void set_projection_operator(IProjectionOperator<vector_type, vector_type>& P) {m_vProjection[0] = &P;}

	// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u);

	//	Prepare for Linear Operartor L
		virtual bool init(ILinearOperator<vector_type, vector_type>& L);

	//	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	//	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	//	Clone
		ILinearIterator<vector_type,vector_type>* clone()
		{
			AssembledMultiGridCycle<TApproximationSpace, TAlgebra>* clone =
				new AssembledMultiGridCycle<TApproximationSpace, TAlgebra>(	*m_pAss, *m_pApproxSpace,
																			m_surfaceLevel, m_baseLevel, m_cycle_type,
																			*(m_vSmoother[0]), m_nu1, m_nu2, *m_pBaseSolver, m_grid_changes);

			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	//	Destructor
		~AssembledMultiGridCycle();

 	protected:
 		// smooth on level l, restrict defect, call lmgc (..., l-1) and interpolate correction
		bool lmgc(size_t lev);

		// smmoth help function: perform smoothing on level l, nu times
		bool smooth(function_type& d, function_type& c, size_t lev, int nu);

		bool init_common(bool nonlinear);
		bool init_smoother_and_base_solver();
		bool allocate_memory();
		bool free_memory();

	protected:
		// operator to invert (surface level)
		operator_type* m_Op;

		IAssemble<dof_distribution_type, algebra_type>* m_pAss;
		approximation_space_type* m_pApproxSpace;

		size_t m_surfaceLevel;
		size_t m_baseLevel;
		int m_cycle_type;

		int m_nu1;
		int m_nu2;

		std::vector<smoother_type*> m_vSmoother;
		base_solver_type* m_pBaseSolver;

		std::vector<operator_type*> m_A;
		std::vector<projection_operator_type*> m_vProjection;
		std::vector<prolongation_operator_type*> m_vProlongation;

		std::vector<function_type*> m_u;
		std::vector<function_type*> m_c;
		std::vector<function_type*> m_d;
		std::vector<function_type*> m_t;

		// true -> allocate new matrices on every prepare
		bool m_grid_changes;
		bool m_allocated;

#ifdef UG_PARALLEL
		// communicator
		pcl::ParallelCommunicator<IndexLayout> m_Com;
#endif
};

}

#include "mg_solver_impl.hpp"


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__ */
