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
#include "lib_algebra/operator/operator_iterator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/operator/operator_interface.h"

// library intern headers
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/operator/linear_operator/assembled_linear_operator.h"

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
class AssembledMultiGridCycle :
	virtual public ILinearIterator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	///	Approximation Space
		typedef TApproximationSpace approximation_space_type;

	///	Function type this preconditioner works on
		typedef typename TApproximationSpace::function_type function_type;

	///	DoFDistribution Type
		typedef typename TApproximationSpace::dof_distribution_type dof_distribution_type;

	/// Implementation type of DoFDistribution
		typedef typename dof_distribution_type::implementation_type dof_distribution_impl_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Level Operator Type
		typedef AssembledLinearOperator<dof_distribution_impl_type, algebra_type> operator_type;

	///	Prolongation Operator
		typedef IProlongationOperator<vector_type, vector_type> prolongation_operator_type;

	///	Projection Operator
		typedef IProjectionOperator<vector_type, vector_type> projection_operator_type;

	///	Base Solver type
		typedef ILinearOperatorInverse<vector_type, vector_type> base_solver_type;

	///	Smoother type
		typedef ILinearIterator<vector_type, vector_type> smoother_type;

	///	Own type
		typedef ILinearIterator<vector_type, vector_type> base_type;

	public:
	/// default Constructor
		AssembledMultiGridCycle() :
			m_pAss(NULL), m_pApproxSpace(NULL),
			m_topLev(0), m_baseLev(0), m_cycleType(1),
			m_numPreSmooth(1), m_numPostSmooth(1), m_pBaseSolver(NULL),
			m_grid_changes(false), m_allocated(false),
			m_bFullRefined(false),
			m_pDebugWriter(NULL), m_dbgIterCnt(0)
		{
			m_vSmoother.resize(1);
			m_vSmoother[0] = NULL;

			m_vProlongation.resize(1);
			m_vProlongation[0] = NULL;

			m_vProjection.resize(1);
			m_vProjection[0] = NULL;
		};

	///////////////////////////////////////////////////////////////////////////
	//	Setup
	///////////////////////////////////////////////////////////////////////////

	/// sets the assembling procedure that is used to compute coarse grid matrices
		void set_discretization(IAssemble<dof_distribution_impl_type, algebra_type>& ass)
			{m_pAss = &ass;}

	///	sets the approximation space that is used to build up the grid hierarchy
		void set_approximation_space(approximation_space_type& approxSpace)
			{m_pApproxSpace = &approxSpace;}

	///	sets the level where exact solving is performed in the mg cycle
		void set_base_level(int baseLevel) {m_baseLev = baseLevel;}

	///	sets the base solver that is used
		void set_base_solver(base_solver_type& baseSolver)
			{m_pBaseSolver = &baseSolver;}

	///	sets the smoother that is used
		void set_smoother(smoother_type& smoother) {m_vSmoother[0] = & smoother;}

	///	sets the cycle type (1 = V-cycle, 2 = W-cycle, ...)
		void set_cycle_type(int type) {m_cycleType = type;}

	///	sets the number of pre-smoothing steps to be performed
		void set_num_presmooth(int num) {m_numPreSmooth = num;}

	///	sets the number of post-smoothing steps to be performed
		void set_num_postsmooth(int num) {m_numPostSmooth = num;}

	///	sets the prolongation operator
		void set_prolongation_operator(IProlongationOperator<vector_type, vector_type>& P)
			{m_vProlongation[0] = &P;}

	///	sets the projection operator
		void set_projection_operator(IProjectionOperator<vector_type, vector_type>& P)
			{m_vProjection[0] = &P;}

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Geometric MultiGrid";}

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(ILinearOperator<vector_type, vector_type>& L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		base_type* clone();

	///	Destructor
		~AssembledMultiGridCycle();

 	protected:
 	/// smooth on level l, restrict defect, call lmgc (..., l-1) and interpolate correction
		bool lmgc(size_t lev);

	/// performs smoothing on level l, nu times
		bool smooth(function_type& d, function_type& c, size_t lev, int nu);

	///	initializes common part
		bool init_common(bool nonlinear);

	///	initializes the smoother and base solver
		bool init_smoother_and_base_solver();

	///	allocates the memory
		bool allocate_memory();

	///	frees the memory
		bool free_memory();

	protected:
	/// operator to invert (surface grid)
		operator_type* m_Op;

	///	assembling routine for coarse grid matrices
		IAssemble<dof_distribution_impl_type, algebra_type>* m_pAss;

	///	approximation space for level and surface grid
		approximation_space_type* m_pApproxSpace;

	///	top level (i.e. highest level in hierarchy. This is the surface level
	///	in case of non-adaptive refinement)
		size_t m_topLev;

	///	base level (where exact inverse is computed)
		size_t m_baseLev;

	///	cylce type (1 = V-cycle, 2 = W-cylcle, ...)
		int m_cycleType;

	///	number of Presmooth steps
		int m_numPreSmooth;

	///	number of Postsmooth steps
		int m_numPostSmooth;

	///	smoothing iterator for every grid level
		std::vector<smoother_type*> m_vSmoother;

	///	base solver for the coarse problem
		base_solver_type* m_pBaseSolver;

	///	coarse grid operator for each grid level
		std::vector<operator_type*> m_A;

	///	projection operator between grid levels
		std::vector<projection_operator_type*> m_vProjection;

	///	prolongation/restriction operator between grid levels
		std::vector<prolongation_operator_type*> m_vProlongation;

		std::vector<function_type*> m_u;
		std::vector<function_type*> m_c;
		std::vector<function_type*> m_d;
		std::vector<function_type*> m_t;

		// true -> allocate new matrices on every prepare
		bool m_grid_changes;
		bool m_allocated;

	///	flag indicating if grid is full refined
		bool m_bFullRefined;

#ifdef UG_PARALLEL
	/// communicator
		pcl::ParallelCommunicator<IndexLayout> m_Com;
#endif

	public:
	///	set debug output
	/**
	 * If a DebugWriter is passed by this method, the multi grid cycle writes
	 * the level/surface vectors and matrices for debug purposes.
	 *
	 * \param[in]	debugWriter		Debug Writer to use
	 */
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	protected:
	///	writes debug output for a level vector
	/**
	 * This method writes the level vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the vector
	 */
		bool write_level_debug(const vector_type& vec, const char* filename, size_t lev);

	///	writes debug output for a level matrix
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		bool write_level_debug(const matrix_type& mat, const char* filename, size_t lev);

	///	writes debug output for a surface vector
	/**
	 * This method writes the surface vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 */
		bool write_surface_debug(const vector_type& vec, const char* filename);

	///	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;

	///	counter for debug, to distinguish the iterations
		int m_dbgIterCnt;
};

} // end namespace ug

// include implementation
#include "mg_solver_impl.hpp"

#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER__ */
